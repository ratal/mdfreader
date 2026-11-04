"""Export functions for MDF data into various file formats (CSV, NetCDF, HDF5, MATLAB, Excel, Parquet, pandas)."""
from __future__ import annotations

import os
from os.path import splitext
from typing import Any
from warnings import warn

import numpy as np
from numpy import array, datetime64, float32, float64, vstack
from numpy.ma import MaskedArray

from .mdf import (
    dataField,
    descriptionField,
    masterField,
    masterTypeField,
    unitField,
)


def _convert_to_matlab_name(channel: str) -> str:
    """Removes disallowed characters for a MATLAB variable name.

    Parameters
    ----------
    channel : str
        Original channel name.

    Returns
    -------
    str
        MATLAB-compatible identifier (up to 63 characters).
    """
    channel_name = channel.replace('[', '_ls_')
    channel_name = channel_name.replace(']', '_rs_')
    channel_name = channel_name.replace('$', '')
    channel_name = channel_name.replace('.', 'p')
    channel_name = channel_name.replace('\\', '_bs_')
    channel_name = channel_name.replace('/', '_fs_')
    channel_name = channel_name.replace('(', '_lp_')
    channel_name = channel_name.replace(')', '_rp_')
    channel_name = channel_name.replace(',', '_c_')
    channel_name = channel_name.replace('@', '_am_')
    channel_name = channel_name.replace(' ', '_')
    channel_name = channel_name.replace(':', '_co_')
    channel_name = channel_name.replace('-', '_hy_')

    allowed_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.')
    channel_name = ''.join(c for c in channel_name if c in allowed_chars)
    return channel_name[:63]


def _convert_to_hdf5_name(channel: str) -> str:
    """Removes path separators for HDF5 dataset name compatibility.

    Parameters
    ----------
    channel : str
        Channel name.

    Returns
    -------
    str
        HDF5-compatible name.
    """
    return channel.replace('/', '_fs_').replace('\\', '_bs_')


def export_csv(mdf: Any, file_name: str | None = None, sampling: float | None = None) -> str | None:
    """Export MDF data into a CSV file.

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target CSV file path. Defaults to mdf.fileName with .csv extension.
    sampling : float, optional
        Resampling interval before exporting.

    Returns
    -------
    str or None
        Path to the generated CSV file.
    """
    if not mdf:
        warn('no data to be exported')
        return None

    import csv

    if sampling is not None:
        mdf.resample(sampling)

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.csv"

    encoding = 'utf8' if getattr(mdf, 'MDFVersionNumber', 400) >= 400 else 'latin-1'

    with open(file_name, 'wt', encoding=encoding, newline='') as f:
        writer = csv.writer(f, dialect=csv.excel, escapechar='\\')
        # Header: channel names
        valid_channels = [
            name for name in mdf
            if mdf.get_channel_data(name).dtype.kind not in ('S', 'U', 'V')
            and mdf.get_channel_data(name).ndim <= 1
        ]
        writer.writerow(valid_channels)
        # Units
        writer.writerow([mdf.get_channel_unit(name) or '' for name in valid_channels])

        temp = []
        ref_len = None
        for name in mdf:
            data = mdf.get_channel_data(name)
            if data.dtype.kind not in ('S', 'U', 'V') and data.ndim == 1:
                if ref_len is None:
                    ref_len = len(data)
                if len(data) == ref_len:
                    temp.append(data.transpose())
        if temp:
            buf = vstack(temp).transpose()
            r, _ = buf.shape
            chunk_size = 500
            for i in range(0, r, chunk_size):
                writer.writerows(buf[i:i + chunk_size].tolist())

    return file_name


def export_netcdf(mdf: Any, file_name: str | None = None, sampling: float | None = None) -> str | None:
    """Export MDF data into a NetCDF 3 file.

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target NetCDF file path. Defaults to mdf.fileName with .nc extension.
    sampling : float, optional
        Resampling interval before exporting.

    Returns
    -------
    str or None
        Path to the generated NetCDF file.
    """
    try:
        from scipy.io import netcdf_file
    except ImportError:
        warn('scipy.io module not found')
        return None

    def clean_name(name: str) -> str:
        allowed_str = ' ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-+_.@'
        return ''.join(c for c in name if c in allowed_str)

    def set_attribute(f: Any, name: str, value: Any) -> None:
        if value is not None and len(value) > 0:
            if isinstance(value, dict) and 'name' in value:
                value = value['name']
            if isinstance(value, bytes):
                value = value.decode('utf-8', 'ignore')
            setattr(f, name, clean_name(str(value)))

    if sampling is not None:
        mdf.resample(sampling)

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.nc"

    f = netcdf_file(file_name, 'w')
    setattr(f, 'Time', mdf.fileMetadata.get('time', 0))
    set_attribute(f, 'Author', mdf.fileMetadata.get('author'))
    set_attribute(f, 'Organization', mdf.fileMetadata.get('organisation'))
    set_attribute(f, 'ProjectName', mdf.fileMetadata.get('project'))
    set_attribute(f, 'Subject', mdf.fileMetadata.get('subject'))
    set_attribute(f, 'Comment', mdf.fileMetadata.get('comment'))

    for master in mdf.masterChannelList:
        channels_in_group = mdf.masterChannelList[master]
        if channels_in_group:
            f.createDimension(master, len(mdf.get_channel_data(channels_in_group[0])))

    var: dict[str, Any] = {}
    master_channel_set = set(mdf.masterChannelList.keys())
    for name in mdf:
        data = mdf.get_channel_data(name)
        if data.dtype == 'float64':
            data_type = 'd'
        elif data.dtype == 'float32':
            data_type = 'f'
        elif data.dtype in {'int8', 'int16', 'uint8', 'uint16'}:
            data_type = 'h'
        elif data.dtype in {'int32', 'uint32'}:
            data_type = 'i'
        elif data.dtype.kind in {'S', 'U'}:
            data_type = 'c'
        else:
            data_type = None
            warn(f'Can not process numpy type {data.dtype} of channel {name}')

        if data_type is not None:
            cleaned_name = clean_name(name)
            if len(mdf.masterChannelList) == 1:
                var[name] = f.createVariable(
                    cleaned_name, data_type, (next(iter(mdf.masterChannelList)),))
            else:
                var[name] = f.createVariable(
                    cleaned_name, data_type, (mdf.get_channel_master(name),))

            set_attribute(var[name], 'title', cleaned_name)
            set_attribute(var[name], 'units', mdf.get_channel_unit(name))
            set_attribute(var[name], 'Description', mdf.get_channel_desc(name))
            if name in master_channel_set:
                set_attribute(var[name], 'Type', 'Master Channel')
                set_attribute(var[name], 'datatype', 'master')
            else:
                set_attribute(var[name], 'Type', 'Data Channel')

    for name in mdf:
        if name in var:
            var[name] = mdf.get_channel_data(name)
    f.close()
    return file_name


def export_hdf5(
    mdf: Any,
    file_name: str | None = None,
    sampling: float | None = None,
    compression: str | None = None,
    compression_opts: int | None = None,
) -> str | None:
    """Export MDF data structure into an HDF5 file.

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target HDF5 file path. Defaults to mdf.fileName with .hdf extension.
    sampling : float, optional
        Resampling interval before exporting.
    compression : str, optional
        Compression filter ('gzip', 'lzf').
    compression_opts : int, optional
        gzip compression level (0-9).

    Returns
    -------
    str or None
        Path to the generated HDF5 file.
    """
    try:
        import h5py
    except ImportError:
        warn('h5py not found')
        return None

    def set_attribute(obj: Any, name: str, value: Any) -> None:
        if value is not None and len(value) > 0:
            try:
                if isinstance(value, dict) and 'name' in value:
                    value = value['name']
                obj.attrs[name] = value
            except Exception:
                pass

    if sampling is not None:
        mdf.resample(sampling)

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.hdf"

    if compression is not None:
        compression = compression.lower()
        if compression not in {'gzip', 'lzf'}:
            compression = None
            compression_opts = None
        elif compression == 'lzf':
            compression_opts = None
        if compression_opts not in range(10):
            compression_opts = None
    else:
        compression_opts = None

    f = h5py.File(file_name, 'w')
    file_group = f.create_group(os.path.basename(file_name))
    file_group.attrs['Time'] = mdf.fileMetadata.get('time', 0)
    set_attribute(file_group, 'Author', mdf.fileMetadata.get('author'))
    set_attribute(file_group, 'Organization', mdf.fileMetadata.get('organisation'))
    set_attribute(file_group, 'ProjectName', mdf.fileMetadata.get('project'))
    set_attribute(file_group, 'Subject', mdf.fileMetadata.get('subject'))
    set_attribute(file_group, 'Comment', mdf.fileMetadata.get('comment'))

    master_type_dict = {
        0: 'None', 1: 'Time', 2: 'Angle', 3: 'Distance', 4: 'Index', None: 'None',
    }

    chunks = True if compression is not None else None

    if len(mdf.masterChannelList) > 1:
        groups: dict[str, int] = {}
        n_groups = 0
        grp: dict[int, Any] = {}
        for channel in mdf:
            channel_data = mdf.get_channel_data(channel)
            master_name = mdf.get_channel_master(channel)
            if masterField in mdf[channel] and master_name not in groups:
                n_groups += 1
                group_name = master_name if master_name else f'{masterField}{n_groups}'
                groups[group_name] = n_groups
                grp[n_groups] = file_group.create_group(group_name)
                set_attribute(grp[n_groups], masterField, master_name)
                set_attribute(grp[n_groups], masterTypeField,
                              master_type_dict.get(mdf.get_channel_master_type(channel), 'None'))
            elif masterField in mdf[channel] and master_name in groups:
                group_name = master_name

            if channel_data.dtype.kind not in ('U', 'O') and channel_data.dtype.itemsize > 0:
                channel_name = _convert_to_hdf5_name(channel)
                dset = grp[groups[group_name]].create_dataset(
                    channel_name,
                    data=channel_data,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                unit = mdf.get_channel_unit(channel)
                if unit:
                    set_attribute(dset, unitField, unit)
                if descriptionField in mdf[channel]:
                    desc = mdf.get_channel_desc(channel)
                    if desc:
                        set_attribute(dset, descriptionField, desc)
    else:
        master_name = next(iter(mdf.masterChannelList)) if mdf.masterChannelList else ''
        set_attribute(file_group, masterField, master_name)
        set_attribute(file_group, masterTypeField,
                      master_type_dict.get(mdf.get_channel_master_type(master_name), 'None'))
        for channel in mdf:
            channel_data = mdf.get_channel_data(channel)
            if channel_data.dtype.kind not in ('U', 'O') and channel_data.dtype.itemsize > 0:
                channel_name = _convert_to_hdf5_name(channel)
                dset = file_group.create_dataset(
                    channel_name,
                    data=channel_data,
                    compression=compression,
                    compression_opts=compression_opts,
                    chunks=chunks,
                )
                unit = mdf.get_channel_unit(channel)
                if unit:
                    set_attribute(dset, unitField, unit)
                if descriptionField in mdf[channel]:
                    desc = mdf.get_channel_desc(channel)
                    if desc:
                        set_attribute(dset, descriptionField, desc)

    f.close()
    return file_name


def export_matlab(mdf: Any, file_name: str | None = None) -> str | None:
    """Export MDF data into a MATLAB .mat file (format 7.3 preferred).

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target .mat file path. Defaults to mdf.fileName with .mat extension.

    Returns
    -------
    str or None
        Path to the generated .mat file.
    """
    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.mat"

    # Fast direct HDF5 MAT 7.3 writer if h5py is installed
    try:
        import h5py
        matlab_class_map = {
            'f8': b'double', 'f4': b'single', 'i8': b'int64', 'u8': b'uint64',
            'i4': b'int32', 'u4': b'uint32', 'i2': b'int16', 'u2': b'uint16',
            'i1': b'int8', 'u1': b'uint8', 'b1': b'logical',
        }
        with h5py.File(file_name, 'w', userblock_size=512) as f:
            g_mcl = f.create_group('masterChannelList')
            g_mcl.attrs['MATLAB_class'] = np.bytes_(b'struct')
            g_masters = {}
            for master in mdf.masterChannelList:
                master_mat_name = _convert_to_matlab_name(master) or master
                gm = g_mcl.create_group(master_mat_name)
                gm.attrs['MATLAB_class'] = np.bytes_(b'struct')
                g_masters[master] = gm

            for channel in mdf:
                data = mdf.get_channel_data(channel)
                if data.dtype.kind not in ('S', 'U', 'V'):
                    channel_name = _convert_to_matlab_name(channel)
                    if not channel_name:
                        warn(f'Could not export {channel}, name is not compatible with Matlab')
                        continue
                    if getattr(data.dtype, 'base', None) is not None and data.dtype.base.name == 'float16':
                        data = data.astype(float32)
                    if isinstance(data, MaskedArray):
                        data = data.filled()
                    data_view = data.reshape(-1, 1) if data.ndim == 1 else data
                    dset = f.create_dataset(channel_name, data=data_view)
                    cls_name = matlab_class_map.get(data.dtype.str[1:], b'double')
                    dset.attrs['MATLAB_class'] = np.bytes_(cls_name)

                    unit = mdf.get_channel_unit(channel)
                    desc = mdf.get_channel_desc(channel)
                    master = mdf.get_channel_master(channel)
                    if unit or desc or master == channel:
                        gm = g_masters.get(master)
                        if gm is not None:
                            gc = gm.create_group(channel_name)
                            gc.attrs['MATLAB_class'] = np.bytes_(b'struct')
                            if desc:
                                ds_desc = gc.create_dataset(
                                    'description',
                                    data=np.array([[ord(c) for c in desc]], dtype='uint16'),
                                )
                                ds_desc.attrs['MATLAB_class'] = np.bytes_(b'char')
                                ds_desc.attrs['MATLAB_int_decode'] = np.int64(2)
                            if unit:
                                ds_unit = gc.create_dataset(
                                    'unit',
                                    data=np.array([[ord(c) for c in unit]], dtype='uint16'),
                                )
                                ds_unit.attrs['MATLAB_class'] = np.bytes_(b'char')
                                ds_unit.attrs['MATLAB_int_decode'] = np.int64(2)
                            if master == channel:
                                mt = mdf.get_channel_master_type(channel)
                                ds_mt = gc.create_dataset(
                                    'masterType',
                                    data=np.array([[mt if mt is not None else 0]], dtype='int64'),
                                )
                                ds_mt.attrs['MATLAB_class'] = np.bytes_(b'int64')

        hdr = bytearray(512)
        text = b'MATLAB 7.3 MAT-file, Platform: POSIX, Created by mdfreader HDF5 schema 0.5 .'
        hdr[:len(text)] = text
        hdr[124:128] = b'\x00\x02IM'
        with open(file_name, 'r+b') as fb:
            fb.seek(0)
            fb.write(hdr)

        return file_name
    except Exception:
        pass

    try:
        from hdf5storage import savemat
    except ImportError:
        warn('hdf5storage module not found')
        try:
            from scipy.io import savemat
        except ImportError:
            warn('scipy module also not found')
            return None

    temp: dict[str, Any] = {'masterChannelList': {}}
    for master in mdf.masterChannelList:
        temp['masterChannelList'][master] = {}

    for channel in mdf:
        data = mdf.get_channel_data(channel)
        if data.dtype.kind not in ('S', 'U', 'V'):
            channel_name = _convert_to_matlab_name(channel)
            if channel_name:
                temp[channel_name] = data
                if getattr(data.dtype, 'base', None) is not None and data.dtype.base.name == 'float16':
                    temp[channel_name] = data.astype(float32)
                if isinstance(data, MaskedArray):
                    temp[channel_name] = data.filled()
                unit = mdf.get_channel_unit(channel)
                desc = mdf.get_channel_desc(channel)
                master = mdf.get_channel_master(channel)
                if unit or desc or master == channel:
                    if master not in temp['masterChannelList']:
                        temp['masterChannelList'][master] = {}
                    temp['masterChannelList'][master][channel_name] = {}
                    if desc:
                        temp['masterChannelList'][master][channel_name]['description'] = desc
                    if unit:
                        temp['masterChannelList'][master][channel_name]['unit'] = unit
                    if master == channel:
                        temp['masterChannelList'][master][channel_name]['masterType'] = (
                            mdf.get_channel_master_type(channel)
                        )
            else:
                warn(f'Could not export {channel}, name is not compatible with Matlab')

    try:
        savemat(file_name, temp, format='7.3', long_field_names=True, oned_as='column',
                structured_numpy_ndarray_as_struct=True)
    except Exception:
        savemat(file_name, temp, long_field_names=True, format='5')

    return file_name


def export_excel(mdf: Any, file_name: str | None = None) -> str | None:
    """Export MDF data into Excel 95-2003 (.xls) file.

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target .xls file path. Defaults to mdf.fileName with .xls extension.

    Returns
    -------
    str or None
        Path to the generated Excel file.
    """
    try:
        import xlwt
    except ImportError:
        warn('xlwt module missing')
        return None

    from math import ceil

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.xls"

    wb = xlwt.Workbook(encoding='utf-8')
    channel_list = list(mdf.keys())
    max_cols = 255
    max_lines = 65535
    workbook_number = int(ceil(len(channel_list) / max_cols))
    too_long_channels: list[str] = []

    for workbook in range(workbook_number):
        ws = wb.add_sheet(f'Sheet{workbook}')
        if workbook == workbook_number - 1:
            column_range = list(range(workbook * max_cols, len(channel_list)))
        else:
            column_range = list(range(workbook * max_cols, (workbook + 1) * max_cols))

        for col in column_range:
            col_idx = col - workbook * max_cols
            ws.write(0, col_idx, channel_list[col])
            ws.write(1, col_idx, mdf.get_channel_unit(channel_list[col]))
            vector = mdf.get_channel_data(channel_list[col])
            is_numeric = vector.dtype.kind not in {'S', 'U'}
            write_count = min(len(vector), max_lines)

            if len(vector) > max_lines:
                too_long_channels.append(channel_list[col])

            for row in range(write_count):
                if is_numeric:
                    ws.row(row + 2).set_cell_number(col_idx, vector[row])
                else:
                    ws.row(row + 2).set_cell_text(col_idx, str(vector[row]))

    wb.save(file_name)
    if too_long_channels:
        warn(f'Following channels were too long to be processed completely, '
             f'maybe you should resample: {too_long_channels}')

    return file_name


def export_xlsx(mdf: Any, file_name: str | None = None) -> str | None:
    """Export MDF data into Excel OpenXML (.xlsx) file.

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target .xlsx file path. Defaults to mdf.fileName with .xlsx extension.

    Returns
    -------
    str or None
        Path to the generated xlsx file.
    """
    try:
        import openpyxl
        from openpyxl.cell.cell import ILLEGAL_CHARACTERS_RE
        from openpyxl.utils.exceptions import IllegalCharacterError
    except ImportError:
        warn('Module openpyxl missing')
        return None

    def _clean_str(val: Any) -> Any:
        if isinstance(val, str):
            return ILLEGAL_CHARACTERS_RE.sub('', val)
        return val

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.xlsx"

    warn('Creating Excel sheet')
    wb = openpyxl.workbook.Workbook()
    ws = wb.active

    n_col = 1
    for master in mdf.masterChannelList:
        channels = mdf.masterChannelList[master]
        for channel in channels:
            data = mdf.get_channel_data(channel)
            if data.ndim <= 1:
                try:
                    ws.cell(row=1, column=n_col).value = _clean_str(channel)
                    ws.cell(row=2, column=n_col).value = _clean_str(mdf.get_channel_unit(channel))
                    if data.dtype.kind in ('i', 'u') or 'f4' in data.dtype.str:
                        for r, cell_data in enumerate(data):
                            ws.cell(row=r + 3, column=n_col).value = float64(cell_data)
                    elif data.dtype.kind != 'V':
                        for r, cell_data in enumerate(data):
                            cell_val = cell_data
                            if isinstance(cell_val, str):
                                cell_val = _clean_str(cell_val)
                            ws.cell(row=r + 3, column=n_col).value = cell_val
                    n_col += 1
                except IllegalCharacterError:
                    warn(f'could not export {channel}')

    warn('Writing file, please wait')
    wb.save(file_name)
    return file_name


def return_pandas_dataframe(mdf: Any, master_channel_name: str) -> Any:
    """Return a pandas DataFrame for a specific raster described by master channel name.

    Parameters
    ----------
    mdf : Mdf
        MDF object.
    master_channel_name : str
        Name of master channel for the desired raster group.

    Returns
    -------
    pandas.DataFrame or None
    """
    try:
        import pandas as pd
    except ImportError:
        warn('Module pandas missing')
        return None

    if master_channel_name not in mdf:
        warn('Master channel name not in mdf')
        return None

    group_channels = mdf.masterChannelList.get(master_channel_name, [])
    if master_channel_name not in group_channels:
        warn(f'no master channel in group {master_channel_name}')

    if master_channel_name in group_channels:
        if mdf.get_channel_master_type(master_channel_name) == 1:
            epoch_ns = int(mdf.fileMetadata.get('time', 0) * 1e9)
            datetime_info = datetime64(epoch_ns, 'ns')
            time_idx = datetime_info + array(
                mdf.get_channel_data(master_channel_name) * 1e9, dtype='timedelta64[ns]'
            )
            df = pd.DataFrame(index=time_idx)
        else:
            df = pd.DataFrame(index=mdf.get_channel_data(master_channel_name))
    else:
        df = pd.DataFrame()

    channel_dict: dict[str, Any] = {}
    for key in group_channels:
        data = mdf.get_channel_data(key)
        if data.dtype.byteorder not in ['=', '|']:
            data = data.byteswap().view(data.dtype.newbyteorder())
        if data.ndim == 1 and data.shape[0] == df.shape[0] and data.dtype.char != 'V':
            channel_dict[key] = data

    return pd.DataFrame(data=channel_dict, index=df.index)


def get_master_channel_datetimes(mdf: Any, master_channel_name: str) -> np.ndarray | None:
    """Return a time master channel as absolute datetime64[ns] array.

    Parameters
    ----------
    mdf : Mdf
        MDF object.
    master_channel_name : str
        Name of a time master channel.

    Returns
    -------
    numpy.ndarray of datetime64[ns], or None.
    """
    if master_channel_name not in mdf or mdf.get_channel_master_type(master_channel_name) != 1:
        return None

    t = mdf.fileMetadata.get('time', 0)
    offset_seconds = mdf.fileMetadata.get('time_offset') or 0
    epoch_ns = int((t + offset_seconds) * 1e9)
    relative_ns = array(mdf.get_channel_data(master_channel_name) * 1e9, dtype='timedelta64[ns]')
    return datetime64(epoch_ns, 'ns') + relative_ns


def to_pandas(mdf: Any, sampling: float | None = None) -> None:
    """Convert MDF data structure in-place into pandas DataFrame(s).

    Parameters
    ----------
    mdf : Mdf
        MDF object to convert.
    sampling : float, optional
        Resampling interval before conversion.
    """
    try:
        import pandas as pd  # noqa: F401
    except ImportError:
        warn('Module pandas missing')
        return

    if sampling is not None:
        mdf.resample(sampling)

    groups = list(mdf.masterChannelList.keys())
    for group in groups:
        mdf[f'{group}_group'] = return_pandas_dataframe(mdf, group)
        for channel in mdf.masterChannelList.get(group, []):
            mdf[channel].pop(dataField, None)
            mdf[channel].pop(masterField, None)

    mdf.masterGroups = [f'{g}_group' for g in groups]
    mdf.masterChannelList = {}
    mdf._pandasframe = True


def export_parquet(mdf: Any, file_name: str | None = None) -> None:
    """Export MDF data into Parquet file(s).

    Parameters
    ----------
    mdf : Mdf
        MDF object to export.
    file_name : str, optional
        Target Parquet file path.
    """
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        warn('pyarrow not installed, cannot export to parquet')
        return

    if file_name is None:
        file_name = f"{splitext(mdf.fileName)[0]}.parquet"

    masters = list(mdf.masterChannelList.keys())
    base, ext = splitext(file_name)
    for i, master_channel_name in enumerate(masters):
        frame = return_pandas_dataframe(mdf, master_channel_name)
        if frame is None:
            continue
        out = file_name if len(masters) == 1 else f'{base}_{i}{ext}'
        frame.to_parquet(out, engine='pyarrow', compression='gzip')

