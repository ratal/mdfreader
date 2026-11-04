# -*- coding: utf-8 -*-
""" Measured Data Format file reader main module

Platform and python version
----------------------------------------
With Unix, macOS and Windows for python 3.9+

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Sun Oct 10 12:57:28 2010

Dependencies
-------------------
- Python >=3.9 <https://www.python.org>
- Numpy >1.14 <https://numpy.org>
- Sympy to convert channels with formula
- bitarray for not byte aligned data parsing (only without the compiled dataRead extension)
- Matplotlib <https://matplotlib.org> and mplcursors for plotting
- scipy for NetCDF
- h5py for the HDF5 export
- xlwt for the excel 95-2003 export
- openpyxl >2.0 for the excel 2007 export
- hdf5storage for the Matlab file conversion
- pyzstd and lz4 for ZStandard / LZ4 compressed data blocks
- zlib to uncompress data block if needed

mdfreader
--------------------------
"""

from __future__ import annotations

from argparse import ArgumentParser
from io import open
from math import ceil
from os import name as osname
from os import remove
from os.path import splitext
from struct import unpack
from typing import Any, Sequence
from warnings import warn

from numpy import (
    all,
    arange,
    array,
    clip,
    datetime64,
    diff,
    empty,
    float32,
    float64,
    interp,
    linspace,
    mean,
    nan,
    searchsorted,
    vstack,
)
from numpy.ma import MaskedArray, empty as ma_empty, masked

from .export import (
    _convert_to_hdf5_name,
    _convert_to_matlab_name,
    export_csv,
    export_excel,
    export_hdf5,
    export_matlab,
    export_netcdf,
    export_parquet,
    export_xlsx,
    get_master_channel_datetimes,
    return_pandas_dataframe,
    to_pandas,
)
from .mdf import (
    CompressedData,
    MdfSkeleton,
    _open_mdf,
    conversionField,
    dataField,
    descriptionField,
    idField,
    invalidChannel,
    invalidPosField,
    masterField,
    masterTypeField,
    unitField,
)
from .mdf3reader import Mdf3
from .mdf4reader import Mdf4
from .mdfinfo3 import Info3, _generate_dummy_mdf3
from .mdfinfo4 import Info4, _generate_dummy_mdf4


class MdfInfo(dict):
    __slots__ = ['fileName', 'fid', 'zipfile',
                 'mdfversion', 'filterChannelNames']
    """ MDFINFO is a class gathering information from block headers in a MDF (Measure Data Format) file.
    Structure is nested dicts. Primary key is Block type, then data group, channel group and channel number.
    Examples of dicts
    - mdfinfo['HDBlock'] header block
    - mdfinfo['DGBlock'][dataGroup] Data Group block
    - mdfinfo['CGBlock'][dataGroup][channelGroup] Channel Group block
    - mdfinfo['CNBlock'][dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    - mdfinfo['CCBlock'][dataGroup][channelGroup][channel] Channel conversion information

    Attributes
    --------------
    fileName : str
        file name
    mdfversion : int
        mdf file version number
    filterChannelNames : bool
        flag to filter long channel names including module names separated by a '.'
    fid
        file identifier
    zipfile
        flag to indicate the mdf4 is packaged in a zip

    Methods
    ------------
    read_info( fileName = None, filterChannelNames=False )
        Reads MDF file and extracts its complete structure
    list_channels( file_name = None )
        Read MDF file blocks and returns a list of contained channels

    Examples
    --------------
    >>> import mdfreader
    >>> FILENAME='toto.dat'
    >>> yop=mdfreader.MdfInfo(FILENAME)
    or if you are just interested to have only list of channels
    >>> yop=mdfreader.MdfInfo() # creates new instance of mdfinfo class
    >>> yop.list_channels(FILENAME) # returns a simple list of channel names
    """

    def __init__(
        self,
        file_name: str | None = None,
        filter_channel_names: bool = False,
        fid: Any = None,
        minimal: int = 0,
    ) -> None:
        """ You can give optionally to constructor a file name that will be parsed

        Parameters
        ----------------
        file_name : str, optional
            file name
        filter_channel_names : bool, optional
            flag to filter long channel names including module names separated by a '.'
        fid : file identifier, optional
        minimal : int, optional
        """
        self.fileName = file_name
        self.filterChannelNames = filter_channel_names
        self.mdfversion = 410
        self.fid = fid
        self.zipfile = False
        if file_name is not None:
            self.read_info(file_name, fid, minimal)

    def read_info(
        self,
        file_name: str | None = None,
        fid: Any = None,
        minimal: int = 0,
    ) -> None:
        """ Reads MDF file and extracts its complete structure

        Parameters
        ----------------
        file_name : str, optional
            file name. If not input, uses fileName attribute
        fid : file identifier, optional
        minimal : int
            0 will load every metadata
            1 will load DG, CG, CN and CC
            2 will load only DG
        """
        if self.fileName is None or file_name is not None:
            self.fileName = file_name

        # Open file
        if self.fid is None or (self.fid is not None and self.fid.closed):
            (self.fid, self.fileName, self.zipfile) = _open_mdf(self.fileName)

        # read Identifier block
        self.fid.seek(28)
        mdf_version_number = unpack('<H', self.fid.read(2))
        self.mdfversion = mdf_version_number[0]
        if self.mdfversion < 400:  # up to version 3.x not compatible with version 4.x
            self.update(
                Info3(None, self.fid, self.filterChannelNames, minimal))
        else:  # MDF version 4.x
            self.update(
                Info4(None, self.fid, self.filterChannelNames, minimal))
            if self.zipfile and fid is None:  # not from mdfreader.read()
                remove(self.fileName)

    def list_channels(self, file_name: str | None = None) -> list[str]:
        """ Read MDF file blocks and returns a list of contained channels

        Parameters
        ----------------
        file_name : string
            file name

        Returns
        -----------
        nameList : list of string
            list of channel names
        """
        if self.fileName is None or file_name is not None:
            self.fileName = file_name
        # Open file
        (self.fid, self.fileName, zipfile) = _open_mdf(self.fileName)
        # read Identifier block
        self.fid.seek(28)
        mdf_version_number = unpack('<H', self.fid.read(2))
        self.mdfversion = mdf_version_number[0]
        if self.mdfversion < 400:  # up to version 3.x not compatible with version 4.x
            channel_name_list = Info3()
            name_list = channel_name_list.list_channels3(
                self.fileName, self.fid)
        else:
            channel_name_list = Info4()
            name_list = channel_name_list.list_channels4(
                self.fileName, self.fid)
            if zipfile:  # not from mdfreader.read()
                remove(self.fileName)
        return name_list

    def _generate_dummy_mdf(self, channel_list: Sequence[str] | None = None) -> Any:
        """ Parse MDF file structure and create a dummy mdf object structure

        Parameters
        ----------------
        channel_list : Sequence of str, optional
            list of channels
        """
        if self.mdfversion < 400:  # up to version 3.x not compatible with version 4.x
            return _generate_dummy_mdf3(self, channel_list)
        else:  # MDF version 4.x
            return _generate_dummy_mdf4(self, channel_list)


class Mdf(MdfSkeleton):
    """ Mdf class

    Attributes
    --------------
    fileName : str
        file name
    MDFVersionNumber : int
        mdf file version number
    masterChannelList : dict
        Represents data structure: a key per master channel with corresponding value containing a list of channels
        One key or master channel represents then a data group having same sampling interval.
    multiProc : bool
        Flag to request channel conversion multi processed for performance improvement.
        One thread per data group.
    fileMetadata : dict
        file metadata with minimum keys : author, organisation, project, subject, comment, time, date

    Methods
    ------------
    read( file_name = None, multi_processed = False, channel_list=None, convert_after_read=True,
    filter_channel_names=False, no_data_loading=False, compression=False)
        reads mdf file version 3.x and 4.x
    write( file_name=None )
        writes simple mdf file
    get_channel_data( channel_name )
        returns channel numpy array
    convert_all_channel()
        converts all channel data according to CCBlock information
    get_channel_unit( channel_name )
        returns channel unit
    plot( channels )
        Plot channels with Matplotlib
    resample( sampling_time = 0.1, master_channel=None )
        Resamples all data groups
    export_to_csv( file_name = None, sampling = 0.1 )
        Exports mdf data into CSV file
    export_to_NetCDF( file_name = None, sampling = None )
        Exports mdf data into netcdf file
    export_to_hdf5( file_name = None, sampling = None )
        Exports mdf class data structure into hdf5 file
    export_to_matlab( file_name = None )
        Exports mdf class data structure into Matlab file
    export_to_excel( file_name = None )
        Exports mdf data into excel 95 to 2003 file
    export_to_xlsx( file_name=None )
        Exports mdf data into excel 2007 and 2010 file
    convert_to_pandas( sampling=None )
        converts mdf data structure into pandas dataframe(s)
    keep_channels( channel_list )
        keeps only list of channels and removes the other channels
    merge_mdf( mdf_class ):
        Merges data of 2 mdf classes

    Notes
    --------
    mdf class is a nested dict.
    Channel name is the primary dict key of mdf class.
    At a higher level, each channel includes the following keys :
        - 'data' : containing vector of data (numpy)
        - 'unit' : unit (string)
        - 'master' : master channel of channel (time, crank angle, etc.)
        - 'description' : Description of channel
        - 'conversion': mdfinfo nested dict for CCBlock.
            Exist if channel not converted, used to convert with getChannelData method

    Examples
    --------------
    >>> import mdfreader
    >>> yop=mdfreader.Mdf('NameOfFile')
    >>> yop.keys() # list channels names
    # list channels grouped by raster or master channel
    >>> yop.masterChannelList
    >>> yop.plot('channelName') or yop.plot({'channel1','channel2'})
    >>> yop.resample(0.1) or yop.resample()
    >>> yop.export_to_csv(sampling=0.01)
    >>> yop.export_to_NetCDF()
    >>> yop.export_to_hdf5()
    >>> yop.export_to_matlab()
    >>> yop.export_to_excel()
    >>> yop.export_to_xlsx()
    >>> yop.convert_to_pandas() # converts data groups into pandas dataframes
    >>> yop.write() # writes mdf file
    # drops all the channels except the one in argument
    >>> yop.keep_channels(['channel1','channel2','channel3'])
    >>> yop.get_channel_data('channelName') # returns channel numpy array
    >>> yop=mdfreader.Mdf()  # create an empty Mdf object
    # add channel in Mdf object
    >>> yop.add_channel(channel_name, data, master_channel, master_type, unit='lumen', description='what you want')
    >>> yop.write('filename') # change version with yop.MDFVersionNumber or specifically use write3/4()
    """

    def __init__(
        self,
        file_name: str | None = None,
        channel_list: Sequence[str] | set[str] | None = None,
        convert_after_read: bool = True,
        filter_channel_names: bool = False,
        no_data_loading: bool = False,
        compression: bool | str | int = False,
        convert_tables: bool = True,
        metadata: int = 2,
        finalization_writing_to_file: bool = False,
        force_file_integrity_check: bool = False,
    ) -> None:
        super().__init__(
            file_name=file_name,
            channel_list=channel_list,
            convert_after_read=convert_after_read,
            filter_channel_names=filter_channel_names,
            no_data_loading=no_data_loading,
            compression=compression,
            convert_tables=convert_tables,
            metadata=metadata,
            finalization_writing_to_file=finalization_writing_to_file,
            force_file_integrity_check=force_file_integrity_check,
        )

    # ── Format-specific reader / writer delegations (MDF3 & MDF4) ────────────
    def read3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3.read3(self, *args, **kwargs)

    def write3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3.write3(self, *args, **kwargs)

    def _get_channel_data3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3._get_channel_data3(self, *args, **kwargs)

    def _convert_all_channel3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3._convert_all_channel3(self, *args, **kwargs)

    def _convert_channel3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3._convert_channel3(self, *args, **kwargs)

    def _convert3(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf3._convert3(self, *args, **kwargs)

    def read4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.read4(self, *args, **kwargs)

    def write4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.write4(self, *args, **kwargs)

    def _write4_non_column(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4._write4_non_column(self, *args, **kwargs)

    def _write4_column(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4._write4_column(self, *args, **kwargs)

    def _get_channel_data4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4._get_channel_data4(self, *args, **kwargs)

    def _convert_all_channel4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4._convert_all_channel4(self, *args, **kwargs)

    def _convert_channel4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4._convert_channel4(self, *args, **kwargs)

    @staticmethod
    def _convert_channel_data4(*args: Any, **kwargs: Any) -> Any:
        return Mdf4._convert_channel_data4(*args, **kwargs)

    def apply_invalid_bit(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.apply_invalid_bit(self, *args, **kwargs)

    def get_invalid_mask(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.get_invalid_mask(self, *args, **kwargs)

    def apply_all_invalid_bit(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.apply_all_invalid_bit(self, *args, **kwargs)

    def get_channel_name4(self, *args: Any, **kwargs: Any) -> Any:
        return Mdf4.get_channel_name4(self, *args, **kwargs)

    def read(
        self,
        file_name: str | None = None,
        multi_processed: bool = False,
        channel_list: Sequence[str] | set[str] | None = None,
        convert_after_read: bool = True,
        filter_channel_names: bool = False,
        no_data_loading: bool = False,
        compression: bool | str | int = False,
        metadata: int = 2,
        finalization_writing_to_file: bool = False,
        force_file_integrity_check: bool = False,
    ) -> None:
        """ reads mdf file version 3.x and 4.x

        Parameters
        ----------------
        file_name : str, optional
            file name

        multi_processed : bool
            flag to activate multiprocessing of channel data conversion.

        channel_list : list of str, optional
            list of channel names to be read.
            If you use channelList, reading might be much slower but it will save you memory.
            Can be used to read big files.

        convert_after_read : bool, optional
            flag to convert channel after read, True by default.
            If you use convertAfterRead by setting it to false, all data from channels will be kept raw,
            no conversion applied. If many float are stored in file, you can gain from 3 to 4 times memory footprint.
            To calculate value from channel, you can then use method getChannelData()

        filter_channel_names : bool, optional
            flag to filter long channel names from its module names separated by '.'

        no_data_loading : bool, optional
            Flag to read only file info but no data to have minimum memory use.

        compression : bool or str, optional
            To compress data in memory using blosc or bcolz, takes cpu time.
            if compression = int(1 to 9), uses bcolz for compression.
            if compression = 'blosc', uses blosc for compression.
            Choice given, efficiency depends of data.

        metadata: int, optional, default = 2
            Reading metadata has impact on performance, especially for mdf 4.x using xml.
            2: minimal metadata reading (mostly channel blocks).
            1: used for noDataLoading.
            0: all metadata reading, including Source Information, Attachment, etc..

        finalization_writing_to_file : bool, optional, False by default
            If file is detected not finalised (id_unfin_flags!=0), file is corrected
            writing in the blocks if set to True, otherwise correction is done only
            to file representation in memory. Valid from version MDF 4.11

        force_file_integrity_check : bool, optional, False by default
            Perform block sizes check for potentially corrupted file without finalization
            flags (id_unfin_flags==0). Combined with finalization_writing_to_file is
            very experimental and risky, correction should be tried in memory first.

        Notes
        --------
        If you keep convertAfterRead to true, you can set attribute mdf.multiProc to activate channel conversion
         in multiprocessing. Gain in reading time can be around 30% if file is big and using a lot of float channels

        Warning
        ------------
        MultiProc use should be avoided when reading several files in a batch, it is not thread safe.
        You should better multi process instances of mdf rather than using multiproc in mdf class
        (see implementation of mdfconverter)
        """
        if self.fileName is None or file_name is not None:
            self.fileName = file_name

        # Open file
        (self.fid, self.fileName, self.zipfile) = _open_mdf(self.fileName)

        # read Identifier block
        self.fid.seek(28)
        mdf_version_number = unpack('<H', self.fid.read(2))
        self.MDFVersionNumber = mdf_version_number[0]

        if self.MDFVersionNumber < 400:  # up to version 3.x not compatible with version 4.x
            if not no_data_loading:
                self.read3(self.fileName, None, multi_processed, channel_list,
                           convert_after_read, filter_channel_names, compression)
            else:  # populate minimum mdf structure
                self._noDataLoading = True
                self.info = Info3(None, fid=self.fid,
                                  filter_channel_names=filter_channel_names, minimal=1)
                (self.masterChannelList, mdf_dict) = _generate_dummy_mdf3(
                    self.info, channel_list)
                self.update(mdf_dict)
        else:  # MDF version 4.x
            if not no_data_loading:
                self.read4(self.fileName, None, multi_processed, channel_list,
                           convert_after_read, filter_channel_names, compression, metadata,
                           finalization_writing_to_file, force_file_integrity_check)
            else:  # populate minimum mdf structure
                self._noDataLoading = True
                self.info = Info4(None, fid=self.fid,
                                  filter_channel_names=filter_channel_names, minimal=1)
                (self.masterChannelList, mdf_dict) = _generate_dummy_mdf4(
                    self.info, channel_list)
                self.update(mdf_dict)

        if not self.fid.closed:  # close file
            self.fid.close()

    def write(
        self,
        file_name: str | None = None,
        compression: bool = False,
        column_oriented: bool = False,
    ) -> None:
        """Writes simple mdf file, same format as originally read, default is 4.x

        Parameters
        ----------------
        file_name : str, optional
            Name of file
            If file name is not input, written file name will be the one read with
            appended '_new' string before extension
        compression : bool
            Flag to store data compressed (from mdf version 4.1)
            If activated, will write in version 4.1 even if original file is in version 3.x
        column_oriented : bool
            Flag to store , column oriented channels

        Notes
        --------
        All channels will be converted, so size might be bigger than original file
        """
        if file_name is None:
            split_name = splitext(self.fileName)
            if split_name[-1] in ('.mfxz', '.MFXZ'):
                split_name[-1] = '.mfx'  # do not resave in compressed file
            file_name = ''.join([split_name[-2], '_New', split_name[-1]])
        # makes sure all channels are converted
        self.convert_all_channels()
        if self.MDFVersionNumber < 400 and not compression:
            self.write3(file_name=file_name)
        else:
            self.write4(file_name=file_name, compression=compression,
                        column_oriented=column_oriented)

    def get_channel_data(self, channel_name: str, raw_data: bool = False) -> Any:
        """Return channel numpy array

        Parameters
        ----------------
        channel_name : str
            channel name
        raw_data: bool
            flag to return non converted data

        Returns
        -----------
        numpy array
            converted, if not already done, data corresponding to channel name

        Notes
        ------
        This method is the safest to get channel data as numpy array from 'data' dict key might contain raw data
        """
        if channel_name in self:
            chan_dict = self[channel_name]
            vec = chan_dict.get(dataField)
            if vec is not None and not self._noDataLoading and not isinstance(vec, CompressedData):
                if raw_data or conversionField not in chan_dict:
                    return vec
        if self.MDFVersionNumber < 400:
            vector = self._get_channel_data3(channel_name, raw_data)
        else:
            vector = self._get_channel_data4(channel_name, raw_data)
        if self._noDataLoading:
            # remove data loaded in object to save memory
            self.set_channel_data(channel_name, None)
        return vector

    def convert_all_channels(self) -> Any:
        """Converts all channels from raw data to converted data according to CCBlock information.
        Converted data will take more memory.
        """
        if self.MDFVersionNumber < 400:
            return self._convert_all_channel3()
        else:
            return self._convert_all_channel4()

    def plot(self, channel_name_list_of_list: str | Sequence[str] | Sequence[Sequence[str]]) -> Any:
        """Plot channels with Matplotlib

        Parameters
        ----------------
        channel_name_list_of_list : str or list of str or list of list of str
            channel name or list of channel names or list of list of channel names
            list of list will create multiplots

        Notes
        ---------
        Channel description and unit will be tentatively displayed with axis labels
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            warn('matplotlib not found')
            return
        try:
            from mplcursors import cursor as datacursor
            cursor_possible = True
        except ImportError:
            try:  # legacy, unmaintained package
                from mpldatacursor import datacursor
                cursor_possible = True
            except ImportError:
                warn('no cursor available, to use cursor, please install mplcursors')
                cursor_possible = False
        if isinstance(channel_name_list_of_list, str):
            channel_name_list_of_list = [
                channel_name_list_of_list]  # converts in list
        for channel_name_list in channel_name_list_of_list:
            fig = plt.subplot(len(channel_name_list_of_list), 1,
                              channel_name_list_of_list.index(channel_name_list) + 1)
            if isinstance(channel_name_list, str):
                channel_name_list = [channel_name_list]  # converts in list
            for channelName in channel_name_list:
                if channelName in self:
                    data = self.get_channel_data(channelName)
                    # if channel not a string
                    if data.dtype.kind not in {'S', 'U'}:
                        # self.fig = plt.figure()  # could be needed
                        # plot using matplotlib the channel versus master channel
                        # Resampled signals or only one master
                        if len(self.masterChannelList) == 1:
                            master_name = next(iter(self.masterChannelList))
                            if not master_name:  # resampled channels, only one time channel probably called 'master'
                                master_name = 'master'
                            master_data = self.get_channel_data(master_name)
                            if master_data is None:  # no master channel
                                master_data = arange(0, len(data), 1)
                            if master_name in self.masterChannelList:  # time channel properly defined
                                plt.plot(master_data, data, label=channelName)
                                plt.xlabel(f'{master_name} [{self.get_channel_unit(master_name)}]')
                            else:  # no time channel found
                                plt.plot(data)
                        else:  # not resampled
                            master_name = self.get_channel_master(channelName)
                            master_data = self.get_channel_data(master_name)
                            if master_data is None:  # no master channel
                                master_data = arange(0, len(data), 1)
                            if master_name in self.masterChannelList:  # master channel is proper channel name
                                plt.plot(master_data, data, label=channelName)
                                plt.xlabel(f'{master_name} [{self.get_channel_unit(master_name)}]')
                            else:
                                plt.plot(data)

                        plt.title(self.get_channel_desc(channelName))
                        if self.get_channel_unit(channelName) == {}:
                            plt.ylabel(channelName)
                        else:
                            plt.ylabel(f'{channelName} [{self.get_channel_unit(channelName)}]')
                else:
                    warn(f'Channel {channelName} not existing')
            plt.grid(True)
            plt.legend(loc="upper left", frameon=False)
            plt.title('')
        if cursor_possible:
            datacursor()
        plt.show()
        return fig

    def plot_all(self) -> None:
        # plot all channels in the object, be careful for test purpose only,
        # can display many many many plots overloading your computer
        for Name in self:
            try:
                self.plot(Name)
            except Exception:
                warn(Name)

    def resample(
        self,
        sampling: float | None = None,
        channel: str | None = None,
        master_channel: str | None = None,
        interpolation_kind: str | None = None,
    ) -> None:
        """ Resamples as much as possible all data groups into one data group having defined
        sampling interval or sharing same defined master channel

        Parameters
        ----------------
        sampling : float, optional
            resampling interval, None by default. If None, will rely on channel or master_channel
             parameters to define reference data group. If both are undefined, picking the first master
        interpolation_kind : str, optional
            interpolation type for floating data from scipy interp1d list (‘linear’, ‘nearest’, ‘zero’,
             ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’)
        ** or | and **
        channel : str, optional
            channel name to be resampled
        ** or | and **
        master_channel : str, optional
            master channel name to be used as reference

        Notes
        --------
        1. resampling will be applied only to master channels that have same type as the one
        given by channel or master_channel parameters (applicable only to mdf4)

        2. resampling will convert all your channels so be careful for big files
        and memory consumption
        """
        if self:  # mdf contains data
            # must make sure all channels are converted
            self.convert_all_channels()

            if channel is not None and master_channel is None:
                master_channel_name = self.get_channel_master(channel)
                master_channel_type = self.get_channel_master_type(
                    master_channel_name)
            elif channel is None and master_channel is None:
                # pick the first master
                master_channel_name = next(iter(self.masterChannelList))
                master_channel_type = self.get_channel_master_type(
                    master_channel_name)
            else:
                master_channel_name = master_channel
                master_channel_type = self.get_channel_master_type(
                    master_channel)
            master_data = self.get_channel_data(master_channel_name)
            if sampling is not None:
                master_data = arange(master_data[0], master_data[-1], sampling)

            if master_channel_name is None or master_channel_name not in self.masterChannelList or\
                    master_channel_name not in self.masterChannelList[master_channel_name]:
                # No master channel in selected group, looking at other groups
                min_master = []
                max_master = []
                length = []
                master_channel_name = 'master'
                for master in self.masterChannelList:
                    if master is not None and master != '' and \
                            master in self and self.masterChannelList[master] and \
                            self.get_channel_master_type(master) == master_channel_type:
                        master_data = self.get_channel_data(master)
                        # consider groups having minimum size
                        if master in self and len(master_data) > 5:
                            min_master.append(master_data[0])
                            max_master.append(master_data[-1])
                            length.append(len(master_data))
                            master_channel_name = master
                if min_master:  # at least 1 datagroup has a master channel to be resampled
                    if sampling is None:
                        master_data = linspace(
                            min(min_master), max(max_master), num=max(length))
                    else:
                        master_data = arange(
                            min(min_master), max(max_master), sampling)
                    self.add_channel(master_channel_name, master_data, master_channel_name,
                                     master_type=self.get_channel_master_type(
                                         master),
                                     unit=self.get_channel_unit(master),
                                     description=self.get_channel_desc(master), conversion=None)
                else:
                    warn(
                        'no master channel existing, considering first channel as master !')
                    current_master = next(iter(self.masterChannelList))
                    # pick the first channel from the first data group to be master, maybe lucky
                    master_channel_name = self.masterChannelList[current_master][0]
                    # changing master to first channel
                    self.masterChannelList[master_channel_name] = self.masterChannelList.pop(
                        current_master)
                    for channel in self:
                        self.set_channel_master(channel, master_channel_name)
                    master_channel_type = self.get_channel_master_type(
                        master_channel_name)
                    master_data = self.get_channel_data(master_channel_name)
                    if sampling is not None:
                        master_data = arange(
                            master_data[0], master_data[-1], sampling)

            # Interpolate channels
            for master in list(self.masterChannelList):
                if self.get_channel_master_type(master) == master_channel_type:
                    self.resample_group(sampling, master, new_master_data=master_data,
                                        interpolation_kind=interpolation_kind)
                    # remove old master channel
                    if not master_channel_name == master:
                        # removing previous master from list
                        self.masterChannelList[master].remove(master)
                        for channel in self.masterChannelList[master]:
                            # assigning new master to resampled channels
                            self.set_channel_master(
                                channel, master_channel_name)
                        # merging channel lists
                        self.masterChannelList[master_channel_name] += self.masterChannelList[master]
                        # removing old master
                        self.masterChannelList.pop(master)
                        self.pop(master)
        else:
            warn('no data to be resampled')

    def resample_group(
        self,
        sampling: float,
        channel: str,
        new_master_data: Any = None,
        interpolation_kind: str | None = None,
    ) -> None:
        """ Resamples one channel along with its dataGroup

        Parameters
        ----------------
        sampling : float
            resampling interval

        channel : str
            channel name to be resampled (could the master channel)

        new_master_data : array, optional
            master channel data to be applied to the group identified by channel

        interpolation_kind : str, optional
            interpolation type for floating data from scipy interp1d list (‘linear’, ‘nearest’, ‘zero’,
             ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’)

        Notes
        --------
        Resampling will convert all channels so be careful for big files
        and memory consumption
        """
        def interp_close_point(x, new_x, kind):
            # interpolates with closest point at 'kind' side
            idx = searchsorted(x, new_x, side=kind)
            idx -= 1
            if len(idx) == 0:
                return idx
            idx = clip(idx, 0, idx[-1])
            return idx

        def interpolate(new_x, x, y, interpolation_kind):
            # select right interpolation method depending of
            # interpolation kind and data type
            if interpolation_kind is None:
                if y.dtype.kind == 'f':
                    return interp(new_x, x, y)
                else:
                    return y[interp_close_point(x, new_x, 'right')]
            else:
                if interpolation_kind == 'previous':
                    return y[interp_close_point(x, new_x, 'left')]
                elif interpolation_kind == 'next':
                    return y[interp_close_point(x, new_x, 'right')]
                elif y.dtype.kind not in ('U', 'S'):
                    from scipy.interpolate import interp1d
                    f = interp1d(x, y, kind=interpolation_kind,
                                 fill_value="extrapolate")
                    return f(new_x)
                else:
                    return y[interp_close_point(x, new_x, 'right')]

        master_channel = self.get_channel_master(channel)
        old_master_data = self.get_channel_data(master_channel)
        if new_master_data is None:
            new_master_data = arange(
                old_master_data[0], old_master_data[-1], sampling)

        cached_idx_right = None
        cached_idx_left = None

        def get_idx(side):
            nonlocal cached_idx_right, cached_idx_left
            if side == 'right':
                if cached_idx_right is None:
                    cached_idx_right = interp_close_point(old_master_data, new_master_data, 'right')
                return cached_idx_right
            else:
                if cached_idx_left is None:
                    cached_idx_left = interp_close_point(old_master_data, new_master_data, 'left')
                return cached_idx_left

        def interpolate_cached(new_x, x, y, interpolation_kind):
            if interpolation_kind is None:
                if y.dtype.kind == 'f':
                    return interp(new_x, x, y)
                else:
                    return y[get_idx('right')]
            elif interpolation_kind == 'previous':
                return y[get_idx('left')]
            elif interpolation_kind == 'next':
                return y[get_idx('right')]
            elif y.dtype.kind not in ('U', 'S'):
                from scipy.interpolate import interp1d
                f = interp1d(x, y, kind=interpolation_kind, fill_value="extrapolate")
                return f(new_x)
            else:
                return y[get_idx('right')]
        for Name in self.masterChannelList[master_channel]:
            # forces list() because masterChannelList is dynamic, channels can be removed
            if Name == master_channel:  # master channel
                self.set_channel_data(Name, new_master_data)
            else:
                channel_data = self.get_channel_data(Name)
                if channel_data.dtype.kind not in ('S', 'U', 'V') and channel_data.ndim == 1:
                    # if channel not array of string
                    try:
                        self.set_channel_data(Name, interpolate_cached(new_master_data, old_master_data,
                                                                       channel_data, interpolation_kind))
                    except Exception:
                        if not all(diff(old_master_data) > 0):
                            warn(f'{Name} has non regularly increasing master channel {master_channel}.\n'
                                 ' Faulty samples will be dropped in related data group')
                            self._clean_uneven_master_data(master_channel)
                            self.set_channel_data(Name, interpolate(new_master_data, old_master_data,
                                                                    channel_data, interpolation_kind))
                        elif old_master_data is not None and len(old_master_data) != len(channel_data):
                            warn(f'{Name} and master channel {master_channel} do not have same length')
                    self.remove_channel_conversion(Name)
                else:  # can not interpolate strings, remove channel containing string
                    self.remove_channel(Name)

    def _clean_uneven_master_data(self, master_channel_name: str) -> None:
        """ clean data group having non evenly increasing master

            Parameters
            ----------------
            master_channel_name : str
                master channel name referring to data group to be cleaned

        """
        # create mask
        mask = diff(self.get_channel_data(master_channel_name)) > 0
        # remove samples from data group defined by mask
        for channel in self.masterChannelList[master_channel_name]:
            data = self.get_channel_data(channel).view(MaskedArray)
            data.mask = mask
            self.set_channel_data(channel, data.compressed())

    def cut(
        self,
        master_channel: str,
        begin: float | None = None,
        end: float | None = None,
    ) -> None:
        """ Cut data

        Parameters
        ----------------
        master_channel : str
            channel to cut data (can be master channel)

        begin : float
            beginning value in master channel from which to start cutting in all channels

        end : float
            ending value in master channel from which to start cutting in all channels

        Notes
        ------
        Only the data groups with same master type as master_channel will be cut (only for mdf4)

        """
        if begin is None and end is None:
            raise Exception(
                'Please input at least one beginning or ending value to cut data')

        master_channel = self.get_channel_master(master_channel)
        master_channel_type = self.get_channel_master_type(master_channel)
        for master in self.masterChannelList:  # for each channel group
            # find corresponding indexes to cut
            master_data = self.get_channel_data(master)
            if master_data is not None and len(master_data) > 0 and \
                    self.get_channel_master_type(master) == master_channel_type:
                # not empty data and same master type
                if begin is not None:
                    start_index = searchsorted(master_data, begin, side='left')
                else:
                    start_index = 0
                if end is not None:
                    end_index = searchsorted(master_data, end, side='right')
                else:
                    end_index = len(master_data)
                if start_index == end_index:
                    # empty array
                    for channel in self.masterChannelList[master]:
                        self.set_channel_data(channel, array([]))
                else:
                    for channel in self.masterChannelList[master]:
                        data = self.get_channel_data(channel)
                        self.set_channel_data(
                            channel, data[start_index: end_index])

    # ── Exporters (delegated to mdfreader.export) ───────────────────────────
    def export_to_csv(self, file_name: str | None = None, sampling: float | None = None) -> str | None:
        """Exports mdf data into CSV file."""
        return export_csv(self, file_name=file_name, sampling=sampling)

    def export_to_NetCDF(self, file_name: str | None = None, sampling: float | None = None) -> str | None:
        """Exports mdf data into NetCDF file."""
        return export_netcdf(self, file_name=file_name, sampling=sampling)

    def export_to_hdf5(
        self,
        file_name: str | None = None,
        sampling: float | None = None,
        compression: str | None = None,
        compression_opts: int | None = None,
    ) -> str | None:
        """Exports mdf class data structure into HDF5 file."""
        return export_hdf5(
            self,
            file_name=file_name,
            sampling=sampling,
            compression=compression,
            compression_opts=compression_opts,
        )

    def export_to_matlab(self, file_name: str | None = None) -> str | None:
        """Export mdf data into Matlab file preferrably in format 7.3."""
        return export_matlab(self, file_name=file_name)

    def export_to_excel(self, file_name: str | None = None) -> str | None:
        """Exports mdf data into excel 95 to 2003 file."""
        return export_excel(self, file_name=file_name)

    def export_to_xlsx(self, file_name: str | None = None) -> str | None:
        """Exports mdf data into excel 2007 and 2010 file."""
        return export_xlsx(self, file_name=file_name)

    def keep_channels(self, channel_set: Sequence[str] | set[str]) -> None:
        """ keeps only a list or set of channels and removes the other channels

        Parameters
        ----------------
        channel_set : list or set of str
            list or set of channel names
        """
        if isinstance(channel_set, (list, tuple)):
            channel_set = set(channel_set)
        remove_channels = []
        master_channel_set = set(self.masterChannelList.keys())
        for channel in self:
            if channel not in channel_set and masterField not in channel \
                    and channel not in master_channel_set:
                # avoid to remove master channels otherwise problems with resample
                remove_channels.append(channel)
        if not len(remove_channels) == 0:
            [self.masterChannelList[self.get_channel_master(channel)].remove(
                channel) for channel in remove_channels]
            [self.pop(channel) for channel in remove_channels]

    def concat_mdf(self, mdf_class: Mdf) -> None:
        """Concatenate data of input mdf class after the current one.

        Parameters
        ----------------
        mdf_class : Mdf
            mdf class instance to be concatenated after self

        Notes
        --------
        It creates union of both channel lists and fills with Nan for unknown sections in channels
        If one channel is not present in both classes, masked array is created
        If invalid bytes are present, masked array are created
        """
        first_class_masters = set(self.masterChannelList.keys())
        second_class_masters = set(mdf_class.masterChannelList.keys())
        union_masters = first_class_masters | second_class_masters
        # find max time of all masters from both classes
        first_masters = {}
        for master_channel_name in first_class_masters:
            type = self.get_channel_master_type(master_channel_name)
            data = self.get_channel_data(master_channel_name)
            if data is None or len(data) == 0:
                continue
            if type not in first_masters:
                first_masters[type] = {}
                first_masters[type]['max'] = data[-1]
            # first_masters[type]['sampling'] = mean(diff(data))  # sampling
            first_masters[type]['max'] = max(
                [data[-1], first_masters[type]['max']])
        second_masters = {}
        for master_channel_name in second_class_masters:
            type = mdf_class.get_channel_master_type(master_channel_name)
            data = mdf_class.get_channel_data(master_channel_name)
            if data is None or len(data) == 0:
                continue
            if type not in second_masters:
                second_masters[type] = {}
                second_masters[type]['max'] = data[-1]
            # second_masters[type]['sampling'] = mean(diff(data))  # sampling
            second_masters[type]['max'] = max(
                [data[-1], second_masters[type]['max']])

        for master_channel_name in union_masters:
            if master_channel_name in first_class_masters and master_channel_name in second_class_masters:
                # same master name in both classes
                _first_data = self.get_channel_data(master_channel_name)
                if _first_data is None:
                    continue
                first_class_length = len(_first_data)
                first_class_channels = set(
                    self.masterChannelList[master_channel_name])
                second_class_channels = set(
                    mdf_class.masterChannelList[master_channel_name])
                if self.MDFVersionNumber >= 400:
                    first_class_group_number = self[master_channel_name].get(idField, [[None]])[0][0]
                    invalid_channel = f'invalid_bytes{first_class_group_number}'
                    if invalid_channel in first_class_channels:
                        # invalid bits present, converting to masked array
                        # but invalid bit channels should be removed, updating channel sets
                        self.apply_all_invalid_bit()
                        first_class_channels.remove(invalid_channel)
                if mdf_class.MDFVersionNumber >= 400:
                    second_class_group_number = mdf_class[master_channel_name].get(idField, [[None]])[0][0]
                    invalid_channel = f'invalid_bytes{second_class_group_number}'
                    if invalid_channel in second_class_channels:
                        mdf_class.apply_all_invalid_bit()
                        second_class_channels.remove(invalid_channel)
                unioned_set = first_class_channels | second_class_channels
                second_class_length = len(
                    mdf_class.get_channel_data(master_channel_name))
                total_length = first_class_length + second_class_length
                # remove master channel
                unioned_set.remove(master_channel_name)
                # merge master channels
                data = self.get_channel_data(master_channel_name)
                mdf_data = mdf_class.get_channel_data(master_channel_name)
                temp = empty(total_length, dtype=data.dtype)
                temp[:first_class_length] = data
                if self.get_channel_master_type(master_channel_name) == 1:
                    # master of type time
                    offset = mean(diff(mdf_data)) if len(mdf_data) > 1 else 0  # sampling
                    offset = (data[-1] + offset) if len(data) > 0 else offset  # offset
                    temp[first_class_length:] = mdf_data + offset
                else:
                    temp[first_class_length:] = mdf_data
                self.set_channel_data(master_channel_name, temp)
                for channel in unioned_set:
                    if channel in first_class_channels and channel in second_class_channels:
                        # channel exists in both classes
                        data = self.get_channel_data(channel)
                        if data is None:
                            continue
                        mdf_data2 = mdf_class.get_channel_data(channel)
                        if mdf_data2 is None:
                            continue
                        if data.ndim == 1:
                            ch_total = len(data) + len(mdf_data2)
                            if isinstance(data, MaskedArray):
                                temp = ma_empty(ch_total, dtype=data.dtype)
                            else:
                                temp = empty(ch_total, dtype=data.dtype)
                            temp[:len(data)] = data
                            temp[len(data):] = mdf_data2
                        else:
                            from numpy import concatenate as np_concatenate
                            temp = np_concatenate([data, mdf_data2], axis=0)
                        self.set_channel_data(channel, temp)
                    elif channel in second_class_channels:
                        # initialise all fields, units, descriptions, etc.
                        self[channel] = mdf_class[channel]
                        # new channel for self from mdfClass
                        data = self.get_channel_data(channel)
                        if data is None or data.ndim != 1:
                            continue  # skip non-1D channels for padding
                        if isinstance(data, MaskedArray):
                            temp = ma_empty(total_length, dtype=data.dtype)
                        else:
                            temp = empty(total_length, dtype=data.dtype)
                        temp[:first_class_length] = nan  # fill with NANs
                        temp[first_class_length:] = data
                        self.set_channel_data(channel, temp)
                    else:  # channel missing in mdfClass
                        data = self.get_channel_data(channel)
                        if data is None or data.ndim != 1:
                            continue  # skip non-1D channels for padding
                        if isinstance(data, MaskedArray):
                            temp = ma_empty(total_length, dtype=data.dtype)
                        else:
                            temp = empty(total_length, dtype=data.dtype)
                        temp[:first_class_length] = data
                        temp[first_class_length:] = nan  # fill with NANs
                        self.set_channel_data(channel, temp)
            elif master_channel_name in first_class_masters:
                first_class_channels = set(
                    self.masterChannelList[master_channel_name])
                if self.MDFVersionNumber >= 400:
                    first_class_group_number = self[master_channel_name].get(idField, [[None]])[0][0]
                    invalid_channel = f'invalid_bytes{first_class_group_number}'
                    if invalid_channel in first_class_channels:
                        # invalid bits present, converting to masked array
                        # but invalid bit channels has been removed, updating channel sets
                        self.apply_all_invalid_bit()
                        first_class_channels.remove(invalid_channel)
                master_type = self.get_channel_master_type(master_channel_name)
                if master_type in second_masters:
                    # same kind of master existing in both classes
                    master_data = self.get_channel_data(master_channel_name)
                    first_class_length = len(master_data)
                    first_master_end = master_data[-1]
                    total_length = first_class_length + 1
                    # extend by 1 sample all channels and convert to masked array
                    for channel in first_class_channels:
                        data = self.get_channel_data(channel)
                        temp = ma_empty(total_length, dtype=data.dtype)
                        temp[:first_class_length] = data.view(MaskedArray)
                        temp[-1] = masked  # mask last sample
                        self.set_channel_data(channel, temp)
                    master_data = self.get_channel_data(master_channel_name)
                    # last master channel value adjusted and becomes valid
                    master_data[-1] = first_master_end + \
                        second_masters[master_type]['max']
                    self.set_channel_data(master_channel_name, master_data)
            else:
                second_class_channels = set(
                    mdf_class.masterChannelList[master_channel_name])
                if mdf_class.MDFVersionNumber >= 400:
                    second_class_group_number = mdf_class[master_channel_name].get(idField, [[None]])[0][0]
                    invalid_channel = f'invalid_bytes{second_class_group_number}'
                    if invalid_channel in second_class_channels:
                        mdf_class.apply_all_invalid_bit()
                        second_class_channels.remove(invalid_channel)
                master_type = mdf_class.get_channel_master_type(
                    master_channel_name)
                if master_type in first_masters:
                    # same kind of master existing in both classes
                    master_data = mdf_class.get_channel_data(
                        master_channel_name)
                    second_class_length = len(master_data)
                    second_master_end = master_data[-1]
                    total_length = second_class_length + 1
                    for channel in second_class_channels:
                        self[channel] = mdf_class[channel]
                        data = self.get_channel_data(channel)
                        temp = ma_empty(total_length, dtype=data.dtype)
                        temp[1:] = data
                        temp[0] = masked  # mask first sample
                        self.set_channel_data(channel, temp)
                    master_data = self.get_channel_data(master_channel_name)
                    # first master channel value adjusted and becomes valid
                    master_data += first_masters[master_type]['max']
                    master_data[0] = 0
                    self.set_channel_data(master_channel_name, master_data)

    def merge_mdf(self, mdf_class: Mdf) -> None:
        """merge data of input mdf class with the current one.

        Parameters
        ----------------
        mdf_class : Mdf
            mdf class instance to be merged with self

        Notes
        --------
        If there are common channel names between the 2 mdf, channels are renamed to make them unique
        """
        # apply eventual invalid bytes, most probably
        if self.MDFVersionNumber >= 400:
            self.apply_all_invalid_bit()
        if mdf_class.MDFVersionNumber >= 400:
            mdf_class.apply_all_invalid_bit()
        # check for same channel names
        first_channels = set(self.keys())
        second_channels = set(mdf_class.keys())
        common_channels = first_channels & second_channels
        if common_channels:
            # Batch rename conflicting channels in mdf_class efficiently
            rename_map = {}
            for channel in common_channels:
                channel_id = mdf_class[channel].get(idField, [[channel]])[0][0]
                new_name = f'{channel}_{channel_id}'
                rename_map[channel] = new_name
                mdf_class[new_name] = mdf_class.pop(channel)

            new_master_list = {}
            for master, ch_list in mdf_class.masterChannelList.items():
                new_master = rename_map.get(master, master)
                updated_ch_list = [rename_map.get(ch, ch) for ch in ch_list]
                new_master_list[new_master] = updated_ch_list
                if new_master != master:
                    for ch in updated_ch_list:
                        mdf_class.set_channel_master(ch, new_master)
            mdf_class.masterChannelList = new_master_list
        # copy the data
        for channel in mdf_class:
            self[channel] = mdf_class[channel]
        # merge the 2 masterChannelList
        masterChannelList = {**self.masterChannelList,
                             **mdf_class.masterChannelList}
        self.masterChannelList = masterChannelList

    def convert_to_pandas(self, sampling: float | None = None) -> None:
        """Converts mdf data structure into pandas dataframe(s).

        Parameters
        ----------
        sampling : float, optional
            Resampling interval before conversion.
        """
        return to_pandas(self, sampling=sampling)

    def return_pandas_dataframe(self, master_channel_name: str) -> Any:
        """Returns a dataframe of a raster described by its master channel name.

        Parameters
        ----------
        master_channel_name : str
            Master channel name, key to a raster to be returned as pandas dataframe.

        Returns
        -------
        pandas.DataFrame or None
        """
        return return_pandas_dataframe(self, master_channel_name=master_channel_name)

    def get_master_channel_datetimes(self, master_channel_name: str) -> Any:
        """Return a time master channel as absolute datetime64[ns] array.

        Parameters
        ----------
        master_channel_name : str
            Name of a master channel with masterType == 1 (Time).

        Returns
        -------
        numpy.ndarray of datetime64[ns], or None.
        """
        return get_master_channel_datetimes(self, master_channel_name=master_channel_name)

    def export_to_parquet(self, file_name: str | None = None) -> None:
        """Exports mdf data into parquet file.

        Parameters
        ----------
        file_name : str, optional
            Target Parquet file path.
        """
        return export_parquet(self, file_name=file_name)


if __name__ == "__main__":
    if osname == 'nt':
        from multiprocessing import freeze_support
        freeze_support()
    parser = ArgumentParser(prog='mdfreader', description='reads mdf file')
    parser.add_argument('--export', dest='export', default=None,
                        choices=['CSV', 'HDF5', 'Matlab', 'Xlsx',
                                 'Excel', 'NetCDF', 'MDF', 'parquet'],
                        help='Export after parsing to defined file type')
    parser.add_argument('--list_channels', dest='list_channels', action='store_true',
                        help='list of channels in file')
    parser.add_argument('fileName', help='mdf file name')
    parser.add_argument('--channelList', dest='channelList', nargs='+', type=str,
                        default=None, help='list of channels to only read')
    parser.add_argument('--plot', dest='plot_channel_list', nargs='+', type=str,
                        default=None, help='plots list of channels')
    parser.add_argument('--noConversion', dest='convertAfterRead', action='store_false',
                        help='Do not convert raw channel data \
            to physical values just after reading. Useful if you have memory concerns')
    parser.add_argument('--filterChannelNames', action='store_true',
                        help='activates channel name filtering; \
            removes modules names separated by a point character')
    parser.add_argument('--noDataLoading', action='store_true',
                        help='Do not load data in memory, creates dummy \
            mdf structure and load data from file when needed')
    parser.add_argument('--compression', action='store_true',
                        help='compress data in memory')

    args = parser.parse_args()

    temp = Mdf(file_name=args.fileName, channel_list=args.channelList,
               convert_after_read=args.convertAfterRead,
               filter_channel_names=args.filterChannelNames,
               no_data_loading=args.noDataLoading,
               compression=args.compression)
    if args.export is not None:
        if args.export == 'CSV':
            temp.export_to_csv()
        elif args.export == 'HDF5':
            temp.export_to_hdf5()
        elif args.export == 'Matlab':
            temp.export_to_matlab()
        elif args.export == 'Xlsx':
            temp.export_to_xlsx()
        elif args.export == 'Excel':
            temp.export_to_excel()
        elif args.export == 'NetCDF':
            temp.export_to_NetCDF()
        elif args.export == 'MDF':
            temp.write()
        elif args.export == 'parquet':
            temp.export_to_parquet()
    if args.plot_channel_list is not None:
        temp.plot(args.plot_channel_list)
    if args.list_channels:
        print(temp.masterChannelList)
