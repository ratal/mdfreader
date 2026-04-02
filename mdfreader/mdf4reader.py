# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 4.x.

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Thu Dec 10 12:57:28 2013


Dependencies
-------------------
- Python >3.4 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>
- bitarray to parse bits in not aligned bytes
- Sympy to convert channels with formula if needed
- zlib to uncompress data block if needed

mdf4reader
--------------------------

"""
from struct import Struct
from struct import pack, unpack as structunpack
from math import pow
from io import open
from os.path import splitext
from multiprocessing import Queue, Process
from sys import byteorder
import re
from collections import defaultdict, OrderedDict
import numpy as np
if np.lib.NumpyVersion(np.__version__) >= '2.0.0b1':
    from numpy.rec import fromstring, fromarrays
else:
    from numpy.core.records import fromstring, fromarrays
from numpy import array, recarray, asarray, empty, where, frombuffer, reshape
from numpy import arange, right_shift, bitwise_and, all, diff, interp, zeros, concatenate, maximum
from numpy import issubdtype, number as numpy_number
from numpy import max as npmax, min as npmin
from numpy.lib.recfunctions import rename_fields
from numpy.ma import MaskedArray
from zlib import decompress
from warnings import simplefilter, warn
from .mdfinfo4 import Info4, IDBlock, HDBlock, DGBlock, \
    CGBlock, CNBlock, FHBlock, CommentBlock, _load_header, DLBlock, \
    DZBlock, HLBlock, CCBlock, DTBlock, CABlock, DVBlock, LDBlock
from .mdf import MdfSkeleton, _open_mdf, invalidChannel, dataField, \
    conversionField, idField, invalidPosField, CompressedData
from .channel import Channel4
try:
    from dataRead import sorted_data_read, unsorted_data_read4, sd_data_read, vd_data_read
    dataRead_available = True
except ImportError:
    warn('dataRead cannot be imported, compile it with Cython', ImportWarning)
    dataRead_available = False

# reads by chunk of 100Mb, can be tuned for best performance
chunk_size_reading = 100000000
_VLSDStruct = Struct('I')

# Module-level frozensets for O(1) block ID membership checks (used in hot paths)
_DATA_BLOCK_IDS = frozenset((b'##DT', b'##DV', b'##RD', b'##DI', b'##RV', b'##RI',
                              '##DT', '##DV', '##RD', '##DI', '##RV', '##RI'))
_DV_BLOCK_IDS = frozenset((b'##DV', '##DV'))
_SD_BLOCK_IDS = frozenset((b'##SD', '##SD'))
_DZ_BLOCK_IDS = frozenset((b'##DZ', '##DZ'))
_DI_BLOCK_IDS = frozenset((b'##DI', '##DI'))
_DT_RD_DV_BLOCK_IDS = frozenset((b'##DT', b'##DV', b'##RD', '##DT', '##RD', '##DV'))


def _data_block(record, info, parent_block, channel_set=None, n_records=None, sorted_flag=True, vlsd=None):
    """ converts raw data into arrays

    Parameters
    ----------------
    record : class
        record class instance describing a channel group record
    parent_block : class
        MDFBlock class containing at least parent block header
    channel_set : set of str, optional
        defines set of channels to only read, can be slow but saves memory,
        for big files
    n_records: int, optional
        number of records to read
    sorted_flag : bool, optional
        flag to know if data block is sorted (only one Channel Group in block)
        or unsorted (several Channel Groups identified by a recordID).
        As unsorted block can contain CG records in random order, block
        is processed iteratively, not in raw like sorted -> much slower reading
    vlsd : bool
        indicate a sd block, compressed (DZ) or not (SD)

    Returns
    ---------
    a recarray containing the channels data

    Notes
    --------
    This function will read DTBlock, RDBlock, DZBlock (compressed),
    RDBlock (VLSD), sorted or unsorted
    """
    if n_records is None and hasattr(record, 'numberOfRecords'):
        n_records = record.numberOfRecords
    # normal data block
    if parent_block['id'] in _DT_RD_DV_BLOCK_IDS:
        if sorted_flag:
            if channel_set is None and not record.hiddenBytes and\
                    record.byte_aligned:  # No channel list and length of records corresponds to C datatypes
                # for debugging purpose
                # print(n_records, record.numpyDataRecordFormat, record.dataRecordName)
                if info['DG'][record.dataGroup]['unique_channel_in_DG'] and parent_block['id'] in _DV_BLOCK_IDS:
                    return frombuffer(parent_block['data'], dtype={'names': record.dataRecordName,
                                                                   'formats': record.numpyDataRecordFormat})
                else:
                    return frombuffer(parent_block['data'], dtype={'names': record.dataRecordName,
                                                                   'formats': record.numpyDataRecordFormat},
                                      count=n_records).view(recarray).copy()
            else:  # record is not byte aligned or channelSet not None
                return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        else:  # unsorted reading
            return _read_unsorted(record, info, parent_block, record[next(iter(record))]['record'].recordIDsize)

    elif parent_block['id'] in _SD_BLOCK_IDS:
        return _read_sd_block(record[record.VLSD[0]].signal_data_type(info), parent_block['data'],
                              parent_block['length'] - 24, n_records, vlsd)

    elif parent_block['id'] in _DZ_BLOCK_IDS:  # zipped data block
        # uncompress data
        parent_block['data'] = DZBlock.decompress_data_block(parent_block['data'], parent_block['dz_zip_type'],
                                                             parent_block['dz_zip_parameter'],
                                                             parent_block['dz_org_data_length'])
        if vlsd is not None:  # VLSD channel
            return _read_sd_block(record[record.VLSD[0]].signal_data_type(info), parent_block['data'],
                                  parent_block['dz_org_data_length'], n_records, vlsd)
        if channel_set is None and sorted_flag:  # reads all blocks if sorted block and no channelSet defined
            if record.byte_aligned and not record.hiddenBytes:
                return frombuffer(parent_block['data'], dtype={'names': record.dataRecordName,
                                                              'formats': record.numpyDataRecordFormat},
                                   count=n_records).view(recarray).copy()
            else:
                return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        elif channel_set is not None and sorted_flag:  # sorted data but channel list requested
            return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        else:  # unsorted reading
            return _read_unsorted(record, info, parent_block, record[next(iter(record))]['record'].recordIDsize)
    elif parent_block['id'] in _DI_BLOCK_IDS:  # Invalid data block
        return frombuffer(parent_block['data'], dtype={'names': [record.invalid_channel.name],
                                                       'formats': [record.invalid_channel.data_format(info)]})


def _apply_unsorted_bit_masking(buf, record, info):
    """Apply bit-shift and masking for non-byte-aligned channels in unsorted data.

    unsorted_data_read4 (and the Python fallback) copy raw bytes without any
    bit-level manipulation.  Channels whose bit_count is not a multiple of 8
    (e.g. 1-bit flags, 15-bit signals) or that have a non-zero bit_offset must
    be right-shifted and masked after the raw bytes have been viewed as the
    target dtype.
    """
    for record_id in record:
        if record[record_id]['record'].Flags & 0b100001:  # VLSD / VLSC - no masking needed
            continue
        for Channel in record[record_id]['record'].values():
            name = Channel.name
            if name not in buf:
                continue
            sig_dt = Channel.signal_data_type(info)
            if sig_dt not in (0, 1, 2, 3):  # only integer types need masking
                continue
            bit_count = Channel.bit_count(info)
            bit_offset = Channel.bit_offset(info)
            n_bytes_aligned = Channel.nBytes_aligned
            if bit_offset == 0 and bit_count == n_bytes_aligned * 8:
                continue  # already byte-aligned and fills all bytes — no masking needed
            temp = buf[name]
            if temp.dtype.kind not in ('u', 'i'):
                continue
            if bit_offset > 0:
                temp = right_shift(temp, bit_offset)
            mask = int((1 << bit_count) - 1)
            temp = bitwise_and(temp, mask)
            if sig_dt in (2, 3):  # signed: sign-extend from bit_count-1
                sign_bit_mask = 1 << (bit_count - 1)
                sign_extend = ((1 << (temp.itemsize * 8 - bit_count)) - 1) << bit_count
                sign_bits = bitwise_and(temp, sign_bit_mask)
                temp_u = temp.view('u{}'.format(temp.itemsize))
                for idx, sign in enumerate(sign_bits):
                    if sign:
                        temp_u[idx] |= sign_extend
                temp = temp_u.view(temp.dtype)
            buf[name] = temp


def _read_unsorted(record, info, parent_block, record_id_size):
    """ reads only the channels using offset functions, channel by channel within unsorted data

    Parameters
    ------------
    record : class
        record class
    info: class
        info class
    parent_block: class
        MDFBlock class containing at least parent block header
    record_id_size : int
        Size of recordId

    Returns
    --------
    buf : array
        data array
    """
    data_block_length = parent_block['length'] - 24
    if dataRead_available:
        try:
            buf = unsorted_data_read4(record, info, bytes(parent_block['data']), record_id_size, data_block_length)
            _apply_unsorted_bit_masking(buf, record, info)
            return buf
        except Exception as e:
            warn('data_read cython module - unsorted_data_read4 function crashed, using python based parsing backup')
    # initialise data structure
    # key is channel name
    buf = {}
    VLSD = {}
    pos_byte_beg = {}
    pos_byte_end = {}
    numpy_format = {}
    # key is record id
    index = {}
    CGrecordLength = {}
    VLSD_flag = {}
    VLSD_CG_name = {}
    VLSD_CG_signal_data_type = {}
    channel_name_set = {}
    for record_id in record:
        if record[record_id]['record'].Flags & 0b100001:  # VLSD (bit 0) or VLSC compact (bit 5)
            VLSD_flag[record_id] = True
            VLSD[record[record_id]['record'].VLSD_CG[record_id]['channelName']] = []
            VLSD_CG_name[record_id] = record[record_id]['record'].VLSD_CG[record_id]['channelName']
            VLSD_CG_signal_data_type[record_id] = record[record_id]['record'].VLSD_CG[record_id]['channel'].signal_data_type(
                info)
        else:
            VLSD_flag[record_id] = False
            for Channel in record[record_id]['record'].values():
                # if not Channel.VLSD_CG_Flag:
                buf[Channel.name] = empty((record[record_id]['record'].numberOfRecords,),
                                          dtype='V{}'.format(Channel.nBytes_aligned))
                numpy_format[Channel.name] = Channel.data_format(info)
                pos_byte_beg[Channel.name] = record_id_size + \
                    Channel.byteOffset
                pos_byte_end[Channel.name] = pos_byte_beg[Channel.name] + \
                    Channel.nBytes_aligned
            index[record_id] = 0
            CGrecordLength[record_id] = record[record_id]['record'].CGrecordLength
            channel_name_set[record_id] = record[record_id]['record'].channelNames.copy(
            )
    position = 0
    record_id_c_format = record[next(iter(record))]['record'].recordIDCFormat
    # read data
    while position < data_block_length:
        (record_id,) = record_id_c_format.unpack(
            parent_block['data'][position:position + record_id_size])
        if record_id not in VLSD_flag:  # unknown record_id, skip remaining
            break
        if not VLSD_flag[record_id]:  # not VLSD CG
            # list of channel classes from channelSet
            for channel_name in channel_name_set[record_id]:
                buf[channel_name][index[record_id]] = \
                    parent_block['data'][position + pos_byte_beg[channel_name]:position + pos_byte_end[channel_name]]
            index[record_id] += 1
            position += CGrecordLength[record_id]
        else:  # VLSD CG
            position += record_id_size
            (VLSDLen,) = _VLSDStruct.unpack(
                parent_block['data'][position:position + 4])  # VLSD length
            position += 4
            temp = parent_block['data'][position:position + VLSDLen - 1]
            signal_data_type = VLSD_CG_signal_data_type[record_id]
            if signal_data_type == 6:
                temp = temp.decode('ISO8859')
            elif signal_data_type == 7:
                temp = temp.decode('utf-8')
            elif signal_data_type == 8:
                temp = temp.decode('<utf-16')
            elif signal_data_type == 9:
                temp = temp.decode('>utf-16')
            VLSD[VLSD_CG_name[record_id]].append(temp)
            position += VLSDLen
    # changing from bytes type to desired type
    if buf:
        for name in buf.keys():
            buf[name] = buf[name].view(dtype=numpy_format[name])
    # convert list to array for VLSD only
    if VLSD:
        for channel_name in VLSD:
            VLSD[channel_name] = array(VLSD[channel_name])
        buf.update(VLSD)
    return buf


def _read_sd_block(signal_data_type, sd_block, sd_block_length, n_records, pointer):
    """ Reads vlsd channel from its SD Block bytes

    Parameters
    ----------------
    signal_data_type : int

    sd_block : bytes
        SD Block bytes

    sd_block_length: uint64
        SD Block data length (header not included)

    n_records: int
        number of records

    pointer: numpy array of uint
        position of records in SD block

    Returns
    -----------
    array
    """
    if dataRead_available:
        try:
            return sd_data_read(signal_data_type, memoryview(sd_block).tobytes(),
                                sd_block_length, n_records)
        except Exception as e:
            warn('data_read cython module - sd_data_read function crashed, using python based parsing backup')
    VLSDLen = maximum(diff(pointer.astype('<i8')) - 4, 0)
    last_len = max(0, int(sd_block_length) - int(pointer[-1]) - 4)
    VLSDLen = concatenate((VLSDLen,
                           array([last_len], dtype='<i8')), axis=0)
    max_len = int(max(VLSDLen))
    if max_len > sd_block_length:  # pointer values are out of range (e.g. variant/union overlap)
        max_len = 0
    if max_len > 0:
        if signal_data_type < 10:
            if signal_data_type == 6:
                channel_format = 'ISO8859'
            elif signal_data_type == 7:
                channel_format = 'utf-8'
            elif signal_data_type == 8:
                channel_format = '<utf-16'
            elif signal_data_type == 9:
                channel_format = '>utf-16'
            else:
                channel_format = 'utf-8'
                warn('signal_data_type should have fixed length')
            output = zeros((n_records,), dtype="U{:d}".format(max_len))
            for index, position in enumerate(pointer):
                position = int(position + 4)
                output[index] = sd_block[position:int(
                    position + VLSDLen[index])].decode(channel_format).rstrip('\x00')
        else:  # byte arrays or mime types
            output = empty((n_records,), dtype="V{:d}".format(max_len))
            for index, position in enumerate(pointer):
                position = int(position + 4)
                output[index] = bytearray(sd_block[position:int(
                    position + VLSDLen[index])]).rjust(max_len,  b'\x00')
        return output
    else:
        warn('VLSD channel is empty')
        return None


def _read_vd_block_data(fid, pointer):
    """Read raw bytes from a VD/DT/DZ/DL/HL block at pointer."""
    if not pointer:
        return b''
    h = _load_header(fid, pointer)
    if not h:
        return b''
    block_id = h['id']
    if block_id in (b'##VD', '##VD', b'##DT', '##DT'):
        # raw data immediately after header (link_count should be 0)
        return fid.read(h['length'] - 24)
    elif block_id in _DZ_BLOCK_IDS:
        dz = DZBlock()
        dz.read_dz(fid)
        raw = fid.read(dz['dz_data_length'])
        return bytes(DZBlock.decompress_data_block(
            raw, dz['dz_zip_type'], dz['dz_zip_parameter'], dz['dz_org_data_length']))
    elif block_id in (b'##HL', '##HL'):
        hl = HLBlock()
        hl.read_hl(fid)
        return _read_vd_block_data(fid, hl['hl_dl_first'])
    elif block_id in frozenset((b'##DL', '##DL')):
        result = bytearray()
        next_ptr = structunpack('<Q', fid.read(8))[0]
        data_ptrs = structunpack('<{}Q'.format(h['link_count'] - 1),
                                 fid.read(8 * (h['link_count'] - 1)))
        for ptr in data_ptrs:
            if ptr:
                result.extend(_read_vd_block_data(fid, ptr))
        while next_ptr:
            h2 = _load_header(fid, next_ptr)
            if not h2:
                break
            next_ptr = structunpack('<Q', fid.read(8))[0]
            data_ptrs2 = structunpack('<{}Q'.format(h2['link_count'] - 1),
                                      fid.read(8 * (h2['link_count'] - 1)))
            for ptr in data_ptrs2:
                if ptr:
                    result.extend(_read_vd_block_data(fid, ptr))
        return bytes(result)
    return b''


def _decode_vlsc_data(vd_bytes, offsets, sizes, signal_data_type):
    """Decode VLSC string data from VD block bytes using offset/size arrays."""
    n_records = len(offsets)
    if n_records == 0:
        return None
    if dataRead_available:
        try:
            return vd_data_read(signal_data_type, vd_bytes, offsets.astype('<u8'),
                                sizes.astype('<u8'), n_records)
        except Exception:
            pass
    max_len = int(max(sizes)) if len(sizes) > 0 else 0
    if max_len == 0:
        return None
    if signal_data_type == 6:
        encoding = 'ISO-8859-1'
    elif signal_data_type == 7:
        encoding = 'utf-8'
    elif signal_data_type == 8:
        encoding = 'utf-16-le'
    elif signal_data_type == 9:
        encoding = 'utf-16-be'
    elif signal_data_type == 17:
        # BOM-per-value: detect encoding from each value's BOM
        output = zeros((n_records,), dtype='U{}'.format(max_len))
        for i in range(n_records):
            size = int(sizes[i])
            if size > 0:
                offset = int(offsets[i])
                raw = vd_bytes[offset:offset + size]
                if len(raw) >= 3 and raw[:3] == b'\xef\xbb\xbf':
                    s = raw[3:].decode('utf-8', errors='replace')
                elif len(raw) >= 2 and raw[:2] == b'\xff\xfe':
                    s = raw[2:].decode('utf-16-le', errors='replace')
                elif len(raw) >= 2 and raw[:2] == b'\xfe\xff':
                    s = raw[2:].decode('utf-16-be', errors='replace')
                else:
                    s = raw.decode('utf-8', errors='replace')
                output[i] = s.rstrip('\x00')
        return output
    else:
        encoding = 'utf-8'
    output = zeros((n_records,), dtype='U{}'.format(max_len))
    for i in range(n_records):
        size = int(sizes[i])
        if size > 0:
            offset = int(offsets[i])
            raw = vd_bytes[offset:offset + size]
            try:
                output[i] = raw.decode(encoding).rstrip('\x00')
            except (UnicodeDecodeError, ValueError):
                output[i] = raw.decode(encoding, errors='replace').rstrip('\x00')
    return output


def _calculate_aligned_offset(current_bit_pos, alignment_offset, cn_alignment):
    """MDF 4.3 Section 6.25.3 alignment calculation.

    Find smallest byte offset >= ceil(current_bit_pos/8) where
    (offset - alignment_offset) mod 2^cn_alignment == 0.

    Returns byte offset (or bit position if cn_alignment == 255 for bit-packed).
    """
    if cn_alignment == 255:  # bit-packed: no byte alignment
        return current_bit_pos
    current_byte_pos = (current_bit_pos + 7) // 8
    if cn_alignment == 0:  # 1-byte aligned (2^0 = 1)
        return current_byte_pos
    alignment_value = 1 << cn_alignment
    if current_byte_pos <= alignment_offset:
        return alignment_offset
    relative = current_byte_pos - alignment_offset
    remainder = relative % alignment_value
    if remainder == 0:
        return current_byte_pos
    return current_byte_pos + (alignment_value - remainder)


def _decode_datastream(data_bytes, ds_info, cn_info_map, cn_name_map, n_records):
    """Decode inline data stream records using DS/CL composition hierarchy.

    Parameters
    ----------
    data_bytes : bytes
        Raw bytes for a single record's data stream blob.
    ds_info : dict
        DSBlock dict (from mdfinfo4). Must have 'ds_cn_composition' pointer.
    cn_info_map : dict
        {pointer: cn_dict} for all channels in the channel group.
    cn_name_map : dict
        {pointer: channel_name} for looking up resolved channel names.
    n_records : int
        Ignored here — caller iterates per-record (this decodes one blob).

    Returns
    -------
    dict
        {channel_name: decoded_value_bytes} for this record.
    """
    result = {}

    def _decode_compo(data, state, compo_ptr):
        """Recursively decode a composition starting at compo_ptr.
        state is [bit_pos, align_offset] — mutable list to allow updates.
        """
        if not compo_ptr or compo_ptr not in cn_info_map:
            return

        cn = cn_info_map[compo_ptr]
        name = cn_name_map.get(compo_ptr)

        cn_flags = cn.get('cn_flags', 0)
        if cn_flags & (1 << 18):  # CN_F_ALIGNMENT_RESET
            state[1] = (state[0] + 7) // 8

        cn_align = cn.get('cn_alignment', 0)
        byte_off = _calculate_aligned_offset(state[0], state[1], cn_align)
        final_byte = byte_off + cn.get('cn_byte_offset', 0)
        bit_count = cn.get('cn_bit_count', 0)

        if bit_count > 0:
            byte_count = (bit_count + 7) // 8
            if name and final_byte + byte_count <= len(data):
                result[name] = data[final_byte:final_byte + byte_count]
            state[0] = (final_byte + byte_count) * 8
        # else: variable-length — skip for now (handled by VLSC path)

        # Follow cn_cn_next chain
        next_ptr = cn.get('cn_cn_next', 0)
        if next_ptr:
            _decode_compo(data, state, next_ptr)

    ds_comp_ptr = ds_info.get('ds_cn_composition', 0)
    if ds_comp_ptr:
        _decode_compo(data_bytes, [0, 0], ds_comp_ptr)
    return result


class Data(dict):
    __slots__ = ['fid', 'pointer_to_data', 'type']
    """ Data class is organizing record classes itself made of channel class.
    This class inherits from dict. Keys are corresponding to channel group recordID
    A Dataclass corresponds to a data block, a dict of record classes (one per channel group)
    Each record class contains a list of channel class representing the structure of channel record.

    Attributes
    --------------
    fid : io.open
        file identifier
    pointer_to_data : int
        position of Data block in mdf file
    type : str
        'sorted' or 'unsorted' data block

    Methods
    ------------
    add_record(record)
        Adds a new record in DATA class dict
    read(channel_set, zip=None)
        Reads data block
    load(record, zip=None, name_list=None)
        Reads sorted data block from record definition
    read_record(recordID, buf, channel_set=None):
        read record from a buffer
    """

    def __init__(self, fid, pointer):
        """ Constructor

        Parameters
        ----------------
        fid : float
            file identifier
        pointer : int
            position of data block in file
        """
        self.fid = fid
        self.pointer_to_data = pointer
        self.type = 'sorted'

    def add_record(self, record):
        """Adds a new record in Data class dict.

        Parameters
        ----------------
        record: class
            channel group definition listing record channel classes
        """
        if record.recordID not in self:  # do not overwrite already existing record in Data, for VLSD_CG
            self[record.recordID] = {}
            self[record.recordID]['record'] = record
            # detects VLSD CG
            for recordID in self[record.recordID]['record'].VLSD_CG:
                if recordID not in self:  # create new recordId
                    self[recordID] = {}
                    self[recordID]['record'] = record
                self[recordID]['record'].VLSD_CG = self[record.recordID]['record'].VLSD_CG
        else:  # VLSD_CG
            record.VLSD_CG = self[record.recordID]['record'].VLSD_CG
            self[record.recordID]['record'] = record

    def read(self, channel_set, info, filename):
        """Reads data block

        Parameters
        ----------------
        channel_set : set of str
            set of channel names
        info : info object
            contains blocks structures
        filename
            name of file ot read
        """
        # checks if file is closed
        if self.fid is None or self.fid.closed:
            self.fid = open(filename, 'rb')
        if len(self) == 1:  # sorted dataGroup
            recordID = next(iter(self))
            record = self[recordID]['record']
            self[recordID]['data'], self[recordID]['invalid_data'] = self.load(record,
                                                                               info, name_list=channel_set,
                                                                               sorted_flag=True)
            if record.VLSD:  # VLSD channels exist
                self[recordID]['VLSD'] = {}
                for cn in record.VLSD:  # VLSD channels
                    if channel_set is None or record[cn].name in channel_set:
                        # all channels
                        temp = Data(self.fid, record[cn].data(info))
                        # vector of each record position
                        pointer = self[recordID]['data'][record[cn].name]
                        temp, invalid = temp.load(record, info, name_list=channel_set, sorted_flag=True,
                                                  vlsd=pointer.view(dtype='<u{}'.format(pointer.dtype.itemsize)))
                        if temp is not None:
                            # change channel name by appending offset
                            self[recordID]['data'] = rename_fields(self[recordID]['data'],
                                                                   {record[cn].name: '{}_offset'.format(record[cn].name)})
                            self[recordID]['VLSD'][record[cn].name] = temp
            if record.DS:  # DS-composition channels (data stream mode, MDF 4.3)
                # Build a map from file pointer to cn_info dict and to resolved name
                dg = record.dataGroup
                cg = record.channelGroup
                cn_by_ptr = {info['CN'][dg][cg][k].get('pointer', 0): info['CN'][dg][cg][k]
                             for k in info['CN'][dg][cg]}
                name_by_ptr = {info['CN'][dg][cg][k].get('pointer', 0): info['CN'][dg][cg][k].get('name', '')
                               for k in info['CN'][dg][cg]}
                # Find the parent VLSD/VLSC channel whose composition block is a DSBlock
                for parent_cn_num in info['CN'][dg][cg]:
                    parent_info = info['CN'][dg][cg][parent_cn_num]
                    if 'DSBlock' not in parent_info:
                        continue
                    ds_block = parent_info['DSBlock']
                    # Read all VLSD blobs for this parent channel
                    vd_data = _read_vd_block_data(self.fid, parent_info.get('cn_data', 0))
                    if not vd_data:
                        continue
                    # VLSD: the fixed record stores the byte offset into vd_data per record
                    parent_name = parent_info.get('name', '')
                    if self[recordID]['data'] is None:
                        continue
                    try:
                        offsets_field = parent_name + '_offset' if (parent_name + '_offset') in self[recordID]['data'].dtype.names else parent_name
                        offsets_arr = self[recordID]['data'][offsets_field]
                        if offsets_arr.dtype.kind in ('S', 'V'):
                            offsets_arr = frombuffer(offsets_arr.tobytes(), dtype='<u8')
                    except (KeyError, ValueError, AttributeError):
                        continue
                    # Also need sizes (from size channel or end-of-blob heuristic)
                    n = len(offsets_arr)
                    sizes_arr = None
                    size_ptr = parent_info.get('cn_cn_size', 0)
                    if size_ptr:
                        for k in info['CN'][dg][cg]:
                            if info['CN'][dg][cg][k].get('pointer') == size_ptr:
                                sname = info['CN'][dg][cg][k].get('name', '')
                                try:
                                    sizes_arr = self[recordID]['data'][sname].astype('<u8')
                                except (KeyError, ValueError):
                                    pass
                                break
                    if sizes_arr is None and n > 0:
                        # Estimate sizes from consecutive offsets
                        sizes_arr = zeros(n, dtype='<u8')
                        for i in range(n - 1):
                            sizes_arr[i] = offsets_arr[i + 1] - offsets_arr[i]
                        sizes_arr[n - 1] = len(vd_data) - offsets_arr[n - 1]
                    if sizes_arr is None:
                        continue
                    # Decode each blob
                    decoded_per_channel = {}  # {name: list_of_values}
                    for i in range(n):
                        off = int(offsets_arr[i])
                        sz = int(sizes_arr[i])
                        if sz <= 0 or off + sz > len(vd_data):
                            # empty blob — append None/empty for each DS channel
                            for ds_cn_num in record.DS:
                                ch_name = record[ds_cn_num].name
                                decoded_per_channel.setdefault(ch_name, []).append(None)
                            continue
                        blob = vd_data[off:off + sz]
                        rec_result = _decode_datastream(blob, ds_block, cn_by_ptr, name_by_ptr, 1)
                        for ds_cn_num in record.DS:
                            ch_name = record[ds_cn_num].name
                            val = rec_result.get(ch_name)
                            decoded_per_channel.setdefault(ch_name, []).append(val)
                    # Store decoded arrays in VLSD dict
                    self[recordID].setdefault('VLSD', {})
                    for ch_name, vals in decoded_per_channel.items():
                        if any(v is not None for v in vals):
                            self[recordID]['VLSD'][ch_name] = vals

            if record.VLSC:  # VLSC channels exist (MDF 4.3)
                self[recordID].setdefault('VLSD', {})
                for cn in record.VLSC:  # VLSC channels
                    if channel_set is None or record[cn].name in channel_set:
                        cn_info = info['CN'][record[cn].dataGroup][record[cn].channelGroup][record[cn].channelNumber]
                        vd_pointer = cn_info.get('cn_data', 0)
                        cn_size_pointer = cn_info.get('cn_cn_size', 0)
                        # find the size channel by matching its CN block pointer
                        size_channel_name = None
                        if cn_size_pointer:
                            for other_cn_num in info['CN'][record[cn].dataGroup][record[cn].channelGroup]:
                                if info['CN'][record[cn].dataGroup][record[cn].channelGroup][other_cn_num].get('pointer') == cn_size_pointer:
                                    size_channel_name = info['CN'][record[cn].dataGroup][record[cn].channelGroup][other_cn_num]['name']
                                    break
                        if size_channel_name is None or self[recordID]['data'] is None:
                            continue
                        try:
                            sizes = self[recordID]['data'][size_channel_name]
                            offsets = self[recordID]['data'][record[cn].name]
                        except (ValueError, IndexError, KeyError):
                            continue
                        # reinterpret byte-string field as uint64 offset
                        if hasattr(offsets, 'view') and offsets.dtype.kind in ('S', 'V'):
                            offsets = frombuffer(offsets.tobytes(), dtype='<u8')
                        if hasattr(sizes, 'view') and sizes.dtype.kind in ('S', 'V'):
                            sizes = frombuffer(sizes.tobytes(), dtype='<u8')
                        vd_bytes = _read_vd_block_data(self.fid, vd_pointer)
                        temp = _decode_vlsc_data(vd_bytes, offsets, sizes,
                                                 record[cn].signal_data_type(info))
                        if temp is not None:
                            self[recordID]['data'] = rename_fields(self[recordID]['data'],
                                                                   {record[cn].name: '{}_offset'.format(record[cn].name)})
                            self[recordID]['VLSD'][record[cn].name] = temp
        else:  # unsorted DataGroup
            self.type = 'unsorted'
            data, invalid = self.load(
                self, info, name_list=channel_set, sorted_flag=False)
            for recordID in self:
                self[recordID]['data'] = {}
                for channel in self[recordID]['record'].values():
                    self[recordID]['data'][channel.name] = data[channel.name]

    def load(self, record, info, name_list=None, sorted_flag=True, vlsd=None):
        """Reads data block from record definition

        Parameters
        ----------------
        record : class
            channel group definition listing record channel classes
        info : class
            contains blocks
        name_list : list of str, optional
            list of channel names
        sorted_flag : bool, optional
            flag to know if data block is sorted (only one Channel Group in block)
            or unsorted (several Channel Groups identified by a recordID).
            As unsorted block can contain CG records in random order, block
            is processed iteratively, not in raw like sorted -> much slower reading
        vlsd : array or None
            indicate a sd block, compressed (DZ) or not (SD)

        Returns
        -----------
        numpy recarray of data
        """
        temps = defaultdict()
        temps['invalid_data'] = None
        # block header
        temps.update(_load_header(self.fid, self.pointer_to_data))
        if temps['id'] in _DV_BLOCK_IDS:
            # to be optimised by using unpack in case of column oriented storage (only one channel)
            temps['data'] = record.read_sorted_record(
                self.fid, info, channel_set=name_list)
        elif temps['id'] in frozenset((b'##DL', b'##LD', '##DL', '##LD')):  # data list block
            if temps['id'] in frozenset((b'##DL', '##DL')):
                temp = DLBlock()
                temp.read_dl(self.fid, temps['link_count'])
            else:
                temp = LDBlock()
                temp.read_ld(self.fid, temps['link_count'])
            temps.update(temp)
            if temps['next']:
                index = 1
            # reads pointers to all data blocks (DT, RD, SD, DZ)
            while temps['next']:
                temp = defaultdict()
                temp.update(_load_header(self.fid, temps['next']))
                (temps['next'], ) = structunpack('<Q', self.fid.read(8))
                temps['list_data'][index] = \
                    structunpack('<{}Q'.
                                 format(temp['link_count'] - 1),
                                 self.fid.read(8 * (temp['link_count'] - 1)))
                index += 1
            if temps['count']:
                # read and concatenate raw blocks
                if vlsd is not None or not sorted_flag:
                    # need to load all blocks as variable length, cannot process block by block
                    data_block = defaultdict()
                    data_block['data'] = bytearray()
                    data_block_length = 0
                    for DL in temps['list_data']:
                        for pointer in temps['list_data'][DL]:
                            # read fist data blocks linked by DLBlock to identify data block type
                            data_block.update(_load_header(self.fid, pointer))
                            if data_block['id'] in frozenset((b'##SD', b'##DT', '##SD', '##DT')):
                                data_block['data'].extend(
                                    self.fid.read(data_block['length'] - 24))
                                data_block_length += data_block['length'] - 24
                            elif data_block['id'] in _DZ_BLOCK_IDS:
                                temp = DZBlock()
                                temp.read_dz(self.fid)
                                data_block.update(temp)
                                data_block['data'].extend(DZBlock.decompress_data_block(
                                    self.fid.read(
                                        data_block['dz_data_length']),
                                    data_block['dz_zip_type'],
                                    data_block['dz_zip_parameter'],
                                    data_block['dz_org_data_length']))
                                if isinstance(data_block['dz_org_block_type'], str):
                                    data_block['id'] = '##{}'.format(
                                        data_block['dz_org_block_type'])
                                else:
                                    data_block['id'] = '##{}'.format(
                                        data_block['dz_org_block_type'].decode('ASCII'))
                                data_block_length += data_block['dz_org_data_length']
                    data_block['length'] = data_block_length + 24
                    temps['data'] = _data_block(record, info, parent_block=data_block, channel_set=name_list,
                                                n_records=None, sorted_flag=sorted_flag, vlsd=vlsd)
                else:
                    temps['data'] = self.read_data_list('list_data', record.CGrecordLength, temps, record, info,
                                                        name_list, sorted_flag, vlsd)
                    if temps['id'] in (b'##LD', '##LD'):
                        try:  # invalid bytes in DIBlock
                            temps['invalid_data'] = self.read_data_list('inval_data',
                                                                        record.invalid_channel.nBytes_aligned,
                                                                        temps, record, info,
                                                                        name_list, sorted_flag, vlsd)
                        except (KeyError, AttributeError):
                            pass
            else:  # empty datalist
                temps['data'] = None
        elif temps['id'] in (b'##HL', '##HL'):  # header list block for DZBlock
            # link section
            temp = HLBlock()
            temp.read_hl(self.fid)
            temps.update(temp)
            self.pointer_to_data = temps['hl_dl_first']
            temps['data'], temps['invalid_data'] = self.load(record, info, name_list=name_list,
                                                             sorted_flag=sorted_flag, vlsd=vlsd)
        elif temps['id'] in frozenset((b'##DT', b'##RD', '##DT', '##RD')):
            if sorted_flag:  # normal sorted data block, direct read
                temps['data'] = record.read_sorted_record(
                    self.fid, info, channel_set=name_list)
            else:  # VLSD_CG
                temps['data'] = self.fid.read(temps['length'] - 24)
                temps['data'] = _data_block(record, info, parent_block=temps, channel_set=name_list, n_records=None,
                                            sorted_flag=sorted_flag)
        elif temps['id'] in _SD_BLOCK_IDS:  # VLSD, not CG
            temps['data'] = self.fid.read(temps['length'] - 24)
            temps['data'] = _data_block(record, info, parent_block=temps, channel_set=name_list, n_records=None,
                                        sorted_flag=sorted_flag, vlsd=vlsd)
        elif temps['id'] in _DZ_BLOCK_IDS:  # zipped data block
            temp = DZBlock()
            temp.read_dz(self.fid)
            temps.update(temp)
            temps['data'] = self.fid.read(temps['dz_data_length'])
            temps['length'] = temps['dz_org_data_length'] + 24
            temps['data'] = _data_block(record, info, parent_block=temps, channel_set=name_list, n_records=None,
                                        sorted_flag=sorted_flag, vlsd=vlsd)
        elif temps['id'] in (b'##GD', '##GD'):  # guard block (MDF 4.3)
            (gd_link,) = structunpack('<q', self.fid.read(8))  # pointer to guarded block
            (gd_version,) = structunpack('<H', self.fid.read(2))  # minimum MDF version
            if gd_version <= 430 and gd_link:
                self.pointer_to_data = gd_link
                temps['data'], temps['invalid_data'] = self.load(record, info,
                                                                  name_list=name_list,
                                                                  sorted_flag=sorted_flag,
                                                                  vlsd=vlsd)
            else:
                temps['data'] = None
        else:
            raise Exception('unknown data block')
        return temps['data'], temps['invalid_data']

    def read_record(self, record_id, info, buf):
        """ read record from a buffer

        Parameters
        ----------------
        record_id : int
            record identifier
        info : class
            contains blocks
        buf : str
            buffer of data from file to be converted to channel raw data
        """
        return self[record_id]['record'].read_record_buf(buf, info)

    def read_data_list(self, field, nBytes, temps, record, info, name_list, sorted_flag, vlsd):
        previous_index = 0
        data = None
        data_block = defaultdict()
        data_block['data'] = bytearray()
        for DL in temps[field]:
            for pointer in temps[field][DL]:
                # read fist data blocks linked by DLBlock to identify data block type
                data_block.update(_load_header(self.fid, pointer))
                if data_block['id'] in _DATA_BLOCK_IDS:
                    data_block['data'].extend(
                        self.fid.read(data_block['length'] - 24))
                elif data_block['id'] in _DZ_BLOCK_IDS:
                    temp = DZBlock()
                    temp.read_dz(self.fid)
                    data_block.update(temp)
                    data_block['data'].extend(DZBlock.decompress_data_block(
                        self.fid.read(data_block['dz_data_length']),
                        data_block['dz_zip_type'],
                        data_block['dz_zip_parameter'],
                        data_block['dz_org_data_length']))
                    if isinstance(data_block['dz_org_block_type'], str):
                        data_block['id'] = '##{}'.format(
                            data_block['dz_org_block_type'])
                    else:
                        data_block['id'] = '##{}'.format(
                            data_block['dz_org_block_type'].decode('ASCII'))
                nrecord_chunk = len(data_block['data']) // nBytes
                nremain = len(data_block['data']) % nBytes
                if nremain:
                    remain = data_block['data'][-nremain:]
                    del data_block['data'][-nremain:]
                if previous_index + nrecord_chunk > record.numberOfRecords:
                    # there could be more data than needed for the expected number of records
                    nrecord_chunk = record.numberOfRecords - previous_index
                tmp = _data_block(record, info, parent_block=data_block, channel_set=name_list,
                                  n_records=nrecord_chunk, sorted_flag=sorted_flag, vlsd=vlsd)
                if tmp is None or not hasattr(tmp, 'dtype'):
                    continue
                if not previous_index:  # initialise recarray
                    data = recarray(record.numberOfRecords, dtype=tmp.dtype)
                data[previous_index: previous_index + nrecord_chunk] = tmp
                previous_index += nrecord_chunk
                if nremain:
                    data_block['data'] = remain
                else:
                    data_block['data'] = bytearray()  # flush
        return data


class Record(dict):
    __slots__ = ['CGrecordLength', 'recordLength', 'numberOfRecords', 'recordID',
                 'recordIDsize', 'recordIDCFormat', 'dataGroup', 'channelGroup',
                 'numpyDataRecordFormat', 'dataRecordName', 'master',
                 'recordToChannelMatching', 'channelNames', 'Flags', 'VLSD_CG',
                 'VLSD', 'VLSC', 'DS', 'MLSD', 'byte_aligned', 'hiddenBytes', 'invalid_channel',
                 'CANOpen', 'unique_channel_in_DG']
    """ Record class listing channel classes. It is representing a channel group

    Attributes
    --------------
    CGrecordLength : int
        length of record corresponding of channel group in Byte CG Block information
    recordLength : int
        length of record as understood by program based on C datatypes
    numberOfRecords : int
        number of records in data block
    recordID : int
        recordID corresponding to channel group
    recordIDsize : int
        size of recordID
    recordIDCFormat : str
        record identifier C format string as understood by fread
    dataGroup : int:
        data group number
    channelGroup : int
        channel group number
    numpyDataRecordFormat : list
        list of numpy (dtype) for each channel
    dataRecordName : list
        list of channel names used for recarray attribute definition
    master : dict
        define name and number of master channel
    recordToChannelMatching : dict
        helps to identify nested bits in byte
    channelNames : set
        channel names to be stored, useful for low memory consumption but slow
    Flags : bool
        channel flags as from specification
    VLSD_CG : dict
        dict of Channel Group VLSD, key being recordID
    VLSD : list of channel classes
        list of channel classes being VLSD
    MLSD : dict
        copy from info['MLSD'] if existing
    byte_aligned : Bool, True by default
        flag for byte aligned record
    hiddenBytes : Bool, False by default
        flag in case of non declared channels in record, forces to use readBitarray
    invalid_channel : Default None
        invalid_byte class if existing in record otherwise None
    CANOpen : str, Default None
        'time' if record contains CANOpen time channel, same for 'date'

    Methods
    ------------
    add_channel(info, channelNumber)
    load_info(info)
    readSortedRecord(fid, pointer, info, channelSet=None)
    generate_chunks()
    read_all_channels_sorted_record(fid)
    read_not_all_channels_sorted_record(fid, info, channelSet)
    readRecordBuf(buf, info, channelSet=None)
    initialise_recarray(info, channel_set, nrecords, dtype=None, channels_indexes=None)
    read_channels_from_bytes(bita, info, channelSet=None, nrecords=None, dtype=None, channels_indexes=None)
    read_channels_from_bytes_fallback(bita, info, channel_set=None, nrecords=None, dtype=None, channels_indexes=None)
    """

    def __init__(self, data_group, channel_group):
        self.CGrecordLength = 0
        self.recordLength = 0
        self.numberOfRecords = 0
        self.recordID = 0
        self.recordIDsize = 0
        self.recordIDCFormat = ''
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.numpyDataRecordFormat = []
        self.dataRecordName = []
        self.master = 'master_{}'.format(data_group)
        self.Flags = 0
        self.VLSD_CG = {}
        self.VLSD = []
        self.VLSC = []
        self.DS = []
        self.MLSD = {}
        self.recordToChannelMatching = {}
        self.channelNames = set()
        self.byte_aligned = True
        self.hiddenBytes = False
        self.invalid_channel = None
        self.CANOpen = None
        self.unique_channel_in_DG = False

    def __str__(self):
        output = list()
        output.append(
            'Record class content print\nTotal number of channels : {}\n'.format(len(self)))
        for chan in self:
            output.append(chan.name)
            output.append('  Type ')
            output.append(chan.type)
            output.append(' DG {} '.format(chan.dataGroup))
            output.append('CG {} '.format(chan.channelGroup))
            output.append('CN {} '.format(chan.channelNumber))
            output.append('VLSD {} \n'.format(chan.VLSD_CG_Flag))
        output.append('CG block record bytes length : {}\nData group number : {}'
                      '\nByte aligned : {}\nHidden bytes : {}\n'.format(self.CGrecordLength, self.dataGroup,
                                                                        self.byte_aligned, self.hiddenBytes))
        if self.master is not None:
            output.append('Master channel : {}\n'.format(self.master))
        output.append('Numpy records format : \n')
        for record in self.numpyDataRecordFormat:
            output.append(' {} '.format(record))
        output.append('\nVLSD_CG {}\n'.format(self.VLSD_CG))
        output.append('\nCG_flags {}\n'.format(self.Flags))
        return ''.join(output)

    def add_channel(self, info, channel_number):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class
        """
        channel = Channel4(self.dataGroup, self.channelGroup, channel_number)
        self.append(channel.set(info))
        self.channelNames.add(self[-1].name)

    def load_info(self, info):
        """ gathers records related from info class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        """
        self.CGrecordLength = info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        self.recordIDsize = info['DG'][self.dataGroup]['dg_rec_id_size']
        if not self.recordIDsize == 0:  # no record ID
            self.dataRecordName.append('RecordID{}'.format(self.channelGroup))
            if self.recordIDsize == 1:
                self.numpyDataRecordFormat.append('uint8')
                self.recordIDCFormat = Struct('B')
                self.recordLength += 1
                self.CGrecordLength += 1
            elif self.recordIDsize == 2:
                self.numpyDataRecordFormat.append('uint16')
                self.recordIDCFormat = Struct('H')
                self.recordLength += 2
                self.CGrecordLength += 2
            elif self.recordIDsize == 4:
                self.numpyDataRecordFormat.append('uint32')
                self.recordIDCFormat = Struct('I')
                self.recordLength += 4
                self.CGrecordLength += 4
            elif self.recordIDsize == 8:
                self.numpyDataRecordFormat.append('uint64')
                self.recordIDCFormat = Struct('L')
                self.recordLength += 8
                self.CGrecordLength += 8
        self.recordID = info['CG'][self.dataGroup][self.channelGroup]['cg_record_id']
        self.numberOfRecords = info['CG'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        self.Flags = info['CG'][self.dataGroup][self.channelGroup]['cg_flags']
        if 'MLSD' in info:
            self.MLSD = info['MLSD']
        self.unique_channel_in_DG = info['DG'][self.dataGroup]['unique_channel_in_DG']
        embedding_channel = None
        prev_chan = Channel4(self.dataGroup, self.channelGroup, 0)
        for channelNumber in sorted(info['CN'][self.dataGroup][self.channelGroup].keys()):
            cn_info = info['CN'][self.dataGroup][self.channelGroup][channelNumber]
            if cn_info.get('cn_flags', 0) & 0x20000:  # CN_F_DATA_STREAM_MODE: not in fixed record
                self.DS.append(channelNumber)
                channel = Channel4(self.dataGroup, self.channelGroup, channelNumber)
                channel.set(info)
                self[channelNumber] = channel
                self.channelNames.add(channel.name)
                continue
            channel = Channel4(
                self.dataGroup, self.channelGroup, channelNumber)
            channel.set(info)
            channel_type = channel.channel_type(info)
            data_format = channel.data_format(info)
            self.master = info['masters'][info['CN'][self.dataGroup]
                                          [self.channelGroup][channelNumber]['masterCG']]['name']
            if channel_type in (0, 2, 4, 5):  # not virtual channel
                signal_data_type = channel.signal_data_type(info)
                if signal_data_type == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        # new object otherwise only modified
                        channel = Channel4(
                            self.dataGroup, self.channelGroup, channelNumber)
                        channel.set_CANOpen(info, name)
                        self[channelNumber] = channel
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(
                            channel.data_format(info))
                    self.recordLength += 7
                    self.CANOpen = 'date'
                    embedding_channel = None
                elif signal_data_type == 14:
                    for name in ('ms', 'days'):
                        channel = Channel4(
                            self.dataGroup, self.channelGroup, channelNumber)
                        channel.set_CANOpen(info, name)
                        self[channelNumber] = channel
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(
                            channel.data_format(info))
                    self.recordLength += 6
                    self.CANOpen = 'time'
                    embedding_channel = None
                else:
                    self[channelNumber] = channel
                    self.channelNames.add(channel.name)
                    # Checking if several channels are embedded in bytes
                    if len(self) > 1:
                        # all channels are already ordered in record based on byte_offset
                        # and bit_offset so just comparing with previous channel
                        channel_pos_bit_end = channel.pos_bit_end(info)
                        prev_chan_byte_offset = prev_chan.byteOffset
                        prev_chan_n_bytes = prev_chan.nBytes_aligned
                        prev_chan_includes_curr_chan = channel.pos_bit_beg >= 8 * prev_chan_byte_offset \
                            and channel_pos_bit_end <= 8 * (prev_chan_byte_offset + prev_chan_n_bytes)
                        if embedding_channel is not None:
                            embedding_channel_includes_curr_chan = \
                                channel_pos_bit_end <= embedding_channel.pos_byte_end(
                                    info) * 8
                        else:
                            embedding_channel_includes_curr_chan = False
                        if channel.byteOffset >= prev_chan_byte_offset and \
                                channel.pos_bit_beg < 8 * (prev_chan_byte_offset +
                                                           prev_chan_n_bytes) < channel_pos_bit_end:
                            # not byte aligned
                            self.byte_aligned = False
                        if embedding_channel is not None and \
                                channel_pos_bit_end > embedding_channel.pos_byte_end(info) * 8:
                            embedding_channel = None
                        if prev_chan_includes_curr_chan or \
                                embedding_channel_includes_curr_chan:  # bit(s) in byte(s)
                            if embedding_channel is None and prev_chan_includes_curr_chan:
                                embedding_channel = prev_chan  # new embedding channel detected
                            if self.recordToChannelMatching:  # not first channel
                                self.recordToChannelMatching[channel.name] = \
                                    self.recordToChannelMatching[prev_chan.name]
                            else:  # first channels
                                self.recordToChannelMatching[channel.name] = channel.name
                                self.numpyDataRecordFormat.append(data_format)
                                self.dataRecordName.append(channel.name)
                                self.recordLength += channel.nBytes_aligned
                    if embedding_channel is None:  # adding bytes
                        self.recordToChannelMatching[channel.name] = channel.name
                        self.numpyDataRecordFormat.append(data_format)
                        self.dataRecordName.append(channel.name)
                        self.recordLength += channel.nBytes_aligned

            elif channel_type in (3, 6):  # virtual channel
                # channel calculated based on record index later in conversion function
                self[channelNumber] = channel
                self.channelNames.add(channel.name)
                self.recordToChannelMatching[channel.name] = channel.name
            elif channel_type == 1:  # VLSD channel
                self.VLSD.append(channelNumber)
                self[channelNumber] = channel
                self.channelNames.add(channel.name)
                if 'VLSD_CG' in info:  # is there VLSD CG
                    for recordID in info['VLSD_CG']:  # look for VLSD CG Channel
                        if info['VLSD_CG'][recordID]['cg_cn'] == (self.channelGroup, channelNumber):
                            self.VLSD_CG[recordID] = info['VLSD_CG'][recordID]
                            self.VLSD_CG[recordID]['channel'] = channel
                            self.VLSD_CG[recordID]['channelName'] = channel.name
                            self[channelNumber].VLSD_CG_Flag = True
                            break
            elif channel_type == 7:  # VLSC channel (MDF 4.3)
                self.VLSC.append(channelNumber)
                self[channelNumber] = channel
                self.channelNames.add(channel.name)
                # for compact VLSC structure, cn_data points to VLSC CG (like VLSD_CG)
                if 'VLSD_CG' in info:
                    for recordID in info['VLSD_CG']:
                        if info['VLSD_CG'][recordID]['cg_cn'] == (self.channelGroup, channelNumber):
                            self.VLSD_CG[recordID] = info['VLSD_CG'][recordID]
                            self.VLSD_CG[recordID]['channel'] = channel
                            self.VLSD_CG[recordID]['channelName'] = channel.name
                            self[channelNumber].VLSD_CG_Flag = True
                            break
            prev_chan = channel

        if info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']:  # invalid bytes existing
            cg_data_bytes = info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
            # invalid bytes at the end of the data record
            invalid_bytes = Channel4(
                self.dataGroup, self.channelGroup, cg_data_bytes)
            invalid_bytes.set_invalid_bytes(info)
            self.invalid_channel = invalid_bytes
            self[cg_data_bytes] = self.invalid_channel
            self.channelNames.add(self.invalid_channel.name)
            self.recordToChannelMatching[self.invalid_channel.name] = \
                self.invalid_channel.name
            if info['DG'][self.dataGroup]['data_block_header']['id'] not in (b'##LD', '##LD'):
                # with LDBlock, invalid bytes not in record but in DIBlock
                self.CGrecordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
                self.recordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
                self.numpyDataRecordFormat.append(
                    self.invalid_channel.data_format(info))
                self.dataRecordName.append(self.invalid_channel.name)

        # check for hidden bytes or VLSD within at least one channel
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True
        # check record length consistency
        elif self.CGrecordLength < self.recordLength:
            # forces to use dataRead instead of numpy records.
            self.byte_aligned = False

    def read_sorted_record(self, fid, info, channel_set=None):
        """ reads record, only one channel group per datagroup

        Parameters
        ----------------
        fid :
            file identifier
        info
            info class
        channel_set : set of str, optional
            set of channel to read

        Returns
        -----------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)

        Notes
        --------
        If channelSet is None, read data using numpy.core.records.fromfile that is rather quick.
        However, in case of large file, you can use channelSet to load only interesting channels or
        only one channel on demand, but be aware it might be much slower.
        """
        if channel_set is None and self.byte_aligned and not self.hiddenBytes:
            if self.unique_channel_in_DG:
                return self.read_unique_channel(fid, info)
            else:
                return self.read_all_channels_sorted_record(fid)
        else:  # reads only some channels from a sorted data block
            if channel_set is None or len(channel_set & self.channelNames) > 0:
                if self.unique_channel_in_DG:
                    return self.read_unique_channel(fid, info)
                else:
                    return self.read_not_all_channels_sorted_record(fid, info, channel_set)

    def generate_chunks(self):
        """ calculate data split

        Returns
        --------
        (n_record_chunk, chunk_size)
        """
        n_chunks = (self.CGrecordLength *
                    self.numberOfRecords) // chunk_size_reading + 1
        chunk_length = (self.CGrecordLength * self.numberOfRecords) // n_chunks
        n_record_chunk = chunk_length // self.CGrecordLength
        chunks = [(n_record_chunk, self.CGrecordLength *
                   n_record_chunk)] * n_chunks
        n_record_chunk = self.numberOfRecords - n_record_chunk * n_chunks
        if n_record_chunk > 0:
            chunks.append(
                (n_record_chunk, self.CGrecordLength * n_record_chunk))
        return chunks

    def read_all_channels_sorted_record(self, fid):
        """ reads all channels from file using numpy fromstring, chunk by chunk

        Parameters
        ------------
        fid :
            file identifier

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        chunks = self.generate_chunks()
        previous_index = 0
        buf = recarray(self.numberOfRecords, dtype={'names': self.dataRecordName,
                                                    'formats': self.numpyDataRecordFormat})  # initialise array
        simplefilter('ignore', FutureWarning)
        for n_record_chunk, chunk_size in chunks:
            buf[previous_index: previous_index + n_record_chunk] = \
                frombuffer(fid.read(chunk_size),
                           dtype={'names': self.dataRecordName,
                                  'formats': self.numpyDataRecordFormat},
                           count=n_record_chunk)
            previous_index += n_record_chunk
        return buf

    def read_unique_channel(self, fid, info):
        """ reads all channels from file using numpy fromstring, chunk by chunk

            Parameters
            ------------
            fid :
                file identifier
            info
                info class

            Returns
            --------
            rec : numpy recarray
                contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        nbytes = info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes'] * \
            info['CG'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        dtype = {'names': self.dataRecordName, 'formats': self.numpyDataRecordFormat}
        if nbytes == 0 or not self.dataRecordName:
            fid.read(nbytes)
            return frombuffer(b'', dtype=dtype)
        return frombuffer(fid.read(nbytes), dtype=dtype)

    def read_not_all_channels_sorted_record(self, fid, info, channel_set):
        """ reads channels from file listed in channelSet

        Parameters
        ------------
        fid :
            file identifier
        info: info class
        channel_set : set of str, optional
            set of channel to read

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        chunks = self.generate_chunks()
        previous_index = 0
        if channel_set is None:
            channel_set = self.channelNames
        if channel_set is not None and not self.master in channel_set:
            channel_set.add(self.master)  # adds master channel
        rec, channels_indexes = self.initialise_recarray(
            info, channel_set, self.numberOfRecords)
        if rec is not None:
            if dataRead_available:
                for n_record_chunk, chunk_size in chunks:
                    rec[previous_index: previous_index + n_record_chunk] = \
                        self.read_channels_from_bytes(fid.read(chunk_size), info, channel_set, n_record_chunk,
                                                      rec.dtype, channels_indexes)
                    previous_index += n_record_chunk
                return rec
            else:
                for n_record_chunk, chunk_size in chunks:
                    rec[previous_index: previous_index + n_record_chunk] = \
                        self.read_channels_from_bytes_fallback(fid.read(chunk_size), info, channel_set, n_record_chunk,
                                                               rec.dtype, channels_indexes)
                    previous_index += n_record_chunk
                return rec
        else:
            return []

    def read_record_buf(self, buf, info, channel_set=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        info: class
            contains blocks structure
        channel_set : set of str, optional
            set of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        temp = {}
        if channel_set is None:
            channel_set = self.channelNames
        for Channel in self.values():  # list of channel classes from channelSet
            if Channel.name in channel_set and not Channel.VLSD_CG_Flag:
                temp[Channel.name] = \
                    Channel.c_format_structure(info).unpack(
                        buf[Channel.pos_byte_beg(info):Channel.pos_byte_end(info)])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def initialise_recarray(self, info, channel_set, n_records, dtype=None, channels_indexes=None):
        """ Initialise recarray

        Parameters
        ------------
        info: info class
        channel_set : set of str, optional
            set of channel to read
        n_records: int
            number of records
        dtype: numpy dtype, optional
        channels_indexes: list of int, optional

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        if dtype is not None and channels_indexes is not None:
            return recarray(n_records, dtype=dtype), channels_indexes
        else:
            if channel_set is None:
                channel_set = self.channelNames
            formats = []
            names = []
            channels_indexes = []
            for chan in self:
                if self[chan].name in channel_set and self[chan].channel_type(info) not in (3, 6) \
                        and chan not in self.DS:
                    # not virtual channel, not data stream channel, and part of channelSet
                    channels_indexes.append(chan)
                    formats.append(self[chan].native_data_format(info))
                    names.append(self[chan].name)
            if formats:
                rec = recarray(n_records, dtype={
                               'names': names, 'formats': formats})
                return rec, channels_indexes
            else:
                return None, []

    def read_channels_from_bytes(self, bit_stream, info, channel_set=None, n_records=None, dtype=None,
                                 channels_indexes=None):
        """ reads stream of record bytes using dataRead module if available otherwise bitarray

        Parameters
        ------------
        bit_stream : stream
            stream of bytes
        info: info class
        channel_set : set of str, optional
            set of channel to read
        n_records: int
            number of records
        dtype: numpy dtype
        channels_indexes: list of int

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        if n_records is None:
            n_records = self.numberOfRecords
        # initialise recarray
        if dtype is None:
            buf, channels_indexes = self.initialise_recarray(
                info, channel_set, n_records, dtype, channels_indexes)
        else:
            buf = recarray(n_records, dtype=dtype)
        if buf is not None:  # at least some channels should be parsed
            if dataRead_available:  # use rather cython compiled code for performance
                bytes_data = bytes(bit_stream)
                for chan in channels_indexes:
                    if self[chan].is_ca_block(info):
                        ca = self[chan].ca_block(info)
                        array_flag = ca['ca_ndim']
                    else:
                        array_flag = 0
                    buf[self[chan].name] = \
                        sorted_data_read(bytes_data, self[chan].bit_count(info),
                                         self[chan].signal_data_type(info),
                                         self[chan].native_data_format(info),
                                         n_records, self.CGrecordLength,
                                         self[chan].bit_offset(
                                             info), self[chan].pos_byte_beg(info),
                                         self[chan].calc_bytes(info, aligned=False), array_flag)
                return buf
            else:
                return self.read_channels_from_bytes_fallback(bit_stream, info, channel_set, n_records, dtype)
        else:
            return []

    def read_channels_from_bytes_fallback(self, bit_stream, info, channel_set=None, n_records=None, dtype=None,
                                          channels_indexes=None):
        """ reads stream of record bytes using bitarray in case no dataRead available

        Parameters
        ------------
        bit_stream : stream
            stream of bytes
        info: info class
        channel_set : set of str, optional
            set of channel to read
        n_records: int
            number of records
        dtype: numpy dtype
        channels_indexes: list of int

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """

        def signed_int(temp, extension):
            """ extend bits of signed data managing two's complement
            """
            extension.setall(False)
            extension_inv = bitarray(extension, endian='little')
            extension_inv.setall(True)
            for i in range(n_records):  # extend data of bytes to match numpy requirement
                sign_bit = temp[i][-1]
                if not sign_bit:  # positive value, extend with 0
                    temp[i].extend(extension)
                else:  # negative value, extend with 1
                    sign_bit = temp[i].pop(-1)
                    temp[i].extend(extension_inv)
                    temp[i].append(sign_bit)
            return temp
        if n_records is None:
            n_records = self.numberOfRecords
        if dtype is None:
            buf, channels_indexes = self.initialise_recarray(
                info, channel_set, n_records, dtype, channels_indexes)
        else:
            buf = recarray(n_records, dtype=dtype)
        if buf is not None:
            # read data
            from bitarray import bitarray
            bit_array = bitarray(endian="little")  # little endian by default
            bit_array.frombytes(bytes(bit_stream))
            record_bit_size = self.CGrecordLength * 8
            for chan in channels_indexes:
                signal_data_type = self[chan].signal_data_type(info)
                n_bytes_estimated = self[chan].nBytes_aligned
                if not self[chan].type in (1, 2):
                    temp = [bit_array[self[chan].pos_bit_beg + record_bit_size * i:
                                      self[chan].pos_bit_end(info) + record_bit_size * i]
                            for i in range(n_records)]
                    n_bytes = len(temp[0].tobytes())
                    if not n_bytes == n_bytes_estimated and \
                            signal_data_type not in (6, 7, 8, 9, 10, 11, 12):  # not Ctype byte length
                        byte = bitarray(
                            8 * (n_bytes_estimated - n_bytes), endian='little')
                        byte.setall(False)
                        if signal_data_type not in (2, 3):  # not signed integer
                            # extend data of bytes to match numpy requirement
                            for i in range(n_records):
                                temp[i].extend(byte)
                        # signed integer (two's complement), keep sign bit and extend with bytes
                        else:
                            temp = signed_int(temp, byte)
                    n_trail_bits = n_bytes_estimated * \
                        8 - self[chan].bit_count(info)
                    if signal_data_type in (2, 3) and \
                            n_bytes == n_bytes_estimated and \
                            n_trail_bits > 0:  # C type byte length but signed integer
                        trail_bits = bitarray(n_trail_bits, endian='little')
                        temp = signed_int(temp, trail_bits)
                else:  # Channel Array
                    temp = [bit_array[self[chan].pos_bit_beg + record_bit_size * i:
                                      self[chan].pos_bit_beg + 8 * n_bytes_estimated + record_bit_size * i]
                            for i in range(n_records)]
                is_ca = self[chan].type in (1, 2)
                if is_ca:
                    # For CA array channels: read raw bytes and interpret via element dtype
                    raw = b''.join(t.tobytes() for t in temp)
                    cn_info = info['CN'][self[chan].dataGroup][self[chan].channelGroup][self[chan].channelNumber]
                    signal_data_type = self[chan].signal_data_type(info)
                    n_bits = cn_info['cn_bit_count'] + cn_info['cn_bit_offset']
                    n_elem_bytes = (n_bits + 7) // 8
                    from mdfreader.channel import data_type_format4 as _dtf4
                    endian, elem_fmt = _dtf4(signal_data_type, n_elem_bytes)
                    elem_dtype = endian + elem_fmt
                    # collect array dimensions from nested CA blocks
                    ca_block = cn_info.get('CABlock')
                    dims = []
                    block = ca_block
                    while block is not None:
                        d = block['ca_dim_size']
                        if hasattr(d, '__iter__'):
                            dims.extend(int(x) for x in d)
                        else:
                            dims.append(int(d))
                        block = block.get('CABlock')
                    try:
                        temp = frombuffer(raw, dtype=elem_dtype).reshape((n_records, *dims))
                    except (ValueError, TypeError):
                        temp = frombuffer(raw, dtype=elem_dtype)
                elif 's' not in self[chan].c_format(info):
                    sdt = self[chan].signal_data_type(info)
                    if sdt in (15, 16):  # complex: use frombuffer with numpy dtype
                        raw = b''.join(t.tobytes() for t in temp)
                        fmt = self[chan].native_data_format(info)
                        endian = '<' if sdt == 15 else '>'
                        temp = frombuffer(raw, dtype=endian + fmt)
                    else:
                        c_structure = self[chan].c_format_structure(info)
                        if ('>' in self[chan].data_format(info) and byteorder == 'little') or \
                           (byteorder == 'big' and '<' in self[chan].data_format(info)):
                            temp = [c_structure.unpack(temp[i].tobytes())[0] for i in range(n_records)]
                            temp = asarray(temp)
                            temp = temp.byteswap().view(temp.dtype.newbyteorder())
                        else:
                            temp = [c_structure.unpack(temp[i].tobytes())[0] for i in range(n_records)]
                            temp = asarray(temp)
                else:
                    temp = [temp[i].tobytes()
                            for i in range(n_records)]
                    temp = asarray(temp)
                buf[self[chan].name] = temp
            return buf
        else:
            return []


class Mdf4(MdfSkeleton):

    """ mdf file reader class from version 4.0 to 4.1.1

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
    convertAfterRead : bool
        flag to convert raw data to physical just after read
    filterChannelNames : bool
        flag to filter long channel names from its module names separated by '.'
    fileMetadata : dict
        file metadata with minimum keys : author, organisation, project, subject, comment, time, date

    Methods
    ------------
    read4( fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True)
        Reads mdf 4.x file data and stores it in dict
    _get_channel_data_4(channelName)
        Returns channel numpy array
    _convert_channel_data_4(channel, channel_name, convert_tables, multiProc=False, Q=None)
        select right conversion and calculates it
    _convert_channel_4(channelName)
        converts specific channel from raw to physical data according to CCBlock information
    _convert_all_channel_4()
        Converts all channels from raw data to converted data according to CCBlock information
    write4(file_name=None, compression=False)
        writes mdf 4.1 file
    apply_invalid_bit(channel_name)
        mask data from invalid bit channel if existing
    get_channel_name_4(name, path)
        returns a list of tuples
    """

    def read4(self, file_name=None, info=None, multi_processed=False, channel_list=None, convert_after_read=True,
              filter_channel_names=False, compression=False, metadata=2, finalization_writing_to_file=False,
              force_file_integrity_check=False):
        """ Reads mdf 4.x file data and stores it in dict

        Parameters
        ----------------
        file_name : str, optional, None by default
            file name

        info : mdfinfo4.info4 class, optional, None by default
            info4 class containing all MDF Blocks

        multi_processed : bool, False by default
            flag to activate multiprocessing of channel data conversion

        channel_list : list of str, optional, None by default
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory.
            Can be used to read big files

        convert_after_read : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false,
            all data from channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .get_channel_data()

        filter_channel_names : bool, optional, False by default
            flag to filter long channel names from its module names separated by '.'

        compression : bool, optional
            flag to activate data compression with blosc

        metadata: int, optional, default = 2
            Reading metadata has impact on performance, especially for mdf 4.x using xml.
            2: minimal metadata reading (mostly channel blocks)
            1: used for noDataLoading
            0: all metadata reading, including Source Information, Attachment, etc..

        finalization_writing_to_file : bool, optional, False by default
            If file is detected not finalised (id_unfin_flags!=0), file is corrected
            writing in the blocks if set to True, otherwise correction is done only
            to file representation in memory.

        force_file_integrity_check : bool, optional, False by default
            Perform block sizes check for potentially corrupted file without finalization
            flags (id_unfin_flags==0). Combined with finalization_writing_to_file is
            very experimental and risky, correction should be tried in memory first.

        """

        self.multiProc = multi_processed

        if self.fileName is None and info is not None:
            self.fileName = info.fileName
        elif file_name is not None and self.fileName is None:
            self.fileName = file_name

        minimal = metadata  # always read minimum info (2), full info (0)

        # set is more efficient for large number of channels (n^2 vs n*log(n)):
        if channel_list is not None:
            channel_set_file = set(channel_list)  # make sure it is a set
            minimal = 1  # reads at least CN to populate ChannelNamesByDG
        else:
            channel_set_file = None

        # Read information block from file
        if info is None:
            if self.info is None:
                info = Info4(self.fileName, None,
                             filter_channel_names=filter_channel_names, minimal=minimal)
            else:
                info = self.info

        # check for eventual unfinalised file
        info = file_finalization(self.MDFVersionNumber, info, self.fid,
                                 finalization_writing_to_file, force_file_integrity_check)

        if info.fid is None or info.fid.closed:
            info.fid = open(self.fileName, 'rb')

        # keep events
        try:
            self.events = info['EV']
        except KeyError:
            pass

        # keep attachments
        try:
            self.attachments = info['AT']
        except KeyError:
            pass

        # reads metadata
        if not self._noDataLoading:
            ttime = info['HD']['hd_start_time_ns'] / 1E9
            if 1 << 1 & info['HD']['hd_time_flags']:
                # timezone and daylight applicable
                ttime = info['HD']['hd_tz_offset_min'] * \
                    60 + info['HD']['hd_dst_offset_min'] * 60

            def returnField(obj, field):
                try:
                    return obj[field]
                except KeyError:
                    return ''
            if 'Comment' in info['HD']:
                Comment = info['HD']['Comment']
                author = returnField(Comment, 'author')
                organisation = returnField(Comment, 'department')
                project = returnField(Comment, 'project')
                subject = returnField(Comment, 'subject')
                comment = returnField(Comment, 'TX')
                self.add_metadata(author=author, organisation=organisation,
                                  project=project, subject=subject,
                                  comment=comment, time=ttime)
            else:
                self.add_metadata(time=ttime)

        data_groups = info['DG']  # parse all data groups
        if self._noDataLoading and channel_list is not None:
            data_groups = [self[channel][idField][0][0]
                           for channel in channel_list]

        for dataGroup in data_groups:
            channel_set = channel_set_file
            if not info['DG'][dataGroup]['dg_data'] == 0 and \
                    (channel_set is None or
                     len(channel_set & info['ChannelNamesByDG'][dataGroup]) > 0):  # there is data block and channel in
                if minimal > 1 and not self._noDataLoading:  # load CG, CN and CC block info
                    info.read_cg_blocks(info.fid, dataGroup,
                                        channel_set, minimal=minimal)
                data_existing_in_data_group = False
                for dg in info['CG'][dataGroup]:
                    if info['CG'][dataGroup][dg]['cg_cycle_count']:
                        data_existing_in_data_group = True  # data existing
                        break
                if data_existing_in_data_group:
                    # Pointer to data block
                    pointer_to_data = info['DG'][dataGroup]['dg_data']

                    if 'dataClass' not in info['DG'][dataGroup]:
                        buf = Data(info.fid, pointer_to_data)
                        for channelGroup in info['CG'][dataGroup]:
                            # create record class
                            temp = Record(dataGroup, channelGroup)
                            # load all info related to record
                            temp.load_info(info)
                            buf.add_record(temp)  # adds record to DATA
                            record_id = info['CG'][dataGroup][channelGroup]['cg_record_id']
                            if temp.master is not None \
                                    and buf[record_id]['record'].channelNames:
                                if channel_set is not None and not self._noDataLoading\
                                        and temp.master not in channel_set:
                                    # adds master channel in channelSet if missing
                                    channel_set.add(temp.master)
                            if channel_set is not None and buf[record_id]['record'].CANOpen:
                                # adds CANOpen channels if existing in not empty channelSet
                                if buf[record_id]['record'].CANOpen == 'time':
                                    channel_set.update(('ms', 'days'))
                                elif buf[record_id]['record'].CANOpen == 'date':
                                    channel_set.update(
                                        ('ms', 'minute', 'hour', 'day', 'month', 'year'))
                        if self._noDataLoading:
                            self.info['DG'][dataGroup]['dataClass'] = buf
                    else:
                        buf = self.info['DG'][dataGroup]['dataClass']

                    # reads raw data from data block with DATA and _data_block classes
                    buf.read(channel_set, info, self.fileName)

                    channel_groups = buf
                    if self._noDataLoading and channel_list is not None:
                        channel_groups = [info['CG'][dataGroup][self[channel][idField][0][1]]['cg_record_id']
                                          for channel in channel_list]

                    # processing data from buf then transfer to self
                    for record_id in channel_groups:  # for each channel group in data block
                        if 'record' in buf[record_id]:
                            master_channel = buf[record_id]['record'].master

                            if self._noDataLoading and channel_list is not None:
                                channels = [
                                    buf[record_id]['record'][self[channel][idField][0][2]] for channel in channel_list]
                            else:
                                channels = list(
                                    buf[record_id]['record'].values())
                            for chan in channels:  # for each channel class
                                if channel_set is None or chan.name in channel_set:
                                    if not chan.type == 4:  # normal channel
                                        # not virtual channel
                                        if chan.channel_type(info) not in (3, 6):
                                            # in case record is used for several channels
                                            if channel_set is None and not buf[record_id]['record'].hiddenBytes \
                                                    and buf[record_id]['record'].byte_aligned:
                                                record_name = buf[record_id]['record'].recordToChannelMatching[chan.name]
                                            else:
                                                record_name = chan.name
                                            try:  # data in channel group
                                                # extract channel vector
                                                temp = buf[record_id]['data'][record_name]
                                            # no sorted data but maybe VLSD data
                                            except (ValueError, IndexError, KeyError):
                                                try:
                                                    temp = buf[record_id]['VLSD'][record_name]
                                                except (KeyError, TypeError):
                                                    temp = None
                                            except Exception:
                                                temp = None
                                        else:  # virtual channel
                                            temp = arange(
                                                buf[record_id]['record'].numberOfRecords)

                                        # Process concatenated bits inside uint8
                                        bit_count = chan.bit_count(info)
                                        if buf[record_id]['record'].byte_aligned \
                                                and not buf[record_id]['record'].hiddenBytes and \
                                                channel_set is None and\
                                                0 < bit_count < 64 and bit_count not in (8, 16, 32) \
                                                and temp is not None\
                                                and temp.dtype.kind not in ('S', 'U'):
                                            # if channel data do not use complete bytes and Ctypes
                                            signal_data_type = chan.signal_data_type(
                                                info)
                                            # integers
                                            if signal_data_type in (0, 1, 2, 3):
                                                bit_offset = chan.bit_offset(
                                                    info)
                                                if bit_offset > 0:
                                                    temp = right_shift(
                                                        temp, bit_offset)
                                                # masks isBitUnit8
                                                mask = int(
                                                    (1 << bit_count) - 1)
                                                temp = bitwise_and(temp, mask)
                                                if signal_data_type in (2, 3):
                                                    # signed integer, moving bit sign of two's complement
                                                    sign_bit_mask = (
                                                        1 << (bit_count - 1))
                                                    sign_extend = (
                                                        (1 << (temp.itemsize * 8 - bit_count)) - 1) << bit_count
                                                    sign_bit = bitwise_and(
                                                        temp, sign_bit_mask)
                                                    temp_u = temp.view('u{}'.format(temp.itemsize))
                                                    for number, sign in enumerate(sign_bit):
                                                        # negative value, sign extend
                                                        if sign:
                                                            temp_u[number] |= sign_extend
                                                    temp = temp_u.view(temp.dtype)
                                            else:  # should not happen
                                                warn('bit count and offset not applied to correct '
                                                     'data type {}'.format(chan.name))

                                        if temp is not None:  # channel contains data
                                            # half-precision complex: (2,)f2 → complex64
                                            if temp.ndim > 1 and temp.shape[-1] == 2 \
                                                    and temp.dtype.kind == 'f' \
                                                    and temp.dtype.itemsize == 2:
                                                signal_data_type = chan.signal_data_type(info)
                                                if signal_data_type == 16:  # BE: byteswap each f2
                                                    temp = temp.byteswap().view(
                                                        temp.dtype.newbyteorder())
                                                temp = (temp[..., 0].astype('f4')
                                                        + 1j * temp[..., 1].astype('f4'))
                                            # string data decoding
                                            if temp.dtype.kind == 'S':
                                                signal_data_type = chan.signal_data_type(
                                                    info)
                                                if signal_data_type == 6:  # string ISO-8859-1 Latin
                                                    encoding = 'latin-1'
                                                elif signal_data_type == 7:  # UTF-8
                                                    encoding = 'UTF-8'
                                                elif signal_data_type == 8:
                                                    encoding = 'UTF-16LE'
                                                elif signal_data_type == 9:  # UTF-16 big endian
                                                    encoding = 'UTF-16BE'
                                                else:
                                                    encoding = None
                                                if encoding is not None:
                                                    temp2 = empty(
                                                        len(temp), dtype='U{}'.format(temp.dtype.str[-1]))
                                                    for t in range(temp.size):
                                                        try:
                                                            temp2[t] = temp[t].decode(
                                                                encoding, 'ignore')
                                                        except Exception:
                                                            warn(
                                                                'Cannot decode channel {}'.format(chan.name))
                                                            temp2[t] = ''
                                                    temp = temp2

                                            # channel creation
                                            self.add_channel(chan.name, temp, master_channel,
                                                             master_type=chan.channel_sync_type(
                                                                 info),
                                                             unit=chan.unit(info), description=chan.desc(info),
                                                             conversion=chan.conversion(info), info=chan.cn_block(info),
                                                             compression=compression,
                                                             identifier=info.unique_id(chan.dataGroup,
                                                                                       chan.channelGroup,
                                                                                       chan.channelNumber))
                                            # sync channel
                                            if chan.channel_type(info) == 4:
                                                # attach stream to be synchronised
                                                self.set_channel_attachment(
                                                    chan.name, chan.attachment(info.fid, info))
                                            if chan.has_invalid_bit(info) and \
                                                    not info['DG'][dataGroup]['unique_channel_in_DG']:
                                                # has invalid bit
                                                self.set_invalid_bit(
                                                    chan.name, chan.invalid_bit(info))
                                                self.set_invalid_channel(
                                                    chan.name, 'invalid_bytes{}'.format(dataGroup))
                                    else:  # invalid bytes channel
                                        if buf[record_id]['invalid_data'] is None:
                                            invalid_data = buf[record_id]['data'].__getattribute__(
                                                chan.name)
                                        else:
                                            invalid_data = buf[record_id]['invalid_data']
                                        if not info['DG'][dataGroup]['unique_channel_in_DG']:
                                            invalid_data = frombuffer(invalid_data.tobytes(),
                                                                      dtype='u1').reshape(len(invalid_data),
                                                                                          invalid_data.dtype.itemsize)
                                            self.add_channel(chan.name, invalid_data, master_channel,
                                                             master_type=0, unit='', description='', info=None,
                                                             compression=compression, identifier=None)
                                        else:
                                            # unique channel in DG, applying easily maskarray
                                            data = self._get_channel_data4(
                                                channels[0].name)
                                            data = data.view(MaskedArray)
                                            data.mask = invalid_data
                                            self.set_channel_data(
                                                channels[0].name, data)
                            buf[record_id].pop('data', None)
                    del buf
                if minimal > 1:
                    # clean CN, CC and CG info to free memory
                    info.clean_dg_info(dataGroup)
        info.fid.close()  # close file

        if convert_after_read and not compression:
            self._noDataLoading = False
            self._convert_all_channel4()
        # print( 'Finished in ' + str( time.clock() - inttime ) , file=stderr)

    def _get_channel_data4(self, channel_name, raw_data=False):
        """Returns channel numpy array

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
            vector = self.get_channel(channel_name)[dataField]
            if vector is None:  # noDataLoading reading argument flag activated
                if self.info.fid is None or (self.info.fid is not None and self.info.fid.closed):
                    (self.info.fid, self.info.fileName,
                     self.info.zipfile) = _open_mdf(self.fileName)
                self.read4(file_name=None, info=None, channel_list=[
                           channel_name], convert_after_read=False)
            if not raw_data:
                return self._convert_channel_data4(self.get_channel(channel_name), channel_name,
                                                   self.convertTables)[channel_name]
            else:
                return self.get_channel(channel_name)[dataField]
        else:
            return None

    @staticmethod
    def _convert_channel_data4(channel, channel_name, convert_tables, multi_processed=False, q=None):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channel : channel class
            channel4 object
        channel_name : str
            name of channel
        convert_tables : bool
            activates computation intensive loops for conversion with tables. Default is False
        multi_processed : bool, default False
            flag to put data in multiprocess queue
        q : Queue class, default None
            Queue used for multiprocessing

        Returns
        -----------
        dict
            returns dict with channelName key containing numpy array converted to physical values according to
             conversion type
        """
        if channel[dataField] is None:
            vector = channel[dataField]
        else:
            if isinstance(channel[dataField], CompressedData):
                # uncompressed blosc data
                vector = channel[dataField].decompression()
            else:
                # to have bcolz uncompressed data
                vector = channel[dataField][:]
        # there is conversion property
        if conversionField in channel and channel[conversionField]['type']:
            text_type = vector.dtype.kind in [
                'S', 'U', 'V']  # channel of string or not ?
            conversion_type = channel[conversionField]['type']
            conversion_parameter = channel[conversionField]['parameters']
            if conversion_type == 1 and not text_type:
                vector = _linear_conversion(
                    vector, conversion_parameter['cc_val'])
            elif conversion_type == 2 and not text_type:
                vector = _rational_conversion(
                    vector, conversion_parameter['cc_val'])
            elif conversion_type == 3 and not text_type:
                vector = _formula_conversion(
                    vector, conversion_parameter['cc_ref']['Comment'])
            elif conversion_type == 4 and not text_type:
                vector = _value_to_value_table_with_interpolation_conversion(
                    vector, conversion_parameter['cc_val'])
            elif conversion_type == 5 and not text_type:
                vector = _value_to_value_table_without_interpolation_conversion(
                    vector, conversion_parameter['cc_val'])
            elif conversion_type == 6 and not text_type and convert_tables:
                vector = _value_range_to_value_table_conversion(
                    vector, conversion_parameter['cc_val'])
            elif conversion_type == 7 and not text_type and convert_tables:
                vector = _value_to_text_conversion(vector, conversion_parameter['cc_val'],
                                                   conversion_parameter['cc_ref'])
            elif conversion_type == 8 and not text_type and convert_tables:
                vector = _value_range_to_text_conversion(vector, conversion_parameter['cc_val'],
                                                         conversion_parameter['cc_ref'])
            elif conversion_type == 9 and text_type and convert_tables:
                vector = _text_to_value_conversion(vector, conversion_parameter['cc_val'],
                                                   conversion_parameter['cc_ref'])
            elif conversion_type == 10 and text_type and convert_tables:
                vector = _text_to_text_conversion(
                    vector, conversion_parameter['cc_ref'])
            elif conversion_type == 11 and text_type and convert_tables:
                vector = _bitfield_text_table_conversion(vector, conversion_parameter['cc_val'],
                                                         conversion_parameter['cc_ref'])
        L = dict()
        L[channel_name] = vector
        if multi_processed:
            q.put(L)
        else:
            return L

    def _convert_channel4(self, channel_name):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channel_name : str
            Name of channel
        """
        self.set_channel_data(
            channel_name, self._get_channel_data4(channel_name))
        self.remove_channel_conversion(channel_name)

    def _convert_all_channel4(self):
        """Converts all channels from raw data to converted data according to CCBlock information
        Converted data will take more memory.
        """

        if self._noDataLoading:  # no data loaded, load everything
            self.read4(self.fileName, convert_after_read=True)
        else:
            if self.multiProc is False:
                [self._convert_channel4(channelName) for channelName in self]
            else:  # multiprocessing
                proc = []
                Q = Queue()
                L = {}
                for channelName in self:
                    channel = self.get_channel(channelName)
                    if 'conversion' in channel:
                        conversion = self.get_channel_conversion(channelName)
                        # more time in multi proc
                        if conversion['type'] in (1, 2):
                            self._convert_channel4(channelName)
                        else:
                            proc.append(Process(target=self._convert_channel_data4,
                                                args=(channel, channelName, self.convertTables, True, Q)))
                            proc[-1].start()
                for p in proc:
                    # concatenate results of processes in dict
                    L.update(Q.get())
                for p in proc:
                    p.join()
                del Q  # free memory
                for channelName in self:
                    if channelName in L:
                        self.set_channel_data(channelName, L[channelName])
                        self.remove_channel_conversion(channelName)

    def write4(self, file_name=None, compression=False, column_oriented=False):
        """Writes simple mdf file

        Parameters
        ----------------
        file_name : str, optional
            Name of file
            If file name is not input, written file name will be the one read with appended '_new' string before extension
        compression : bool
            flag to store data compressed
        column_oriented : bool
            flag to store data in columns, faster reading channel by channel and not jumping in records

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        # Starts first to write ID and header
        if file_name is None:
            split_name = splitext(self.fileName)
            if split_name[-1] in ('.mfxz', '.MFXZ'):
                split_name[-1] = '.mfx'  # do not resave in compressed file
            file_name = ''.join([split_name[-2], '_New', split_name[-1]])
        fid = open(file_name, 'wb')  # buffering should automatically be set
        # IDBLock writing
        temp = IDBlock()
        if column_oriented:
            temp['id_vers'] = b'4.20    '
            temp['id_ver'] = 420
        else:
            temp['id_vers'] = b'4.11    '
            temp['id_ver'] = 411
        temp.write(fid)

        blocks = OrderedDict()
        pointer = 64
        # Header Block
        blocks['HD'] = HDBlock()
        blocks['HD']['time'] = self.fileMetadata['time']
        blocks['HD']['block_start'] = pointer
        pointer += 104

        # Header Block comments
        blocks['HD']['MD'] = pointer
        blocks['HD_comment'] = CommentBlock()
        blocks['HD_comment']['block_start'] = pointer
        blocks['HD_comment'].load(self.fileMetadata, 'HD')
        pointer = blocks['HD_comment']['block_start'] + \
            blocks['HD_comment']['block_length']

        # file history block
        blocks['FH'] = FHBlock()
        blocks['HD']['FH'] = pointer
        blocks['FH']['block_start'] = pointer
        pointer = blocks['FH']['block_start'] + 56

        # File History comment
        blocks['FH']['MD'] = pointer
        blocks['FH_comment'] = CommentBlock()
        blocks['FH_comment']['block_start'] = pointer
        blocks['FH_comment'].load(self.fileMetadata, 'FH')
        pointer = blocks['FH_comment']['block_start'] + \
            blocks['FH_comment']['block_length']

        # write DG block
        blocks['HD']['DG'] = pointer  # first DG

        # write all files header blocks
        for block in blocks.values():
            block.write(fid)
        if self.masterChannelList:  # some channels exist
            if column_oriented:
                for dataGroup, masterChannel in enumerate(self.masterChannelList):
                    # write master channel
                    master_channel_data = self._get_channel_data4(
                        masterChannel)
                    if master_channel_data is not None:
                        cg_cycle_count = len(master_channel_data)
                    elif self._get_channel_data4(self.masterChannelList[masterChannel][0]).shape[0] == 1:
                        # classification
                        cg_cycle_count = 1
                    elif master_channel_data not in self.masterChannelList[masterChannel]:
                        cg_cycle_count = len(self._get_channel_data4(
                            self.masterChannelList[masterChannel][0]))
                        warn('no master channel in data group {}'.format(dataGroup))
                    else:
                        # no data in data group, skip
                        warn('no data in data group {0} with master channel {1}'.format(
                            dataGroup, masterChannel))
                        continue

                    cg_master_pointer = None
                    for channel in self.masterChannelList[masterChannel]:
                        # not invalid bytes channel
                        if channel.find('invalid_bytes') == -1:
                            # skip zero-size channels (e.g. DS composition containers with V0 data)
                            ch_data = self._get_channel_data4(channel)
                            if ch_data is not None and hasattr(ch_data, 'dtype') and ch_data.dtype.itemsize == 0:
                                continue
                            # write other channels
                            if channel == masterChannel:
                                cg_master_pointer = pointer + 64  # CG after DG at pointer
                                dg, pointer = self._write4_column(fid, pointer, cg_cycle_count, masterChannel,
                                                                  True, compression)
                            else:
                                if cg_master_pointer is not None:
                                    dg, pointer = self._write4_column(fid, pointer, cg_cycle_count, channel,
                                                                      False, compression,
                                                                      cg_master_pointer=cg_master_pointer)
                                else:  # no master
                                    dg, pointer = self._write4_column(fid, pointer, cg_cycle_count, channel,
                                                                      False, compression)

                fid.seek(dg['block_start'] + 24)
                fid.write(pack('Q', 0))  # last DG pointer is null
            else:
                self._write4_non_column(fid, pointer, compression)
        fid.close()

    def _write4_non_column(self, fid, pointer, compression=False):
        """Writes simple mdf 4.1 file with sorted data

        Parameters
        ----------------
        fid
            file identifier
        compression : bool
            flag to store data compressed

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        for dataGroup, masterChannel in enumerate(self.masterChannelList):
            # writes dataGroup Block
            dg = DGBlock()
            dg['block_start'] = pointer
            pointer = dg['block_start'] + 64
            dg['CG'] = pointer  # First CG link

            blocks = OrderedDict()  # initialise blocks for this datagroup
            # write CGBlock
            blocks['CG'] = CGBlock()
            blocks['CG']['master_channel_flag'] = False
            blocks['CG']['length'] = 104
            blocks['CG']['block_start'] = pointer
            pointer = blocks['CG']['block_start'] + blocks['CG']['length']
            blocks['CG']['CN'] = pointer  # First CN link
            blocks['CG']['cg_inval_bytes'] = 0  # no invalidation bytes

            master_channel_data = self._get_channel_data4(masterChannel)
            if master_channel_data is not None:
                cg_cycle_count = len(master_channel_data)
            # classification
            elif self._get_channel_data4(self.masterChannelList[masterChannel][0]).shape[0] == 1:
                cg_cycle_count = 1
            elif master_channel_data not in self.masterChannelList[masterChannel]:
                cg_cycle_count = len(self._get_channel_data4(
                    self.masterChannelList[masterChannel][0]))
                warn('no master channel in data group {}'.format(dataGroup))
            else:
                # no data in data group, skip
                warn('no data in data group {0} with master channel {1}'.format(
                    dataGroup, masterChannel))
                continue
            blocks['CG']['cg_cycle_count'] = cg_cycle_count

            # write channels
            record_byte_offset = 0
            cn_flag = 0
            n_records = 0
            data_list = ()
            last_channel = 0
            previous_n_channel = 0
            for n_channel, channel in enumerate(self.masterChannelList[masterChannel]):
                data = self.get_channel_data(channel)
                # no interest to write invalid bytes as channel, should be processed if needed before writing
                if channel.find('invalid_bytes') == -1 and data is not None and len(data) > 0 \
                        and data.dtype.itemsize > 0:
                    byte_count = data.dtype.itemsize
                    blocks[n_channel] = CNBlock()
                    blocks[n_channel]['cn_byte_offset'] = record_byte_offset

                    last_channel = n_channel
                    data_ndim = data.ndim - 1
                    if not data_ndim:
                        data_list = data_list + (data, )
                        record_byte_offset += byte_count
                    else:  # data contains arrays
                        data_dim_size = data.shape
                        if not cg_cycle_count == data_dim_size[0]:
                            warn(
                                'Array length do not match number of cycled in CG block')
                        data_dim_size = data_dim_size[1:]
                        SNd = 0
                        PNd = 1
                        for x in data_dim_size:
                            SNd += x
                            PNd *= x
                        flattened = reshape(data, (cg_cycle_count, PNd))
                        for i in range(PNd):
                            data_list = data_list + (flattened[:, i],)
                        record_byte_offset += byte_count * PNd

                    if issubdtype(data.dtype, numpy_number):  # is numeric
                        blocks[n_channel]['cn_val_range_min'] = npmin(data)
                        blocks[n_channel]['cn_val_range_max'] = npmax(data)
                        # only Bit 3: Limit range valid flag
                        blocks[n_channel]['cn_flags'] = 8
                    else:
                        blocks[n_channel]['cn_val_range_min'] = 0
                        blocks[n_channel]['cn_val_range_max'] = 0
                        blocks[n_channel]['cn_flags'] = 0
                    if masterChannel is not channel:
                        blocks[n_channel]['cn_type'] = 0
                        blocks[n_channel]['cn_sync_type'] = 0
                    else:
                        blocks[n_channel]['cn_type'] = 2  # master channel
                        blocks[n_channel]['cn_sync_type'] = self.get_channel_master_type(
                            channel)
                        n_records = len(data)

                    cn_numpy_kind = data.dtype.kind
                    if cn_numpy_kind in ('u', 'b'):
                        data_type = 0  # LE
                    elif cn_numpy_kind == 'i':
                        data_type = 2  # LE
                    elif cn_numpy_kind == 'f':
                        data_type = 4  # LE
                    elif cn_numpy_kind == 'c':
                        data_type = 15  # complex LE
                    elif cn_numpy_kind == 'S':
                        data_type = 6
                    elif cn_numpy_kind == 'U':
                        data_type = 7  # UTF-8
                    elif cn_numpy_kind == 'V':
                        data_type = 10  # bytes
                    else:
                        warn('{} {} {}'.format(channel, data.dtype, cn_numpy_kind))
                        raise Exception('Not recognized dtype')
                    blocks[n_channel]['cn_data_type'] = data_type
                    # always byte aligned
                    blocks[n_channel]['cn_bit_offset'] = 0

                    blocks[n_channel]['cn_bit_count'] = byte_count * 8
                    blocks[n_channel]['block_start'] = pointer
                    pointer = blocks[n_channel]['block_start'] + 160

                    # arrays handling
                    if not data_ndim:
                        blocks[n_channel]['Composition'] = 0
                    else:
                        # pointer to CABlock
                        blocks[n_channel]['Composition'] = pointer
                        # creates CABlock
                        ca = ''.join([channel, '_CA'])
                        blocks[ca] = CABlock()
                        blocks[ca]['block_start'] = pointer
                        blocks[ca]['ndim'] = data_ndim
                        blocks[ca]['ndim_size'] = data_dim_size
                        blocks[ca].load(data.itemsize)
                        pointer = blocks[ca]['block_start'] + \
                            blocks[ca]['block_length']

                    if cn_flag:
                        # Next DG
                        blocks[previous_n_channel]['CN'] = blocks[n_channel]['block_start']
                        blocks[n_channel]['CN'] = 0  # initialise 'CN' key
                    else:
                        cn_flag = blocks[n_channel]['block_start']
                        # creates first CN link, null for the moment
                        blocks[n_channel]['CN'] = 0
                    previous_n_channel = n_channel

                    # write channel name
                    blocks[n_channel]['TX'] = pointer
                    blocks[channel] = CommentBlock()
                    blocks[channel]['block_start'] = pointer
                    blocks[channel].load(channel, 'TX')
                    pointer = blocks[channel]['block_start'] + \
                        blocks[channel]['block_length']

                    # write channel unit
                    unit = self.get_channel_unit(channel)
                    if unit is not None and len(unit) > 0:
                        blocks[n_channel]['Unit'] = pointer
                        unit_name = u'{}{}{}'.format(channel, '_U_', n_channel)
                        blocks[unit_name] = CommentBlock()
                        blocks[unit_name]['block_start'] = pointer
                        blocks[unit_name].load(unit, 'TX')
                        pointer = blocks[unit_name]['block_start'] + \
                            blocks[unit_name]['block_length']
                    else:
                        blocks[n_channel]['Unit'] = 0

                    # write channel description
                    desc = self.get_channel_desc(channel)
                    if desc is not None and len(desc) > 0:
                        blocks[n_channel]['Comment'] = pointer
                        desc_name = '{}{}{}'.format(channel, '_C_', n_channel)
                        blocks[desc_name] = CommentBlock()
                        blocks[desc_name]['block_start'] = pointer
                        blocks[desc_name].load(desc, 'TX')
                        pointer = blocks[desc_name]['block_start'] + \
                            blocks[desc_name]['block_length']
                    else:
                        blocks[n_channel]['Comment'] = 0

            if n_records == 0 and masterChannel is not self.masterChannelList[masterChannel]:
                # No master channel in channel group
                if data_list:
                    n_records = len(data_list[0])

            if last_channel in blocks:
                blocks[last_channel]['CN'] = 0  # last CN link is null
            # writes size of record in CG
            blocks['CG']['cg_data_bytes'] = record_byte_offset

            # data pointer in data group
            dg['data'] = pointer
            if compression:
                data = HLBlock()
                data.load(record_byte_offset, n_records, pointer)
                dg_start_position = fid.tell()
                dg['DG'] = 0
            else:
                data = DTBlock()
                data.load(record_byte_offset, n_records, pointer)
                dg['DG'] = data['end_position']

            dg.write(fid)  # write DG block

            # writes all blocks (CG, CN, TX for unit and description) before writing data block
            for block in blocks.values():
                block.write(fid)

            # data block writing
            if not data_list:
                pointer = data.write(fid, b'')
            else:
                pointer = data.write(fid, fromarrays(data_list).tobytes(order='F'))
            if compression:
                # next DG position is not predictable due to DZ Blocks unknown length
                fid.seek(dg_start_position + 24)
                fid.write(pack('<Q', pointer))
                fid.seek(pointer)

        fid.seek(dg['block_start'] + 24)
        fid.write(pack('Q', 0))  # last DG pointer is null

    def _write4_column(self, fid, pointer, cg_cycle_count, channel, master_channel_flag,
                       compression=False, cg_master_pointer=None):
        """Writes simple mdf 4.2 file with column oriented channels
        Mostly efficient for reading but size might be bigger than original file

        Parameters
        ----------------
        fid
            file identifier
        cg_cycle_count : int
            number of records for the channel
        channel : string
            channel name
        master_channel_flag : boolean
            flag set if it is a master channel
        compression : bool
            flag to store data compressed
        cg_master_pointer : int
            point to other CGBlock containing master channel

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """
        record_byte_offset = 0

        # writes dataGroup Block
        dg = DGBlock()
        dg['block_start'] = pointer
        pointer = dg['block_start'] + 64
        dg['CG'] = pointer  # First CG link

        blocks = OrderedDict()  # initialise blocks for this datagroup
        # write CGBlock
        blocks['CG'] = CGBlock()
        blocks['CG']['cg_cg_master'] = cg_master_pointer
        if cg_master_pointer is None:  # master channel writing
            blocks['CG']['length'] = 104
        else:
            blocks['CG']['length'] = 112
        blocks['CG']['block_start'] = pointer
        pointer = blocks['CG']['block_start'] + blocks['CG']['length']
        blocks['CG']['CN'] = pointer  # First CN link
        blocks['CG']['cg_cycle_count'] = cg_cycle_count

        # write channels
        data = self.get_channel_data(channel)
        if self.get_invalid_channel(channel) or isinstance(data, MaskedArray):
            invalid_channel = True
            # as column oriented, only one byte for DIBlock
            blocks['CG']['cg_inval_bytes'] = 1
        else:
            invalid_channel = False
            blocks['CG']['cg_inval_bytes'] = 0
        # no interest to write invalid bytes as channel, should be processed if needed before writing
        if data is not None and len(data) > 0 and data.dtype.itemsize > 0:
            byte_count = data.dtype.itemsize
            data_ndim = data.ndim - 1
            if data_ndim:  # data contains arrays
                data_dim_size = data.shape
                if not cg_cycle_count == data_dim_size[0]:
                    warn('Array length do not match number of cycled in CG block')
                data_dim_size = data_dim_size[1:]
                SNd = 0
                PNd = 1
                for x in data_dim_size:
                    SNd += x
                    PNd *= x
                data = reshape(data, (cg_cycle_count, PNd))
                record_byte_offset += PNd * byte_count
            else:
                record_byte_offset += byte_count

            blocks['CN'] = CNBlock()
            blocks['CN']['cn_flags'] = 0
            if issubdtype(data.dtype, numpy_number) and not invalid_channel:  # is numeric
                blocks['CN']['cn_val_range_min'] = npmin(data)
                blocks['CN']['cn_val_range_max'] = npmax(data)
                # only Bit 3: Limit range valid flag
                blocks['CN']['cn_flags'] |= (1 << 3)
            else:
                blocks['CN']['cn_val_range_min'] = 0
                blocks['CN']['cn_val_range_max'] = 0
            if invalid_channel:
                blocks['CN']['cn_flags'] |= (1 << 1)

            if not master_channel_flag:
                blocks['CN']['cn_type'] = 0
                blocks['CN']['cn_sync_type'] = 0
            else:
                blocks['CN']['cn_type'] = 2  # master channel
                blocks['CN']['cn_sync_type'] = self.get_channel_master_type(
                    channel)

            cn_numpy_kind = data.dtype.kind
            if cn_numpy_kind in ('u', 'b'):
                data_type = 0  # LE
            elif cn_numpy_kind == 'i':
                data_type = 2  # LE
            elif cn_numpy_kind == 'f':
                data_type = 4  # LE
            elif cn_numpy_kind == 'c':
                data_type = 15  # complex LE
            elif cn_numpy_kind == 'S':
                data_type = 6
            elif cn_numpy_kind == 'U':
                data_type = 7  # UTF-8
            elif cn_numpy_kind == 'V':
                data_type = 10  # bytes
            else:
                warn('{} {} {}'.format(channel, data.dtype, cn_numpy_kind))
                raise Exception('Not recognized dtype')
            blocks['CN']['cn_data_type'] = data_type
            blocks['CN']['cn_bit_offset'] = 0  # always byte aligned
            blocks['CN']['cn_byte_offset'] = 0  # only one channel
            blocks['CN']['cn_bit_count'] = byte_count * 8
            blocks['CN']['block_start'] = pointer
            pointer = blocks['CN']['block_start'] + 160
            # creates first CN link, null for the moment
            blocks['CN']['CN'] = 0

            blocks['CN']['Composition'] = 0
            # arrays handling
            if data_ndim:
                blocks['CN']['Composition'] = pointer  # pointer to CABlock
                # creates CABlock
                ca = ''.join([channel, '_CA'])
                blocks[ca] = CABlock()
                blocks[ca]['block_start'] = pointer
                blocks[ca]['ndim'] = data_ndim
                blocks[ca]['ndim_size'] = data_dim_size
                blocks[ca].load(data.itemsize)
                pointer = blocks[ca]['block_start'] + \
                    blocks[ca]['block_length']

            # write channel name
            blocks['CN']['TX'] = pointer
            blocks[channel] = CommentBlock()
            blocks[channel]['block_start'] = pointer
            blocks[channel].load(channel, 'TX')
            pointer = blocks[channel]['block_start'] + \
                blocks[channel]['block_length']

            # write channel unit
            unit = self.get_channel_unit(channel)
            if unit is not None and len(unit) > 0:
                blocks['CN']['Unit'] = pointer
                unit_name = u'{}{}{}'.format(channel, '_U_', 'CN')
                blocks[unit_name] = CommentBlock()
                blocks[unit_name]['block_start'] = pointer
                blocks[unit_name].load(unit, 'TX')
                pointer = blocks[unit_name]['block_start'] + \
                    blocks[unit_name]['block_length']
            else:
                blocks['CN']['Unit'] = 0

            # write channel description
            desc = self.get_channel_desc(channel)
            if desc is not None and len(desc) > 0:
                blocks['CN']['Comment'] = pointer
                desc_name = '{}{}{}'.format(channel, '_C_', 'CN')
                blocks[desc_name] = CommentBlock()
                blocks[desc_name]['block_start'] = pointer
                blocks[desc_name].load(desc, 'TX')
                pointer = blocks[desc_name]['block_start'] + \
                    blocks[desc_name]['block_length']
            else:
                blocks['CN']['Comment'] = 0

            # writes size of record in CG
            blocks['CG']['cg_data_bytes'] = record_byte_offset

            # data pointer in data group
            dg['data'] = pointer
            if compression or invalid_channel:
                data_blocks = LDBlock()
                data_blocks.load(record_byte_offset, cg_cycle_count, pointer,
                                 invalid_bytes=blocks['CG']['cg_inval_bytes'], column_oriented_flag=True)
                dg_start_position = fid.tell()
                dg['DG'] = 0
            else:
                data_blocks = DVBlock()
                data_blocks.load(record_byte_offset, cg_cycle_count, pointer)
                dg['DG'] = data_blocks['end_position']

            dg.write(fid)  # write DG block

            # writes all blocks (CG, CN, TX for unit and description) before writing data block
            for block in blocks.values():
                block.write(fid)

            # data block writing
            if not invalid_channel:
                inval_data = None
            else:  # must be 1 byte length, boolean
                if isinstance(data, MaskedArray):
                    inval_data = data.mask.tobytes()
                else:
                    inval_data = self.get_invalid_mask(
                        channel).astype(bool).tobytes()
            if not invalid_channel and not compression:
                pointer = data_blocks.write(fid, data.tobytes())
            else:
                pointer = data_blocks.write(fid, data.tobytes(),
                                            invalid_data=inval_data, compression_flag=compression)

            if compression or invalid_channel:
                # next DG position is not predictable due to DZ Blocks unknown length
                fid.seek(dg_start_position + 24)
                fid.write(pack('<Q', pointer))
                fid.seek(pointer)

        return dg, pointer

    def apply_invalid_bit(self, channel_name):
        """Mask data of channel based on its invalid bit definition if there is

        Parameters
        ----------------
        channel_name : str
            Name of channel
        """
        try:
            data = self._get_channel_data4(channel_name)
            data = data.view(MaskedArray)
            data.mask = self.get_invalid_mask(channel_name)
            self.set_channel_data(channel_name, data)
            self._remove_channel_field(channel_name, invalidPosField)
            self._remove_channel_field(channel_name, invalidChannel)
        except KeyError:
            pass
            # warn('no invalid data found for channel ')

    def get_invalid_mask(self, channel_name):
        invalid_channel = self._get_channel_data4(
            self.get_invalid_channel(channel_name))
        invalid_bit_pos = self.get_invalid_bit(channel_name)
        if isinstance(invalid_bit_pos, int):  # invalid bit existing
            invalid_byte = invalid_bit_pos >> 3
            return bitwise_and(invalid_channel[:, invalid_byte], invalid_bit_pos & 0x07)
        else:
            return None

    def apply_all_invalid_bit(self):
        """Mask data of all channels based on its invalid bit definition if there is

        """
        for master_channel in self.masterChannelList:
            group_channels = set(self.masterChannelList[master_channel])
            group_number = self[master_channel][idField][0][0]
            invalid_channel = 'invalid_bytes{}'.format(group_number)
            if invalid_channel in group_channels:
                # invalid bytes channel present in this data group
                for channel_name in self.masterChannelList[master_channel]:
                    self.apply_invalid_bit(channel_name)
                # remove invalid bytes channel, redundant
                self.masterChannelList[master_channel].remove(invalid_channel)
                self.pop(invalid_channel)

    def get_channel_name4(self, name, path):
        """finds mdf channel name from name and path

        Parameters
        ----------------
        name : str
            channel name
        path: str
            source path or name, or channel group name, source name or path

        Returns
        -----------
        list of tuples (channel_name, (ndg, ncg, ncn))

        """
        output = []
        for channel_name in self:
            try:
                (ndg, ncg, ncn), (cn, cs, cp), (gn,
                                                gs, gp) = self[channel_name]['id']
                if name == cn and path in (cs, cp) or path in (gn, gs, gp):
                    output.append((channel_name, (ndg, ncg, ncn)))
            except KeyError:  # most probably a invalid bit channel
                pass
        return output


def _linear_conversion(vector, cc_val):
    """ apply linear conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    p1 = cc_val[0]
    p2 = cc_val[1]
    if p2 == 1.0 and p1 in (0.0, -0.0):
        return vector  # keeps dtype probably more compact than float64
    else:
        return vector * p2 + p1


def _rational_conversion(vector, cc_val):
    """ apply rational conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    p1 = cc_val[0]
    p2 = cc_val[1]
    p3 = cc_val[2]
    p4 = cc_val[3]
    p5 = cc_val[4]
    p6 = cc_val[5]
    return (p1 * vector * vector + p2 * vector + p3) / (p4 * vector * vector + p5 * vector + p6)


def _formula_conversion(vector, formula):
    """ apply formula conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    formula : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    try:
        from sympy import lambdify, symbols
    except ImportError:
        warn('Please install sympy to convert channel ')
        return vector
    X = symbols('X')
    expr = lambdify(X, formula, modules='numpy', dummify=False)
    return expr(vector)


def _value_to_value_table_without_interpolation_conversion(vector, cc_val):
    """ apply value to value table without interpolation conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = 2 * int(len(cc_val) / 2)
    val_array = array([cc_val[i] for i in range(val_count)])
    int_val = val_array[::2]
    phys_val = val_array[1::2]
    if all(diff(int_val) > 0):
        try:
            from scipy import interpolate
            # MDF4 spec 6.17.7: below first key → value[0], above last key → value[n-1]
            f = interpolate.interp1d(
                int_val, phys_val, kind='nearest', bounds_error=False,
                fill_value=(phys_val[0], phys_val[-1]))
            return f(vector)
        except ImportError:
            warn('Please install scipy to convert channel')
    else:
        warn('X values for interpolation of channel are not increasing')


def _value_to_value_table_with_interpolation_conversion(vector, cc_val):
    """ apply value to value table with interpolation conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = 2 * int(len(cc_val) / 2)
    val_array = array([cc_val[i] for i in range(val_count)])
    int_val = val_array[::2]
    phys_val = val_array[1::2]
    if all(diff(int_val) > 0):
        return interp(vector, int_val, phys_val)  # with interpolation
    else:
        warn('X values for interpolation of channel are not increasing')


def _value_range_to_value_table_conversion(vector, cc_val):
    """ apply value range to value table conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = int(len(cc_val) / 3)
    # cc_val layout: [lo_0,hi_0,val_0, lo_1,hi_1,val_1, ..., lo_{n-1},hi_{n-1},val_{n-1}, default]
    # so cc_val has 3*n+1 elements; key_min[n] holds the default value
    key_min = [cc_val[i] for i in range(0, 3 * val_count + 1, 3)]
    key_max = [cc_val[i] for i in range(1, 3 * val_count + 1, 3)]
    value = [cc_val[i] for i in range(2, 3 * val_count + 1, 3)]
    default_val = key_min[val_count]
    # MDF4 spec 6.17.8: integer data uses [lo, hi] (inclusive upper);
    # float data uses [lo, hi) (exclusive upper).
    import numpy as np_inner
    is_integer = np_inner.issubdtype(vector.dtype, np_inner.integer)
    for l_index in range(len(vector)):
        result = default_val
        x = vector[l_index]
        for i in range(val_count):
            if is_integer:
                if key_min[i] <= x <= key_max[i]:
                    result = value[i]
                    break
            else:
                if key_min[i] <= x < key_max[i]:
                    result = value[i]
                    break
        vector[l_index] = result
    return vector


def _value_to_text_conversion(vector, cc_val, cc_ref):
    """ apply value to text conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    str_refs = [ref for ref in cc_ref if isinstance(ref, str)]
    maxlen = max(len(ref) for ref in str_refs) if str_refs else 256
    # initialize empty array with max length
    if len(vector) == 0 or vector.ndim > 1:
        return vector  # can't apply value-to-text to empty or multi-dim data
    temp = empty(len(vector), dtype='U{}'.format(maxlen))
    # checks for scaling
    try:
        from sympy import lambdify, symbols
    except ImportError:
        warn('Please install sympy to convert channel ')
        return vector
    X = symbols('X')
    for ref in range(len(cc_ref)):
        if isinstance(cc_ref[ref], CCBlock):
            if cc_ref[ref]['cc_type'] == 3:
                # formula to be applied
                cc_ref[ref] = lambdify(X, cc_ref[ref]['cc_ref']['Comment'],
                                       modules='numpy', dummify=False)
            elif cc_ref[ref]['cc_type'] == 1:  # linear conversion
                cc_ref[ref] = lambdify(X, '{0}* X + {1}'.format(cc_ref[ref]['cc_val'][1], cc_ref[ref]['cc_val'][0]),
                                       modules='numpy', dummify=False)
            else:
                warn('To implement missing conversion, please ask')
        elif not isinstance(cc_ref[ref], str):  # identity, non conversion
            cc_ref[ref] = lambdify(X, 'X', modules='numpy', dummify=False)
    # look up for first value in vector
    key_index = where(vector[0] == cc_val)[0]
    if not len(key_index) == 0:  # value corresponding in cc_val
        temp[0] = cc_ref[key_index[0]]
    else:  # default value
        if callable(cc_ref[-1]):
            temp[0] = cc_ref[-1](vector[0])
        else:
            temp[0] = cc_ref[-1]
    for lindex in range(1, len(vector)):
        # same value as before, no need to look further
        if vector[lindex] == vector[lindex - 1]:
            temp[lindex] = temp[lindex - 1]
        else:  # value changed from previous step
            key_index = where(vector[lindex] == cc_val)[0]
            if not len(key_index) == 0:  # found match
                temp[lindex] = cc_ref[key_index[0]]
            else:  # default
                if callable(cc_ref[-1]):
                    temp[lindex] = cc_ref[-1](vector[lindex])
                else:
                    temp[lindex] = cc_ref[-1]
    return asarray(temp)


def _value_range_to_text_conversion(vector, cc_val, cc_ref):
    """ apply value range to text conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = int(len(cc_val) / 2)
    key_min = [cc_val[i] for i in range(0, 2 * val_count, 2)]
    key_max = [cc_val[i] for i in range(1, 2 * val_count, 2)]
    # checks for scaling
    try:
        from sympy import lambdify, symbols
    except ImportError:
        warn('Please install sympy to convert channel ')
        return vector
    X = symbols('X')
    for ref in range(len(cc_ref)):
        if isinstance(cc_ref[ref], CCBlock):
            if cc_ref[ref]['cc_type'] == 3:
                # formula to be applied
                cc_ref[ref] = lambdify(X, cc_ref[ref]['cc_ref']['Comment'],
                                       modules='numpy', dummify=False)
            elif cc_ref[ref]['cc_type'] == 1:  # linear conversion
                cc_ref[ref] = lambdify(X, '{0} * X + {1}'.format(cc_val[1], cc_val[0]),
                                       modules='numpy', dummify=False)
            else:  # identity, no conversion
                cc_ref[ref] = lambdify(X, '1 * X',
                                       modules='numpy', dummify=False)
            # Otherwise a string
    # look up in range keys
    temp = []
    for value in vector:
        key_index = val_count  # default index if not found
        for i in range(val_count):
            if key_min[i] <= value <= key_max[i]:
                key_index = i
                break
        if callable(cc_ref[key_index]):
            # TXBlock string
            temp.append(cc_ref[key_index](value))
        else:  # scale to be applied
            temp.append(cc_ref[key_index])
    return asarray(temp)


def _text_to_value_conversion(vector, cc_val, cc_ref):
    """ apply text to value conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    ref_count = len(cc_ref)
    temp = []
    for l_index in range(len(vector)):
        key_index = ref_count  # default index if not found
        for i in range(ref_count):
            if vector[l_index] == cc_ref[i]:
                key_index = i
                break
        temp.append(cc_val[key_index])
    return asarray(temp)


def _text_to_text_conversion(vector, cc_ref):
    """ apply text to text conversion to data

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    ref_count = len(cc_ref) - 2
    for l_index in range(len(vector)):
        key_index = ref_count + 1  # default index if not found
        for i in range(0, ref_count, 2):
            if vector[l_index] == cc_ref[i]:
                key_index = i
                break
        vector[l_index] = cc_ref[key_index]
    return vector


def _bitfield_text_table_conversion(vector, cc_val, cc_ref):
    """ apply bitfield to raw data and convert to string

    Parameters
    ----------------
    vector : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to string
    """
    ref_count = len(cc_ref)
    temp = []
    for l_index in range(len(vector)):
        assembled_string = ''
        for i in range(ref_count):
            bitmask = bitwise_and(vector[l_index], cc_val[i])
            if cc_ref[i]:  # not NIL link
                if cc_ref[i]['cc_type'] == 7:
                    assembled_string += _value_to_text_conversion(bitmask,
                                                                  cc_ref[i]['cc_val'],
                                                                  cc_ref[i]['cc_ref'])
                if cc_ref[i]['cc_type'] == 8:
                    assembled_string += _value_range_to_text_conversion(bitmask,
                                                                        cc_ref[i]['cc_val'],
                                                                        cc_ref[i]['cc_ref'])
        temp.append(assembled_string)
    return asarray(temp)


def file_finalization(version, info, fid, finalization_writing_to_file,
                      force_file_integrity_check):
    """ check for eventual unfinalised flag and try to adjust info content

        Parameters
        ----------------
        version : float
            version number
        info : class
            info mdf class
        fid:
            file handler
        finalization_writing : bool, optional, False by default
            If file is detected not finalised (id_unfin_flags!=0), file is corrected
            writing in the blocks if set to True, otherwise correction is done only
            to file representation in memory.
        force_file_integrity_check : bool, optional, False by default
            Perform block sizes check for potentially corrupted file without finalization
            flags (id_unfin_flags==0). Combined with finalization_writing_to_file is
            very experimental and risky, correction should be tried in memory first.

        Returns
        -----------
        info class
    """
    unfinalized_flags = info['ID']['id_unfi_flags']
    if force_file_integrity_check:
        unfinalized_flags = 21  # checks only for what is implemented
        version = 420
    if version >= 411 and unfinalized_flags > 0:
        # looking for all the block in file
        pattern = re.compile(
            rb"(?P<block>##(D[GVTZIL]|AT|C[AGHNC]|EV|FH|HL|LD|MD|R[DVI]|S[IRD]|TX))",
            re.DOTALL | re.MULTILINE,
        )

        fid.seek(0)

        addresses = []
        blocks = {}
        addresses_btype = OrderedDict()

        for match in re.finditer(pattern, fid.read()):
            btype = match.group('block')
            start = match.start()

            addresses_btype[start] = btype
            btype_addresses = blocks.setdefault(btype, [])
            btype_addresses.append(start)
            addresses.append(start)

        fid.seek(0, 2)
        eof_position = fid.tell()

        if not info['CG']:
            # minimal = 1 was used.
            for dataGroup in info['DG']:
                info.read_cg_blocks(fid, dataGroup)

        # treatment of unfinalised file
        if unfinalized_flags & (1 << 4):
            warn('Checking last DL block')
            # check if all DL block are actually full and present
            for dg in info['CG']:
                data_block_starting_position = info['DG'][dg]['dg_data']
                if data_block_starting_position == 0:
                    break  # no data block
                block_type = addresses_btype.get(data_block_starting_position, b'')
                if not block_type:
                    continue  # unknown block type, skip integrity check for this DG

                if block_type in frozenset((b'##DL', b'##HL')):
                    dl_position = data_block_starting_position
                    if block_type == b'##HL':
                        hl = HLBlock()
                        fid.seek(data_block_starting_position + 24)
                        hl.read_hl(fid)
                        dl_position = hl['hl_dl_first']

                    dl = DLBlock()
                    while True:
                        dl.update(_load_header(fid, dl_position))
                        dl.read_dl(fid, dl['link_count'])
                        for pointer in dl['list_data'][0]:
                            dt = _load_header(fid, pointer)
                        if not dl['next']:
                            # last DL block for this DG Block
                            break
                        dl_position = dl['next']

                    faulty_index = [index for index, f in enumerate(
                        dl['list_data'][0]) if f == 0]
                    if faulty_index:
                        # some data blocks were not written — filter them out (O(n))
                        n_faulty = len(faulty_index)
                        faulty_set = set(faulty_index)
                        dl['link_count'] -= n_faulty
                        dl['count'] -= n_faulty
                        dl['list_data'][0] = [v for i, v in enumerate(dl['list_data'][0])
                                              if i not in faulty_set]
                        dl['length'] -= n_faulty * 8
                        if not dl['flags'] & 0b1:  # offset list
                            dl['offset'] = [v for i, v in enumerate(dl['offset'])
                                            if i not in faulty_set]
                            dl['link_count'] -= n_faulty
                            dl['length'] -= n_faulty * 8
                        if dl['flags'] & 0b10:  # time values
                            dl['time_values'] = [v for i, v in enumerate(dl['time_values'])
                                                 if i not in faulty_set]
                            dl['link_count'] -= n_faulty
                            dl['length'] -= n_faulty * 8
                        if dl['flags'] & 0b100:  # angle values
                            dl['angle_values'] = [v for i, v in enumerate(dl['angle_values'])
                                                  if i not in faulty_set]
                            dl['link_count'] -= n_faulty
                            dl['length'] -= n_faulty * 8
                        if dl['flags'] & 0b1000:  # distance values
                            dl['distance_values'] = [v for i, v in enumerate(dl['distance_values'])
                                                     if i not in faulty_set]
                            dl['link_count'] -= n_faulty
                            dl['length'] -= n_faulty * 8
                        if finalization_writing_to_file:
                            dl.write_dl(fid, data_block_starting_position)

        if unfinalized_flags & (1 << 2) and not force_file_integrity_check:
            warn('Checking length for last DT block')

        if unfinalized_flags & 1 and not force_file_integrity_check:
            warn('Checking cycle counters for CG/CA blocks')

        if (unfinalized_flags & 1) | (unfinalized_flags & (1 << 2)):
            for dg in info['CG']:
                data_block_starting_position = info['DG'][dg]['dg_data']

                if len(info['CG'][dg]) == 1 and data_block_starting_position:
                    # not unsorted or empty data block
                    cg = info['DG'][dg]['dg_cg_first']
                    cg_cycle_count = info['CG'][dg][cg]['cg_cycle_count']
                    record_size = info['DG'][dg]['dg_rec_id_size'] + \
                        info['CG'][dg][cg]['cg_data_bytes'] + \
                        info['CG'][dg][cg]['cg_invalid_bytes']
                    expected_data_size = cg_cycle_count * record_size
                    block_type = addresses_btype[data_block_starting_position]

                    if block_type == (b'##DT', b'##DV'):
                        try:
                            data_block_ending_position = addresses[addresses.index(
                                data_block_starting_position) + 1]
                        except IndexError:  # last block at end of file
                            data_block_ending_position = eof_position
                        data_size_in_file = data_block_ending_position \
                            - data_block_starting_position - 24
                        dt = _load_header(fid, data_block_starting_position)
                        if expected_data_size > data_size_in_file > 0:
                            warn(f'DT Block at position {data_block_starting_position} '
                                 f'has a length of {data_size_in_file} '
                                 f'while it is expected to be {expected_data_size}')
                            corrected_data_size = data_size_in_file - data_size_in_file % 8
                            info['CG'][dg][cg]['cg_cycle_count'] = corrected_data_size // record_size
                            if finalization_writing_to_file:
                                # write in cg
                                fid.seek(info['CG'][dg][cg]['pointer'] + 88)
                                fid.write(
                                    bytes(info['CG'][dg][cg]['cg_cycle_count']))
                                # write in dt
                                fid.seek(dt['pointer'] + 8)
                                fid.write(bytes(corrected_data_size + 24))

                    elif block_type == b'##DL':
                        dl = DLBlock()
                        data_size_in_file = 0
                        dl_position = data_block_starting_position
                        data_block_size_in_file = 24
                        while True:
                            dl.update(_load_header(fid, dl_position))
                            dl.read_dl(fid, dl['link_count'])
                            for pointer in dl['list_data'][0]:
                                dt = _load_header(fid, pointer)
                                data_block_size_in_file = dt['length']
                                data_size_in_file += data_block_size_in_file - 24
                            if not dl['next']:
                                # last DL block for this DG Block
                                break
                            dl_position = dl['next']
                        data_size_in_file -= data_block_size_in_file - 24
                        if data_size_in_file > 0:
                            try:
                                dt_ending_position = addresses[addresses.index(
                                    pointer) + 1]
                            except IndexError:
                                dt_ending_position = eof_position
                            data_size_in_file += dt_ending_position - pointer - 24

                            if expected_data_size > data_size_in_file:
                                warn(f'Data Blocks have a length of {data_size_in_file} '
                                     f'while it is expected to be {expected_data_size}')
                                corrected_block_data_size = data_block_size_in_file \
                                    - data_block_size_in_file % 8
                                corrected_data_size = data_size_in_file \
                                    - data_block_size_in_file \
                                    + corrected_block_data_size
                                info['CG'][dg][cg]['cg_cycle_count'] = corrected_data_size // record_size
                                if finalization_writing_to_file:
                                    # write in cg
                                    fid.seek(info['CG'][dg][cg]
                                             ['pointer'] + 88)
                                    fid.write(
                                        bytes(info['CG'][dg][cg]['cg_cycle_count']))
                                    # write in dt
                                    fid.seek(dt['pointer'] + 8)
                                    fid.write(
                                        bytes(corrected_block_data_size + 24))

                    elif block_type in (b'##HL', b'##DZ'):
                        dz = DZBlock()
                        data_size_in_file = 0
                        data_block = bytearray()
                        last_block_is_DT = False
                        if block_type == b'##HL':
                            hl = HLBlock()
                            fid.seek(data_block_starting_position + 24)
                            hl.read_hl(fid)
                            dl = DLBlock()
                            dl_position = hl['hl_dl_first']
                            while True:
                                dl.update(_load_header(fid, dl_position))
                                dl.read_dl(fid, dl['link_count'])
                                dl['list_data'][0] = list(dl['list_data'][0])
                                last_block_pointer = dl['list_data'][0].pop(-1)
                                for pointer in dl['list_data'][0]:
                                    dz.update(_load_header(fid, pointer))
                                    if dz['id'] == b'##DZ':
                                        dz.read_dz(fid)
                                        data_size_in_file += dz['dz_org_data_length']
                                    else:  # could be DT
                                        header_size = dz['link_count'] * 8 + 24
                                        fid.seek(pointer + header_size)
                                        data_size_in_file += dz['length'] - 24
                                if not dl['next']:
                                    # last DL block for this DG Block
                                    break
                                dl_position = dl['next']
                            dz.update(_load_header(fid, last_block_pointer))
                            if dz['id'] == b'##DZ':
                                dz.read_dz(fid)
                                data_block_size_in_file = dz['dz_org_data_length']
                            else:  # last block could be a small DT
                                last_block_is_DT = True
                                data_block_size_in_file = dz['length'] - 24
                        else:  # single DZ block
                            dz.update(_load_header(
                                fid, data_block_starting_position))
                            dz.read_dz(fid)
                            data_block_size_in_file = dz['dz_org_data_length']

                        try:
                            dz_ending_position = addresses[addresses.index(
                                dz['pointer']) + 1]
                        except IndexError:
                            dz_ending_position = eof_position

                        if last_block_is_DT:  # last block is DT
                            data_size_in_file += dz_ending_position - \
                                dz['pointer'] - 24
                        else:  # block to be uncompressed to know its length
                            dz_data_section_length = dz_ending_position - \
                                dz['pointer'] - 48
                            if dz_data_section_length >= dz['dz_data_length']:
                                data_size_in_file += dz['dz_org_data_length']
                            else:  # block shorter than expected by block description
                                data_block = fid.read(dz_data_section_length)
                                data_size_in_file += len(
                                    decompress(data_block))

                        if expected_data_size > data_size_in_file:
                            warn(f'Data Blocks have a length of {data_size_in_file} '
                                 f'while it is expected to be {expected_data_size}')
                            corrected_block_data_size = data_block_size_in_file \
                                - data_block_size_in_file % 8
                            corrected_data_size = data_size_in_file \
                                - data_block_size_in_file \
                                + corrected_block_data_size
                            info['CG'][dg][cg]['cg_cycle_count'] = corrected_data_size // record_size
                            if finalization_writing_to_file:
                                # write in cg
                                fid.seek(info['CG'][dg][cg]['pointer'] + 88)
                                fid.write(
                                    bytes(info['CG'][dg][cg]['cg_cycle_count']))
                                if last_block_is_DT:
                                    fid.seek(dz['pointer'] + 8)
                                    fid.write(
                                        bytes(corrected_block_data_size + 24))
                                else:
                                    # write in dz
                                    fid.seek(dz['pointer'] + 32)
                                    fid.write(bytes(corrected_block_data_size))
                                    fid.seek(dz['pointer'] + 40)
                                    fid.write(bytes(dz_data_section_length))

        if unfinalized_flags & (1 << 5):
            warn('Update of cg_data_bytes and cg_inval_bytes in VLSD '
                 'CG block required')
            # No need to correct
        if unfinalized_flags & (1 << 6):
            warn('Update of offset values for VLSD channel required '
                 'in case a VLSD CG block is used')
            # offset values not used to parse VLSD, no need
        if unfinalized_flags & (1 << 1):
            warn('Update of cycle counters for SR blocks required')
            # SRBlocks not parsed
        if unfinalized_flags & (1 << 3):
            warn('Update of length for last RD block required')
            # RDBlocks not parsed
    return info
