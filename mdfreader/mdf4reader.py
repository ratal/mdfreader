# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 4.x.

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Thu Dec 10 12:57:28 2013


Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>
- bitarray to parse bits in not aligned bytes
- Sympy to convert channels with formula if needed
- zlib to uncompress data block if needed


Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.2+


mdf4reader
--------------------------

"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from __future__ import print_function
from struct import Struct
from struct import pack, unpack as structunpack
from math import pow
from io import open  # for python 3 and 2 consistency
from os.path import splitext
from time import gmtime, strftime
from multiprocessing import Queue, Process
from sys import version_info, byteorder
from collections import defaultdict, OrderedDict
from numpy.core.records import fromstring, fromarrays
from numpy import array, recarray, asarray, empty, where, frombuffer, reshape
from numpy import arange, right_shift, bitwise_and, all, diff, interp
from numpy import issubdtype, number as numpy_number
from numpy import max as npmax, min as npmin
from numpy.lib.recfunctions import rec_append_fields, rename_fields
from numpy.ma import MaskedArray
from warnings import simplefilter, warn
from .mdfinfo4 import Info4, IDBlock, HDBlock, DGBlock, \
    CGBlock, CNBlock, FHBlock, CommentBlock, _load_header, DLBlock, \
    DZBlock, HLBlock, CCBlock, DTBlock, CABlock
from .mdf import MdfSkeleton, _open_mdf, invalidChannel, dataField, \
    conversionField, idField, invalidPosField, CompressedData
from .channel import Channel4
try:
    from dataRead import data_read
    dataRead_available = True
except ImportError:
    warn('dataRead cannot be imported, compile it with Cython', ImportWarning)
    dataRead_available = False

PythonVersion = version_info
PythonVersion = PythonVersion[0]

chunk_size_reading = 100000000  # reads by chunk of 100Mb, can be tuned for best performance


def _data_block(record, info, parent_block, channel_set=None, n_records=None, sorted_flag=True, vlsd=False):
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
    if parent_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal data block
        if sorted_flag:
            if channel_set is None and not record.hiddenBytes and\
                    record.byte_aligned:  # No channel list and length of records corresponds to C datatypes
                # for debugging purpose
                # print(n_records, record.numpyDataRecordFormat, record.dataRecordName)
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=n_records)
            else:  # record is not byte aligned or channelSet not None
                return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        else:  # unsorted reading
            return _read_unsorted(record, info, parent_block, channel_set)

    elif parent_block['id'] in ('##SD', b'##SD'):
        return _read_sd_block(record[record.VLSD[0]].signal_data_type(info), parent_block['data'],
                              parent_block['length'] - 24)

    elif parent_block['id'] in ('##DZ', b'##DZ'):  # zipped data block
        # uncompress data
        parent_block['data'] = DZBlock.decompress_data_block(parent_block['data'], parent_block['dz_zip_type'],
                                                             parent_block['dz_zip_parameter'],
                                                             parent_block['dz_org_data_length'])
        if vlsd:  # VLSD channel
            return _read_sd_block(record[record.VLSD[0]].signal_data_type(info), parent_block['data'],
                                  parent_block['dz_org_data_length'])
        if channel_set is None and sorted_flag:  # reads all blocks if sorted block and no channelSet defined
            if record.byte_aligned and not record.hiddenBytes:
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=n_records)
            else:
                return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        elif channel_set is not None and sorted_flag:  # sorted data but channel list requested
            return record.read_channels_from_bytes(parent_block['data'], info, channel_set, n_records)
        else:  # unsorted reading
            return _read_unsorted(record, info, parent_block, channel_set)


def _read_unsorted(record, info, parent_block, channel_set=None):
    """ reads only the channels using offset functions, channel by channel within unsorted data

    Parameters
    ------------
    record : class
        record class
    info: info class
    channel_set : set of str, optional
        set of channel to read
    parent_block:
        number of records
    channel_set: set
    set of channel names to be parsed

    Returns
    --------
    buf : array
        data array

    Notes
    ------
    channel_set is not used yet
    """
    buf = defaultdict(list)
    position = 0
    record_id_c_format = record[list(record.keys())[0]]['record'].recordIDCFormat
    record_id_size = record[list(record.keys())[0]]['record'].recordIDsize
    VLSDStruct = Struct('I')
    # initialise data structure
    for record_id in record:
        for channelName in record[record_id]['record'].dataRecordName:
            buf[channelName] = []
    # read data
    while position < len(parent_block['data']):
        record_id = record_id_c_format.unpack(parent_block['data'][position:position + record_id_size])[0]
        if not record[record_id]['record'].Flags & 0b1:  # not VLSD CG)
            temp = record.read_record(record_id, info, parent_block['data'][position:position + record[record_id][
                'record'].CGrecordLength + 1])
            position += record[record_id]['record'].CGrecordLength
            for channelName in temp:
                buf[channelName].append(temp[channelName])
        else:  # VLSD CG
            position += record_id_size
            VLSDLen = VLSDStruct.unpack(parent_block['data'][position:position + 4])[0]  # VLSD length
            position += 4
            temp = parent_block['data'][position:position + VLSDLen - 1]
            signal_data_type = record[record_id]['record'].VLSD_CG[record_id]['channel'].signal_data_type(info)
            if signal_data_type == 6:
                temp = temp.decode('ISO8859')
            elif signal_data_type == 7:
                temp = temp.decode('utf-8')
            elif signal_data_type == 8:
                temp = temp.decode('<utf-16')
            elif signal_data_type == 9:
                temp = temp.decode('>utf-16')
            buf[record[record_id]['record'].VLSD_CG[record_id]['channelName']].append(temp)
            position += VLSDLen
    # convert list to array
    for chan in buf:
        buf[chan] = array(buf[chan])
    return buf


def _read_sd_block(signal_data_type, sd_block, sd_block_length):
    """ Reads vlsd channel from its SD Block bytes

    Parameters
    ----------------
    signal_data_type : int

    sd_block : bytes
    SD Block bytes

    sd_block_length: int
    SD Block data length (header not included)

    Returns
    -----------
    array
    """
    pointer = 0
    buf = []
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
        while pointer < sd_block_length:
            VLSDLen = structunpack('I', sd_block[pointer:pointer + 4])[0]  # length of data
            pointer += 4
            buf.append(sd_block[pointer:pointer + VLSDLen].decode(channel_format).rstrip('\x00'))
            pointer += VLSDLen
        buf = _equalize_string_length(buf)
    else:  # byte arrays or mime types
        while pointer < sd_block_length:
            VLSDLen = structunpack('I', sd_block[pointer:pointer + 4])[0]  # length of data
            pointer += 4
            buf.append(sd_block[pointer:pointer + VLSDLen])
            pointer += VLSDLen
        # buf = _equalize_byte_length(buf)
    return array(buf)


def _equalize_string_length(buf):
    """ Makes all strings in a list having same length by appending spaces strings.

    Parameters
    ----------------
    buf : list of str

    Returns
    -----------
    list of str elements all having same length
    """
    max_len = len(max(buf, key=len))
    for i, element in enumerate(buf):  # resize string to same length, numpy constrain
        buf[i] = ''.join([element, ' ' * (max_len - len(element))])
    return buf


def _equalize_byte_length(buf):
    """ Makes all bytes in a list having same length by appending null bytes.

    Parameters
    ----------------
    buf : list of bytes

    Returns
    -----------
    list of bytes all having same length
    """
    max_len = len(max(buf, key=len))
    for i, element in enumerate(buf):  # resize string to same length, numpy constrain
        buf[i] = bytearray(element).extend(b'0' * (max_len - len(element)))
    return buf


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
        self[record.recordID] = {}
        self[record.recordID]['record'] = record
        # detects VLSD CG
        for recordID in self[record.recordID]['record'].VLSD_CG:
            self[recordID]['record'].VLSD_CG = self[record.recordID]['record'].VLSD_CG

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
            recordID = list(self.keys())[0]
            record = self[recordID]['record']
            self[recordID]['data'] = self.load(record, info, name_list=channel_set, sorted_flag=True)
            for cn in record.VLSD:  # VLSD channels
                if channel_set is None or record[cn].name in channel_set:
                    temp = Data(self.fid, record[cn].data(info))  # all channels
                    temp = temp.load(record, info, name_list=channel_set, sorted_flag=True, vlsd=True)
                    self[recordID]['data'] = rename_fields(self[recordID]['data'],
                                                           {record[cn].name: '{}_offset'.format(record[cn].name)})
                    self[recordID]['data'] = rec_append_fields(self[recordID]['data'],
                                                               record[cn].name, temp)
        else:  # unsorted DataGroup
            self.type = 'unsorted'
            data = self.load(self, info, name_list=channel_set, sorted_flag=False)
            for recordID in self:
                self[recordID]['data'] = {}
                for channel in self[recordID]['record']:
                    self[recordID]['data'][channel.name] = data[channel.name]

    def load(self, record, info, name_list=None, sorted_flag=True, vlsd=False):
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
        vlsd : bool
            indicate a sd block, compressed (DZ) or not (SD)

        Returns
        -----------
        numpy recarray of data
        """
        temps = defaultdict()
        # block header
        temps.update(_load_header(self.fid, self.pointer_to_data))
        if temps['id'] in ('##DL', b'##DL'):  # data list block
            temp = DLBlock()
            temp.read(self.fid, temps['link_count'])
            temps.update(temp)
            if temps['dl_dl_next']:
                index = 1
            while temps['dl_dl_next']:  # reads pointers to all data blocks (DT, RD, SD, DZ)
                temp = defaultdict()
                temp.update(_load_header(self.fid, temps['dl_dl_next']))
                temps['dl_dl_next'] = structunpack('<Q', self.fid.read(8))[0]
                temps['dl_data'][index] = \
                    structunpack('<{}Q'.
                                 format(temp['link_count'] - 1),
                                 self.fid.read(8 * (temp['link_count'] - 1)))
                index += 1
            if temps['dl_count']:
                # read and concatenate raw blocks
                previous_index = 0
                data_block = defaultdict()
                data_block['data'] = bytearray()
                if vlsd:  # need to load all blocks as variable length, cannot process block by block
                    for DL in temps['dl_data']:
                        for pointer in temps['dl_data'][DL]:
                            # read fist data blocks linked by DLBlock to identify data block type
                            data_block.update(_load_header(self.fid, pointer))
                            if data_block['id'] in ('##SD', b'##SD'):
                                data_block['data'].extend(self.fid.read(data_block['length'] - 24))
                            elif data_block['id'] in ('##DZ', b'##DZ'):
                                temp = DZBlock()
                                temp.read(self.fid)
                                data_block.update(temp)
                                data_block['data'].extend(DZBlock.decompress_data_block(
                                    self.fid.read(data_block['dz_data_length']),
                                    data_block['dz_zip_type'],
                                    data_block['dz_zip_parameter'],
                                    data_block['dz_org_data_length']))
                                if isinstance(data_block['dz_org_block_type'], str):
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'])
                                else:
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'].decode('ASCII'))
                    temps['data'] = _data_block(record, info, parent_block=data_block, channel_set=name_list,
                                                n_records=None, sorted_flag=sorted_flag, vlsd=vlsd)
                else:
                    for DL in temps['dl_data']:
                        for pointer in temps['dl_data'][DL]:
                            # read fist data blocks linked by DLBlock to identify data block type
                            data_block.update(_load_header(self.fid, pointer))
                            if data_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):
                                data_block['data'].extend(self.fid.read(data_block['length'] - 24))
                            elif data_block['id'] in ('##DZ', b'##DZ'):
                                temp = DZBlock()
                                temp.read(self.fid)
                                data_block.update(temp)
                                data_block['data'].extend(DZBlock.decompress_data_block(
                                                          self.fid.read(data_block['dz_data_length']),
                                                          data_block['dz_zip_type'],
                                                          data_block['dz_zip_parameter'],
                                                          data_block['dz_org_data_length']))
                                if isinstance(data_block['dz_org_block_type'], str):
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'])
                                else:
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'].decode('ASCII'))
                            nrecord_chunk = len(data_block['data']) // record.CGrecordLength
                            nremain = len(data_block['data']) % record.CGrecordLength
                            if nremain:
                                remain = data_block['data'][-nremain:]
                                del data_block['data'][-nremain:]
                            if previous_index + nrecord_chunk > record.numberOfRecords:
                                # there could be more data than needed for the expected number of records
                                nrecord_chunk = record.numberOfRecords - previous_index
                            tmp = _data_block(record, info, parent_block=data_block, channel_set=name_list,
                                              n_records=nrecord_chunk, sorted_flag=sorted_flag, vlsd=vlsd)
                            if 'data' not in temps:  # initialise recarray
                                temps['data'] = recarray(record.numberOfRecords, dtype=tmp.dtype)
                            temps['data'][previous_index: previous_index + nrecord_chunk] = tmp
                            previous_index += nrecord_chunk
                            if nremain:
                                data_block['data'] = remain
                            else:
                                data_block['data'] = bytearray()  # flush
            else:  # empty datalist
                temps['data'] = None
        elif temps['id'] in ('##HL', b'##HL'):  # header list block for DZBlock
            # link section
            temp = HLBlock()
            temp.read(self.fid)
            temps.update(temp)
            self.pointer_to_data = temps['hl_dl_first']
            temps['data'] = self.load(record, info, name_list=name_list, sorted_flag=sorted_flag, vlsd=vlsd)
        elif temps['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal sorted data block, direct read
            temps['data'] = record.read_sorted_record(self.fid, info, channel_set=name_list)
        elif temps['id'] in ('##SD', b'##SD'):  # VLSD
            temps['data'] = self.fid.read(temps['length'] - 24)
            temps['data'] = _data_block(record, info, parent_block=temps, channel_set=name_list, n_records=None,
                                        sorted_flag=sorted_flag)
        elif temps['id'] in ('##DZ', b'##DZ'):  # zipped data block
            temp = DZBlock()
            temp.read(self.fid)
            temps.update(temp)
            temps['data'] = self.fid.read(temps['dz_data_length'])
            temps['data'] = _data_block(record, info, parent_block=temps, channel_set=name_list, n_records=None,
                                        sorted_flag=sorted_flag, vlsd=vlsd)
        else:
            raise Exception('unknown data block')
        return temps['data']

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


class Record(list):
    __slots__ = ['CGrecordLength', 'recordLength', 'numberOfRecords', 'recordID',
                 'recordIDsize', 'recordIDCFormat', 'dataGroup', 'channelGroup',
                 'numpyDataRecordFormat', 'dataRecordName', 'master',
                 'recordToChannelMatching', 'channelNames', 'Flags', 'VLSD_CG',
                 'VLSD', 'MLSD', 'byte_aligned', 'hiddenBytes', 'invalid_channel',
                 'CANOpen']
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
        self.master = {}
        self.master['name'] = 'master_{}'.format(data_group)
        self.master['number'] = None
        self.Flags = 0
        self.VLSD_CG = {}
        self.VLSD = []
        self.MLSD = {}
        self.recordToChannelMatching = {}
        self.channelNames = set()
        self.byte_aligned = True
        self.hiddenBytes = False
        self.invalid_channel = None
        self.CANOpen = None

    def __str__(self):
        output = list()
        output.append('Record class content print\nTotal number of channels : {}\n'.format(len(self)))
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
        if self.master['name'] is not None:
            output.append('Master channel : {}\n'.format(self.master['name']))
        output.append('Numpy records format : \n')
        for record in self.numpyDataRecordFormat:
            output.append(' {} '.format(record))
        output.append('\nVLSD_CG {}\n'.format(self.VLSD_CG))
        return ''.join(output)

    def add_channel(self, info, channel_number):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class
        """
        channel = Channel4()
        self.append(channel.set(info, self.dataGroup, self.channelGroup, channel_number))
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
                self.recordIDCFormat = 'H'
                self.recordLength += 2
                self.CGrecordLength += 2
            elif self.recordIDsize == 3:
                self.numpyDataRecordFormat.append('uint32')
                self.recordIDCFormat = 'I'
                self.recordLength += 3
                self.CGrecordLength += 3
            elif self.recordIDsize == 4:
                self.numpyDataRecordFormat.append('uint64')
                self.recordIDCFormat = 'L'
                self.recordLength += 4
                self.CGrecordLength += 4
        self.recordID = info['CG'][self.dataGroup][self.channelGroup]['cg_record_id']
        self.numberOfRecords = info['CG'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        self.Flags = info['CG'][self.dataGroup][self.channelGroup]['cg_flags']
        if 'MLSD' in info:
            self.MLSD = info['MLSD']
        embedding_channel = None
        for channelNumber in info['CN'][self.dataGroup][self.channelGroup]:
            channel = Channel4()
            channel.set(info, self.dataGroup, self.channelGroup, channelNumber)
            channel_type = channel.channel_type(info)
            data_format = channel.data_format(info)
            if channel_type in (2, 3):  # master channel found
                if self.master['number'] is None or channel.channel_sync_type(info) == 1:
                    # new master channel found
                    # or more than 1 master channel, priority to time channel
                    self.master['number'] = channelNumber
                    self.master['name'] = channel.name
            if channel_type in (0, 1, 2, 4, 5):  # not virtual channel
                signal_data_type = channel.signal_data_type(info)
                if signal_data_type == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        channel = Channel4()  # new object otherwise only modified
                        channel.set_CANOpen(info, self.dataGroup, self.channelGroup, channelNumber, name)
                        self.append(channel)
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(channel.data_format(info))
                    self.recordLength += 7
                    self.CANOpen = 'date'
                    embedding_channel = None
                elif signal_data_type == 14:
                    for name in ('ms', 'days'):
                        channel = Channel4()
                        channel.set_CANOpen(info, self.dataGroup, self.channelGroup, channelNumber, name)
                        self.append(channel)
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(channel.data_format(info))
                    self.recordLength += 6
                    self.CANOpen = 'time'
                    embedding_channel = None
                else:
                    self.append(channel)
                    self.channelNames.add(channel.name)
                    # Checking if several channels are embedded in bytes
                    if len(self) > 1:
                        # all channels are already ordered in record based on byte_offset
                        # and bit_offset so just comparing with previous channel
                        channel_pos_bit_end = channel.pos_bit_end(info)
                        channel_pos_bit_beg = channel.pos_bit_beg(info)
                        prev_chan = self[-2]
                        prev_chan_byte_offset = prev_chan.byteOffset
                        prev_chan_n_bytes = prev_chan.nBytes_aligned
                        prev_chan_includes_curr_chan = channel_pos_bit_beg >= 8 * prev_chan_byte_offset \
                            and channel_pos_bit_end <= 8 * (prev_chan_byte_offset + prev_chan_n_bytes)
                        if embedding_channel is not None:
                            embedding_channel_includes_curr_chan = \
                                channel_pos_bit_end <= embedding_channel.pos_byte_end(info) * 8
                        else:
                            embedding_channel_includes_curr_chan = False
                        if channel.byteOffset >= prev_chan_byte_offset and \
                                channel_pos_bit_beg < 8 * (prev_chan_byte_offset +
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
                    if 'VLSD_CG' in info:  # is there VLSD CG
                        for recordID in info['VLSD_CG']:  # look for VLSD CG Channel
                            if info['VLSD_CG'][recordID]['cg_cn'] == (self.channelGroup, channelNumber):
                                self.VLSD_CG[recordID] = info['VLSD_CG'][recordID]
                                self.VLSD_CG[recordID]['channel'] = channel
                                self.VLSD_CG[recordID]['channelName'] = channel.name
                                self[-1].VLSD_CG_Flag = True
                                break
                    if channel_type == 1:  # VLSD channel
                        self.VLSD.append(channel.channelNumber)
            elif channel_type in (3, 6):  # virtual channel
                self.append(channel)  # channel calculated based on record index later in conversion function
                self.channelNames.add(channel.name)
                self.recordToChannelMatching[channel.name] = channel.name

        if info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']:  # invalid bytes existing
            self.CGrecordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            self.recordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            invalid_bytes = Channel4()
            invalid_bytes.set_invalid_bytes(info, self.dataGroup, self.channelGroup, channelNumber + 1)
            self.invalid_channel = invalid_bytes
            self.append(self.invalid_channel)
            self.channelNames.add(self.invalid_channel.name)
            self.recordToChannelMatching[self.invalid_channel.name] = \
                self.invalid_channel.name
            self.numpyDataRecordFormat.append(self.invalid_channel.data_format(info))
            self.dataRecordName.append(self.invalid_channel.name)
        # check for hidden bytes
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True
        # check record length consistency
        elif self.CGrecordLength < self.recordLength:
            self.byte_aligned = False  # forces to use dataRead instead of numpy records.

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
            return self.read_all_channels_sorted_record(fid)
        else:  # reads only some channels from a sorted data block
            if channel_set is None or len(channel_set & self.channelNames) > 0:
                return self.read_not_all_channels_sorted_record(fid, info, channel_set)

    def generate_chunks(self):
        """ calculate data split

        Returns
        --------
        (n_record_chunk, chunk_size)
        """
        n_chunks = (self.CGrecordLength * self.numberOfRecords) // chunk_size_reading + 1
        chunk_length = (self.CGrecordLength * self.numberOfRecords) // n_chunks
        n_record_chunk = chunk_length // self.CGrecordLength
        chunks = [(n_record_chunk, self.CGrecordLength * n_record_chunk)] * n_chunks
        n_record_chunk = self.numberOfRecords - n_record_chunk * n_chunks
        if n_record_chunk > 0:
            chunks.append((n_record_chunk, self.CGrecordLength * n_record_chunk))
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
                fromstring(fid.read(chunk_size),
                           dtype={'names': self.dataRecordName,
                                  'formats': self.numpyDataRecordFormat},
                           shape=n_record_chunk)
            previous_index += n_record_chunk
        return buf

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
        if channel_set is not None and not self.master['name'] in channel_set:
            channel_set.add(self.master['name'])  # adds master channel
        rec, channels_indexes = self.initialise_recarray(info, channel_set, self.numberOfRecords)
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
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channel_set and not Channel.VLSD_CG_Flag:
                temp[Channel.name] = \
                    Channel.c_format_structure(info).unpack(buf[Channel.pos_byte_beg(info):Channel.pos_byte_end(info)])[0]
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
            for chan in range(len(self)):
                if self[chan].name in channel_set and self[chan].channel_type(info) not in (3, 6):
                    # not virtual channel and part of channelSet
                    channels_indexes.append(chan)
                    formats.append(self[chan].native_data_format(info))
                    names.append(self[chan].name)
            if formats:
                rec = recarray(n_records, dtype={'names': names, 'formats': formats})
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
            buf, channels_indexes = self.initialise_recarray(info, channel_set, n_records, dtype, channels_indexes)
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
                        data_read(bytes_data, self[chan].bit_count(info),
                                  self[chan].signal_data_type(info),
                                  self[chan].native_data_format(info),
                                  n_records, self.CGrecordLength,
                                  self[chan].bit_offset(info), self[chan].pos_byte_beg(info),
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
            buf, channels_indexes = self.initialise_recarray(info, channel_set, n_records, dtype, channels_indexes)
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
                    temp = [bit_array[self[chan].pos_bit_beg(info) + record_bit_size * i:
                            self[chan].pos_bit_end(info) + record_bit_size * i]
                            for i in range(n_records)]
                    n_bytes = len(temp[0].tobytes())
                    if not n_bytes == n_bytes_estimated and \
                            signal_data_type not in (6, 7, 8, 9, 10, 11, 12):  # not Ctype byte length
                        byte = bitarray(8 * (n_bytes_estimated - n_bytes), endian='little')
                        byte.setall(False)
                        if signal_data_type not in (2, 3):  # not signed integer
                            for i in range(n_records):  # extend data of bytes to match numpy requirement
                                temp[i].extend(byte)
                        else:  # signed integer (two's complement), keep sign bit and extend with bytes
                            temp = signed_int(temp, byte)
                    n_trail_bits = n_bytes_estimated*8 - self[chan].bit_count(info)
                    if signal_data_type in (2, 3) and \
                            n_bytes == n_bytes_estimated and \
                            n_trail_bits > 0:  # C type byte length but signed integer
                        trail_bits = bitarray(n_trail_bits, endian='little')
                        temp = signed_int(temp, trail_bits)
                else:  # Channel Array
                    temp = [bit_array[self[chan].pos_bit_beg(info) + record_bit_size * i:
                                      self[chan].pos_bit_beg(info) + 8 * n_bytes_estimated + record_bit_size * i]
                            for i in range(n_records)]
                if 's' not in self[chan].c_format(info):
                    c_structure = self[chan].c_format_structure(info)
                    if ('>' in self[chan].data_format(info) and byteorder == 'little') or \
                       (byteorder == 'big' and '<' in self[chan].data_format(info)):
                        temp = [c_structure.unpack(temp[i].tobytes())[0]
                                for i in range(n_records)]
                        temp = asarray(temp).byteswap().newbyteorder()
                    else:
                        temp = [c_structure.unpack(temp[i].tobytes())[0]
                                for i in range(n_records)]
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
              compression=False, metadata=2):
        """ Reads mdf 4.x file data and stores it in dict

        Parameters
        ----------------
        file_name : str, optional
            file name

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        multi_processed : bool
            flag to activate multiprocessing of channel data conversion

        channel_list : list of str, optional
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory.
            Can be used to read big files

        convert_after_read : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false,
            all data from channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .get_channel_data()

        compression : bool, optional
            flag to activate data compression with blosc

        metadata: int, optional, default = 2
            Reading metadata has impact on performance, especially for mdf 4.x using xml.
            2: minimal metadata reading (mostly channel blocks)
            1: used for noDataLoading
            0: all metadata reading, including Source Information, Attachment, etc..

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
                info = Info4(self.fileName, None, minimal=minimal)
            else:
                info = self.info

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
            file_date_time = gmtime(info['HD']['hd_start_time_ns'] / 1000000000)
            ddate = strftime('%Y-%m-%d', file_date_time)
            ttime = strftime('%H:%M:%S', file_date_time)

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
                                  project=project, subject=subject, comment=comment,
                                  date=ddate, time=ttime)
            else:
                self.add_metadata(date=ddate, time=ttime)

        data_groups = info['DG']  # parse all data groups
        if self._noDataLoading and channel_list is not None:
            data_groups = [self[channel][idField][0] for channel in channel_list]

        for dataGroup in data_groups:
            channel_set = channel_set_file
            if not info['DG'][dataGroup]['dg_data'] == 0 and \
                    (channel_set is None or
                     len(channel_set & info['ChannelNamesByDG'][dataGroup]) > 0):  # there is data block and channel in
                if minimal > 1 and not self._noDataLoading:  # load CG, CN and CC block info
                    info.read_cg_block(info.fid, dataGroup, channel_set, minimal=minimal)
                if info['CG'][dataGroup][0]['cg_cycle_count']:  # data exists
                    # Pointer to data block
                    pointer_to_data = info['DG'][dataGroup]['dg_data']

                    if 'dataClass' not in info['DG'][dataGroup]:
                        buf = Data(info.fid, pointer_to_data)
                        for channelGroup in info['CG'][dataGroup]:
                            temp = Record(dataGroup, channelGroup)  # create record class
                            temp.load_info(info)  # load all info related to record
                            buf.add_record(temp)  # adds record to DATA
                            record_id = info['CG'][dataGroup][channelGroup]['cg_record_id']
                            if temp.master['name'] is not None \
                                    and buf[record_id]['record'].channelNames:
                                if channel_set is not None and not self._noDataLoading\
                                        and temp.master['name'] not in channel_set:
                                    channel_set.add(temp.master['name'])  # adds master channel in channelSet if missing
                            if channel_set is not None and buf[record_id]['record'].CANOpen:
                                # adds CANOpen channels if existing in not empty channelSet
                                if buf[record_id]['record'].CANOpen == 'time':
                                    channel_set.update(('ms', 'days'))
                                elif buf[record_id]['record'].CANOpen == 'date':
                                    channel_set.update(('ms', 'minute', 'hour', 'day', 'month', 'year'))
                        if self._noDataLoading:
                            self.info['DG'][dataGroup]['dataClass'] = buf
                    else:
                        buf = self.info['DG'][dataGroup]['dataClass']

                    # reads raw data from data block with DATA and _data_block classes
                    buf.read(channel_set, info, self.fileName)

                    channel_groups = buf
                    if self._noDataLoading and channel_list is not None:
                        channel_groups = [info['CG'][dataGroup][self[channel][idField][1]]['cg_record_id']
                                          for channel in channel_list]

                    # processing data from buf then transfer to self
                    for record_id in channel_groups:  # for each channel group in data block
                        if 'record' in buf[record_id]:
                            master_channel = buf[record_id]['record'].master['name']

                            channels = buf[record_id]['record']
                            if self._noDataLoading and channel_list is not None:
                                channels = [channels[self[channel][idField][2]] for channel in channel_list]
                            for chan in channels:  # for each channel class
                                if channel_set is None or chan.name in channel_set:
                                    if not chan.type == 4:  # normal channel
                                        if chan.channel_type(info) not in (3, 6):  # not virtual channel
                                            # in case record is used for several channels
                                            if channel_set is None and not buf[record_id]['record'].hiddenBytes \
                                                    and buf[record_id]['record'].byte_aligned:
                                                record_name = buf[record_id]['record'].recordToChannelMatching[chan.name]
                                            else:
                                                record_name = chan.name
                                            if 'data' in buf[record_id] and \
                                                    buf[record_id]['data'] is not None:  # no data in channel group
                                                temp = buf[record_id]['data'][record_name]  # extract channel vector
                                            else:
                                                temp = None
                                        else:  # virtual channel
                                            temp = arange(buf[record_id]['record'].numberOfRecords)

                                        # Process concatenated bits inside uint8
                                        bit_count = chan.bit_count(info)
                                        if buf[record_id]['record'].byte_aligned \
                                                and not buf[record_id]['record'].hiddenBytes and \
                                                channel_set is None and\
                                                0 < bit_count < 64 and bit_count not in (8, 16, 32) \
                                                and temp is not None\
                                                and temp.dtype.kind not in ('S', 'U'):
                                            # if channel data do not use complete bytes and Ctypes
                                            signal_data_type = chan.signal_data_type(info)
                                            if signal_data_type in (0, 1, 2, 3):  # integers
                                                bit_offset = chan.bit_offset(info)
                                                if bit_offset > 0:
                                                    temp = right_shift(temp, bit_offset)
                                                mask = int(pow(2, bit_count) - 1)  # masks isBitUnit8
                                                temp = bitwise_and(temp, mask)
                                                if signal_data_type in (2, 3):
                                                    # signed integer, moving bit sign of two's complement
                                                    sign_bit_mask = (1 << (bit_count - 1))
                                                    sign_extend = ((1 << (temp.itemsize * 8 - bit_count)) - 1) << bit_count
                                                    sign_bit = bitwise_and(temp, sign_bit_mask)
                                                    for number, sign in enumerate(sign_bit):
                                                        # negative value, sign extend
                                                        if sign:
                                                            temp[number] |= sign_extend
                                            else:  # should not happen
                                                warn('bit count and offset not applied to correct '
                                                     'data type {}'.format(chan.name))
                                        else:  # data using full bytes
                                            pass

                                        if temp is not None:  # channel contains data
                                            # string data decoding
                                            if temp.dtype.kind == 'S':
                                                signal_data_type = chan.signal_data_type(info)
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
                                                    temp2 = empty(len(temp), dtype='U{}'.format(temp.dtype.str[-1]))
                                                    for t in range(temp.size):
                                                        try:
                                                            temp2[t] = temp[t].decode(encoding, 'ignore')
                                                        except:
                                                            warn('Cannot decode channel {}'.format(chan.name))
                                                            temp2[t] = ''
                                                    temp = temp2

                                            # channel creation
                                            self.add_channel(chan.name, temp, master_channel,
                                                             master_type=chan.channel_sync_type(info),
                                                             unit=chan.unit(info), description=chan.desc(info),
                                                             conversion=chan.conversion(info), info=chan.cn_block(info),
                                                             compression=compression,
                                                             identifier=info.unique_id(chan.dataGroup,
                                                                                       chan.channelGroup,
                                                                                       chan.channelNumber))
                                            if chan.channel_type(info) == 4:  # sync channel
                                                # attach stream to be synchronised
                                                self.set_channel_attachment(chan.name, chan.attachment(info.fid, info))
                                            if chan.has_invalid_bit(info):  # has invalid bit
                                                self.set_invalid_bit(chan.name, chan.invalid_bit(info))
                                                self.set_invalid_channel(chan.name, 'invalid_bytes{}'.format(dataGroup))
                                    else:
                                        # invalid bytes channel
                                        data = buf[record_id]['data'].__getattribute__(chan.name)
                                        data = frombuffer(data.tobytes(), dtype='u1').reshape(len(data),
                                                                                              data.dtype.itemsize)
                                        self.add_channel(chan.name, data, master_channel, master_type=0, unit='',
                                                         description='', info=None, compression=compression,
                                                         identifier=None)
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
                    (self.info.fid, self.info.fileName, self.info.zipfile) = _open_mdf(self.fileName)
                self.read4(file_name=None, info=None, channel_list=[channel_name], convert_after_read=False)
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
                vector = channel[dataField].decompression()  # uncompressed blosc data
            else:
                vector = channel[dataField][:]  # to have bcolz uncompressed data
        if conversionField in channel and channel[conversionField]['type']:  # there is conversion property
            text_type = vector.dtype.kind in ['S', 'U', 'V']  # channel of string or not ?
            conversion_type = channel[conversionField]['type']
            conversion_parameter = channel[conversionField]['parameters']
            if conversion_type == 1 and not text_type:
                vector = _linear_conversion(vector, conversion_parameter['cc_val'])
            elif conversion_type == 2 and not text_type:
                vector = _rational_conversion(vector, conversion_parameter['cc_val'])
            elif conversion_type == 3 and not text_type:
                vector = _formula_conversion(vector, conversion_parameter['cc_ref']['Comment'])
            elif conversion_type == 4 and not text_type:
                vector = _value_to_value_table_with_interpolation_conversion(vector, conversion_parameter['cc_val'])
            elif conversion_type == 5 and not text_type:
                vector = _value_to_value_table_without_interpolation_conversion(vector, conversion_parameter['cc_val'])
            elif conversion_type == 6 and not text_type and convert_tables:
                vector = _value_range_to_value_table_conversion(vector, conversion_parameter['cc_val'])
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
                vector = _text_to_text_conversion(vector, conversion_parameter['cc_ref'])
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
        self.set_channel_data(channel_name, self._get_channel_data4(channel_name))
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
                        if conversion['type'] in (1, 2):  # more time in multi proc
                            self._convert_channel4(channelName)
                        else:
                            proc.append(Process(target=self._convert_channel_data4,
                                                args=(channel, channelName, self.convertTables, True, Q)))
                            proc[-1].start()
                for p in proc:
                    L.update(Q.get())  # concatenate results of processes in dict
                for p in proc:
                    p.join()
                del Q  # free memory
                for channelName in self:
                    if channelName in L:
                        self.set_channel_data(channelName, L[channelName])
                        self.remove_channel_conversion(channelName)

    def write4(self, file_name=None, compression=False):
        """Writes simple mdf 4.1 file

        Parameters
        ----------------
        file_name : str, optional
            Name of file
            If file name is not input, written file name will be the one read with appended '_new' string before extension
        compression : bool
            flag to store data compressed

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
        temp.write(fid)

        blocks = OrderedDict()
        pointer = 64
        # Header Block
        blocks['HD'] = HDBlock()
        blocks['HD']['block_start'] = pointer
        pointer += 104

        # Header Block comments
        blocks['HD']['MD'] = pointer
        blocks['HD_comment'] = CommentBlock()
        blocks['HD_comment']['block_start'] = pointer
        blocks['HD_comment'].load(self.fileMetadata, 'HD')
        pointer = blocks['HD_comment']['block_start'] + blocks['HD_comment']['block_length']

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
        pointer = blocks['FH_comment']['block_start'] + blocks['FH_comment']['block_length']

        # write DG block
        blocks['HD']['DG'] = pointer  # first DG

        # write all files header blocks
        for block in blocks.values():
            block.write(fid)
        if self.masterChannelList:  # some channels exist
            for dataGroup, masterChannel in enumerate(self.masterChannelList):
                # writes dataGroup Block
                dg = DGBlock()
                dg['block_start'] = pointer
                pointer = dg['block_start'] + 64
                dg['CG'] = pointer  # First CG link

                blocks = OrderedDict()  # initialise blocks for this datagroup
                # write CGBlock
                blocks['CG'] = CGBlock()
                blocks['CG']['block_start'] = pointer
                pointer = blocks['CG']['block_start'] + 104
                blocks['CG']['CN'] = pointer  # First CN link

                master_channel_data = self._get_channel_data4(masterChannel)
                if master_channel_data is not None:
                    cg_cycle_count = len(master_channel_data)
                elif self._get_channel_data4(self.masterChannelList[masterChannel][0]).shape[0] == 1:  # classification
                    cg_cycle_count = 1
                elif master_channel_data not in self.masterChannelList[masterChannel]:
                    cg_cycle_count = len(self._get_channel_data4(self.masterChannelList[masterChannel][0]))
                    warn('no master channel in data group {}'.format(dataGroup))
                else:
                    # no data in data group, skip
                    warn('no data in data group {0} with master channel {1}'.format(dataGroup, masterChannel))
                    continue
                blocks['CG']['cg_cycle_count'] = cg_cycle_count

                # write channels
                record_byte_offset = 0
                cn_flag = 0
                number_of_channel = 0
                n_records = 0
                data_list = ()
                last_channel = 0
                for n_channel, channel in enumerate(self.masterChannelList[masterChannel]):
                    data = self.get_channel_data(channel)
                    # no interest to write invalid bytes as channel, should be processed if needed before writing
                    if channel.find('invalid_bytes') == -1 and data is not None and len(data) > 0:
                        last_channel = n_channel
                        data_ndim = data.ndim - 1
                        if not data_ndim:
                            data_list = data_list + (data, )
                        else:  # data contains arrays
                            data_dim_size = data.shape
                            if not cg_cycle_count == data_dim_size[0]:
                                warn('Array length do not match number of cycled in CG block')
                            data_dim_size = data_dim_size[1:]
                            SNd = 0
                            PNd = 1
                            for x in data_dim_size:
                                SNd += x
                                PNd *= x
                            flattened = reshape(data, (cg_cycle_count, PNd))
                            for i in range(PNd):
                                data_list = data_list + (flattened[:, i],)
                        number_of_channel += 1
                        blocks[n_channel] = CNBlock()
                        if issubdtype(data.dtype, numpy_number):  # is numeric
                            blocks[n_channel]['cn_val_range_min'] = npmin(data)
                            blocks[n_channel]['cn_val_range_max'] = npmax(data)
                            blocks[n_channel]['cn_flags'] = 8  # only Bit 3: Limit range valid flag
                        else:
                            blocks[n_channel]['cn_val_range_min'] = 0
                            blocks[n_channel]['cn_val_range_max'] = 0
                            blocks[n_channel]['cn_flags'] = 0
                        if masterChannel is not channel:
                            blocks[n_channel]['cn_type'] = 0
                            blocks[n_channel]['cn_sync_type'] = 0
                        else:
                            blocks[n_channel]['cn_type'] = 2  # master channel
                            blocks[n_channel]['cn_sync_type'] = 1  # default is time channel
                            n_records = len(data)

                        cn_numpy_kind = data.dtype.kind
                        if cn_numpy_kind in ('u', 'b'):
                            data_type = 0  # LE
                        elif cn_numpy_kind == 'i':
                            data_type = 2  # LE
                        elif cn_numpy_kind == 'f':
                            data_type = 4  # LE
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
                        blocks[n_channel]['cn_bit_offset'] = 0  # always byte aligned
                        blocks[n_channel]['cn_byte_offset'] = record_byte_offset
                        byte_count = data.dtype.itemsize
                        record_byte_offset += byte_count
                        blocks[n_channel]['cn_bit_count'] = byte_count * 8
                        blocks[n_channel]['block_start'] = pointer
                        pointer = blocks[n_channel]['block_start'] + 160

                        # arrays handling
                        if not data_ndim:
                            blocks[n_channel]['Composition'] = 0
                        else:
                            blocks[n_channel]['Composition'] = pointer  # pointer to CABlock
                            # creates CABlock
                            ca = ''.join([channel, '_CA'])
                            blocks[ca] = CABlock()
                            blocks[ca]['block_start'] = pointer
                            blocks[ca]['ndim'] = data_ndim
                            blocks[ca]['ndim_size'] = data_dim_size
                            blocks[ca].load(data.itemsize)
                            pointer = blocks[ca]['block_start'] + blocks[ca]['block_length']

                        if cn_flag:
                            # Next DG
                            blocks[n_channel-1]['CN'] = blocks[n_channel]['block_start']
                            blocks[n_channel]['CN'] = 0  # initialise 'CN' key
                        else:
                            cn_flag = blocks[n_channel]['block_start']
                            blocks[n_channel]['CN'] = 0  # creates first CN link, null for the moment

                        # write channel name
                        blocks[n_channel]['TX'] = pointer
                        blocks[channel] = CommentBlock()
                        blocks[channel]['block_start'] = pointer
                        blocks[channel].load(channel, 'TX')
                        pointer = blocks[channel]['block_start'] + blocks[channel]['block_length']

                        # write channel unit
                        unit = self.get_channel_unit(channel)
                        if unit is not None and len(unit) > 0:
                            blocks[n_channel]['Unit'] = pointer
                            unit_name = '{}{}{}'.format(channel, '_U_', n_channel)
                            blocks[unit_name] = CommentBlock()
                            blocks[unit_name]['block_start'] = pointer
                            blocks[unit_name].load(unit, 'TX')
                            pointer = blocks[unit_name]['block_start'] + blocks[unit_name]['block_length']
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
                            pointer = blocks[desc_name]['block_start'] + blocks[desc_name]['block_length']
                        else:
                            blocks[n_channel]['Comment'] = 0

                if n_records == 0 and masterChannel is not self.masterChannelList[masterChannel]:
                    # No master channel in channel group
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
                pointer = data.write(fid, fromarrays(data_list).tobytes(order='F'))
                if compression:
                    # next DG position is not predictable due to DZ Blocks unknown length
                    fid.seek(dg_start_position + 24)
                    fid.write(pack('<Q', pointer))
                    fid.seek(pointer)

            fid.seek(dg['block_start'] + 24)
            fid.write(pack('Q', 0))  # last DG pointer is null

        fid.close()

    def apply_invalid_bit(self, channel_name):
        """Mask data of channel based on its invalid bit definition if there is

        Parameters
        ----------------
        channel_name : str
            Name of channel
        """
        try:
            invalid_bit_pos = self.get_invalid_bit(channel_name)
            if isinstance(invalid_bit_pos, int):  # invalid bit existing
                mask = self._get_channel_data4(self.get_invalid_channel(channel_name))
                data = self._get_channel_data4(channel_name)
                data = data.view(MaskedArray)
                invalid_byte = invalid_bit_pos >> 3
                data.mask = bitwise_and(mask[:, invalid_byte], invalid_bit_pos & 0x07)
                self.set_channel_data(channel_name, data)
                self._remove_channel_field(channel_name, invalidPosField)
                self._remove_channel_field(channel_name, invalidChannel)
        except KeyError:
            pass
            # warn('no invalid data found for channel ')

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
                (ndg, ncg, ncn), (cn, cs, cp), (gn, gs, gp) = self[channel_name]['id']
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
    except:
        warn('Please install sympy to convert channel ')
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
    int_val = [cc_val[i] for i in range(0, val_count, 2)]
    phys_val = [cc_val[i] for i in range(1, val_count, 2)]
    if all(diff(int_val) > 0):
        try:
            from scipy import interpolate
            f = interpolate.interp1d(int_val, phys_val, kind='nearest', bounds_error=False)  # nearest
            return f(vector)  # fill with Nan out of bounds while should be bounds
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
    int_val = [cc_val[i] for i in range(0, val_count, 2)]
    phys_val = [cc_val[i] for i in range(1, val_count, 2)]
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
    key_min = [cc_val[i] for i in range(0, 3 * val_count + 1, 3)]
    key_max = [cc_val[i] for i in range(1, 3 * val_count + 1, 3)]
    value = [cc_val[i] for i in range(2, 3 * val_count + 1, 3)]
    # look up in range keys
    for l_index in range(len(vector)):
        key_index = 0  # default index if not found
        for i in range(val_count):
            if key_min[i] < vector[l_index] < key_max[i]:
                key_index = i
                break
        vector[l_index] = value[key_index]
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
    maxlen = max([len(ref) for ref in cc_ref])
    temp = empty(len(vector), dtype='U{}'.format(maxlen))  # initialize empty array with max length
    # checks for scaling
    try:
        from sympy import lambdify, symbols
    except:
        warn('Please install sympy to convert channel ')
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
        elif (PythonVersion > 3 and not isinstance(cc_ref[ref], str)) or \
                (PythonVersion < 3 and not isinstance(cc_ref[ref], unicode)):  # identity, non conversion
            cc_ref[ref] = lambdify(X, 'X', modules='numpy', dummify=False)
    key_index = where(vector[0] == cc_val)[0]  # look up for first value in vector
    if not len(key_index) == 0:  # value corresponding in cc_val
        temp[0] = cc_ref[key_index[0]]
    else:  # default value
        if callable(cc_ref[-1]):
            temp[0] = cc_ref[-1](vector[0])
        else:
            temp[0] = cc_ref[-1]
    for lindex in range(1, len(vector)):
        if vector[lindex] == vector[lindex - 1]:  # same value as before, no need to look further
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
    except:
        warn('Please install sympy to convert channel ')
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
