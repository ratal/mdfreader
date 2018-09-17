# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 4.x.

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

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

mdf4reader module
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
from numpy import arange, right_shift, bitwise_and, bitwise_or, all, diff, interp
from numpy import issubdtype, number as numpy_number
from numpy import max as npmax, min as npmin
from numpy.lib.recfunctions import rec_append_fields, rename_fields
from numpy.ma import MaskedArray
from warnings import simplefilter, warn
from .mdfinfo4 import info4, IDBlock, HDBlock, DGBlock, \
    CGBlock, CNBlock, FHBlock, CommentBlock, _loadHeader, DLBlock, \
    DZBlock, HLBlock, CCBlock, DTBlock, CABlock
from .mdf import mdf_skeleton, _open_MDF, invalidChannel, dataField, \
    conversionField, idField, invalidPosField, compressed_data
from .channel import channel4
try:
    from dataRead import dataRead
    dataRead_available = True
except ImportError:
    warn('dataRead cannot be imported, compile it with cython', ImportWarning)
    dataRead_available = False

PythonVersion = version_info
PythonVersion = PythonVersion[0]

chunk_size_reading = 100000000  # reads by chunk of 100Mb, can be tuned for best performance


def DATABlock(record, info, parent_block, channelSet=None, nrecords=None, sortedFlag=True, vlsd=False):
    """ DATABlock converts raw data into arrays

    Parameters
    ----------------
    record : class
        record class instance describing a channel group record
    parent_block : class
        MDFBlock class containing at least parent block header
    channelSet : set of str, optional
        defines set of channels to only read, can be slow but saves memory,
        for big files
    nrecords: int, optional
        number of records to read
    sortedFlag : bool, optional
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
    if nrecords is None and hasattr(record, 'numberOfRecords'):
        nrecords = record.numberOfRecords
    if parent_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal data block
        if sortedFlag:
            if channelSet is None and not record.hiddenBytes and\
                    record.byte_aligned:  # No channel list and length of records corresponds to C datatypes
                # for debugging purpose
                # print(nrecords, record.numpyDataRecordFormat, record.dataRecordName)
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=nrecords)
            else:  # record is not byte aligned or channelSet not None
                return record.read_channels_from_bytes(parent_block['data'], info, channelSet, nrecords)
        else:  # unsorted reading
            return readUnsorted(record, info, parent_block, channelSet)

    elif parent_block['id'] in ('##SD', b'##SD'):
        return read_sdblock(record[record.VLSD[0]].signalDataType(info),
                            parent_block['data'], parent_block['length'] - 24)

    elif parent_block['id'] in ('##DZ', b'##DZ'):  # zipped data block
        # uncompress data
        parent_block['data'] = DZBlock.decompress_datablock(parent_block['data'], parent_block['dz_zip_type'],
                parent_block['dz_zip_parameter'], parent_block['dz_org_data_length'])
        if vlsd:  # VLSD channel
            return read_sdblock(record[record.VLSD[0]].signalDataType(info),
                                parent_block['data'], parent_block['dz_org_data_length'])
        if channelSet is None and sortedFlag:  # reads all blocks if sorted block and no channelSet defined
            if record.byte_aligned and not record.hiddenBytes:
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=nrecords)
            else:
                return record.read_channels_from_bytes(parent_block['data'], info, channelSet, nrecords)
        elif channelSet is not None and sortedFlag:  # sorted data but channel list requested
            return record.read_channels_from_bytes(parent_block['data'], info, channelSet, nrecords)
        else:  # unsorted reading
            return readUnsorted(record, info, parent_block, channelSet)


def readUnsorted(record, info, parent_block, channelSet=None):
    # reads only the channels using offset functions, channel by channel.
    buf = defaultdict(list)
    position = 0
    recordIdCFormat = record[list(record.keys())[0]]['record'].recordIDCFormat
    recordIDsize = record[list(record.keys())[0]]['record'].recordIDsize
    VLSDStruct = Struct('I')
    # initialise data structure
    for recordID in record:
        for channelName in record[recordID]['record'].dataRecordName:
            buf[channelName] = []  # empty(record.numberOfRecords,dtype=record[recordID]['record'].dataFormat)
            # index[channelName]=0
    # read data
    while position < len(parent_block['data']):
        recordID = recordIdCFormat.unpack(parent_block['data'][position:position + recordIDsize])[0]
        if not record[recordID]['record'].Flags & 0b1:  # not VLSD CG)
            temp = record.readRecord(recordID, info, parent_block['data'][position:position + record[recordID]['record'].CGrecordLength + 1], channelSet)
            position += record[recordID]['record'].CGrecordLength
            for channelName in temp:
                buf[channelName].append(temp[channelName])
        else:  # VLSD CG
            position += recordIDsize
            VLSDLen = VLSDStruct.unpack(parent_block['data'][position:position + 4])[0]  # VLSD length
            position += 4
            temp = parent_block['data'][position:position + VLSDLen - 1]
            signalDataType = record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType(info)
            if signalDataType == 6:
                temp = temp.decode('ISO8859')
            elif signalDataType == 7:
                temp = temp.decode('utf-8')
            elif signalDataType == 8:
                temp = temp.decode('<utf-16')
            elif signalDataType == 9:
                temp = temp.decode('>utf-16')
            buf[record[recordID]['record'].VLSD_CG[recordID]['channelName']].append(temp)
            position += VLSDLen
    # convert list to array
    for chan in buf:
        buf[chan] = array(buf[chan])
    return buf


def read_sdblock(signal_data_type, sdblock, sdblock_length):
    """ Reads vlsd channel from its SD Block bytes

        Parameters
        ----------------
        signal_data_type : int

        sdblock : bytes
        SD Block bytes

        sdblock_length: int
        SD Block data length (header not included)

        Returns
        -----------
        array
    """
    if signal_data_type == 6:
        channel_format = 'ISO8859'
    elif signal_data_type == 7:
        channel_format = 'utf-8'
    elif signal_data_type == 8:
        channel_format = '<utf-16'
    elif signal_data_type == 9:
        channel_format = '>utf-16'
    pointer = 0
    buf = []

    while pointer < sdblock_length:
        VLSDLen = structunpack('I', sdblock[pointer:pointer + 4])[0]  # length of data
        pointer += 4
        buf.append(sdblock[pointer:pointer + VLSDLen].decode(channel_format).rstrip('\x00'))
        pointer += VLSDLen
    buf = equalizeStringLength(buf)
    return array(buf)


def equalizeStringLength(buf):
    """ Makes all strings in a list having same length by appending spaces strings.

    Parameters
    ----------------
    buf : list of str

    Returns
    -----------
    list of str elements all having same length
    """
    maxlen = len(max(buf, key=len))
    for i in range(len(buf)):  # resize string to same length, numpy constrain
        buf[i] = ''.join([buf[i], ' ' * (maxlen - len(buf[i]))])
    return buf


class DATA(dict):
    __slots__ = ['fid', 'pointerTodata', 'type']
    """ DATA class is organizing record classes itself made of channel class.
    This class inherits from dict. Keys are corresponding to channel group recordID
    A DATA class corresponds to a data block, a dict of record classes (one per channel group)
    Each record class contains a list of channel class representing the structure of channel record.

    Attributes
    --------------
    fid : io.open
        file identifier
    pointerToData : int
        position of Data block in mdf file
    type : str
        'sorted' or 'unsorted' data block

    Methods
    ------------
    addRecord(record)
        Adds a new record in DATA class dict
    read(channelSet, zip=None)
        Reads data block
    load(record, zip=None, nameList=None)
        Reads sorted data block from record definition
    readRecord(recordID, buf, channelSet=None):
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
        self.pointerTodata = pointer
        self.type = 'sorted'

    def addRecord(self, record):
        """Adds a new record in DATA class dict.

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        """
        self[record.recordID] = {}
        self[record.recordID]['record'] = record
        # detects VLSD CG
        for recordID in self[record.recordID]['record'].VLSD_CG:
            self[recordID]['record'].VLSD_CG = self[record.recordID]['record'].VLSD_CG

    def read(self, channelSet, info, filename):
        """Reads data block

        Parameters
        ----------------
        channelSet : set of str
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
            self[recordID]['data'] = self.load(record, info, nameList=channelSet, sortedFlag=True)
            for cn in record.VLSD:  # VLSD channels
                if channelSet is None or record[cn].name in channelSet:
                    temp = DATA(self.fid, record[cn].data(info))  # all channels
                    temp = temp.load(record, info, nameList=channelSet, sortedFlag=True, vlsd=True)
                    self[recordID]['data'] = rename_fields(self[recordID]['data'],
                                                           {record[cn].name: '{}_offset'.format(record[cn].name)})
                    self[recordID]['data'] = rec_append_fields(self[recordID]['data'],
                                                               record[cn].name, temp)
        else:  # unsorted DataGroup
            self.type = 'unsorted'
            data = self.load(self, info, nameList=channelSet, sortedFlag=False)
            for recordID in self:
                self[recordID]['data'] = {}
                for channel in self[recordID]['record']:
                    self[recordID]['data'][channel.name] = data[channel.name]

    def load(self, record, info, nameList=None, sortedFlag=True, vlsd=False):
        """Reads data block from record definition

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        info class
            contains blocks
        nameList : list of str, optional
            list of channel names
        sortedFlag : bool, optional
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
        temps.update(_loadHeader(self.fid, self.pointerTodata))
        if temps['id'] in ('##DL', b'##DL'):  # data list block
            temp = DLBlock()
            temp.read(self.fid, temps['link_count'])
            temps.update(temp)
            if temps['dl_dl_next']:
                index = 1
            while temps['dl_dl_next']:  # reads pointers to all data blocks (DT, RD, SD, DZ)
                temp = defaultdict()
                temp.update(_loadHeader(self.fid, temps['dl_dl_next']))
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
                            data_block.update(_loadHeader(self.fid, pointer))
                            if data_block['id'] in ('##SD', b'##SD'):
                                data_block['data'].extend(self.fid.read(data_block['length'] - 24))
                            elif data_block['id'] in ('##DZ', b'##DZ'):
                                temp = DZBlock()
                                temp.read(self.fid)
                                data_block.update(temp)
                                data_block['data'].extend(DZBlock.decompress_datablock(
                                    self.fid.read(data_block['dz_data_length']),
                                    data_block['dz_zip_type'],
                                    data_block['dz_zip_parameter'],
                                    data_block['dz_org_data_length']))
                                if isinstance(data_block['dz_org_block_type'], str):
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'])
                                else:
                                    data_block['id'] = '##{}'.format(data_block['dz_org_block_type'].decode('ASCII'))
                    temps['data'] = DATABlock(record, info, parent_block=data_block,
                                              channelSet=nameList, nrecords=None,
                                              sortedFlag=sortedFlag, vlsd=vlsd)
                else:
                    for DL in temps['dl_data']:
                        for pointer in temps['dl_data'][DL]:
                            # read fist data blocks linked by DLBlock to identify data block type
                            data_block.update(_loadHeader(self.fid, pointer))
                            if data_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):
                                data_block['data'].extend(self.fid.read(data_block['length'] - 24))
                            elif data_block['id'] in ('##DZ', b'##DZ'):
                                temp = DZBlock()
                                temp.read(self.fid)
                                data_block.update(temp)
                                data_block['data'].extend(DZBlock.decompress_datablock(
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
                            tmp = DATABlock(record, info, parent_block=data_block,
                                            channelSet=nameList, nrecords=nrecord_chunk,
                                            sortedFlag=sortedFlag, vlsd=vlsd)
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
            self.pointerTodata = temps['hl_dl_first']
            temps['data'] = self.load(record, info, nameList=nameList, sortedFlag=sortedFlag, vlsd=vlsd)
        elif temps['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal sorted data block, direct read
            temps['data'] = record.readSortedRecord(self.fid, info, channelSet=nameList)
        elif temps['id'] in ('##SD', b'##SD'):  # VLSD
            temps['data'] = self.fid.read(temps['length'] - 24)
            temps['data'] = DATABlock(record, info, parent_block=temps, channelSet=nameList,
                                      nrecords=None, sortedFlag=sortedFlag)
        elif temps['id'] in ('##DZ', b'##DZ'):  # zipped data block
            temp = DZBlock()
            temp.read(self.fid)
            temps.update(temp)
            temps['data'] = self.fid.read(temps['dz_data_length'])
            temps['data'] = DATABlock(record, info, parent_block=temps,
                                      channelSet=nameList, nrecords=None,
                                      sortedFlag=sortedFlag, vlsd=vlsd)
        else:
            raise Exception('unknown data block')
        return temps['data']

    def readRecord(self, recordID, info, buf, channelSet=None):
        """ read record from a buffer

        Parameters
        ----------------
        recordID : int
            record identifier
        info class
            contains blocks
        buf : str
            buffer of data from file to be converted to channel raw data
        channelSet : set of str
            setof channel names to be read
        """
        return self[recordID]['record'].readRecordBuf(buf, info, channelSet)


class record(list):
    __slots__ = ['CGrecordLength', 'recordLength', 'numberOfRecords', 'recordID',
                 'recordIDsize', 'recordIDCFormat', 'dataGroup', 'channelGroup',
                 'numpyDataRecordFormat', 'dataRecordName', 'master',
                 'recordToChannelMatching', 'channelNames', 'Flags', 'VLSD_CG',
                 'VLSD', 'MLSD', 'byte_aligned', 'hiddenBytes', 'invalid_channel',
                 'CANOpen']
    """ record class listing channel classes. It is representing a channel group

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
    addChannel(info, channelNumber)
    loadInfo(info)
    readSortedRecord(fid, pointer, info, channelSet=None)
    readRecordBuf(buf, info, channelSet=None)
    read_channels_from_bytes(bita, info, channelSet=None)
    """

    def __init__(self, dataGroup, channelGroup):
        self.CGrecordLength = 0
        self.recordLength = 0
        self.numberOfRecords = 0
        self.recordID = 0
        self.recordIDsize = 0
        self.recordIDCFormat = ''
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.numpyDataRecordFormat = []
        self.dataRecordName = []
        self.master = {}
        self.master['name'] = 'master_{}'.format(dataGroup)
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
        output.append('CG block record bytes length : {}\nDatagroup number : {}'
                      '\nByte aligned : {}\nHidden bytes : {}\n'.format(self.CGrecordLength, self.dataGroup,
                                                                        self.byte_aligned, self.hiddenBytes))
        if self.master['name'] is not None:
            output.append('Master channel : {}\n'.format(self.master['name']))
        output.append('Numpy records format : \n')
        for record in self.numpyDataRecordFormat:
            output.append(' {} '.format(record))
        output.append('\nVLSD_CG {}\n'.format(self.VLSD_CG))
        return ''.join(output)

    def addChannel(self, info, channelNumber):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        """
        Channel = channel4()
        self.append(Channel.set(info, self.dataGroup, self.channelGroup, channelNumber))
        self.channelNames.add(self[-1].name)

    def loadInfo(self, info):
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
            Channel = channel4()
            Channel.set(info, self.dataGroup, self.channelGroup, channelNumber)
            channelType = Channel.channelType(info)
            dataFormat = Channel.dataFormat(info)
            if channelType in (2, 3):  # master channel found
                if self.master['number'] is None or Channel.channelSyncType(info) == 1:
                    # new master channel found
                    # or more than 1 master channel, priority to time channel
                    self.master['number'] = channelNumber
                    self.master['name'] = Channel.name
            if channelType in (0, 1, 2, 4, 5):  # not virtual channel
                signalDataType = Channel.signalDataType(info)
                if signalDataType == 13:
                    for name in ('ms', 'minute', 'hour', 'day', 'month', 'year'):
                        Channel = channel4()  # new object otherwise only modified
                        Channel.setCANOpen(info, self.dataGroup, self.channelGroup, channelNumber, name)
                        self.append(Channel)
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(Channel.dataFormat(info))
                    self.recordLength += 7
                    self.CANOpen = 'date'
                    embedding_channel = None
                elif signalDataType == 14:
                    for name in ('ms', 'days'):
                        Channel = channel4()
                        Channel.setCANOpen(info, self.dataGroup, self.channelGroup, channelNumber, name)
                        self.append(Channel)
                        self.channelNames.add(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(Channel.dataFormat(info))
                    self.recordLength += 6
                    self.CANOpen = 'time'
                    embedding_channel = None
                else:
                    self.append(Channel)
                    self.channelNames.add(Channel.name)
                    # Checking if several channels are embedded in bytes
                    if len(self) > 1:
                        # all channels are already ordered in record based on byte_offset
                        # and bit_offset so just comparing with previous channel
                        Channel_posBitEnd = Channel.posBitEnd(info)
                        Channel_posBitBeg = Channel.posBitBeg(info)
                        prev_chan = self[-2]
                        prev_chan_byteOffset = prev_chan.byteOffset(info)
                        prev_chan_nBytes = prev_chan.nBytes(info)
                        prev_chan_includes_curr_chan = Channel_posBitBeg >= 8 * prev_chan_byteOffset \
                                and Channel_posBitEnd <= 8 * (prev_chan_byteOffset + prev_chan_nBytes)
                        if embedding_channel is not None:
                            embedding_channel_includes_curr_chan = \
                                Channel_posBitEnd <= embedding_channel.posByteEnd(info) * 8
                        else:
                            embedding_channel_includes_curr_chan = False
                        if Channel.byteOffset(info) >= prev_chan_byteOffset and \
                                Channel_posBitBeg < 8 * (prev_chan_byteOffset + prev_chan_nBytes) < Channel_posBitEnd:
                            # not byte aligned
                            self.byte_aligned = False
                        if embedding_channel is not None and \
                                Channel_posBitEnd > embedding_channel.posByteEnd(info) * 8:
                            embedding_channel = None
                        if prev_chan_includes_curr_chan or \
                                embedding_channel_includes_curr_chan:  # bit(s) in byte(s)
                            if embedding_channel is None and prev_chan_includes_curr_chan:
                                embedding_channel = prev_chan  # new embedding channel detected
                            if self.recordToChannelMatching: # not first channel
                                self.recordToChannelMatching[Channel.name] = \
                                    self.recordToChannelMatching[prev_chan.name]
                            else: # first channels
                                self.recordToChannelMatching[Channel.name] = Channel.name
                                self.numpyDataRecordFormat.append(dataFormat)
                                self.dataRecordName.append(Channel.name)
                                self.recordLength += Channel.nBytes(info)
                    if embedding_channel is None:  # adding bytes
                        self.recordToChannelMatching[Channel.name] = Channel.name
                        self.numpyDataRecordFormat.append(dataFormat)
                        self.dataRecordName.append(Channel.name)
                        self.recordLength += Channel.nBytes(info)
                    if 'VLSD_CG' in info:  # is there VLSD CG
                        for recordID in info['VLSD_CG']:  # look for VLSD CG Channel
                            if info['VLSD_CG'][recordID]['cg_cn'] == (self.channelGroup, channelNumber):
                                self.VLSD_CG[recordID] = info['VLSD_CG'][recordID]
                                self.VLSD_CG[recordID]['channel'] = Channel
                                self.VLSD_CG[recordID]['channelName'] = Channel.name
                                self[-1].VLSD_CG_Flag = True
                                break
                    if channelType == 1:  # VLSD channel
                        self.VLSD.append(Channel.channelNumber)
            elif channelType in (3, 6):  # virtual channel
                self.append(Channel)  # channel calculated based on record index later in conversion function
                self.channelNames.add(Channel.name)
                self.recordToChannelMatching[Channel.name] = Channel.name

        if info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']:  # invalid bytes existing
            self.CGrecordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            self.recordLength += info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            invalid_bytes = channel4()
            invalid_bytes.setInvalidBytes(info, self.dataGroup, self.channelGroup, channelNumber + 1)
            self.invalid_channel = invalid_bytes
            self.append(self.invalid_channel)
            self.channelNames.add(self.invalid_channel.name)
            self.recordToChannelMatching[self.invalid_channel.name] = \
                self.invalid_channel.name
            self.numpyDataRecordFormat.append(self.invalid_channel.dataFormat(info))
            self.dataRecordName.append(self.invalid_channel.name)
        # check for hidden bytes
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True
        # check record length consistency
        elif self.CGrecordLength < self.recordLength:
            self.byte_aligned = False  # forces to use dataRead instead of numpy records.

    def readSortedRecord(self, fid, info, channelSet=None):
        """ reads record, only one channel group per datagroup

        Parameters
        ----------------
        fid :
            file identifier
        pointer
            position in file of data block beginning
        channelSet : set of str, optional
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
        if channelSet is None and self.byte_aligned and not self.hiddenBytes:
            return self.read_all_channels_sorted_record(fid)
        else:  # reads only some channels from a sorted data block
            if channelSet is None or len(channelSet & self.channelNames) > 0:
                return self.read_not_all_channels_sorted_record(fid, info, channelSet)

    def generate_chunks(self):
        """ calculate data split

        Returns
        --------
        (nrecord_chunk, chunk_size)
        """
        nchunks = (self.CGrecordLength * self.numberOfRecords) // chunk_size_reading + 1
        chunk_length = (self.CGrecordLength * self.numberOfRecords) // nchunks
        nrecord_chunk = chunk_length // self.CGrecordLength
        chunks = [(nrecord_chunk, self.CGrecordLength * nrecord_chunk)] * nchunks
        nrecord_chunk = self.numberOfRecords - nrecord_chunk * nchunks
        if nrecord_chunk > 0:
            chunks.append((nrecord_chunk, self.CGrecordLength * nrecord_chunk))
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
        for nrecord_chunk, chunk_size in chunks:
            buf[previous_index: previous_index + nrecord_chunk] = \
                fromstring(fid.read(chunk_size),
                           dtype={'names': self.dataRecordName,
                                  'formats': self.numpyDataRecordFormat},
                           shape=nrecord_chunk)
            previous_index += nrecord_chunk
        return buf

    def read_not_all_channels_sorted_record(self, fid, info, channelSet):
        """ reads channels from file listed in channelSet

        Parameters
        ------------
        fid :
            file identifier
        info: info class
        channelSet : set of str, optional
            set of channel to read

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        chunks = self.generate_chunks()
        previous_index = 0
        if channelSet is None:
            channelSet = self.channelNames
        if channelSet is not None and not self.master['name'] in channelSet:
            channelSet.add(self.master['name'])  # adds master channel
        rec, channels_indexes = self.initialise_recarray(info, channelSet, self.numberOfRecords)
        if rec is not None:
            if dataRead_available:
                for nrecord_chunk, chunk_size in chunks:
                    rec[previous_index: previous_index + nrecord_chunk] = \
                        self.read_channels_from_bytes(fid.read(chunk_size),
                                                      info, channelSet, nrecord_chunk,
                                                      rec.dtype, channels_indexes)
                    previous_index += nrecord_chunk
                return rec
            else:
                for nrecord_chunk, chunk_size in chunks:
                    rec[previous_index: previous_index + nrecord_chunk] = \
                        self.read_channels_from_bytes_fallback(fid.read(chunk_size),
                                                               info, channelSet, nrecord_chunk,
                                                               rec.dtype, channels_indexes)
                    previous_index += nrecord_chunk
                return rec
        else:
            return []

    def readRecordBuf(self, buf, info, channelSet=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        info class
            contains blocks structure
        channelSet : set of str, optional
            set of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        temp = {}
        if channelSet is None:
            channelSet = self.channelNames
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channelSet and not Channel.VLSD_CG_Flag:
                temp[Channel.name] = \
                    Channel.CFormat(info).unpack(buf[Channel.posByteBeg(info):Channel.posByteEnd(info)])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def initialise_recarray(self, info, channelSet, nrecords, dtype=None, channels_indexes=None):
        """ Initialise recarray

        Parameters
        ------------
        info: info class
        channelSet : set of str, optional
            set of channel to read
        nrecords: int
            number of records
        dtype: numpy dtype, optional
        channels_indexes: list of int, optional

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        if dtype is not None and channels_indexes is not None:
            return recarray(nrecords, dtype=dtype), channels_indexes
        else:
            if channelSet is None:
                channelSet = self.channelNames
            formats = []
            names = []
            channels_indexes = []
            for chan in range(len(self)):
                if self[chan].name in channelSet and self[chan].channelType(info) not in (3, 6):
                    # not virtual channel and part of channelSet
                    channels_indexes.append(chan)
                    formats.append(self[chan].nativedataFormat(info))
                    names.append(self[chan].name)
            if formats:
                rec = recarray(nrecords, dtype={'names': names, 'formats': formats})
                return rec, channels_indexes
            else:
                return None, []

    def read_channels_from_bytes(self, bita, info, channelSet=None, nrecords=None,
                                 dtype=None, channels_indexes=None):
        """ reads stream of record bytes using dataRead module if available otherwise bitarray
        
        Parameters
        ------------
        bita : stream
            stream of bytes
        info: info class
        channelSet : set of str, optional
            set of channel to read
        nrecords: int
            number of records
        dtype: numpy dtype
        channels_indexes: list of int
        
        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        if nrecords is None:
            nrecords = self.numberOfRecords
        # initialise recarray
        if dtype is None:
            buf, channels_indexes = self.initialise_recarray(info, channelSet, nrecords, dtype, channels_indexes)
        else:
            buf = recarray(nrecords, dtype=dtype)
        if buf is not None:  # at least some channels should be parsed
            if dataRead_available:  # use rather cython compiled code for performance
                bytesdata = bytes(bita)
                for chan in channels_indexes:
                    if self[chan].isCABlock(info):
                        ca = self[chan].CABlock(info)
                        array_flag = ca['ca_ndim']
                    else:
                        array_flag = 0
                    buf[self[chan].name] = \
                        dataRead(bytesdata, self[chan].bitCount(info),
                                 self[chan].signalDataType(info),
                                 self[chan].nativedataFormat(info),
                                 nrecords, self.CGrecordLength,
                                 self[chan].bitOffset(info), self[chan].posByteBeg(info),
                                 self[chan].posByteEnd(info), array_flag)
                return buf
            else:
                return self.read_channels_from_bytes_fallback(bita, info, channelSet, nrecords, dtype)
        else:
            return []

    def read_channels_from_bytes_fallback(self, bita, info, channelSet=None, nrecords=None,
                                          dtype=None, channels_indexes=None):
        """ reads stream of record bytes using bitarray in case no dataRead available

        Parameters
        ------------
        bita : stream
            stream of bytes
        info: info class
        channelSet : set of str, optional
            set of channel to read
        nrecords: int
            number of records
        dtype: numpy dtype
        channels_indexes: list of int

        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """

        def signedInt(temp, extension):
            """ extend bits of signed data managing two's complement
            """
            extension.setall(False)
            extensionInv = bitarray(extension, endian='little')
            extensionInv.setall(True)
            for i in range(nrecords):  # extend data of bytes to match numpy requirement
                signBit = temp[i][-1]
                if not signBit:  # positive value, extend with 0
                    temp[i].extend(extension)
                else:  # negative value, extend with 1
                    signBit = temp[i].pop(-1)
                    temp[i].extend(extensionInv)
                    temp[i].append(signBit)
            return temp
        if nrecords is None:
            nrecords = self.numberOfRecords
        if dtype is None:
            buf, channels_indexes = self.initialise_recarray(info, channelSet, nrecords, dtype, channels_indexes)
        else:
            buf = recarray(nrecords, dtype=dtype)
        if buf is not None:
            # read data
            from bitarray import bitarray
            B = bitarray(endian="little")  # little endian by default
            B.frombytes(bytes(bita))
            record_bit_size = self.CGrecordLength * 8
            for chan in channels_indexes:
                signalDataType = self[chan].signalDataType(info)
                nBytes = self[chan].nBytes(info)
                if not self[chan].type in ('CA', 'NestCA'):
                    temp = [B[self[chan].posBitBeg(info) + record_bit_size * i:
                            self[chan].posBitEnd(info) + record_bit_size * i]
                            for i in range(nrecords)]
                    nbytes = len(temp[0].tobytes())
                    if not nbytes == nBytes and \
                            signalDataType not in (6, 7, 8, 9, 10, 11, 12):  # not Ctype byte length
                        byte = bitarray(8 * (nBytes - nbytes), endian='little')
                        byte.setall(False)
                        if signalDataType not in (2, 3):  # not signed integer
                            for i in range(nrecords):  # extend data of bytes to match numpy requirement
                                temp[i].extend(byte)
                        else:  # signed integer (two's complement), keep sign bit and extend with bytes
                            temp = signedInt(temp, byte)
                    nTrailBits = nBytes*8 - self[chan].bitCount(info)
                    if signalDataType in (2, 3) and \
                            nbytes == nBytes and \
                            nTrailBits > 0:  # Ctype byte length but signed integer
                        trailBits = bitarray(nTrailBits, endian='little')
                        temp = signedInt(temp, trailBits)
                else:  # Channel Array
                    temp = [B[self[chan].posBitBeg(info) + record_bit_size * i:
                              self[chan].posBitBeg(info) + 8 * nBytes + record_bit_size * i]
                            for i in range(nrecords)]
                if 's' not in self[chan].Format(info):
                    CFormat = self[chan].CFormat(info)
                    if ('>' in self[chan].dataFormat(info) and byteorder == 'little') or \
                       (byteorder == 'big' and '<' in self[chan].dataFormat(info)):
                        temp = [CFormat.unpack(temp[i].tobytes())[0]
                                for i in range(nrecords)]
                        temp = asarray(temp).byteswap().newbyteorder()
                    else:
                        temp = [CFormat.unpack(temp[i].tobytes())[0]
                                for i in range(nrecords)]
                        temp = asarray(temp)
                else:
                    temp = [temp[i].tobytes()
                            for i in range(nrecords)]
                    temp = asarray(temp)
                buf[self[chan].name] = temp
            return buf
        else:
            return []


class mdf4(mdf_skeleton):

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
    file_metadata : dict
        file metadata with minimum keys : author, organisation, project, subject, comment, time, date

    Methods
    ------------
    read4( fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True)
        Reads mdf 4.x file data and stores it in dict
    _getChannelData4(channelName)
        Returns channel numpy array
    _convertChannel4(channelName)
        converts specific channel from raw to physical data according to CCBlock information
    _convertAllChannel4()
        Converts all channels from raw data to converted data according to CCBlock information
    """

    def read4(self, fileName=None, info=None, multiProc=False, channelList=None,
              convertAfterRead=True, filterChannelNames=False, compression=False):
        """ Reads mdf 4.x file data and stores it in dict

        Parameters
        ----------------
        fileName : str, optional
            file name

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        multiProc : bool
            flag to activate multiprocessing of channel data conversion

        channelList : list of str, optional
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory. Can be used to read big files

        convertAfterRead : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false, all data from channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .getChannelData()

        compression : bool, optional
            falg to activate data compression with blosc
        """

        self.multiProc = multiProc

        if self.fileName is None and info is not None:
            self.fileName = info.fileName
        elif fileName is not None and self.fileName is None:
            self.fileName = fileName

        minimal = 2  # always read minimum info

        # set is more efficient for large number of channels (n^2 vs n*log(n)):
        if channelList is not None:
            channelSetFile = set(channelList)  # make sure it is a set
            minimal = 1  # reads at least CN to populate ChannelNamesByDG
        else:
            channelSetFile = None

        # Read information block from file
        if info is None:
            if self.info is None:
                info = info4(self.fileName, None, minimal=minimal)
            else:
                info = self.info

        if info.fid is None or info.fid.closed:
            info.fid = open(self.fileName, 'rb')

        # reads metadata
        if not self._noDataLoading:
            fileDateTime = gmtime(info['HD']['hd_start_time_ns'] / 1000000000)
            ddate = strftime('%Y-%m-%d', fileDateTime)
            ttime = strftime('%H:%M:%S', fileDateTime)

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
        if self._noDataLoading and channelList is not None:
            data_groups = [self[channel][idField][0] for channel in channelList]

        for dataGroup in data_groups:
            channelSet = channelSetFile
            if not info['DG'][dataGroup]['dg_data'] == 0 and \
                    (channelSet is None or
                     len(channelSet & info['ChannelNamesByDG'][dataGroup]) > 0):  # there is data block and channel in
                if minimal > 1 and not self._noDataLoading:  # load CG, CN and CC block info
                    info.readCGBlock(info.fid, dataGroup, channelSet, minimal=minimal)
                if info['CG'][dataGroup][0]['cg_cycle_count']:  # data exists
                    # Pointer to data block
                    pointerToData = info['DG'][dataGroup]['dg_data']

                    if 'dataClass' not in info['DG'][dataGroup]:
                        buf = DATA(info.fid, pointerToData)
                        for channelGroup in info['CG'][dataGroup]:
                            temp = record(dataGroup, channelGroup)  # create record class
                            temp.loadInfo(info)  # load all info related to record
                            buf.addRecord(temp)  # adds record to DATA
                            recordID = info['CG'][dataGroup][channelGroup]['cg_record_id']
                            if temp.master['name'] is not None \
                                    and buf[recordID]['record'].channelNames:
                                if channelSet is not None and not self._noDataLoading\
                                        and temp.master['name'] not in channelSet:
                                    channelSet.add(temp.master['name'])  # adds master channel in channelSet if missing
                            if channelSet is not None and buf[recordID]['record'].CANOpen:
                                # adds CANOpen channels if existing in not empty channelSet
                                if buf[recordID]['record'].CANOpen == 'time':
                                    channelSet.update(('ms', 'days'))
                                elif buf[recordID]['record'].CANOpen == 'date':
                                    channelSet.update(('ms', 'minute', 'hour', 'day', 'month', 'year'))
                        if self._noDataLoading:
                            self.info['DG'][dataGroup]['dataClass'] = buf
                    else:
                        buf = self.info['DG'][dataGroup]['dataClass']

                    # reads raw data from data block with DATA and DATABlock classes
                    buf.read(channelSet, info, self.fileName)

                    channel_groups = buf
                    if self._noDataLoading and channelList is not None:
                        channel_groups = [info['CG'][dataGroup][self[channel][idField][1]]['cg_record_id']
                                          for channel in channelList]

                    # processing data from buf then transfer to self
                    for recordID in channel_groups:  # for each channel group in data block
                        if 'record' in buf[recordID]:
                            master_channel = buf[recordID]['record'].master['name']

                            channels = buf[recordID]['record']
                            if self._noDataLoading and channelList is not None:
                                channels = [channels[self[channel][idField][2]] for channel in channelList]
                            for chan in channels:  # for each channel class
                                if channelSet is None or chan.name in channelSet:
                                    if not chan.type == 'Inv':  # normal channel
                                        if chan.channelType(info) not in (3, 6):  # not virtual channel
                                            # in case record is used for several channels
                                            if channelSet is None and not buf[recordID]['record'].hiddenBytes \
                                                    and buf[recordID]['record'].byte_aligned:
                                                recordName = buf[recordID]['record'].recordToChannelMatching[chan.name]
                                            else:
                                                recordName = chan.name
                                            if 'data' in buf[recordID] and \
                                                    buf[recordID]['data'] is not None:  # no data in channel group
                                                temp = buf[recordID]['data'][recordName]  # extract channel vector
                                            else:
                                                temp = None
                                        else:  # virtual channel
                                            temp = arange(buf[recordID]['record'].numberOfRecords)

                                        # Process concatenated bits inside uint8
                                        bitCount = chan.bitCount(info)
                                        if buf[recordID]['record'].byte_aligned \
                                                and not buf[recordID]['record'].hiddenBytes and \
                                                channelSet is None and\
                                                0 < bitCount < 64 and bitCount not in (8, 16, 32) \
                                                and temp is not None\
                                                and temp.dtype.kind not in ('S', 'U'):
                                            # if channel data do not use complete bytes and Ctypes
                                            signal_data_type = chan.signalDataType(info)
                                            if signal_data_type in (0, 1, 2, 3):  # integers
                                                bitOffset = chan.bitOffset(info)
                                                if bitOffset > 0:
                                                    temp = right_shift(temp, bitOffset)
                                                mask = int(pow(2, bitCount) - 1)  # masks isBitUnit8
                                                temp = bitwise_and(temp, mask)
                                                if signal_data_type in (2, 3):
                                                    # signed integer, moving bit sign of two's complement
                                                    signBitMask = (1 << (bitCount - 1))
                                                    signExtend = ((1 << (temp.itemsize * 8 - bitCount)) - 1) << bitCount
                                                    signBit = bitwise_and(temp, signBitMask)
                                                    for number, sign in enumerate(signBit):
                                                        # negative value, sign extend
                                                        if sign:
                                                            temp[number] |= signExtend
                                            else:  # should not happen
                                                warn('bit count and offset not applied to correct '
                                                     'data type {}'.format(chan.name))
                                        else:  # data using full bytes
                                            pass

                                        if temp is not None:  # channel contains data
                                            # string data decoding
                                            if temp.dtype.kind == 'S':
                                                signalDataType = chan.signalDataType(info)
                                                if signalDataType == 6:  # string ISO-8859-1 Latin
                                                    encoding = 'latin-1'
                                                elif signalDataType == 7:  # UTF-8
                                                    encoding = 'UTF-8'
                                                elif signalDataType == 8:
                                                    encoding = 'UTF-16LE'
                                                elif signalDataType == 9:  # UTF-16 big endian
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
                                            self.add_channel(dataGroup, chan.name, temp,
                                                             master_channel,
                                                             master_type=chan.channelSyncType(info),
                                                             unit=chan.unit(info),
                                                             description=chan.desc(info),
                                                             conversion=chan.conversion(info),
                                                             info=chan.CNBlock(info),
                                                             compression=compression)
                                            if chan.channelType(info) == 4:  # sync channel
                                                # attach stream to be synchronised
                                                self.setChannelAttachment(chan.name, chan.attachment(info.fid, info))
                                            if chan.has_invalid_bit(info):  # has invalid bit
                                                self.setInvalidBit(chan.name, chan.invalid_bit(info))
                                                self.setInvalidChannel(chan.name, 'invalid_bytes{}'.format(dataGroup))
                                    else:
                                        # invalid bytes channel
                                        data = buf[recordID]['data'].__getattribute__(chan.name)
                                        data = frombuffer(data.tobytes(), dtype='u1').reshape(len(data),
                                                                                              data.dtype.itemsize)
                                        self.add_channel(dataGroup, chan.name,
                                                         data,
                                                         master_channel,
                                                         master_type=0,
                                                         unit='',
                                                         description='',
                                                         info=None,
                                                         compression=compression)
                            buf[recordID].pop('data', None)
                    del buf
                if minimal > 1:
                    # clean CN, CC and CG info to free memory
                    info.cleanDGinfo(dataGroup)
        info.fid.close()  # close file

        if convertAfterRead and not compression:
            self._noDataLoading = False
            self._convertAllChannel4()
        # print( 'Finished in ' + str( time.clock() - inttime ) , file=stderr)

    def _getChannelData4(self, channelName, raw_data=False):
        """Returns channel numpy array

        Parameters
        ----------------
        channelName : str
            channel name
        raw_data: bool
            flag to return non converted data

        Returns:
        -----------
        numpy array
            converted, if not already done, data corresponding to channel name

        Notes
        ------
        This method is the safest to get channel data as numpy array from 'data' dict key might contain raw data
        """
        if channelName in self:
            vect = self.getChannel(channelName)[dataField]
            if vect is None:  # noDataLoading reading argument flag activated
                if self.info.fid is None or (self.info.fid is not None and self.info.fid.closed):
                    (self.info.fid, self.info.fileName, self.info.zipfile) = _open_MDF(self.fileName)
                self.read4(fileName=None, info=None, channelList=[channelName], convertAfterRead=False)
            if not raw_data:
                return self._convertChannelData4(self.getChannel(channelName),
                        channelName, self.convert_tables)[channelName]
            else:
                return self.getChannel(channelName)[dataField]
        else:
            return None

    def _convertChannelData4(self, channel, channelName, convert_tables, multiProc=False, Q=None):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channel : channel class
            channel4 object
        channelName : str
            name of channel
        convert_tables : bool
            activates computation intensive loops for conversion with tables. Default is False
        multiProc : bool, default False
            flag to put data in multiprocess queue
        Q : Queue class, default None
            Queue used for multiprocessing

        Returns
        -----------
        dict
            returns dict with channelName key containing numpy array converted to physical values according to conversion type
        """
        if channel[dataField] is None:
            vect = channel[dataField]
        else:
            if isinstance(channel[dataField], compressed_data):
                vect = channel[dataField].decompression()  # uncompressed blosc data
            else:
                vect = channel[dataField][:]  # to have bcolz uncompressed data
        if conversionField in channel and channel[conversionField]['type']:  # there is conversion property
            text_type = vect.dtype.kind in ['S', 'U', 'V']  # channel of string or not ?
            conversion_type = channel[conversionField]['type']
            conversion_parameter = channel[conversionField]['parameters']
            if conversion_type == 1 and not text_type:
                vect = linearConv(vect, conversion_parameter['cc_val'])
            elif conversion_type == 2 and not text_type:
                vect = rationalConv(vect, conversion_parameter['cc_val'])
            elif conversion_type == 3 and not text_type:
                vect = formulaConv(vect, conversion_parameter['cc_ref']['Comment'])
            elif conversion_type == 4 and not text_type:
                vect = valueToValueTableWInterpConv(vect, conversion_parameter['cc_val'])
            elif conversion_type == 5 and not text_type:
                vect = valueToValueTableWOInterpConv(vect, conversion_parameter['cc_val'])
            elif conversion_type == 6 and not text_type and convert_tables:
                vect = valueRangeToValueTableConv(vect, conversion_parameter['cc_val'])
            elif conversion_type == 7 and not text_type and convert_tables:
                vect = valueToTextConv(vect, conversion_parameter['cc_val'], conversion_parameter['cc_ref'])
            elif conversion_type == 8 and not text_type and convert_tables:
                vect = valueRangeToTextConv(vect, conversion_parameter['cc_val'], conversion_parameter['cc_ref'])
            elif conversion_type == 9 and text_type and convert_tables:
                vect = textToValueConv(vect, conversion_parameter['cc_val'], conversion_parameter['cc_ref'])
            elif conversion_type == 10 and text_type and convert_tables:
                vect = textToTextConv(vect, conversion_parameter['cc_ref'])
        L = {}
        L[channelName] = vect
        if multiProc:
            Q.put(L)
        else:
            return L

    def _convertChannel4(self, channelName):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channelName : str
            Name of channel
        """
        self.setChannelData(channelName, self._getChannelData4(channelName))
        self.remove_channel_conversion(channelName)

    def _convertAllChannel4(self):
        """Converts all channels from raw data to converted data according to CCBlock information
        Converted data will take more memory.
        """

        if self._noDataLoading:  # no data loaded, load everything
            self.read4(self.fileName, convertAfterRead=True)
        else:
            if self.multiProc is False:
                [self._convertChannel4(channelName) for channelName in self]
            else:  # multiprocessing
                proc = []
                Q = Queue()
                L = {}
                for channelName in self:
                    channel = self.getChannel(channelName)
                    if 'conversion' in channel:
                        conversion = self.getChannelConversion(channelName)
                        if conversion['type'] in (1, 2):  # more time in multi proc
                            self._convertChannel4(channelName)
                        else:
                            proc.append(Process(target=self._convertChannelData4,
                                    args=(channel, channelName, self.convert_tables, True, Q)))
                            proc[-1].start()
                for p in proc:
                    L.update(Q.get())  # concatenate results of processes in dict
                for p in proc:
                    p.join()
                del Q  # free memory
                for channelName in self:
                    if channelName in L:
                        self.setChannelData(channelName, L[channelName])
                        self.remove_channel_conversion(channelName)

    def write4(self, fileName=None, compression=False):
        """Writes simple mdf 4.1 file

        Parameters
        ----------------
        fileName : str, optional
            Name of file
            If file name is not input, written file name will be the one read with appended '_new' string before extension
        compression : bool
            flag to store data compressed

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        # Starts first to write ID and header
        if fileName is None:
            splitName = splitext(self.fileName)
            if splitName[-1] in ('.mfxz', '.MFXZ'):
                splitName[-1] = '.mfx'  # do not resave in compressed file
            fileName = ''.join([splitName[-2], '_New', splitName[-1]])
        fid = open(fileName, 'wb')  # buffering should automatically be set
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
        blocks['HD_comment'].load(self.file_metadata, 'HD')
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
        blocks['FH_comment'].load(self.file_metadata, 'FH')
        pointer = blocks['FH_comment']['block_start'] + blocks['FH_comment']['block_length']

        # write DG block
        blocks['HD']['DG'] = pointer  # first DG

        # write all files header blocks
        for block in blocks.values():
            block.write(fid)
        if self.masterChannelList:  # some channels exist
            for dataGroup, masterChannel in enumerate(self.masterChannelList):
                # writes dataGroup Block
                DG = DGBlock()
                DG['block_start'] = pointer
                pointer = DG['block_start'] + 64
                DG['CG'] = pointer  # First CG link

                blocks = OrderedDict()  # initialise blocks for this datagroup
                # write CGBlock
                blocks['CG'] = CGBlock()
                blocks['CG']['block_start'] = pointer
                pointer = blocks['CG']['block_start'] + 104
                blocks['CG']['CN'] = pointer  # First CN link

                master_channel_data = self._getChannelData4(masterChannel)
                if master_channel_data is not None:
                    cg_cycle_count = len(master_channel_data)
                elif self._getChannelData4(self.masterChannelList[masterChannel][0]).shape[0] == 1:  # classification
                    cg_cycle_count = 1
                elif master_channel_data not in self.masterChannelList[masterChannel]:
                    cg_cycle_count = len(self._getChannelData4(self.masterChannelList[masterChannel][0]))
                    warn('no master channel in datagroup {}'.format(dataGroup))
                else:
                    # no data in datagroup, skip
                    warn('no data in datagroup {0} with master channel {1}'.format(dataGroup, masterChannel))
                    continue
                blocks['CG']['cg_cycle_count'] = cg_cycle_count

                # write channels
                record_byte_offset = 0
                CN_flag = 0
                number_of_channel = 0
                nRecords = 0
                dataList = ()
                last_channel = 0
                for nchannel, channel in enumerate(self.masterChannelList[masterChannel]):
                    data = self.getChannelData(channel)
                    # no interest to write invalid bytes as channel, should be processed if needed before writing
                    if channel.find('invalid_bytes') == -1 and data is not None and len(data) > 0:
                        last_channel = nchannel
                        data_ndim = data.ndim - 1
                        if not data_ndim:
                            dataList = dataList + (data, )
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
                                dataList = dataList + (flattened[:, i],)
                        number_of_channel += 1
                        blocks[nchannel] = CNBlock()
                        if issubdtype(data.dtype, numpy_number):  # is numeric
                            blocks[nchannel]['cn_val_range_min'] = npmin(data)
                            blocks[nchannel]['cn_val_range_max'] = npmax(data)
                            blocks[nchannel]['cn_flags'] = 8  # only Bit 3: Limit range valid flag
                        else:
                            blocks[nchannel]['cn_val_range_min'] = 0
                            blocks[nchannel]['cn_val_range_max'] = 0
                            blocks[nchannel]['cn_flags'] = 0
                        if masterChannel is not channel:
                            blocks[nchannel]['cn_type'] = 0
                            blocks[nchannel]['cn_sync_type'] = 0
                        else:
                            blocks[nchannel]['cn_type'] = 2  # master channel
                            blocks[nchannel]['cn_sync_type'] = 1  # default is time channel
                            nRecords = len(data)

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
                        blocks[nchannel]['cn_data_type'] = data_type
                        blocks[nchannel]['cn_bit_offset'] = 0  # always byte aligned
                        blocks[nchannel]['cn_byte_offset'] = record_byte_offset
                        byte_count = data.dtype.itemsize
                        record_byte_offset += byte_count
                        blocks[nchannel]['cn_bit_count'] = byte_count * 8
                        blocks[nchannel]['block_start'] = pointer
                        pointer = blocks[nchannel]['block_start'] + 160

                        # arrays handling
                        if not data_ndim:
                            blocks[nchannel]['Composition'] = 0
                        else:
                            blocks[nchannel]['Composition'] = pointer  # pointer to CABlock
                            # creates CABlock
                            ca = ''.join([channel, '_CA'])
                            blocks[ca] = CABlock()
                            blocks[ca]['block_start'] = pointer
                            blocks[ca]['ndim'] = data_ndim
                            blocks[ca]['ndim_size'] = data_dim_size
                            blocks[ca].load(data.itemsize)
                            pointer = blocks[ca]['block_start'] + blocks[ca]['block_length']

                        if CN_flag:
                            # Next DG
                            blocks[nchannel-1]['CN'] = blocks[nchannel]['block_start']
                            blocks[nchannel]['CN'] = 0  # initialise 'CN' key
                        else:
                            CN_flag = blocks[nchannel]['block_start']
                            blocks[nchannel]['CN'] = 0  # creates first CN link, null for the moment

                        # write channel name
                        blocks[nchannel]['TX'] = pointer
                        blocks[channel] = CommentBlock()
                        blocks[channel]['block_start'] = pointer
                        blocks[channel].load(channel, 'TX')
                        pointer = blocks[channel]['block_start'] + blocks[channel]['block_length']

                        # write channel unit
                        unit = self.getChannelUnit(channel)
                        if unit is not None and len(unit) > 0:
                            blocks[nchannel]['Unit'] = pointer
                            unit_name = ''.join([channel, '_U'])
                            blocks[unit_name] = CommentBlock()
                            blocks[unit_name]['block_start'] = pointer
                            blocks[unit_name].load(unit, 'TX')
                            pointer = blocks[unit_name]['block_start'] + blocks[unit_name]['block_length']
                        else:
                            blocks[nchannel]['Unit'] = 0

                        # write channel description
                        desc = self.getChannelDesc(channel)
                        if desc is not None and len(desc) > 0:
                            blocks[nchannel]['Comment'] = pointer
                            desc_name = ''.join([channel, '_C'])
                            blocks[desc_name] = CommentBlock()
                            blocks[desc_name]['block_start'] = pointer
                            blocks[desc_name].load(desc, 'TX')
                            pointer = blocks[desc_name]['block_start'] + blocks[desc_name]['block_length']
                        else:
                            blocks[nchannel]['Comment'] = 0

                if nRecords == 0 and masterChannel is not self.masterChannelList[masterChannel]:
                    # No master channel in channel group
                    nRecords = len(dataList[0])

                if last_channel in blocks:
                    blocks[last_channel]['CN'] = 0  # last CN link is null
                # writes size of record in CG
                blocks['CG']['cg_data_bytes'] = record_byte_offset

                # data pointer in datagroup
                DG['data'] = pointer
                if compression:
                    data = HLBlock()
                    data.load(record_byte_offset, nRecords, pointer)
                    DG_start_position = fid.tell()
                    DG['DG'] = 0
                else:
                    data = DTBlock()
                    data.load(record_byte_offset, nRecords, pointer)
                    DG['DG'] = data['end_position']

                DG.write(fid)  # write DG block

                # writes all blocks (CG, CN, TX for unit and description) before writing data block
                for block in blocks.values():
                    block.write(fid)

                # data block writing
                pointer = data.write(fid, fromarrays(dataList).tobytes(order='F'))
                if compression:
                    # next DG position is not predictable due to DZ Blocks unknown length
                    fid.seek(DG_start_position + 24)
                    fid.write(pack('<Q', pointer))
                    fid.seek(pointer)

            fid.seek(DG['block_start'] + 24)
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
            invalid_bit_pos = self.getInvalidBit(channel_name)
            if isinstance(invalid_bit_pos, int):  # invalid bit existing
                mask = self._getChannelData4(self.getInvalidChannel(channel_name))
                data = self._getChannelData4(channel_name)
                data = data.view(MaskedArray)
                invalid_byte = invalid_bit_pos >> 3
                data.mask = bitwise_and(mask[:, invalid_byte], invalid_bit_pos & 0x07)
                self.setChannelData(channel_name, data)
                self._remove_channel_field(channel_name, invalidPosField)
                self._remove_channel_field(channel_name, invalidChannel)
        except KeyError:
            pass
            # warn('no invalid data found for channel ')


def linearConv(vect, cc_val):
    """ apply linear conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    P1 = cc_val[0]
    P2 = cc_val[1]
    if P2 == 1.0 and P1 in (0.0, -0.0):
        return vect  # keeps dtype probably more compact than float64
    else:
        return vect * P2 + P1


def rationalConv(vect, cc_val):
    """ apply rational conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    P1 = cc_val[0]
    P2 = cc_val[1]
    P3 = cc_val[2]
    P4 = cc_val[3]
    P5 = cc_val[4]
    P6 = cc_val[5]
    return (P1 * vect * vect + P2 * vect + P3) / (P4 * vect * vect + P5 * vect + P6)


def formulaConv(vect, formula):
    """ apply formula conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

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
    return expr(vect)


def valueToValueTableWOInterpConv(vect, cc_val):
    """ apply value to value table without interpolation conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = 2 * int(len(cc_val) / 2)
    intVal = [cc_val[i] for i in range(0, val_count, 2)]
    physVal = [cc_val[i] for i in range(1, val_count, 2)]
    if all(diff(intVal) > 0):
        try:
            from scipy import interpolate
        except:
            raise ImportError('Please install scipy to convert channel')
        f = interpolate.interp1d(intVal, physVal, kind='nearest', bounds_error=False)  # nearest
        return f(vect)  # fill with Nan out of bounds while should be bounds
    else:
        warn('X values for interpolation of channel are not increasing')


def valueToValueTableWInterpConv(vect, cc_val):
    """ apply value to value table with interpolation conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    val_count = 2 * int(len(cc_val) / 2)
    intVal = [cc_val[i] for i in range(0, val_count, 2)]
    physVal = [cc_val[i] for i in range(1, val_count, 2)]
    if all(diff(intVal) > 0):
        return interp(vect, intVal, physVal)  # with interpolation
    else:
        warn('X values for interpolation of channel are not increasing')


def valueRangeToValueTableConv(vect, cc_val):
    """ apply value range to value table conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
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
    for Lindex in range(len(vect)):
        key_index = 0  # default index if not found
        for i in range(val_count):
            if key_min[i] < vect[Lindex] < key_max[i]:
                key_index = i
                break
        vect[Lindex] = value[key_index]
    return vect


def valueToTextConv(vect, cc_val, cc_ref):
    """ apply value to text conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    maxlen = max([len(ref) for ref in cc_ref])
    temp = empty(len(vect), dtype='U{}'.format(maxlen))  # initialize empty array with max length
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
                cc_ref[ref] = lambdify(X, cc_ref[ref]['cc_ref']['Comment']
                                       , modules='numpy', dummify=False)
            elif cc_ref[ref]['cc_type'] == 1: # linear conversion
                cc_ref[ref] = lambdify(X, '{0}* X + {1}'.format(cc_ref[ref]['cc_val'][1], cc_ref[ref]['cc_val'][0]) 
                                       , modules='numpy', dummify=False)
            else:
                warn('To implement missing conversion, please ask')
        elif (PythonVersion > 3 and not isinstance(cc_ref[ref], str)) or \
                (PythonVersion < 3 and not isinstance(cc_ref[ref], unicode)):  # identity, non conversion
            cc_ref[ref] = lambdify(X, 'X'
                                    , modules='numpy', dummify=False)
    key_index = where(vect[0] == cc_val)[0]  # look up for first value in vect
    if not len(key_index) == 0:  # value corresponding in cc_val
        temp[0] = cc_ref[key_index[0]]
    else:  # default value
        if callable(cc_ref[-1]):
            temp[0] = cc_ref[-1](vect[0])
        else:
            temp[0] = cc_ref[-1]
    for Lindex in range(1, len(vect)):
        if vect[Lindex] == vect[Lindex - 1]:  # same value as before, no need to look further
            temp[Lindex] = temp[Lindex - 1]
        else:  # value changed from previous step
            key_index = where(vect[Lindex] == cc_val)[0]
            if not len(key_index) == 0:  # found match
                temp[Lindex] = cc_ref[key_index[0]]
            else:  # default
                if callable(cc_ref[-1]):
                    temp[Lindex] = cc_ref[-1](vect[Lindex])
                else:
                    temp[Lindex] = cc_ref[-1]
    return asarray(temp)


def valueRangeToTextConv(vect, cc_val, cc_ref):
    """ apply value range to text conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
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
                cc_ref[ref] = lambdify(X, cc_ref[ref]['cc_ref']['Comment']
                                       , modules='numpy', dummify=False)
            elif cc_ref[ref]['cc_type'] == 1:  # linear conversion
                cc_ref[ref] = lambdify(X, '{0} * X + {1}'.format(cc_val[1], cc_val[0])
                                       , modules='numpy', dummify=False)
            else:  # identity, no conversion
                cc_ref[ref] = lambdify(X, '1 * X'
                                       , modules='numpy', dummify=False)
            # Otherwise a string
    # look up in range keys
    temp = []
    for value in vect:
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


def textToValueConv(vect, cc_val, cc_ref):
    """ apply text to value conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_val : cc_val from mdfinfo4.info4 conversion block ('CCBlock') dict
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    ref_count = len(cc_ref)
    temp = []
    for Lindex in range(len(vect)):
        key_index = ref_count  # default index if not found
        for i in range(ref_count):
            if vect[Lindex] == cc_ref[i]:
                key_index = i
                break
        temp.append(cc_val[key_index])
    return asarray(temp)


def textToTextConv(vect, cc_ref):
    """ apply text to text conversion to data

    Parameters
    ----------------
    vect : numpy 1D array
        raw data to be converted to physical value
    cc_ref : cc_ref from mdfinfo4.info4 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    ref_count = len(cc_ref) - 2
    for Lindex in range(len(vect)):
        key_index = ref_count + 1  # default index if not found
        for i in range(0, ref_count, 2):
            if vect[Lindex] == cc_ref[i]:
                key_index = i
                break
        vect[Lindex] = cc_ref[key_index]
    return vect
