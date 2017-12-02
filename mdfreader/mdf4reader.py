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
from __future__ import print_function
from struct import Struct
from struct import unpack as structunpack
from math import pow
from io import open  # for python 3 and 2 consistency
from time import gmtime, strftime
from multiprocessing import Queue, Process
from sys import version_info, exc_info, byteorder, stderr, path
from os.path import dirname, abspath
from collections import defaultdict
from numpy.core.records import fromstring, fromfile, fromarrays
from numpy import array, recarray, append, asarray, empty, where
from numpy import arange, right_shift, bitwise_and, all, diff, interp
from numpy import max as npmax, min as npmin
from numpy.lib.recfunctions import rec_append_fields, rename_fields

_root = dirname(abspath(__file__))
path.append(_root)
from mdfinfo4 import info4, ATBlock, IDBlock, HDBlock, DGBlock, \
    CGBlock, CNBlock, FHBlock, CommentBlock, _loadHeader, DLBlock, \
    DZBlock, HLBlock, CCBlock, _writePointer, _writeHeader
from mdf import mdf_skeleton, _open_MDF, _bits_to_bytes, \
    dataField, conversionField, compressed_data
from channel import channel4

PythonVersion = version_info
PythonVersion = PythonVersion[0]

chunk_size_reading = 100000000  # reads by chunk of 100Mb, can be tuned for best performance


def DATABlock(record, info, parent_block, channelSet=None, sortedFlag=True, vlsd=False):
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
    if parent_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal data block
        if sortedFlag:
            if channelSet is None and not record.hiddenBytes and\
                    record.byte_aligned:  # No channel list and length of records corresponds to C datatypes
                # for debugging purpose
                # print(record.numberOfRecords, record.numpyDataRecordFormat, record.dataRecordName)
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=record.numberOfRecords)
            else:  # record is not byte aligned or channelSet not None
                return record.readBitarray(parent_block['data'], info, channelSet)
        else:  # unsorted reading
            # not so much tested, missing example file
            return readUnsorted(record, info, parent_block, channelSet, sortedFlag)

    elif parent_block['id'] in ('##SD', b'##SD'):
        return read_sdblock(record[record.VLSD[0]].signalDataType(info),
                            parent_block['data'], parent_block['length'] - 24)

    elif parent_block['id'] in ('##DZ', b'##DZ'):  # zipped data block
        # uncompress data
        parent_block['data'] = decompress_datablock(parent_block['data'], parent_block['dz_zip_type'],
                parent_block['dz_zip_parameter'], parent_block['dz_org_data_length'])
        if vlsd:  # VLSD channel
            return read_sdblock(record[record.VLSD[0]].signalDataType(info),
                                parent_block['data'], parent_block['dz_org_data_length'])
        if channelSet is None and sortedFlag:  # reads all blocks if sorted block and no channelSet defined
            if record.byte_aligned and not record.hiddenBytes:
                record.numberOfRecords = parent_block['dz_org_data_length'] // record.CGrecordLength
                return fromstring(parent_block['data'], dtype={'names': record.dataRecordName,
                                                               'formats': record.numpyDataRecordFormat},
                                  shape=record.numberOfRecords)
            else:
                return record.readBitarray(parent_block['data'], info, channelSet)
        elif channelSet is not None and sortedFlag:  # sorted data but channel list requested
            return record.readBitarray(parent_block['data'], info, channelSet)
        else:  # unsorted reading
            return readUnsorted(record, info, parent_block, channelSet, sortedFlag)


def readUnsorted(record, info, parent_block, channelSet=None, sortedFlag=True):
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


def decompress_datablock(block, zip_type, zip_parameter, org_data_length):
    """ decompress datablock.

    Parameters
    --------------
    block : bytes
        raw data compressed
    zip_type : int
        0 for non transposed, 1 for transposed data
    zip_parameter : int
        first dimension of matrix to be transposed
    org_data_length : int
        uncompressed data length

    Returns
    ---------
    uncompressed raw data
    """
    from zlib import decompress
    if PythonVersion > 3:  # decompress data
        block = decompress(block)
    else:
        block = decompress(bytes(block))  # decompress data
    if zip_type == 1:  # data bytes transposed
        M = org_data_length // zip_parameter
        temp = array(memoryview(block[:M * zip_parameter]))
        tail = array(memoryview(block[M * zip_parameter:]))
        temp = temp.reshape(zip_parameter, M).T.ravel()
        if len(tail) > 0:
            temp = append(temp, tail)
        block = temp.tostring()
    return block


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

    def read(self, channelSet, info, zip=None):
        """Reads data block

        Parameters
        ----------------
        channelSet : set of str
            set of channel names
        info : info object
            contains blocks structures
        zip : bool, optional
            flag to track if data block is compressed
        """
        if len(self) == 1:  # sorted dataGroup
            recordID = list(self.keys())[0]
            record = self[recordID]['record']
            self[recordID]['data'] = self.load(record, info, zip=None, nameList=channelSet, sortedFlag=True)
            for cn in record.VLSD:  # VLSD channels
                if channelSet is None or record[cn].name in channelSet:
                    temp = DATA(self.fid, record[cn].data(info))  # all channels
                    temp = temp.load(record, info, zip=None, nameList=channelSet, sortedFlag=True, vlsd=True)
                    self[recordID]['data'] = rename_fields(self[recordID]['data'],
                                                           {record[cn].name: '{}_offset'.format(record[cn].name)})
                    self[recordID]['data'] = rec_append_fields(self[recordID]['data'],
                                                               record[cn].name, temp)
        else:  # unsorted DataGroup
            self.type = 'unsorted'
            data = self.load(self, info, zip=None, nameList=channelSet, sortedFlag=False)
            for recordID in self:
                self[recordID]['data'] = {}
                for channel in self[recordID]['record']:
                    self[recordID]['data'][channel.name] = data[channel.name]

    def load(self, record, info, zip=None, nameList=None, sortedFlag=True, vlsd=False):
        """Reads data block from record definition

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        zip : bool, optional
            flag to track if data block is compressed
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
            temps.update(DLBlock(self.fid, temps['link_count']))
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
                buf = bytearray()
                for DL in temps['dl_data']:
                    for pointer in temps['dl_data'][DL]:
                        # read fist data blocks linked by DLBlock to identify data block type
                        data_block = defaultdict()
                        data_block.update(_loadHeader(self.fid, pointer))
                        if data_block['id'] in ('##DT', '##RD', b'##DT',
                                                b'##RD', '##SD', b'##SD'):
                            buf.extend(self.fid.read(
                                       data_block['length'] - 24))
                        elif data_block['id'] in ('##DZ', b'##DZ'):
                            data_block.update(DZBlock(self.fid))
                            data_block['data'] = \
                                decompress_datablock(
                                    self.fid.read(data_block['dz_data_length']),
                                    data_block['dz_zip_type'],
                                    data_block['dz_zip_parameter'],
                                    data_block['dz_org_data_length'])
                            buf.extend(data_block['data'])
                            data_block['id'] = '##DT'  # do not uncompress in DATABlock function
                data_block['data'] = buf
                temps['data'] = DATABlock(record, info, parent_block=data_block,
                                          channelSet=nameList, sortedFlag=sortedFlag, vlsd=vlsd)
            else:  # empty datalist
                temps['data'] = None
        elif temps['id'] in ('##HL', b'##HL'):  # header list block for DZBlock
            # link section
            temps.update(HLBlock(self.fid))
            self.pointerTodata = temps['hl_dl_first']
            temps['data'] = self.load(record, info, zip=temps['hl_zip_type'], nameList=nameList, sortedFlag=sortedFlag)
        elif temps['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal sorted data block, direct read
            temps['data'] = record.readSortedRecord(self.fid, self.pointerTodata, info, channelSet=nameList)
        elif temps['id'] in ('##SD', b'##SD'):  # VLSD
            temps['data'] = self.fid.read(temps['length'] - 24)
            temps['data'] = DATABlock(record, info, parent_block=temps, channelSet=nameList, sortedFlag=sortedFlag)
        elif temps['id'] in ('##DZ', b'##DZ'):  # zipped data block
            temps.update(DZBlock(self.fid))
            temps['data'] = self.fid.read(temps['dz_data_length'])
            temps['data'] = DATABlock(record, info, parent_block=temps,
                                      channelSet=nameList, sortedFlag=sortedFlag, vlsd=vlsd)
        else:
            raise Exception('unknown data block')
        return temps['data']

    def readRecord(self, recordID, info, buf, channelSet=None):
        """ read record from a buffer

        Parameters
        ----------------
        recordID : int
            record identifier
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
    readBitarray(bita, info, channelSet=None)
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
        output = 'Record class content print\nTotal number of channels : {}\n'.format(len(self))
        for chan in self:
            output.join([chan.name, '  Type ', chan.type, '\n',
                        'DG {} '.format(chan.dataGroup),
                        'CG {} '.format(chan.channelGroup),
                        'CN {} '.format(chan.channelNumber),
                        'VLSD {} '.format(chan.VLSD_CG_Flag)])
        output.append(''.join(['CG block record bytes length : {}\n'.format(self.CGrecordLength),
                               'Datagroup number : {}\n'.format(self.dataGroup),
                               'Byte aligned : {}\n'.format(self.byte_aligned),
                               'Hidden bytes : {}\n'.format(self.hiddenBytes)]))
        if self.master['name'] is not None:
            output.append(''.join(['Master channel : ', self.master['name'], '\n']))
        output.append('Numpy records format : \n')
        for record in self.numpyDataRecordFormat:
            output.append('{}\n'.format(record))
        output.append('VLSD_CG {}'.format(self.VLSD_CG))
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
        for channelNumber in range(len(info['CN'][self.dataGroup][self.channelGroup])):
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
                                Channel_posBitBeg < 8 * (prev_chan_byteOffset + prev_chan_nBytes) and \
                                Channel_posBitEnd > 8 * (prev_chan_byteOffset + prev_chan_nBytes):  # not byte aligned
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

    def readSortedRecord(self, fid, pointer, info, channelSet=None):
        """ reads record, only one channel group per datagroup

        Parameters
        ----------------
        fid : float
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
        if channelSet is None:  # reads all, quickest but memory consuming
            if self.byte_aligned and not self.hiddenBytes:
                # print(self) # for debugging purpose
                return fromfile(fid, dtype={'names': self.dataRecordName,
                                            'formats': self.numpyDataRecordFormat},
                                shape=self.numberOfRecords)
            else:
                return self.readBitarray(fid.read(self.CGrecordLength * self.numberOfRecords), info, channelSet)
        else:  # reads only some channels from a sorted data block
            if len(channelSet & self.channelNames) > 0:  # are channelSet in this dataGroup
                try:
                    return self.readBitarray(fid.read(self.CGrecordLength * self.numberOfRecords), channelSet)
                except:  # still memory efficient but takes time
                    print('Unexpected error:', exc_info(), file=stderr)
                    print('dataRead crashed, back to python data reading', file=stderr)
                    # check if master channel is in the list
                    if not self.master['name'] in channelSet:
                        channelSet.add(self.master['name'])  # adds master channel
                    rec = {}
                    recChan = []
                    dataRecordName = []
                    numpyDataRecordFormat = []
                    for channel in channelSet:  # initialise data structure
                        rec[channel] = 0
                    for channel in self:  # list of channel classes from channelSet
                        if channel.name in channelSet:
                            recChan.append(channel)
                            dataRecordName.append(channel.name)
                            numpyDataRecordFormat.append(channel.dataFormat(info))
                    rec = recarray((self.numberOfRecords, ), dtype={'names': dataRecordName,
                                                                 'formats': numpyDataRecordFormat})
                    recordLength = self.recordIDsize + self.CGrecordLength
                    for r in range(self.numberOfRecords):  # for each record,
                        buf = fid.read(recordLength)
                        for channel in recChan:
                            rec[channel.name][r] = \
                                channel.CFormat(info).unpack(buf[channel.posByteBeg(info):channel.posByteEnd(info)])[0]
                    return rec.view(recarray)

    def readRecordBuf(self, buf, info, channelSet=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
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
                    Channel.CFormat(info).unpack(buf[Channel.posByteBeg(info):\
                                                     Channel.posByteEnd(info)])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def readBitarray(self, bita, info, channelSet=None):
        """ reads stream of record bytes using bitarray module needed for not byte aligned data
        
        Parameters
        ------------
        bita : stream
            stream of bytes
        channelSet : set of str, optional
            set of channel to read
        
        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        # initialise data structure
        if channelSet is None:
            channelSet = self.channelNames
        try:  # use rather cython compiled code for performance
            formats = []
            names = []
            for chan in range(len(self)):
                if self[chan].name in channelSet:
                    # dataRead will take care of byte order, switch to native
                    formats.append(self[chan].nativedataFormat(info))
                    names.append(self[chan].name)
            if formats:  # at least some channels should be parsed
                buf = recarray(self.numberOfRecords, formats=formats, names=names)
                from dataRead import dataRead
                bytesdata = bytes(bita)
                for chan in range(len(self)):
                    if self[chan].name in channelSet:
                        buf[self[chan].name] = \
                            dataRead(bytesdata, self[chan].bitCount(info),
                                     self[chan].signalDataType(info),
                                     self[chan].nativedataFormat(info),
                                     self.numberOfRecords, self.CGrecordLength,
                                     self[chan].bitOffset(info), self[chan].posByteBeg(info),
                                     self[chan].posByteEnd(info))
                return buf
        except:
            print('Unexpected error:', exc_info(), file=stderr)
            print('dataRead crashed, back to python data reading', file=stderr)
            from bitarray import bitarray
            B = bitarray(endian="little")  # little endian by default
            B.frombytes(bytes(bita))
            def signedInt(temp, extension):
                """ extend bits of signed data managing two's complement
                """
                extension.setall(False)
                extensionInv = bitarray(extension, endian='little')
                extensionInv.setall(True)
                for i in range(self.numberOfRecords):  # extend data of bytes to match numpy requirement
                    signBit = temp[i][-1]
                    if not signBit:  # positive value, extend with 0
                        temp[i].extend(extension)
                    else:  # negative value, extend with 1
                        signBit = temp[i].pop(-1)
                        temp[i].extend(extensionInv)
                        temp[i].append(signBit)
                return temp
            # read data
            record_bit_size = self.CGrecordLength * 8
            for chan in range(len(self)):
                if self[chan].name in channelSet:
                    temp = [B[self[chan].posBitBeg(info) + record_bit_size * i:\
                            self[chan].posBitEnd(info) + record_bit_size * i]\
                                for i in range(self.numberOfRecords)]
                    nbytes = len(temp[0].tobytes())
                    signalDataType = self[chan].signalDataType(info)
                    nBytes = self[chan].nBytes(info)
                    if not nbytes == nBytes and \
                            signalDataType not in (6, 7, 8, 9, 10, 11, 12): # not Ctype byte length
                        byte = bitarray(8 * (nBytes - nbytes), endian='little')
                        byte.setall(False)
                        if signalDataType not in (2, 3):  # not signed integer
                            for i in range(self.numberOfRecords):  # extend data of bytes to match numpy requirement
                                temp[i].extend(byte)
                        else:  # signed integer (two's complement), keep sign bit and extend with bytes
                            temp = signedInt(temp, byte)
                    nTrailBits = nBytes*8 - self[chan].bitCount(info)
                    if signalDataType in (2, 3) and \
                            nbytes == nBytes and \
                            nTrailBits > 0:  # Ctype byte length but signed integer
                        trailBits = bitarray(nTrailBits, endian='little')
                        temp = signedInt(temp, trailBits)
                    if 's' not in self[chan].Format(info):
                        CFormat = self[chan].CFormat(info)
                        if ('>' in self[chan].dataFormat(info) and byteorder == 'little') or \
                           (byteorder == 'big' and '<' in self[chan].dataFormat(info)):
                            temp = [CFormat.unpack(temp[i].tobytes())[0]
                                    for i in range(self.numberOfRecords)]
                            temp = asarray(temp).byteswap().newbyteorder()
                        else: 
                            temp = [CFormat.unpack(temp[i].tobytes())[0]
                                    for i in range(self.numberOfRecords)]
                            temp = asarray(temp)
                    else:
                        temp = [temp[i].tobytes()
                                for i in range(self.numberOfRecords)]
                        temp = asarray(temp)
                    buf[self[chan].name] = temp
            return buf


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

    def read4(self, fileName=None, info=None, multiProc=False, channelList=None, \
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
        fileDateTime = gmtime(info['HD']['hd_start_time_ns'] / 1000000000)
        ddate = strftime('%Y-%m-%d', fileDateTime)
        ttime = strftime('%H:%M:%S', fileDateTime)
        def returnField(obj, field):
            if field in obj:
                return obj[field]
            else:
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

        for dataGroup in info['DG']:
            channelSet = channelSetFile
            if not info['DG'][dataGroup]['dg_data'] == 0 and \
                    (channelSet is None or
                     len(channelSet & info['ChannelNamesByDG'][dataGroup]) > 0):  # there is data block and channel in
                if minimal > 1:  # load CG, CN and CC block info
                    info.readCGBlock(info.fid, dataGroup, channelSet, minimal=minimal)
                if info['CG'][dataGroup][0]['cg_cycle_count']:  # data exists
                    # Pointer to data block
                    pointerToData = info['DG'][dataGroup]['dg_data']
                    buf = DATA(info.fid, pointerToData)

                    for channelGroup in info['CG'][dataGroup]:
                        temp = record(dataGroup, channelGroup)  # create record class
                        temp.loadInfo(info)  # load all info related to record
                        buf.addRecord(temp)  # adds record to DATA
                        recordID = info['CG'][dataGroup][channelGroup]['cg_record_id']
                        if temp.master['name'] is not None \
                                and buf[recordID]['record'].channelNames:
                            if channelSet is not None and temp.master['name'] not in channelSet and not self._noDataLoading:
                                channelSet.add(temp.master['name'])  # adds master channel in channelSet if missing
                        if self._noDataLoading and channelSet is not None and \
                                len(channelSet & buf[recordID]['record'].channelNames) > 0:
                            channelSet = None  # will load complete datagroup
                        if channelSet is not None and buf[recordID]['record'].CANOpen:
                            # adds CANOpen channels if existing in not empty channelSet
                            if buf[recordID]['record'].CANOpen == 'time':
                                channelSet.update(('ms', 'days'))
                            elif buf[recordID]['record'].CANOpen == 'date':
                                channelSet.update(('ms', 'minute', 'hour', 'day', 'month', 'year'))

                    buf.read(channelSet, info)  # reads raw data from data block with DATA and DATABlock classes

                    # processing data from buf then transfer to self
                    for recordID in list(buf.keys()): # for each record in data block
                        if 'record' in buf[recordID]:
                            master_channel = buf[recordID]['record'].master['name']
                            if master_channel in self and self[master_channel][dataField] is not None:
                                master_channel = ''.join([master_channel, '_{}'.format(dataGroup)])
                            for chan in buf[recordID]['record']: # for each channel class
                                if (channelSet is None or chan.name in channelSet) \
                                        and not chan.type == 'Inv': # normal channel
                                    if chan.channelType(info) not in (3, 6):  # not virtual channel
                                        # in case record is used for several channels
                                        recordName = buf[recordID]['record'].recordToChannelMatching[chan.name]
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
                                            and not buf[recordID]['record'].hiddenBytes and\
                                            0 < bitCount < 64 and bitCount not in (8, 16, 32) \
                                            and temp is not None\
                                            and temp.dtype.kind not in ('S', 'U'):
                                        # if channel data do not use complete bytes and Ctypes
                                        if chan.signalDataType(info) in (0, 1, 2, 3):  # integers
                                            if chan.bitOffset(info) > 0:
                                                temp = right_shift(temp, chan.bitOffset(info))
                                            mask = int(pow(2, bitCount) - 1)  # masks isBitUnit8
                                            temp = bitwise_and(temp, mask)
                                        else:  # should not happen
                                            print('bit count and offset not applied to correct data type ',  chan.name, file=stderr)
                                    else:  # data using full bytes
                                        pass

                                    if temp is not None: # channel contains data
                                        # string data decoding
                                        if temp.dtype.kind == 'S':
                                            signalDataType = chan.signalDataType(info)
                                            if signalDataType == 6: # string ISO-8859-1 Latin
                                                encoding = 'latin-1'
                                            elif signalDataType == 7: # UTF-8
                                                encoding = 'UTF-8'
                                            elif signalDataType == 8: # UTF-16 low endian
                                                encoding = 'UTF-16LE'
                                            elif signalDataType == 9: # UTF-16 big endian
                                                encoding = 'UTF-16BE'
                                            else:
                                                encoding = None
                                            if encoding is not None:
                                                temp2 = empty(len(temp), dtype='U{}'.format(temp.dtype.str[-1]))
                                                for t in range(temp.size):
                                                    try:
                                                        temp2[t] = temp[t].decode(encoding, 'ignore')
                                                    except:
                                                        print('Cannot decode channel ' + chan.name, file=stderr)
                                                        temp2[t] = ''
                                                temp = temp2

                                        # channel creation
                                        self.add_channel(dataGroup, chan.name, temp, \
                                            master_channel, \
                                            master_type=chan.channelSyncType(info), \
                                            unit=chan.unit(info), \
                                            description=chan.desc(info), \
                                            conversion=chan.conversion(info), \
                                            info=chan.CNBlock(info), \
                                            compression=compression)
                                        if chan.channelType(info) == 4:  # sync channel
                                            # attach stream to be synchronised
                                            self.setChannelAttachment(chan.name, chan.attachment(info.fid, info))

                                elif chan.type == 'Inv' and \
                                        (channelSet is None or chan.name in channelSet):
                                    # invalid bytes, no bits processing
                                    self.add_channel(dataGroup, chan.name,
                                            buf[recordID]['data'].__getattribute__(chan.name),
                                            master_channel,
                                            master_type=0,
                                            unit='',
                                            description='',
                                            info=None,
                                            compression=compression)
                        del buf[recordID]
                    del buf
                if minimal > 1:
                    # clean CN, CC and CG info to free memory
                    info.cleanDGinfo(dataGroup)
        info.fid.close()  # close file

        if convertAfterRead and not compression:
            self._noDataLoading = False
            self._convertAllChannel4()
        # print( 'Finished in ' + str( time.clock() - inttime ) , file=stderr)

    def _getChannelData4(self, channelName):
        """Returns channel numpy array

        Parameters
        ----------------
        channelName : str
            channel name

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
            if vect is None: # noDataLoading reading argument flag activated
                if self.info.fid is None or (self.info.fid is not None and self.info.fid.closed):
                    (self.info.fid, self.info.fileName, self.info.zipfile) = _open_MDF(self.fileName)
                self.read4(fileName=None, info=self.info, channelList=[channelName], convertAfterRead=False)
            return self._convertChannelData4(self.getChannel(channelName),
                        channelName, self.convert_tables)[channelName]
        else:
            return None

    def _convertChannelData4(self, channel, channelName, convert_tables, multiProc=False, Q=None):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channelName : dict
            channel dict containing keys like 'data', 'unit', 'comment' and potentially 'conversion' dict
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
        if conversionField in channel:  # there is conversion property
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
                            proc.append(Process(target=convertChannelData4, \
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

    def write4(self, fileName=None):
        """Writes simple mdf 4.1 file

        Parameters
        ----------------
        fileName : str, optional
            Name of file
            If file name is not input, written file name will be the one read with appended '_new' string before extension

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        pointers = {}  # records pointers of blocks when writing

        # Starts first to write ID and header
        fid = open(fileName, 'wb')  # buffering should automatically be set
        # IDBLock writing
        temp = IDBlock()
        temp.write(fid)

        # Header Block
        temp = HDBlock()
        pointers.update(temp.write(fid))

        # Header Block comments
        _writePointer(fid, pointers['HD']['MD'], fid.tell())
        temp = CommentBlock()
        temp.write(fid, self.file_metadata, 'HD')

        # file history block
        temp = FHBlock()
        _writePointer(fid, pointers['HD']['FH'], fid.tell()) 
        pointers.update(temp.write(fid))
        # File History comment
        _writePointer(fid, pointers['FH']['MD'], fid.tell())
        temp = CommentBlock()
        temp.write(fid, self.file_metadata, 'FH')

        # write DG block
        _writePointer(fid, pointers['HD']['DG'], fid.tell()) # first DG

        pointers['DG'] = {}
        pointers['CG'] = {}
        pointers['CN'] = {}
        DG_flag = 0
        ndataGroup = len(self.masterChannelList)
        for dataGroup in range(ndataGroup):
            masterChannel = list(self.masterChannelList.keys())[dataGroup]
            # writes dataGroup Block
            temp = DGBlock()
            DG_pointer = fid.tell()
            pointers['DG'][dataGroup] = temp.write(fid)
            if DG_flag:
                _writePointer(fid, pointers['DG'][dataGroup-1]['DG'], pointers['DG'][dataGroup]['block_start'])  # Next DG
            DG_flag = pointers['DG'][dataGroup]['block_start']
            _writePointer(fid, pointers['DG'][dataGroup]['CG'], fid.tell())  # First CG

            # write CGBlock
            temp = CGBlock()
            cg_cycle_count = len(self._getChannelData4(masterChannel))
            cg_data_bytes = 0
            pointers['CG'][dataGroup] = temp.write(fid, cg_cycle_count, cg_data_bytes)

            # write channels
            record_byte_offset = 0
            CN_flag = 0
            nChannel = len(self.masterChannelList[masterChannel])
            number_of_channel = 0
            nRecords = 0
            dataList = ()
            dataTypeList = ''
            for nchannel in range(nChannel):
                channel = self.masterChannelList[masterChannel][nchannel]
                data = self.getChannelData(channel)
                # no interest to write invalid bytes as channel
                if channel.find('invalid_bytes') == -1 and len(data) > 0:
                    dataList = dataList + (data, )
                    number_of_channel +=1
                    temp = CNBlock()
                    temp['cn_val_range_min'] = npmin(data)
                    temp['cn_val_range_max'] = npmax(data)
                    temp['cn_flags'] = 16 # only Bit 4: Limit range valid flag
                    if masterChannel is not channel:
                        temp['cn_type'] = 0
                        temp['cn_sync_type'] = 0
                    else:
                        temp['cn_type'] = 2  # master channel
                        temp['cn_sync_type'] = 1  # default is time channel
                        nRecords = len(data)

                    cn_numpy_dtype = data.dtype
                    cn_numpy_kind = data.dtype.kind
                    if cn_numpy_dtype in ('uint8', 'uint16', 'uint32', 'uint64', 'bool'):
                        data_type = 0  # LE
                    elif cn_numpy_dtype in ('int8', 'int16', 'int32', 'int64'):
                        data_type = 2  # LE
                    elif cn_numpy_dtype in ('float32', 'float64'):
                        data_type = 4  # LE
                    elif cn_numpy_kind == 'S':
                        data_type = 6
                    elif cn_numpy_kind == 'U':
                        data_type = 7  # UTF-8
                    elif cn_numpy_kind == 'V':
                        data_type = 10  # bytes
                    else:
                        print(cn_numpy_dtype, cn_numpy_kind, file=stderr)
                        raise Exception('Not recognized dtype')
                    temp['cn_data_type'] = data_type
                    temp['cn_bit_offset'] = 0 # always byte aligned
                    if cn_numpy_dtype in ('float64', 'int64', 'uint64'):
                        bit_count = 64
                        byte_count = 8
                    elif cn_numpy_dtype in ('float32', 'int32', 'uint32'):
                        bit_count = 32
                        byte_count = 4
                    elif cn_numpy_dtype in ('uint16', 'int16'):
                        bit_count = 16
                        byte_count = 2
                    elif cn_numpy_dtype in ('uint8', 'int8', 'bool'):
                        bit_count = 8
                        byte_count = 1
                    elif cn_numpy_kind in ('S', 'U', 'V'):
                        byte_count = len(data[0])
                        # if string, not considered
                        bit_count = 8 * byte_count
                    else:
                        bit_count = 8
                        byte_count = 1
                    temp['cn_byte_offset'] = record_byte_offset
                    record_byte_offset += byte_count
                    temp['cn_bit_count'] = bit_count
                    if data.dtype.kind not in ('S', 'U', 'V'):
                        dataTypeList = ''.join([dataTypeList, data.dtype.char])
                    else:
                        dataTypeList = ''.join([dataTypeList, '{}s'.format(data.dtype.itemsize)])
                    # write channel block
                    pointers['CN'][nchannel] = temp.write(fid)
                    if CN_flag:
                        _writePointer(fid, pointers['CN'][nchannel-1]['CN'], pointers['CN'][nchannel]['block_start'])  # Next DG
                    else:
                        CN_flag = pointers['CN'][nchannel]['block_start']
                        _writePointer(fid, pointers['CG'][dataGroup]['CN'], CN_flag)  # first CN block pointer in CG
                    # write channel name
                    _writePointer(fid, pointers['CN'][nchannel]['TX'], fid.tell())
                    temp = CommentBlock()
                    temp.write(fid, channel, 'TX')

                    # write channel unit
                    unit = self.getChannelUnit(channel)
                    if unit is not None and len(unit) > 0:
                        temp = CommentBlock()
                        block_start = temp.write(fid, unit, 'TX')
                        _writePointer(fid, pointers['CN'][nchannel]['Unit'], block_start)

                    # write channel description
                    desc = self.getChannelDesc(channel)
                    if desc is not None and len(desc) > 0:
                        temp = CommentBlock()
                        block_start = temp.write(fid, desc, 'TX')
                        _writePointer(fid, pointers['CN'][nchannel]['Comment'], block_start)

            # writes size of record in CG
            _writePointer(fid, pointers['CG'][dataGroup]['cg_data_bytes'], record_byte_offset)

            # data writing
            # write data pointer in datagroup
            DTposition = _writeHeader(fid, b'##DT', 24 + record_byte_offset * nRecords, 0)
            _writePointer(fid, pointers['DG'][dataGroup]['data'], DTposition)
            # dumps data vector from numpy
            fid.write(fromarrays(dataList).tobytes(order='F'))

        # print(pointers, file=stderr)
        fid.close()


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
        print('Please install sympy to convert channel ', file=stderr)
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
        print('X values for interpolation of channel are not increasing', file=stderr)


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
        print('X values for interpolation of channel are not increasing', file=stderr)


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
    maxlen = max([len(str(ref)) for ref in cc_ref])
    temp = empty(len(vect), dtype='U{}'.format(maxlen))  # initialize empty array with max length
    # checks for scaling
    try:
        from sympy import lambdify, symbols
    except:
        print('Please install sympy to convert channel ', file=stderr)
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
                print('implement missing conversion')
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
        print('Please install sympy to convert channel ', file=stderr)
    X = symbols('X')
    for ref in range(len(cc_ref)):
        if isinstance(cc_ref[ref], CCBlock):
            if cc_ref[ref]['cc_type'] == 3:
                # formula to be applied
                cc_ref[ref] = lambdify(X, cc_ref[ref]['cc_ref']['Comment']
                                       , modules='numpy', dummify=False)
            elif cc_ref[ref]['cc_type'] == 1: # linear conversion
                cc_ref[ref] = lambdify(X, '{0}* X + {1}'.format(cc_val[1], cc_val[0])
                                       , modules='numpy', dummify=False)
            else:  # identity, non conversion
                cc_ref[ref] = lambdify(X, 'X'
                                       , modules='numpy', dummify=False)
            # Otherwise a string
    # look up in range keys
    temp = []
    for Lindex in range(len(vect)):
        key_index = val_count  # default index if not found
        for i in range(val_count):
            if key_min[i] < vect[Lindex] < key_max[i]:
                key_index = i
                break
        if callable(cc_ref[key_index]):
            # TXBlock string
            temp.append(cc_ref[key_index](vect[i]))
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
