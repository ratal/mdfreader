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
    Python version currently running, needed for compatibility of both python 2.6+ and 3.2+

mdf4reader module
--------------------------

"""
from numpy.core.records import fromstring, fromfile
from numpy import array, recarray, append, asarray, empty, zeros, dtype, where
from numpy import arange, right_shift, bitwise_and, all, diff, interp
from struct import Struct, pack
from struct import unpack as structunpack
from math import pow
from sys import version_info, exc_info, byteorder
from io import open  # for python 3 and 2 consistency
#try:
from mdfreader.mdfinfo4 import info4, MDFBlock, ATBlock, IDBlock, HDBlock, DGBlock, CGBlock, CNBlock, MDFBlock, FHBlock, CommentBlock
from mdfreader.mdf import mdf_skeleton
# except:
#     from mdfinfo4 import info4, MDFBlock, ATBlock, IDBlock, HDBlock, DGBlock, CGBlock, CNBlock, MDFBlock, FHBlock, CommentBlock
#     from mdf import mdf_skeleton
from time import gmtime, strftime
from multiprocessing import Queue, Process
PythonVersion = version_info
PythonVersion = PythonVersion[0]

LINK = '<Q'
REAL = '<d'
BOOL = '<h'
UINT8 = '<B'
BYTE = '<c'
INT16 = '<h'
UINT16 = '<H'
UINT32 = '<I'
INT32 = '<i'
UINT64 = '<Q'
INT64 = '<q'
CHAR = '<c'


def DATABlock(record, parent_block, channelList=None, sortedFlag=True):
    """ DATABlock converts raw data into arrays

    Parameters
    ----------------
    record : class
        record class instance describing a channel group record
    parent_block : class
        MDFBlock class containing at least parent block header
    channelList : list of str, optional
        defines list of channels to only read, can be slow but saves memory, for big files
    sortedFlag : bool, optional
        flag to know if data block is sorted (only one Channel Group in block)
        or unsorted (several Channel Groups identified by a recordID). As unsorted block can
        contain CG records in random order, block is processed iteratively,
        not in raw like sorted -> much slower reading

    Returns
    ---------
    a recarray containing the channels data

    Notes
    --------
    This function will read DTBlock, RDBlock, DZBlock (compressed), RDBlock (VLSD), sorted or unsorted
    """

    if parent_block['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal data block
        if sortedFlag:
            if channelList is None and not record.hiddenBytes and\
                    record.byte_aligned: # No channel list and length of records corresponds to C datatypes
                return fromstring(parent_block['data'], dtype=record.numpyDataRecordFormat, shape=record.numberOfRecords, names=record.dataRecordName)
            else: # record is not byte aligned or channelList not None
                return record.readBitarray(parent_block['data'], channelList)
        else:  # unsorted reading
            print('not implemented yet unsorted data block reading')  # to be implemented if needed, missing example file

    elif parent_block['id'] in ('##SD', b'##SD'):
        if record.signalDataType == 6:
            format = 'ISO8859'
        elif record.signalDataType == 7:
            format = 'utf-8'
        elif record.signalDataType == 8:
            format = '<utf-16'
        elif record.signalDataType == 9:
            format = '>utf-16'
        pointer = 0
        buf = []  # buf=array(record.numberOfRecords,dtype='s'+str(record.maxLengthVLSDRecord))

        while pointer < parent_block['length'] - 24:
            VLSDLen = structunpack('I', parent_block['data'][pointer:pointer + 4])[0]  # length of data
            pointer += 4
            buf.append(parent_block['data'][pointer:pointer + VLSDLen].decode(format).rstrip('\x00'))  # to be improved, removing append
            # buf[nElement]=fid.read(VLSDLen).decode(format).rstrip('\x00')
            pointer += VLSDLen
        buf = equalizeStringLength(buf)
        return array(buf)

    elif parent_block['id'] in ('##DZ', b'##DZ'):  # zipped data block
        # uncompress data
        parent_block['data'] = decompress_datablock(parent_block['data'], parent_block['dz_zip_type'],
                parent_block['dz_zip_parameter'], parent_block['dz_org_data_length'])
        if channelList is None and sortedFlag:  # reads all blocks if sorted block and no channelList defined
            if record.byte_aligned and not record.hiddenBytes:
                record.numberOfRecords = parent_block['dz_org_data_length'] // record.CGrecordLength
                return fromstring(parent_block['data'], dtype=record.numpyDataRecordFormat, shape=record.numberOfRecords, names=record.dataRecordName)
            else:
                return record.readBitarray(parent_block['data'], channelList)
        elif channelList is not None and sortedFlag:  # sorted data but channel list requested
            return record.readBitarray(parent_block['data'], channelList)
        else:  # unsorted reading
            # reads only the channels using offset functions, channel by channel.
            buf = {}
            position = 0
            recordIdCFormat = record[list(record.keys())[0]]['record'].recordIDCFormat
            recordIDsize = record[list(record.keys())[0]]['record'].recordIDsize
            VLSDStruct = Struct('I')
            # initialise data structure
            for recordID in record:
                for channelName in record[recordID]['record'].channelNames:
                    buf[channelName] = []  # empty(record.numberOfRecords,dtype=record[recordID]['record'].dataFormat)
                    # index[channelName]=0
            # read data
            while position < len(parent_block['data']):
                recordID = recordIdCFormat.unpack(parent_block['data'][position:position + recordIDsize])[0]
                if not record[recordID]['record'].Flags & 0b1:  # not VLSD CG)
                    temp = record.readRecord(recordID, parent_block['data'][position:position + record[recordID]['record'].CGrecordLength + 1], channelList)
                    position += record[recordID]['record'].CGrecordLength
                    for channelName in temp:
                        buf[channelName].append(temp[channelName])  # to remove append
                        # buf[channelName][index[recordID]]=temp[channelName]
                        # index[recordID] += 1
                else:  # VLSD CG
                    position += recordIDsize
                    VLSDLen = VLSDStruct.unpack(parent_block['data'][position:position + 4])[0]  # VLSD length
                    position += 4
                    temp = parent_block['data'][position:position + VLSDLen - 1]
                    if record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType == 6:
                        temp = temp.decode('ISO8859')
                    elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType == 7:
                        temp = temp.decode('utf-8')
                    elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType == 8:
                        temp = temp.decode('<utf-16')
                    elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType == 9:
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
        buf[i] += ' ' * (maxlen - len(buf[i]))
    return buf

def append_field(rec, name, arr, numpy_dtype=None):
    """ append new field in a recarray

    Parameters
    ----------------
    rec : numpy recarray
    name : str
        name of field to be appended
    arr : numpy array to be appended
    numpy_dtype : numpy dtype, optional
        apply same dtype as arr by default but can be modified

    Returns
    -----------
    numpy recarray appended
    """
    arr = asarray(arr)
    if numpy_dtype is None:
        numpy_dtype = arr.dtype
    newdtype = dtype(rec.dtype.descr + [(name, numpy_dtype)])
    newrec = empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec

def change_field_name(arr, old_name, new_name):
    """ modifies name of field in a recarray

    Parameters
    ----------------
    arr : numpy recarray
    old_name : str
        old field
    new_name : str
        new field

    Returns
    -----------
    numpy recarray with modified field name
    """
    names = list(arr.dtype.names)
    for n in range(len(names)):
        if names[n] == old_name:
            break
    names[n] = new_name
    names = tuple(names)
    arr.dtype.names = names
    return arr

def bits_to_bytes(nBits):
    """ Converts number of bits into number of aligned bytes
    
    Parameters
    -------------
    nBits : int
        number of bits
        
    Returns
    ----------
    number of equivalent bytes
    """
    if nBits <= 8:
        nBytes = 1
    elif nBits <= 16:
        nBytes = 2
    elif nBits <= 32:
        nBytes = 4
    elif nBits <= 64:
        nBytes = 8
    else:
        nBytes = nBits // 8
        if not nBits %8  == 0:
            nBytes += 1
    return nBytes
    
class DATA(dict):

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
    read(channelList, zip=None)
        Reads data block
    load(record, zip=None, nameList=None)
        Reads sorted data block from record definition
    readRecord(recordID, buf, channelList=None):
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

    def read(self, channelList, zip=None):
        """Reads data block

        Parameters
        ----------------
        channelList : list of str
            list of channel names
        zip : bool, optional
            flag to track if data block is compressed
        """
        if len(self) == 1:  # sorted dataGroup
            recordID = list(self.keys())[0]
            record = self[recordID]['record']
            self[recordID]['data'] = self.load(record, zip=None, nameList=channelList, sortedFlag=True)
            for cn in record.VLSD:  # VLSD channels
                if channelList is None or record[cn].name in channelList:
                    temp = DATA(self.fid, record[cn].data)  # all channels
                    rec = self[recordID]['data']  # recarray from data block
                    # determine maximum length of values in VLSD for array dtype
                    # record[cn].maxLengthVLSDRecord=max(diff(rec[convertName(record[cn].name)])-4)
                    temp = temp.load(record[cn], zip=None, nameList=channelList, sortedFlag=True)
                    rec = change_field_name(rec, convertName(record[cn].name), convertName(record[cn].name) + '_offset')
                    rec = append_field(rec, convertName(record[cn].name), temp)
                    self[recordID]['data'] = rec.view(recarray)
        else:  # unsorted DataGroup
            self.type = 'unsorted'
            data = self.load(self, zip=None, nameList=channelList, sortedFlag=False)
            for recordID in list(self.keys()):
                self[recordID]['data'] = {}
                for channel in self[recordID]['record']:
                    self[recordID]['data'][convertName(channel.name)] = data[channel.name]

    def load(self, record, zip=None, nameList=None, sortedFlag=True):  # reads sorted data
        """Reads data block from record definition

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        zip : bool, optional
            flag to track if data block is compressed
        nameList : list of str, optional
            list of channel names

        Returns
        -----------
        numpy recarray of data
        """
        temps = MDFBlock()
        # block header
        temps.loadHeader(self.fid, self.pointerTodata)
        if temps['id'] in ('##DL', b'##DL'):  # data list block
            # link section
            temps['dl_dl_next'] = temps.mdfblockread(self.fid, LINK, 1)
            temps['dl_data'] = {}
            temps['dl_data'][0] = [temps.mdfblockread(self.fid, LINK, 1) for Link in range(temps['link_count'] - 1)]
            # data section
            temps['dl_flags'] = temps.mdfblockread(self.fid, UINT8, 1)
            temps['dl_reserved'] = temps.mdfblockreadBYTE(self.fid, 3)
            temps['dl_count'] = temps.mdfblockread(self.fid, UINT32, 1)
            if temps['dl_flags']:  # equal length datalist
                temps['dl_equal_length'] = temps.mdfblockread(self.fid, UINT64, 1)
            else:  # datalist defined by byte offset
                temps['dl_offset'] = temps.mdfblockread(self.fid, UINT64, temps['dl_count'])
            if temps['dl_dl_next']:
                index = 1
            while temps['dl_dl_next']:  # reads pointers to all data blocks (DT, RD, SD, DZ)
                temp = MDFBlock()
                temp.loadHeader(self.fid, temps['dl_dl_next'])
                temps['dl_dl_next'] = temp.mdfblockread(self.fid, LINK, 1)
                temps['dl_data'][index] = [temp.mdfblockread(self.fid, LINK, 1) for Link in range(temp['link_count'] - 1)]
                index += 1
            if temps['dl_count']:
                # read and concatenate raw blocks
                buf = bytearray()
                for DL in temps['dl_data']:
                    for pointer in temps['dl_data'][DL]:
                        # read fist data blocks linked by DLBlock to identify data block type
                        data_block = MDFBlock()
                        data_block.loadHeader(self.fid, pointer)
                        if data_block['id'] in ('##DT', '##RD', b'##DT', b'##RD', '##SD', b'##SD'):
                            buf.extend(self.fid.read(data_block['length'] - 24))
                        elif data_block['id'] in ('##DZ', b'##DZ'):
                            data_block['dz_org_block_type'] = data_block.mdfblockreadCHAR(self.fid, 2)
                            data_block['dz_zip_type'] = data_block.mdfblockread(self.fid, UINT8, 1)
                            data_block['dz_reserved'] = data_block.mdfblockreadBYTE(self.fid, 1)
                            data_block['dz_zip_parameter'] = data_block.mdfblockread(self.fid, UINT32, 1)
                            data_block['dz_org_data_length'] = data_block.mdfblockread(self.fid, UINT64, 1)
                            data_block['dz_data_length'] = data_block.mdfblockread(self.fid, UINT64, 1)
                            data_block['data'] = decompress_datablock(self.fid.read(data_block['dz_data_length']),
                                    data_block['dz_zip_type'],
                                    data_block['dz_zip_parameter'], data_block['dz_org_data_length'])
                            buf.extend(data_block['data'])
                            data_block['id'] = '##DT' # do not uncompress in DATABlock function
                data_block['data'] = buf
                temps['data'] = DATABlock(record, parent_block=data_block, channelList=nameList, sortedFlag=sortedFlag)
            else:  # empty datalist
                temps['data'] = None
        elif temps['id'] in ('##HL', b'##HL'):  # header list block for DZBlock
            # link section
            temps['hl_dl_first'] = temps.mdfblockread(self.fid, LINK, 1)
            # data section
            temps['hl_flags'] = temps.mdfblockread(self.fid, UINT16, 1)
            temps['hl_zip_type'] = temps.mdfblockread(self.fid, UINT8, 1)
            temps['hl_reserved'] = temps.mdfblockreadBYTE(self.fid, 5)
            self.pointerTodata = temps['hl_dl_first']
            temps['data'] = self.load(record, zip=temps['hl_zip_type'], nameList=nameList, sortedFlag=sortedFlag)
        elif temps['id'] in ('##DT', '##RD', b'##DT', b'##RD'):  # normal sorted data block, direct read
            temps['data'] = record.readSortedRecord(self.fid, self.pointerTodata, channelList=nameList)
        elif temps['id'] in ('##SD', b'##SD'):  # VLSD
            temps['data'] = self.fid.read(temps['length'] - 24)
            temps['data'] = DATABlock(record, parent_block=temps, channelList=nameList, sortedFlag=sortedFlag)
        elif temps['id'] in ('##DZ', b'##DZ'):  # zipped data block
            temps['dz_org_block_type'] = temps.mdfblockreadCHAR(self.fid, 2)
            temps['dz_zip_type'] = temps.mdfblockread(self.fid, UINT8, 1)
            temps['dz_reserved'] = temps.mdfblockreadBYTE(self.fid, 1)
            temps['dz_zip_parameter'] = temps.mdfblockread(self.fid, UINT32, 1)
            temps['dz_org_data_length'] = temps.mdfblockread(self.fid, UINT64, 1)
            temps['dz_data_length'] = temps.mdfblockread(self.fid, UINT64, 1)
            temps['data'] = self.fid.read(temps['dz_data_length'])
            temps['data'] = DATABlock(record, parent_block=temps, channelList=nameList, sortedFlag=sortedFlag)
        else:
            raise Exception('unknown data block')
        return temps['data']

    def readRecord(self, recordID, buf, channelList=None):
        """ read record from a buffer

        Parameters
        ----------------
        recordID : int
            record identifier
        buf : str
            buffer of data from file to be converted to channel raw data
        channelList : list of str
            list of channel names to be read
        """
        return self[recordID]['record'].readRecordBuf(buf, channelList)

class channel():

    """ channel class gathers all about channel structure in a record

    Attributes
    --------------
    name : str
        Name of channel
    unit : str, default empty string
        channel unit
    desc : str
        channel description
    type : str
        channel type. Can be 'standard', 'CANOpen' or 'InvalidBytes'
    conversion : info class
        conversion dictionnary
    channelNumber : int
        channel number corresponding to mdfinfo3.info3 class
    channelGroup : int
        channel group number corresponding to mdfinfo3.info3 class
    dataGroup : int
        data group number corresponding to mdfinfo3.info3 class
    signalDataType : int
        signal type according to specification
    bitCount : int
        number of bits used to store channel record
    nBytes : int
        number of bytes (1 byte = 8 bits) taken by channel record
    dataFormat : str
        numpy dtype as string
    Format :
        C format understood by fread
    CFormat : struct class instance
        struct instance to convert from C Format
    byteOffset : int
        position of channel record in complete record in bytes
    bitOffset : int
        bit position of channel value inside byte in case of channel having bit count below 8
    RecordFormat : list of str
        dtype format used for numpy.core.records functions ((name,name_title),str_stype)
    channelType : int
        channel type ; 0 fixed length data, 1 VLSD, 2 master, 3 virtual master,
        4 sync, 5 MLSD, 6 virtual data
    channelSyncType : int
        channel synchronisation type ; 0 None, 1 Time, 2 Angle, 3 Distance, 4 Index
    posByteBeg : int
        start position in number of byte of channel record in complete record
    posByteEnd : int
        end position in number of byte of channel record in complete record
    posBitBeg : int
        start position in number of bit of channel record in complete record
    posBitEnd : int
        end position in number of bit of channel record in complete record
    VLSD_CG_Flag : bool
        flag when Channel Group VLSD is used
    data : int
        pointer to data block linked to a channel (VLSD, MLSD)

    Methods
    ------------
    __init__()
        constructor
    __str__()
        to print class attributes
    attachment(fid, info)
        in case of sync channel attached
    set(info, dataGroup, channelGroup, channelNumber, recordIDsize)
        standard channel initialisation
    setCANOpen(info, dataGroup, channelGroup, channelNumber, recordIDsize, name)
        CANOpen channel initialisation
    setInvalidBytes(info, dataGroup, channelGroup, recordIDsize, byte_aligned)
        Invalid Bytes channel initialisation
    """

    def __init__(self):
        """ channel class constructor
        """
        self.name = ''
        type = 'standard'
        self.channelNumber = 0
        self.dataGroup = 0
        self.channelGroup = 0
        self.signalDataType = 0
        self.channelSyncType = 0
        self.bitCount = 0
        self.nBytes = 0
        self.little_endian = True
        self.dataFormat = ''
        self.CABlock = 0
        self.channelType = 0
        self.data = 0
        self.maxLengthVLSDRecord = 0
        self.RecordFormat = (('', ''), '')
        self.Format = ''
        self.bitOffset = 0
        self.byteOffset = 0
        self.posByteBeg = 0
        self.posByteEnd = 0
        self.posBitBeg = 0
        self.posBitEnd = 0
        self.VLSD_CG_Flag = False
        self.unit = ''
        self.desc = ''
        self.conversion = None

    def __str__(self):
        output = str(self.channelNumber) + ' '
        output += str(self.RecordFormat) + ' '
        output += str(self.Format) + ' '
        output += 'ChannelType ' + str(self.channelType) + ' '
        output += 'DataType ' + str(self.signalDataType) + ' '
        output += 'bitOffset ' + str(self.bitOffset) + ' '
        output += 'ByteBeg ' + str(self.posByteBeg) + ' '
        output += 'ByteEnd ' + str(self.posByteEnd) + ' '
        output += 'BitBeg ' + str(self.posBitBeg) + ' '
        output += 'BitEnd ' + str(self.posBitEnd)
        output += 'unit ' + self.unit
        output += 'description ' + self.desc
        return output
    
    def attachment(self, fid, info): # in case of sync channel attached to channel
        return ATBlock(fid, info['CNBlock'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data'])

    def set(self, info, dataGroup, channelGroup, channelNumber, recordIDsize):
        """ standard record channel initialisation

        Parameters
        ------------
        info : mdfinfo4.info4 class
        dataGroup : int
            data group number in mdfinfo4.info4 class
        channelGroup : int
            channel group number in mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        recordIDsize : int
            size of record ID in Bytes
        """
        self.type = 'standard'
        self.name = info['CNBlock'][dataGroup][channelGroup][channelNumber]['name']
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data_type']
        self.channelSyncType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_sync_type']
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_count']
        self.nBytes = bits_to_bytes(self.bitCount)
        if self.signalDataType in (1, 3, 5, 9):  # endianness
            self.little_endian = False
        else:
            self.little_endian = True
        if info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_composition'] and \
                'CABlock' in info['CNBlock'][dataGroup][channelGroup][channelNumber]:  # channel array
            self.dataFormat = arrayformat4(self.signalDataType, self.bitCount)
            self.CABlock = info['CNBlock'][dataGroup][channelGroup][channelNumber]['CABlock']
            self.nBytes *= self.CABlock['PNd']  # calculates total array size in bytes
            if self.CABlock['ca_ndim'] > 1:
                array_desc = str(tuple(self.CABlock['ca_dim_size']))
            else:
                array_desc = str(self.CABlock['ca_dim_size'])
            self.dataFormat = array_desc + self.dataFormat
        else:  # not channel array
            self.dataFormat = arrayformat4(self.signalDataType, self.bitCount)
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_type']
        if self.channelType in (1, 4, 5):  # if VSLD or sync or max length channel
            self.data = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data']
            if self.channelType == 1:  # if VSLD
                self.dataFormat = arrayformat4(0, self.bitCount)
                self.maxLengthVLSDRecord = 0  # initialises max length of SDBlock elements to 0 for later calculation
        if self.channelType in (2, 3):  # master channel, add datagroup number to avoid overwriting same sames like time
            self.name += str(dataGroup)
        self.RecordFormat = ((self.name, convertName(self.name)), self.dataFormat)
        if self.signalDataType not in (13, 14):  # processed by several channels
            if not self.channelType == 1:  # if not VSLD
                self.Format = datatypeformat4(self.signalDataType, self.bitCount)
            else:  # VLSD
                self.Format = datatypeformat4(0, self.bitCount)
            self.CFormat = Struct(self.Format)
        self.bitOffset = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_offset']
        self.byteOffset = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_byte_offset']
        self.posByteBeg = recordIDsize + self.byteOffset
        self.posByteEnd = self.posByteBeg + self.nBytes
        self.posBitBeg = self.posByteBeg * 8 + self.bitOffset
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.VLSD_CG_Flag = False
        if 'unit' in list(info['CCBlock'][dataGroup][channelGroup][channelNumber].keys()):
            self.unit = info['CCBlock'][dataGroup][channelGroup][channelNumber]['unit']
        elif 'unit' in list(info['CNBlock'][dataGroup][channelGroup][channelNumber].keys()):
            self.unit = info['CNBlock'][dataGroup][channelGroup][channelNumber]['unit']
        else:
            self.unit = ''
        if 'Comment' in self.unit:
            self.unit = self.unit['Comment']
        if 'Comment' in list(info['CNBlock'][dataGroup][channelGroup][channelNumber].keys()):
            self.desc = info['CNBlock'][dataGroup][channelGroup][channelNumber]['Comment']
            if self.desc is not None and 'description' in self.desc:
                self.desc = self.desc['description']
            if self.desc is not None and 'name' in self.desc:
                self.desc = self.desc['name']
        else:
            self.desc = ''
        self.conversion = info['CCBlock'][dataGroup][channelGroup][channelNumber]

    def setCANOpen(self, info, dataGroup, channelGroup, channelNumber, recordIDsize, name):
        """ CANOpen channel intialisation

        Parameters
        ------------
        info : mdfinfo4.info4 class
        dataGroup : int
            data group number in mdfinfo4.info4 class
        channelGroup : int
            channel group number in mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        recordIDsize : int
            size of record ID in Bytes
        name : str
            name of channel. Should be in ('ms', 'day', 'days', 'hour', 'month', 'min', 'year')
        """
        self.type = 'CANOpen'
        self.name = name
        self.unit = name
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data_type']
        self.channelSyncType = 0
        if self.name == 'ms':
            offset = 0
            self.desc = 'ms'
            if self.signalDataType == 13:
                self.bitCount = 16
                self.dataFormat = '<u2'
            else:
                self.bitCount = 32
                self.dataFormat = '<u4'
        elif self.name == 'days':
            self.bitCount = 16
            self.dataFormat = '<u2'
            offset = 4
            self.desc = 'days'
        else:
            self.bitCount = 8
            self.dataFormat = '<u1'
            if self.name == 'min':
                offset = 2
                self.desc = 'min'
            elif self.name == 'hour':
                offset = 3
                self.desc = 'hour'
            elif self.name == 'day':
                offset = 4
                self.desc = 'day'
            elif self.name == 'month':
                offset = 5
                self.desc = 'month'
            elif self.name == 'year':
                offset = 6
                self.desc = 'year'
            else:
                raise('CANopen type not understood')
        self.nBytes = bits_to_bytes(self.bitCount)
        self.little_endian = True
        self.channelType = 0
        self.maxLengthVLSDRecord = 0  # initialises max length of SDBlock elements to 0 for later calculation
        self.RecordFormat = ((self.name, convertName(self.name)), self.dataFormat)
        self.bitOffset = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_offset'] + offset*8
        self.byteOffset = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_byte_offset'] + offset
        self.posByteBeg = recordIDsize + self.byteOffset
        self.posByteEnd = self.posByteBeg + self.nBytes
        self.posBitBeg = self.posByteBeg * 8 + self.bitOffset
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.signalDataType = 0

    def setInvalidBytes(self, info, dataGroup, channelGroup, recordIDsize, byte_aligned=True):
        """ invalid_bytes channel initialisation

        Parameters
        ----------
        info : mdfinfo4.info4 class
        dataGroup : int
            data group number in mdfinfo4.info4 class
        channelGroup : int
            channel group number in mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        recordIDsize : int
            size of record ID in Bytes
        byte_aligned : Bool
            Flag for byte alignement
        """
        self.type = 'InvalidBytes'
        self.name = u'invalid_bytes' + str(dataGroup)
        self.unit = ''
        self.desc = ''
        self.conversion = ''
        if byte_aligned:
            self.signalDataType = 10  # byte array
        else:
            self.signalDataType = 0  # uint LE
        self.little_endian = True  # default to little endian
        self.nBytes = info['CGBlock'][dataGroup][channelGroup]['cg_invalid_bytes']
        self.bitCount = self.nBytes * 8
        self.channelType = 0  # fixed length data
        self.dataFormat = str(self.nBytes) + 'V'
        self.RecordFormat = ((self.name, convertName(self.name)), self.dataFormat)
        self.Format = str(self.nBytes) + 's'
        self.CFormat = Struct(self.Format)
        self.bitOffset = 0
        self.byteOffset = info['CGBlock'][dataGroup][channelGroup]['cg_data_bytes']
        self.posByteBeg = recordIDsize + self.byteOffset
        self.posByteEnd = self.posByteBeg + self.nBytes
        self.posBitBeg = self.posByteBeg * 8 + self.bitOffset
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.VLSD_CG_Flag = False
        self.invalid_bit = {}
        for channelNumber in info['CNBlock'][dataGroup][channelGroup]:
            name = info['CNBlock'][dataGroup][channelGroup][channelNumber]['name']
            self.invalid_bit[name] = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_invalid_bit_pos']
        self.invalid_bytes = None

    def validity_channel(self, channelName):
        """ extract channel validity bits

        Parameters
        ----------
        channelName : str
            channel name

        Return
        -------
        Numpy vector of validity

        """
        if self.type == 'InvalidBytes':
            return bitwise_and(right_shift(self.invalid_bytes, self.invalid_bit[channelName]), 1)
        else:
            raise('asking for invalid byte array but channel is not invalid byte type')

class record(list):

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
    channelNames : list
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
    readSortedRecord(fid, pointer, channelList=None)
    readRecordBuf(buf, channelList=None)
    readBitarray(bita, channelList=None)
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
        self.master['name'] = 'master_' + str(dataGroup)
        self.master['number'] = None
        self.Flags = 0
        self.VLSD_CG = {}
        self.VLSD = []
        self.MLSD = {}
        self.recordToChannelMatching = {}
        self.channelNames = []
        self.byte_aligned = True
        self.hiddenBytes = False
        self.invalid_channel = None
        self.CANOpen = None

    def __str__(self):
        output = 'Channels : ' + str(len(self)) + '\n'
        for chan in self:
            output += chan.name + '\n'
            output += str(chan.channelNumber) + ' '
            output += str(chan.RecordFormat) + ' '
            output += str(chan.Format) + ' '
            output += 'ChannelType ' + str(chan.channelType) + ' '
            output += 'DataType ' + str(chan.signalDataType) + '\n'
        output += 'Datagroup number : ' + str(self.dataGroup) + '\n'
        output += 'Byte aligned : ' + str(self.byte_aligned) + '\n'
        if self.master['name'] is not None:
            output += 'Master channel : ' + self.master['name'] + '\n'
        output += 'Numpy records format : \n'
        for record in self.numpyDataRecordFormat:
            output += str(record) + '\n'
        output += 'VLSD_CG' + str(self.VLSD_CG)
        return output

    def addChannel(self, info, channelNumber):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        """
        Channel = channel()
        self.append(Channel.set(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize))
        self.channelNames.append(self[-1].name)

    def loadInfo(self, info):
        """ gathers records related from info class

        Parameters
        ----------------
        info : mdfinfo4.info4 class
        """
        self.CGrecordLength = info['CGBlock'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        self.recordIDsize = info['DGBlock'][self.dataGroup]['dg_rec_id_size']
        if not self.recordIDsize == 0:  # no record ID
            self.dataRecordName.append('RecordID' + str(self.channelGroup))
            format = (self.dataRecordName[-1], convertName(self.dataRecordName[-1]))
            if self.recordIDsize == 1:
                self.numpyDataRecordFormat.append((format, 'uint8'))
                self.recordIDCFormat = Struct('B')
                self.recordLength += 1
                self.CGrecordLength += 1
            elif self.recordIDsize == 2:
                self.numpyDataRecordFormat.append((format, 'uint16'))
                self.recordIDCFormat = 'H'
                self.recordLength += 2
                self.CGrecordLength += 2
            elif self.recordIDsize == 3:
                self.numpyDataRecordFormat.append((format, 'uint32'))
                self.recordIDCFormat = 'I'
                self.recordLength += 3
                self.CGrecordLength += 3
            elif self.recordIDsize == 4:
                self.numpyDataRecordFormat.append((format, 'uint64'))
                self.recordIDCFormat = 'L'
                self.recordLength += 4
                self.CGrecordLength += 4
        self.recordID = info['CGBlock'][self.dataGroup][self.channelGroup]['cg_record_id']
        self.numberOfRecords = info['CGBlock'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        self.Flags = info['CGBlock'][self.dataGroup][self.channelGroup]['cg_flags']
        if 'MLSD' in info:
            self.MLSD = info['MLSD']
        for channelNumber in list(info['CNBlock'][self.dataGroup][self.channelGroup].keys()):
            Channel = channel()
            Channel.set(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize)
            if Channel.channelType in (2, 3):  # master channel found
                if self.master['number'] is None or Channel.channelSyncType==1:  # new master channel found
                    # or more than 1 master channel, priority to time channel
                    self.master['number'] = channelNumber
                    self.master['name'] = Channel.name
            if Channel.channelType in (0, 1, 2, 4, 5):  # not virtual channel
                if Channel.signalDataType == 13:
                    for name in ('ms', 'min', 'hour', 'day', 'month', 'year'):
                        Channel = channel() # new object otherwise only modified
                        Channel.setCANOpen(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize, name)
                        self.append(Channel)
                        self.channelNames.append(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(Channel.RecordFormat)
                    self.recordLength += 7
                    self.CANOpen = 'date'
                elif Channel.signalDataType == 14:
                    for name in ('ms', 'days'):
                        Channel = channel()
                        Channel.setCANOpen(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize, name)
                        self.append(Channel)
                        self.channelNames.append(name)
                        self.dataRecordName.append(name)
                        self.recordToChannelMatching[name] = name
                        self.numpyDataRecordFormat.append(Channel.RecordFormat)
                    self.recordLength += 6
                    self.CANOpen = 'time'
                else:
                    self.append(Channel)
                    self.channelNames.append(Channel.name)
                    # Checking if several channels are embedded in bytes
                    embedded_bytes = False
                    if len(self) > 1:
                        # all channels are already ordered in record based on byte_offset and bit_offset
                        # so just comparing with previous channel
                        if Channel.byteOffset >= self[-2].byteOffset and \
                                Channel.posBitBeg < 8 * (self[-2].byteOffset + self[-2].nBytes) and \
                                Channel.posBitEnd > 8 * (self[-2].byteOffset + self[-2].nBytes):  # not byte aligned
                            self.byte_aligned = False
                        if Channel.posBitBeg >= 8 * self[-2].byteOffset \
                                and Channel.posBitEnd <= 8 * (self[-2].byteOffset + self[-2].nBytes):  # bit(s) in byte(s)
                            embedded_bytes = True
                            if self.recordToChannelMatching: # not first channel
                                self.recordToChannelMatching[Channel.name] = self.recordToChannelMatching[self[-2].name]
                            else: # first channel
                                self.recordToChannelMatching[Channel.name] = Channel.name
                                self.numpyDataRecordFormat.append(Channel.RecordFormat)
                                self.dataRecordName.append(Channel.name)
                                self.recordLength += Channel.nBytes
                    if not embedded_bytes:  # adding bytes
                        self.recordToChannelMatching[Channel.name] = Channel.name
                        self.numpyDataRecordFormat.append(Channel.RecordFormat)
                        self.dataRecordName.append(Channel.name)
                        self.recordLength += Channel.nBytes
                    if 'VLSD_CG' in info:  # is there VLSD CG
                        for recordID in info['VLSD_CG']:  # look for VLSD CG Channel
                            if info['VLSD_CG'][recordID]['cg_cn'] == (self.channelGroup, channelNumber):
                                self.VLSD_CG[recordID] = info['VLSD_CG'][recordID]
                                self.VLSD_CG[recordID]['channel'] = Channel
                                self.VLSD_CG[recordID]['channelName'] = Channel.name
                                self[-1].VLSD_CG_Flag = True
                                break
                    if Channel.channelType == 1:  # VLSD channel
                        self.VLSD.append(Channel.channelNumber)
            elif Channel.channelType in (3, 6):  # master virtual channel
                self.append(Channel)  # channel calculated based on record index later in conversion function
                self.channelNames.append(Channel.name)

        if info['CGBlock'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']:  # invalid bytes existing
            self.CGrecordLength += info['CGBlock'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            self.recordLength += info['CGBlock'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
            invalid_bytes = channel()
            invalid_bytes.setInvalidBytes(info, self.dataGroup, self.channelGroup, self.recordIDsize, self.byte_aligned)
            self.invalid_channel = invalid_bytes
            self.append(self.invalid_channel)
            self.channelNames.append(self.invalid_channel.name)
            self.recordToChannelMatching[self.invalid_channel.name] = self.invalid_channel.name
            self.numpyDataRecordFormat.append(self.invalid_channel.RecordFormat)
            self.dataRecordName.append(self.invalid_channel.name)
        # check for hidden bytes
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True

    def readSortedRecord(self, fid, pointer, channelList=None):
        """ reads record, only one channel group per datagroup

        Parameters
        ----------------
        fid : float
            file identifier
        pointer
            position in file of data block beginning
        channelList : list of str, optional
            list of channel to read

        Returns
        -----------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)

        Notes
        --------
        If channelList is None, read data using numpy.core.records.fromfile that is rather quick.
        However, in case of large file, you can use channelList to load only interesting channels or
        only one channel on demand, but be aware it might be much slower.
        """
        if channelList is None:  # reads all, quickest but memory consuming
            if self.byte_aligned:
                #print(self.numpyDataRecordFormat)
                return fromfile(fid, dtype=self.numpyDataRecordFormat, shape=self.numberOfRecords, names=self.dataRecordName)
            else:
                return self.readBitarray(fid.read(self.CGrecordLength * self.numberOfRecords), channelList)
        else:  # reads only some channels from a sorted data block
            if len(list(set(channelList) & set(self.channelNames))) > 0:  # are channelList in this dataGroup
                try:
                    return self.readBitarray(fid.read(self.CGrecordLength * self.numberOfRecords), channelList)
                except:  # still memory efficient but takes time
                    print('Unexpected error:', exc_info())
                    print('dataRead crashed, back to python data reading')
                    # check if master channel is in the list
                    if not self.master['name'] in channelList:
                        channelList.append(self.master['name'])  # adds master channel
                    rec = {}
                    recChan = []
                    numpyDataRecordFormat = []
                    for channel in channelList:  # initialise data structure
                        rec[channel] = 0
                    for channel in self:  # list of channel classes from channelList
                        if channel.name in channelList:
                            recChan.append(channel)
                            numpyDataRecordFormat.append(channel.RecordFormat)
                    rec = zeros((self.numberOfRecords, ), dtype=numpyDataRecordFormat)
                    recordLength = self.recordIDsize + self.CGrecordLength
                    for r in range(self.numberOfRecords):  # for each record,
                        buf = fid.read(recordLength)
                        for channel in recChan:
                            rec[channel.name][r] = channel.CFormat.unpack(buf[channel.posByteBeg:channel.posByteEnd])[0]
                    return rec.view(recarray)

    def readRecordBuf(self, buf, channelList=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        channelList : list of str, optional
            list of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        temp = {}
        if channelList is None:
            channelList = self.channelNames
        for Channel in self:  # list of channel classes from channelList
            if Channel.name in channelList and not Channel.VLSD_CG_Flag:
                temp[Channel.name] = Channel.CFormat.unpack(buf[Channel.posByteBeg:Channel.posByteEnd])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def readBitarray(self, bita, channelList=None):
        """ reads stream of record bytes using bitarray module needed for not byte aligned data
        
        Parameters
        ------------
        bita : stream
            stream of bytes
        channelList : List of str, optional
            list of channel to read
        
        Returns
        --------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
        """
        # initialise data structure
        if channelList is None:
            channelList = self.channelNames
        format = []
        for channel in self:
            if channel.name in channelList:
                format.append(channel.RecordFormat)
        if format:  # at least some channels should be parsed
            buf = recarray(self.numberOfRecords, format)
            try: # use rather cython compiled code for performance
                from .dataRead import dataRead
                for chan in range(len(self)):
                    if self[chan].name in channelList:
                        buf[self[chan].name] = dataRead(bytes(bita), self[chan].bitCount, \
                                self[chan].signalDataType, self[chan].RecordFormat[1], \
                                self.numberOfRecords, self.CGrecordLength, \
                                self[chan].bitOffset, self[chan].posByteBeg, \
                                self[chan].posByteEnd)
                        # dataRead already took care of byte order, switch to native
                        if (buf[self[chan].name].dtype.byteorder == '>' and byteorder == 'little') or \
                                buf[self[chan].name].dtype.byteorder == '<' and byteorder == 'big':
                            buf[self[chan].name] = buf[self[chan].name].newbyteorder()
                return buf
            except:
                print('Unexpected error:', exc_info())
                print('dataRead crashed, back to python data reading')
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
                    if self[chan].name in channelList:
                        temp = [B[self[chan].posBitBeg + record_bit_size * i:\
                                self[chan].posBitEnd + record_bit_size * i]\
                                 for i in range(self.numberOfRecords)]
                        nbytes = len(temp[0].tobytes())
                        if not nbytes == self[chan].nBytes and \
                                self[chan].signalDataType not in (6, 7, 8, 9, 10, 11, 12): # not Ctype byte length
                            byte = bitarray(8 * (self[chan].nBytes - nbytes), endian='little')
                            byte.setall(False)
                            if self[chan].signalDataType not in (2, 3):  # not signed integer
                                for i in range(self.numberOfRecords):  # extend data of bytes to match numpy requirement
                                    temp[i].extend(byte)
                            else:  # signed integer (two's complement), keep sign bit and extend with bytes
                                temp = signedInt(temp, byte)
                        nTrailBits = self[chan].nBytes*8 - self[chan].bitCount
                        if self[chan].signalDataType in (2, 3) and \
                                nbytes == self[chan].nBytes and \
                                nTrailBits > 0:  # Ctype byte length but signed integer
                            trailBits = bitarray(nTrailBits, endian='little')
                            temp = signedInt(temp, trailBits)
                        if 's' not in self[chan].Format:
                            temp = [self[chan].CFormat.unpack(temp[i].tobytes())[0] \
                                    for i in range(self.numberOfRecords)]
                        else:
                            temp = [temp[i].tobytes() \
                                    for i in range(self.numberOfRecords)]
                        buf[self[chan].name] = asarray(temp)
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

    def read4(self, fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True, filterChannelNames=False):
        """ Reads mdf 4.x file data and stores it in dict

        Parameters
        ----------------
        fileName : str, optional
            file name

        info : mdfinfo4.info4 class
            info3 class containing all MDF Blocks

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
        """
        self.multiProc = multiProc

        if self.fileName is None and info is not None:
            self.fileName = info.fileName
        elif fileName is not None:
            self.fileName = fileName

        # inttime = time.clock()
        # Read information block from file
        if info is None:
            info = info4(self.fileName, None)

        # reads metadata
        fileDateTime = gmtime(info['HDBlock']['hd_start_time_ns'] / 1000000000)
        ddate = strftime('%Y-%m-%d', fileDateTime)
        ttime = strftime('%H:%M:%S', fileDateTime)
        def returnField(obj, field):
            if field in obj:
                return obj[field]
            else:
                return ''
        if 'Comment' in info['HDBlock']:
            Comment =  info['HDBlock']['Comment']
            author = returnField(Comment, 'author')
            organisation = returnField(Comment, 'department')
            project = returnField(Comment,'project')
            subject = returnField(Comment,'subject')
            comment = returnField(Comment,'TX')
            self.add_metadata(author=author, organisation=organisation, \
                    project=project, subject=subject, comment=comment, \
                    date=ddate, time=ttime)
        else:
            self.add_metadata(date=ddate, time=ttime)
        try:
            fid = open(self.fileName, 'rb')
        except IOError:
            raise Exception('Can not find file ' + self.fileName)

        for dataGroup in info['DGBlock'].keys():
            if not info['DGBlock'][dataGroup]['dg_data'] == 0 and \
                    info['CGBlock'][dataGroup][0]['cg_cycle_count']:  # data exists
                # Pointer to data block
                pointerToData = info['DGBlock'][dataGroup]['dg_data']
                buf = DATA(fid, pointerToData)

                for channelGroup in info['CGBlock'][dataGroup].keys():
                    temp = record(dataGroup, channelGroup)  # create record class
                    temp.loadInfo(info)  # load all info related to record
                    buf.addRecord(temp)  # adds record to DATA
                    recordID = info['CGBlock'][dataGroup][channelGroup]['cg_record_id']
                    if temp.master['name'] is not None \
                            and buf[recordID]['record'].channelNames:
                        self.masterChannelList[temp.master['name']] = []
                        if channelList is not None and temp.master['name'] not in channelList:
                            channelList.append(temp.master['name'])  # adds master channel in channelList if missing
                    if channelList is not None and buf[recordID]['record'].CANOpen: # adds CANOpen channels if existing in not empty channelList
                        if buf[recordID]['record'].CANOpen == 'time':
                            channelList.extend(['ms', 'days'])
                        elif buf[recordID]['record'].CANOpen == 'date':
                            channelList.extend(['ms', 'min', 'hour', 'day', 'month', 'year'])

                buf.read(channelList)  # reads raw data from data block with DATA and DATABlock classes

                # processing data from buf then transfer to self
                for recordID in list(buf.keys()): # for each record in data block
                    if 'record' in buf[recordID]:
                        master_channel = buf[recordID]['record'].master['name']
                        if master_channel in self.keys():
                            master_channel += '_' + str(dataGroup)
                        for chan in buf[recordID]['record']: # for each channel class
                            if (channelList is None or chan.name in channelList) \
                                    and not chan.type == 'InvalidBytes': # normal channel
                                if chan.channelType not in (3, 6):  # not virtual channel
                                    recordName = buf[recordID]['record'].recordToChannelMatching[chan.name]  # in case record is used for several channels
                                    if 'data' in buf[recordID] and \
                                            buf[recordID]['data'] is not None: # no data in channel group
                                        if isinstance(buf[recordID]['data'],recarray):
                                            temp = buf[recordID]['data'].__getattribute__(convertName(recordName))  # extract channel vector
                                        else:
                                            temp = buf[recordID]['data'][convertName(recordName)] 
                                    else:
                                        temp = None
                                else:  # virtual channel
                                    temp = arange(buf[recordID]['record'].numberOfRecords)

                                # Process concatenated bits inside uint8
                                if buf[recordID]['record'].byte_aligned and \
                                        0 < chan.bitCount < 64 and chan.bitCount not in (8, 16, 32) \
                                        and temp is not None\
                                        and temp.dtype.kind not in ('S', 'U'):  # if channel data do not use complete bytes and Ctypes
                                    if chan.signalDataType in (0, 1, 2, 3):  # integers
                                        if chan.bitOffset > 0:
                                            temp = right_shift(temp, chan.bitOffset)
                                        mask = int(pow(2, chan.bitCount) - 1)  # masks isBitUnit8
                                        temp = bitwise_and(temp, mask)
                                    else:  # should not happen
                                        print('bit count and offset not applied to correct data type ',  chan.name)
                                else:  # data using full bytes
                                    pass

                                if temp is not None: # channel contains data
                                    # string data decoding
                                    if temp.dtype.kind != 'U':
                                        for t in range(temp.size):
                                            try:
                                                if chan.signalDataType == 6: # string ISO-8859-1 Latin
                                                        temp[t] = temp[t].decode('latin-1', 'ignore')
                                                elif chan.signalDataType == 7: # UTF-8
                                                    temp[t] = temp[t].decode('UTF-8', 'ignore')
                                                elif chan.signalDataType == 8: # UTF-16 low endian
                                                    temp[t] = temp[t].decode('UTF-16LE', 'ignore')
                                                elif chan.signalDataType == 9: # UTF-16 big endian
                                                    temp[t] = temp[t].decode('UTF-16BE', 'ignore')
                                            except:
                                                    print('Cannot decode channel ' + chan.name)
                                                    temp[t] = ''
                                    # channel creation
                                    self.add_channel(dataGroup, chan.name, temp, \
                                        master_channel, \
                                        master_type=chan.channelSyncType, \
                                        unit=chan.unit, \
                                        description=chan.desc, \
                                        conversion=chan.conversion)
                                    if chan.channelType == 4:  # sync channel
                                        # attach stream to be synchronised
                                        self.setChannelAttachment(chan.name, chan.attachment(fid, info))

                            elif chan.type == 'InvalidBytes' and \
                                    (channelList is None or chan.name in channelList):  # invalid bytes, no bits processing
                                self.add_channel(dataGroup, chan.name, \
                                        buf[recordID]['data'].__getattribute__(convertName(chan.name)), \
                                        master_channel, \
                                        master_type=0, \
                                        unit='', \
                                        description='')
                    del buf[recordID]
        fid.close()  # close file

        if convertAfterRead:
            self._convertAllChannel4()
        # print( 'Finished in ' + str( time.clock() - inttime ) )

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
            return convertChannelData4(self.getChannel(channelName), \
                    channelName, self.convert_tables)[channelName]
        else:
            raise KeyError(channelName + ' channel not in dictionary')
            

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
        MDFBlock.writePointer(fid, pointers['HD']['MD'], fid.tell())
        temp = CommentBlock()
        temp.write(fid, self.file_metadata, 'HD')

        # file history block
        temp = FHBlock()
        MDFBlock.writePointer(fid, pointers['HD']['FH'], fid.tell()) 
        pointers.update(temp.write(fid))
        # File History comment
        MDFBlock.writePointer(fid, pointers['FH']['MD'], fid.tell())
        temp = CommentBlock()
        temp.write(fid, self.file_metadata, 'FH')

        # write DG block
        MDFBlock.writePointer(fid, pointers['HD']['DG'], fid.tell()) # first DG

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
                MDFBlock.writePointer(fid, pointers['DG'][dataGroup-1]['DG'], pointers['DG'][dataGroup]['block_start'])  # Next DG
            DG_flag = pointers['DG'][dataGroup]['block_start']
            MDFBlock.writePointer(fid, pointers['DG'][dataGroup]['CG'], fid.tell())  # First CG

            # write CGBlock
            temp = CGBlock()
            cg_cycle_count = len(self.getChannelData(masterChannel))
            cg_data_bytes = 0
            pointers['CG'][dataGroup] = temp.write(fid, cg_cycle_count, cg_data_bytes)
            MDFBlock.writePointer(fid, pointers['CG'][dataGroup]['CN'], fid.tell())  # first channel block pointer

            # write channels
            record_byte_offset = 0
            CN_flag = 0
            nChannel = len(self.masterChannelList[masterChannel])
            for nchannel in range(nChannel):
                channel = self.masterChannelList[masterChannel][nchannel]
                temp = CNBlock()
                if masterChannel is not channel:
                    temp['cn_type'] = 0
                    temp['cn_sync_type'] = 0
                else:
                    temp['cn_type'] = 2  # master channel
                    temp['cn_sync_type'] = 1  # default is time channel

                data = self.getChannelData(channel)
                cn_numpy_dtype = data.dtype
                cn_numpy_kind = data.dtype.kind
                if cn_numpy_dtype in ('uint8', 'uint16', 'uint32', 'uint64'):
                    data_type = 0  # LE
                elif cn_numpy_dtype in ('int8', 'int16', 'int32', 'int64'):
                    data_type = 2  # LE
                elif cn_numpy_dtype in ('float32', 'float64'):
                    data_type = 4  # LE
                elif cn_numpy_kind == 'S':
                    data_type = 6
                elif cn_numpy_kind == 'U':
                    data_type = 7  # UTF-8
                else:
                    raise Exception('Not recognized dtype')
                    return cn_numpy_dtype
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
                elif cn_numpy_dtype in ('uint8', 'int8'):
                    bit_count = 8
                    byte_count = 1
                else:
                    bit_count = 8  # if string, not considered
                    byte_count = 1
                temp['cn_byte_offset'] = record_byte_offset
                record_byte_offset += byte_count
                temp['cn_bit_count'] = bit_count
                pointers['CN'][nchannel] = temp.write(fid)
                if CN_flag:
                    MDFBlock.writePointer(fid, pointers['CN'][nchannel-1]['CN'], pointers['CN'][nchannel]['block_start'])  # Next DG
                CN_flag = pointers['CN'][nchannel]['block_start']
                # write channel name
                MDFBlock.writePointer(fid, pointers['CN'][nchannel]['TX'], fid.tell())
                temp = CommentBlock()
                temp.write(fid, channel, 'TX')

                # write channel unit
                unit = self.getChannelUnit(channel)
                if unit is not None and not len(unit) == 0:
                    temp = CommentBlock()
                    block_start = temp.write(fid, unit, 'TX')
                    MDFBlock.writePointer(fid, pointers['CN'][nchannel]['Unit'], block_start)

                # Conversion blocks writing

            MDFBlock.writePointer(fid, pointers['CG'][dataGroup]['cg_data_bytes'], record_byte_offset)  # writes size of record in CG

            # data writing

        # print(pointers)
        fid.close()


def convertName(channelName):
    """ Adds '_title' to channel name for numpy.core.records methods."""
    if PythonVersion < 3:
        channelIdentifier = channelName.encode('utf-8') + '_title'
    else:
        channelIdentifier = str(channelName) + '_title'
    return channelIdentifier


def arrayformat4(signalDataType, numberOfBits):
    """ function returning numpy style string from channel data type and number of bits

    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBits : int
        number of bits taken by channel data in a record

    Returns
    -----------
    dataType : str
        numpy dtype format used by numpy.core.records to read channel raw data
    """

    if signalDataType in (0, 1):  # unsigned
        if numberOfBits <= 8:
            dataType = 'u1'
        elif numberOfBits <= 16:
            dataType = 'u2'
        elif numberOfBits <= 32:
            dataType = 'u4'
        elif numberOfBits <= 64:
            dataType = 'u8'
        else:
            dataType = str(bits_to_bytes(numberOfBits) // 8) + 'V'

    elif signalDataType in (2, 3):  # signed int
        if numberOfBits <= 8:
            dataType = 'i1'
        elif numberOfBits <= 16:
            dataType = 'i2'
        elif numberOfBits <= 32:
            dataType = 'i4'
        elif numberOfBits <= 64:
            dataType = 'i8'
        else:
            print('Unsupported number of bits for signed int ' + str(numberOfBits))

    elif signalDataType in (4, 5):  # floating point
        if numberOfBits == 32:
            dataType = 'f4'
        elif numberOfBits == 64:
            dataType = 'f8'
        else:
            print('Unsupported number of bit for floating point ' + str(numberOfBits))

    elif signalDataType == 6:  # string ISO-8859-1 Latin
        dataType = 'S' + str(numberOfBits // 8)
    elif signalDataType == 7:  # UTF-8
        dataType = 'S' + str(numberOfBits // 8)
    elif signalDataType == 8:  # UTF-16 low endian
        dataType = 'S' + str(numberOfBits // 8)
    elif signalDataType == 9:  # UTF-16 big endian
        dataType = 'S' + str(numberOfBits // 8)
    elif signalDataType == 10:  # bytes array
        dataType = 'V' + str(numberOfBits // 8)
    elif signalDataType in (11, 12):  # MIME sample or MIME stream
        dataType = 'V' + str(int(numberOfBits / 8))
    elif signalDataType in (13, 14):  # CANOpen date or time
        dataType = None
    else:
        print('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits)

    # deal with byte order
    if signalDataType in (0, 2, 4, 8):  # low endian
        dataType = '<' + dataType
    elif signalDataType in (1, 3, 5, 9):  # big endian
        dataType = '>' + dataType

    return dataType


def datatypeformat4(signalDataType, numberOfBits):
    """ function returning C format string from channel data type and number of bits

    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBits : int
        number of bits taken by channel data in a record

    Returns
    -----------
    dataType : str
        C format used by fread to read channel raw data
    """

    if signalDataType in (0, 1):  # unsigned int
        if numberOfBits <= 8:
            dataType = 'B'
        elif numberOfBits <= 16:
            dataType = 'H'
        elif numberOfBits <= 32:
            dataType = 'I'
        elif numberOfBits <= 64:
            dataType = 'Q'
        else:
            dataType = str(bits_to_bytes(numberOfBits) // 8) + 's'

    elif signalDataType in (2, 3):  # signed int
        if numberOfBits <= 8:
            dataType = 'b'
        elif numberOfBits <= 16:
            dataType = 'h'
        elif numberOfBits <= 32:
            dataType = 'i'
        elif numberOfBits <= 64:
            dataType = 'q'
        else:
            print(('Unsupported number of bits for signed int ' + str(signalDataType)))

    elif signalDataType in (4, 5):  # floating point
        if numberOfBits == 32:
            dataType = 'f'
        elif numberOfBits == 64:
            dataType = 'd'
        else:
            print(('Unsupported number of bit for floating point ' + str(signalDataType)))

    elif signalDataType in (6, 7, 8, 9, 10, 11, 12):  # string/bytes
        dataType = str(numberOfBits // 8) + 's'
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits))
    # deal with byte order
    if signalDataType in (0, 2, 4, 8):  # low endian
        dataType = '<' + dataType
    elif signalDataType in (1, 3, 5, 9):  # big endian
        dataType = '>' + dataType

    return dataType


def convertChannelData4(channel, channelName, convert_tables, multiProc=False, Q=None):
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
    vect = channel['data']
    if 'conversion' in channel:  # there is conversion property
        text_type = channel['data'].dtype.kind in ['S', 'U', 'V']  # channel of string or not ?
        conversion_type = channel['conversion']['type']
        conversion_parameter = channel['conversion']['parameters']
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
        print('Please install sympy to convert channel ')
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
        print('X values for interpolation of channel are not increasing')


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
        print('X values for interpolation of channel are not increasing')


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
    val_count = int(len(cc_val))
    cc_val = list(cc_val.values())  # converts dict to lists
    cc_ref = list(cc_ref.values())
    maxlen = len(max(cc_ref, key=len))
    temp = empty(len(vect), dtype='S' + str(maxlen))  # initialize empty array with max length
    key_index = where(vect[0] == cc_val)[0]  # look up for first value in vect
    if not len(key_index) == 0:  # value corresponding in cc_val
        temp[0] = cc_ref[key_index[0]]
    else:  # default value
        temp[0] = cc_ref[val_count]
    for Lindex in set(range(1, len(vect))):
        if vect[Lindex] == vect[Lindex - 1]:  # same value as before, no need to look further
            temp[Lindex] = temp[Lindex - 1]
        else:  # value changed from previous step
            key_index = where(vect[Lindex] == cc_val)[0]
            if not len(key_index) == 0:  # found match
                temp[Lindex] = cc_ref[key_index[0]]
            else:  # default
                temp[Lindex] = cc_ref[val_count]
    return temp


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
    # look up in range keys
    temp = []
    for Lindex in range(len(vect)):
        key_index = val_count  # default index if not found
        for i in range(val_count):
            if key_min[i] < vect[Lindex] < key_max[i]:
                key_index = i
                break
        temp.append(cc_ref[key_index])
    return temp


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
    return temp


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