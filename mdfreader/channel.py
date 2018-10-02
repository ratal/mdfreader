# -*- coding: utf-8 -*-
""" Measured Data Format file reader module.

Platform and python version
----------------------------------------
With Unix and Windows for python 2.7 and 3.4+

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Wed Oct 04 21:13:28 2017

Dependencies
-------------------
- Python >2.6, >3.4 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.4+

channel module
--------------------------

"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from struct import Struct
from warnings import warn
from .mdfinfo4 import ATBlock
from .mdf import _bits_to_bytes, _convertName

CAN_open_offset = {'ms': 0, 'days': 4, 'minute': 2, 'hour': 3, 'day': 4, 'month': 5, 'year': 6}


class channel4(object):
    __slots__ = ['channelNumber', 'channelGroup', 'dataGroup',
                 'type', 'name', 'VLSD_CG_Flag', 'nBytes', 'byteOffset']
    """ channel class gathers all about channel structure in a record

    Attributes
    --------------
    name : str
        Name of channel
    type : int
        channel type. Can be 'standard' : 0, 'Channel Array' : 1, 'Nested Channel Array' : 2,
        'CAN' : 3 or 'Invalid' : 4
    channelNumber : int
        channel number corresponding to mdfinfo4.info4 class
    channelGroup : int
        channel group number corresponding to mdfinfo4.info4 class
    dataGroup : int
        data group number corresponding to mdfinfo4.info4 class
    VLSD_CG_Flag : bool
        flag when Channel Group VLSD is used
    nBytes : int
        number of bytes taken by channel at each sampling

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
    setCANOpen(info, dataGroup, channelGroup, channelNumber,
                recordIDsize, name)
        CANOpen channel initialisation
    setInvalidBytes(info, dataGroup, channelGroup, recordIDsize, byte_aligned)
        Invalid Bytes channel initialisation
    recAttributeName : str
        Name of channel compatible with python attribute name conventions
    unit : str, default empty string
        channel unit
    desc : str
        channel description
    conversion : info class
        conversion dictionnary
    CNBlock : info class
        Channel Block info class
    signalDataType : int
        signal type according to specification
    bitCount : int
        number of bits used to store channel record
    calc_bytes : int
        number of bytes (1 byte = 8 bits) taken by channel record
    little_endian : Bool
        flag to inform of channel data endian
    dataFormat : str
        numpy dtype as string
    Format :
        C format understood by fread
    CFormat : struct class instance
        struct instance to convert from C Format
    byteOffset : int
        position of channel record in complete record in bytes
    bitOffset : int
        bit position of channel value inside byte in case of channel
        having bit count below 8
    RecordFormat : nested tuple of str
        dtype format used for numpy.core.records functions
        ((name_title,name),str_stype)
    nativeRecordFormat : nested tuple of str
        same as RecordFormat but using recAttributeName instead of name
    channelType : int
        channel type ; 0 fixed length data, 1 VLSD, 2 master, 3 virtual master,
        4 sync, 5 MLSD, 6 virtual data
    channelSyncType : int
        channel synchronisation type ; 0 None, 1 Time, 2 Angle,
        3 Distance, 4 Index
    posByteBeg : int
        start position in number of byte of channel record in complete record
    posByteEnd : int
        end position in number of byte of channel record in complete record
    posBitBeg : int
        start position in number of bit of channel record in complete record
    posBitEnd : int
        end position in number of bit of channel record in complete record
    maxLengthVLSDRecord :

    CABlock : CABlock class
        contains CABLock
    data : int
        pointer to data block linked to a channel (VLSD, MLSD)
    invalid_bit : dict
        dict of invalid bit channels data
    invalid_bytes : bytes
        byte containing invalid bit for each channel
    bit_masking_needed : boolean
        False if channel needs bit masking
    """

    def __init__(self):
        """ channel class constructor

        Attributes
        --------------

        name : str
            Name of channel
        type : int, default 0
            Can be 'standard' : 0, 'Channel Array' : 1, 'Nested Channel Array' : 2,
            'CAN' : 3 or 'Invalid' : 4
        channelNumber : int, default 0
            channel number corresponding to mdfinfo4.info4 class
        channelGroup : int, default 0
            channel group number corresponding to mdfinfo4.info4 class
        dataGroup : int, default 0
            data group number corresponding to mdfinfo4.info4 class
        VLSD_CG_Flag : bool, default False
            flag when Channel Group VLSD is used
        """
        self.name = ''
        self.type = 0
        self.channelNumber = 0
        self.dataGroup = 0
        self.channelGroup = 0
        self.VLSD_CG_Flag = False
        self.nBytes = 0
        self.byteOffset = 0

    def __str__(self):
        """ channel object attributes print

        """
        return '{0} {1} {2} {3} {4} {5}'.format(self.name,
                                                self.type,
                                                self.dataGroup,
                                                self.channelGroup,
                                                self.channelNumber,
                                                self.VLSD_CG_Flag)

    def attachment(self, fid, info):
        """ In case of sync channel attached to channel

        Parameters
        ----------------

        fid : class
            file identifier

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        ATBlock class from mdfinfo4 module
        """
        try:
            return ATBlock(fid, info['CN'][self.dataGroup][self.channelGroup]\
                                [self.channelNumber]['cn_data'])
        except KeyError:
            print('No Attachment block for this channel')

    def data(self, info):
        """ returns data block pointer for VLSD, MLD or sync channels"""
        try:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data']
        except KeyError:
            return None

    def CNBlock(self, info):
        """ channel block

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        CNBlock class from mdfinfo4 module
        """
        try:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]
        except KeyError:
            return None

    def signalDataType(self, info, byte_aligned=True):
        """ extract signal data type from info4 class

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks
        byte_aligned : bool
            flag activated if channel is part of a record byte aligned

        Returns
        -----------
        integer corresponding to channel data type
        0 unsigned integer little endian
        1 unsigned integer big endian
        2 signed integer little endian
        3 signed integer big endian
        4 float little endian
        5 float big endian
        6 string latin
        7 string utf-8
        9 string utf-16
        10 byte array
        11 mime sample
        12 mime stream
        13 CANopen date
        14 CANopen time
        """
        if not self.type == 4:  # Invalid bit channel
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_data_type']
        else:
            if byte_aligned:
                return 10  # byte array
            else:
                return 0  # uint LE

    def bitCount(self, info):
        """ calculates channel number of bits

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer corresponding to channel number of bits
        """
        if self.type in (0, 1, 2):  # standard or array channel
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count']
        elif self.type == 3:  # CAN channel
            if self.name == 'ms':
                if self.signalDataType(info) == 13:
                    return 16
                else:
                    return 32
            elif self.name == 'days':
                return 16
            else:
                return 8
        elif self.type == 4:  # Invalid bit channel
            return self.nBytes * 8
        else:
            warn('Not found channel type')

    def channelSyncType(self, info):
        """ Extracts channel sync type from info4

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer corresponding to channel sync type
        0 no sync, normal data
        1 time
        2 angle
        3 distance
        4 index
        """
        try:
            return info['CN'][self.dataGroup][self.channelGroup]\
                    [self.channelNumber]['cn_sync_type']
        except KeyError:
            return 0  # in case of invalid bytes channel

    def CABlock(self, info):
        """ Extracts channel CA Block from info4

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        CABlock object from mdfinfo4 module
        """
        try:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
        except KeyError:
            return None

    def isCABlock(self, info):
        try:
            info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            return True
        except KeyError:
            return False

    def recordIDsize(self, info):
        """ Extracts record id size from info4

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer describing record id size
        0 no record id used
        1 uint8
        2 uint16
        4 uint32
        8 uint64
        """
        return info['DG'][self.dataGroup]['dg_rec_id_size']

    def channelType(self, info):
        """ Extracts channel type from info4

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer describing channel type
        0 normal channel
        1 variable length
        2 master channel
        3 virtual master channel
        4 sync channel
        5 max length data
        6 virtual data channel
        """
        try:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_type']
        except KeyError:
            return 0  # in case of invalid bytes channel

    def isnumeric(self, info):
        """ check this is numeric channel from data type

            Parameters
            ----------------

            info : mdfinfo4.info4 class
                info4 class containing all MDF Blocks

            Returns
            -----------
            boolean, true if numeric channel, otherwise false
        """
        if info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'] < 6:
            return True
        else:
            return False

    def calc_bytes(self, info):
        """ calculates channel bytes number

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        number of bytes integer
        """
        if not self.type == 4:  # not channel containing invalid bit
            nBytes = _bits_to_bytes(info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count']
                                    + info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]
                                    ['cn_bit_offset'], self.isnumeric(info))
            if self.type in (1, 2):  # array channel
                nBytes *= self.CABlock(info)['PNd']
                Block = self.CABlock(info)
                while 'CABlock' in Block:  # nested array
                    Block = Block['CABlock']
                    nBytes *= Block['PNd']
            if self.type == 3:  # CAN channel
                if self.name == 'ms':
                    if info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'] == 13:
                        nBytes = 2
                    else:
                        nBytes = 4
                elif self.name == 'days':
                    nBytes = 2
                else:
                    nBytes = 1
            return nBytes
        else:
            return info['CG'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']

    def little_endian(self, info):
        """ check if channel is little endian

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        boolean
        """
        if info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]\
                ['cn_data_type'] in (1, 3, 5, 9):  # endianness
            return False
        else:
            return True

    def recAttributeName(self, info):
        """ clean up channel name from unauthorised characters

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        channel name compliant to python attributes names (for recarray)
        """
        return _convertName(self.name)

    def numpy_format(self, info):
        """ channel numpy.core.records data format

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        endian, dataType : string data format
        """
        endian = ''
        if self.type == 4:  # Invalid bit channel
            dataformat = '{}V'.format(self.nBytes)
        elif info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_composition'] and \
                'CABlock' in info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]:  # channel array
            CABlock = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            endian, dataformat = arrayformat4(
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'],
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count'] // 8)
            # calculates total array size in bytes
            array_desc = CABlock['ca_dim_size']
            Block = CABlock
            while 'CABlock' in Block:  # nested array
                Block = Block['CABlock']
                if isinstance(array_desc, list):
                    array_desc.append(Block['ca_dim_size'])
                else:
                    array_desc = [array_desc, Block['ca_dim_size']]
            if isinstance(array_desc, list):
                array_desc = str(tuple(array_desc))
            else:
                array_desc = str(array_desc)
            dataformat = array_desc + dataformat
        elif self.type == 3:  # CAN channel
            if self.name == 'ms':
                if self.signalDataType(info) == 13:
                    dataformat = 'u2'
                else:
                    dataformat = 'u4'
            elif self.name == 'days':
                dataformat = 'u2'
            else:
                dataformat = 'u1'
            endian = '<'
        else:  # not channel array
            endian, dataformat = arrayformat4(
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'],
                self.nBytes)
        return endian, dataformat

    def dataFormat(self, info):
        """ channel numpy.core.records data format

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data format
        """
        endian, dataformat = self.numpy_format(info)
        return ''.join([endian, dataformat])

    def nativedataFormat(self, info):
        endian, nativedataFormat = self.numpy_format(info)
        return nativedataFormat

    def Format(self, info):
        """ channel data C format

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data C format
        """
        signalDataType = self.signalDataType(info)
        if self.type == 0:  # standard channel
            if signalDataType not in (13, 14):
                if not self.channelType(info) == 1:  # if not VSLD
                    endian, dataType = datatypeformat4(signalDataType, self.nBytes)
                else:  # VLSD
                    endian, dataType = datatypeformat4(0, self.nBytes)
                return '{}{}'.format(endian, dataType)
        elif self.type in (1, 2):  # array channel
            CA = self.CABlock(info)
            nBytes = _bits_to_bytes(info['CN'][self.dataGroup][self.channelGroup]
                                    [self.channelNumber]['cn_bit_count'] + info['CN'][self.dataGroup][self.channelGroup]
                                    [self.channelNumber]['cn_bit_offset'], self.isnumeric(info))
            endian, dataType = datatypeformat4(signalDataType, nBytes)
            return '{}{}{}'.format(endian, CA['PNd'], dataType)
        elif self.type == 3:  # CAN channel
            if self.name == 'ms':
                if signalDataType == 13:
                    return 'H'
                else:
                    return 'I'
            elif self.name == 'days':
                return 'H'
            else:
                return 'B'
        elif self.type == 4:  # Invalid bit channel
            return '{}s'.format(self.nBytes)
        else:
            warn('Not found channel type')

    def CFormat(self, info):
        """ channel data C format struct object

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data C format struct object
        """
        return Struct(self.Format(info))

    def CANOpenOffset(self):
        """ CANopen channel bytes offset

        Returns
        -----------
        integer, channel bytes offset
        """
        try:
            return CAN_open_offset[self.name]
        except KeyError:
            warn('CANopen type not understood')
            return None

    def bitOffset(self, info):
        """ channel data bit offset in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bit offset
        """
        if self.type in (0, 1, 2):  # standard or channel array
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_offset']
        elif self.type == 3:  # CAN channel
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_offset'] \
                   + self.CANOpenOffset() * 8
        elif self.type == 4:  # Invalid bit channel
            return 0
        else:
            warn('Not found channel type')

    def calc_byteOffset(self, info):
        """ channel data bytes offset in record (without record id)

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes offset
        """
        if self.type in (0, 1, 2):  # standard or channel array
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_byte_offset']
        elif self.type == 3:  # CAN channel
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_byte_offset'] \
                   + self.CANOpenOffset()
        elif self.type == 4:  # Invalid bit channel
            return info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        else:
            warn('Not found channel type')

    def posByteBeg(self, info):
        """ channel data bytes starting position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes starting position
        """
        return self.recordIDsize(info) + self.byteOffset

    def posByteEnd(self, info):
        """ channel data bytes ending position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes ending position
        """
        return self.posByteBeg(info) + self.nBytes

    def posBitBeg(self, info):
        """ channel data bit starting position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bit starting position
        """
        return self.posByteBeg(info) * 8 + self.bitOffset(info)

    def posBitEnd(self, info):
        """ channel data bit ending position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bit ending position
        """
        return self.posBitBeg(info) + self.bitCount(info)

    def unit(self, info):
        """ channel unit

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        channel unit string
        """
        if self.channelNumber not in info['CC'][self.dataGroup][self.channelGroup]:
            return ''
        if 'unit' in info['CC'][self.dataGroup][self.channelGroup][self.channelNumber]:
            unit = info['CC'][self.dataGroup][self.channelGroup][self.channelNumber]['unit']
        elif 'unit' in info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]:
            unit = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['unit']
        else:
            unit = ''
        if 'Comment' in unit:
            unit = unit['Comment']
        return unit

    def desc(self, info):
        """ channel description

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        channel description string
        """
        if not self.type == 3:  # CAN channel
            if self.channelNumber in info['CN'][self.dataGroup][self.channelGroup]:
                if 'Comment' in info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]:
                    desc = info['CN'][self.dataGroup][self.channelGroup]\
                                [self.channelNumber]['Comment']
                    if (desc is not None) and isinstance(desc, dict):
                        if 'description' in desc:
                            desc = desc['description']
                        elif 'name' in desc:
                            desc = desc['name']
                else:
                    desc = ''
                return desc
            else:
                return 'Invalid Bytes DGgroup {0} CGgroup {1}'.format(self.dataGroup,
                                                                      self.channelGroup)
        else:
            return self.name

    def conversion(self, info):
        """ channel conversion CCBlock

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        CCBlock
        """
        try:
            return info['CC'][self.dataGroup][self.channelGroup][self.channelNumber]
        except KeyError:
            return None

    def set(self, info, dataGroup, channelGroup, channelNumber):
        """ channel initialisation

        Parameters
        ------------

        info : mdfinfo4.info4 class
        dataGroup : int
            data group number in mdfinfo4.info4 class
        channelGroup : int
            channel group number in mdfinfo4.info4 class
        channelNumber : int
            channel number in mdfinfo4.info4 class
        """
        self.name = info['CN'][dataGroup][channelGroup][channelNumber]['name']
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.type = 0
        if info['CN'][dataGroup][channelGroup][channelNumber]['cn_composition'] and \
                'CABlock' in info['CN'][dataGroup][channelGroup][channelNumber]:
            # channel array
            self.type = 1
            Block = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            if 'CABlock' in Block:  # nested array
                self.type = 2
        self.nBytes = self.calc_bytes(info)
        self.byteOffset = self.calc_byteOffset(info)

    def setCANOpen(self, info, dataGroup, channelGroup, channelNumber, name):
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
        name : str
            name of channel. Should be in ('ms', 'day', 'days', 'hour',
            'month', 'minute', 'year')
        """
        self.type = 3
        self.name = name
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.nBytes = self.calc_bytes(info)
        self.byteOffset = self.calc_byteOffset(info)

    def setInvalidBytes(self, info, dataGroup, channelGroup, channelNumber):
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
        """
        self.type = 4
        self.name = 'invalid_bytes{}'.format(dataGroup)
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.nBytes = self.calc_bytes(info)
        self.byteOffset = self.calc_byteOffset(info)

    def invalid_bit(self, info):
        """ extracts from info4 the channels valid bits positions

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        channel valid bit position
        """
        return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_invalid_bit_pos']

    def has_invalid_bit(self, info):
        flags = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_flags']
        if flags & 0b10:
            return True
        else:
            return False
    
    def changeChannelName(self, channelGroup):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channelGroup : int
            channelGroup bumber
        """
        self.name = '{0}_{1}'.format(self.name, channelGroup)

    def bit_masking_needed(self, info):
        """ Valid if bit masking need

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        boolean True if channel needs bit masking, otherwise False

        """
        if not self.nBytes == self.bitCount(info) / 8:
            self.bit_masking_needed = True
        else:
            self.bit_masking_needed = False


def arrayformat4(signalDataType, numberOfBytes):
    """ function returning numpy style string from channel data type and number of bits

    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    endian, dataType : str
        numpy dtype format used by numpy.core.records to read channel raw data
    """

    if signalDataType == 0:  # unsigned, low endian
        if numberOfBytes == 1:
            dataType = 'u1'
        elif numberOfBytes == 2:
            dataType = 'u2'
        elif numberOfBytes <= 4:
            dataType = 'u4'
        elif numberOfBytes <= 8:
            dataType = 'u8'
        else:
            dataType = '{}V'.format(numberOfBytes)
        endian = '<'

    elif signalDataType == 2:  # signed int, low endian
        if numberOfBytes == 1:
            dataType = 'i1'
        elif numberOfBytes == 2:
            dataType = 'i2'
        elif numberOfBytes <= 4:
            dataType = 'i4'
        elif numberOfBytes <= 8:
            dataType = 'i8'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(numberOfBytes))
        endian = '<'

    elif signalDataType == 4:  # floating point, low endian
        if numberOfBytes == 4:
            dataType = 'f4'
        elif numberOfBytes == 8:
            dataType = 'f8'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(numberOfBytes))
        endian = '<'

    elif signalDataType == 1:  # unsigned, big endian
        if numberOfBytes == 1:
            dataType = 'u1'
        elif numberOfBytes == 2:
            dataType = 'u2'
        elif numberOfBytes <= 4:
            dataType = 'u4'
        elif numberOfBytes <= 8:
            dataType = 'u8'
        else:
            dataType = '{}V'.format(numberOfBytes)
        endian = '>'

    elif signalDataType == 3:  # signed int, big endian
        if numberOfBytes == 1:
            dataType = 'i1'
        elif numberOfBytes == 2:
            dataType = 'i2'
        elif numberOfBytes <= 4:
            dataType = 'i4'
        elif numberOfBytes <= 8:
            dataType = 'i8'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(numberOfBytes))
        endian = '>'

    elif signalDataType == 5:  # floating point, big endian
        if numberOfBytes == 4:
            dataType = 'f4'
        elif numberOfBytes == 8:
            dataType = 'f8'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(numberOfBytes))
        endian = '>'

    elif signalDataType == 6:  # string ISO-8859-1 Latin
        dataType = 'S{}'.format(numberOfBytes)
        endian = ''
    elif signalDataType == 7:  # UTF-8
        dataType = 'S{}'.format(numberOfBytes)
        endian = ''
    elif signalDataType == 8:  # UTF-16 low endian
        dataType = 'S{}'.format(numberOfBytes)
        endian = '<'
    elif signalDataType == 9:  # UTF-16 big endian
        dataType = 'S{}'.format(numberOfBytes)
        endian = '>'
    elif signalDataType == 10:  # bytes array
        dataType = 'V{}'.format(numberOfBytes)
        endian = ''
    elif signalDataType in (11, 12):  # MIME sample or MIME stream
        dataType = 'V{}'.format(numberOfBytes)
        endian = ''
    elif signalDataType in (13, 14):  # CANOpen date or time
        dataType = ''
        endian = ''
    else:
        warn('Unsupported Signal Data Type {} {}'.format(signalDataType, numberOfBytes))

    return endian, dataType


def datatypeformat4(signalDataType, numberOfBytes):
    """ function returning C format string from channel data type and number of bits

    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    dataType : str
        C format used by fread to read channel raw data
    """

    if signalDataType == 0:  # unsigned int
        if numberOfBytes == 1:
            dataType = 'B'
        elif numberOfBytes == 2:
            dataType = 'H'
        elif numberOfBytes <= 4:
            dataType = 'I'
        elif numberOfBytes <= 8:
            dataType = 'Q'
        else:
            dataType = '{}s'.format(numberOfBytes)
        endian = '<'

    elif signalDataType == 1:  # unsigned int
        if numberOfBytes == 1:
            dataType = 'B'
        elif numberOfBytes == 2:
            dataType = 'H'
        elif numberOfBytes <= 4:
            dataType = 'I'
        elif numberOfBytes <= 8:
            dataType = 'Q'
        else:
            dataType = '{}s'.format(numberOfBytes)
        endian = '>'

    elif signalDataType == 2:  # signed int
        if numberOfBytes == 1:
            dataType = 'b'
        elif numberOfBytes == 2:
            dataType = 'h'
        elif numberOfBytes <= 4:
            dataType = 'i'
        elif numberOfBytes <= 8:
            dataType = 'q'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(signalDataType))
        endian = '<'

    elif signalDataType == 3:  # signed int
        if numberOfBytes == 1:
            dataType = 'b'
        elif numberOfBytes == 2:
            dataType = 'h'
        elif numberOfBytes <= 4:
            dataType = 'i'
        elif numberOfBytes <= 8:
            dataType = 'q'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(signalDataType))
        endian = '>'

    elif signalDataType == 4:  # floating point
        if numberOfBytes == 4:
            dataType = 'f'
        elif numberOfBytes == 8:
            dataType = 'd'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(signalDataType))
        endian = '<'

    elif signalDataType == 5:  # floating point
        if numberOfBytes == 4:
            dataType = 'f'
        elif numberOfBytes == 8:
            dataType = 'd'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(signalDataType))
        endian = '>'

    elif signalDataType in (6, 7, 10, 11, 12):  # string/bytes
        dataType = '{}s'.format(numberOfBytes)
        endian = ''

    elif signalDataType == 8:  # UTF16 string/bytes
        dataType = '{}s'.format(numberOfBytes)
        endian = '<'

    elif signalDataType == 9:  # UTF16 string/bytes
        dataType = '{}s'.format(numberOfBytes)
        endian = '>'

    else:
        warn('Unsupported Signal Data Type {} {}'.format(signalDataType, numberOfBytes))

    return endian, dataType


class Channel3:

    """ Channel class gathers all about channel structure in a record

    Attributes
    --------------
    name : str
        Name of channel
    unit : str, default empty string
        channel unit
    desc : str
        channel description
    conversion : info class
        conversion dictionnary
    channelNumber : int
        channel number corresponding to mdfinfo3.info3 class
    signalDataType : int
        signal type according to specification
    bitCount : int
        number of bits used to store channel record
    nBytes : int
        number of bytes (1 byte = 8 bits) taken by channel record
    dataFormat : str
        numpy dtype as string
    CFormat : struct class instance
        struct instance to convert from C Format
    byteOffset : int
        position of channel record in complete record in bytes
    bitOffset : int
        bit position of channel value inside byte in case of channel
        having bit count below 8
    recAttributeName : str
        channel name compliant to a valid python identifier
        (recarray attribute)
    RecordFormat : list of str
        dtype format used for numpy.core.records functions
        ((name_title,name),str_stype)
    channelType : int
        channel type
    posByteBeg : int
        start position in number of bit of channel record in complete record
    posByteEnd : int
        end position in number of bit of channel record in complete record
    bit_masking_needed : bool, default false
        True if bit masking needed after data read

    Methods
    ------------
    __init__(info, dataGroup, channelGroup, channelNumber, recordIDnumber)
        constructor
    __str__()
        to print class attributes
    """

    def __init__(self, info, dataGroup, channelGroup,
                 channelNumber, recordIDnumber):
        """ Channel class constructor

        Parameters
        ------------
        info : mdfinfo3.info3 class
        dataGroup : int
            data group number in mdfinfo3.info3 class
        channelGroup : int
            channel group number in mdfinfo3.info3 class
        channelNumber : int
            channel number in mdfinfo3.info3 class
        recordIDnumber : int
            Number of record IDs, each one Byte
        """
        self.name = info['CNBlock'][dataGroup][channelGroup][channelNumber]['signalName']
        self.channelNumber = channelNumber
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['signalDataType']
        if not self.signalDataType in (7, 8):
            numeric = True
        else:
            numeric = False
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfBits']
        ByteOrder = info['IDBlock']['ByteOrder']
        self.posBitBeg = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfTheFirstBits']
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.byteOffset = self.posBitBeg // 8
        self.bitOffset = self.posBitBeg % 8
        self.nBytes = _bits_to_bytes(self.bitCount + self.bitOffset, numeric)
        (self.dataFormat, self.nativedataFormat) = \
            _arrayformat3(self.signalDataType, self.nBytes, ByteOrder)
        self.CFormat = Struct(_datatypeformat3(self.signalDataType, self.nBytes, ByteOrder))
        self.embedding_channel_bitOffset = self.bitOffset  # for channel containing other channels
        self.posByteBeg = recordIDnumber + self.byteOffset
        self.posByteEnd = recordIDnumber + self.byteOffset + self.nBytes
        if not self.nBytes == self.bitCount / 8:
            self.bit_masking_needed = True
        else:
            self.bit_masking_needed = False
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['channelType']
        if 'physicalUnit' in info['CCBlock'][dataGroup][channelGroup][channelNumber]:
            self.unit = info['CCBlock'][dataGroup][channelGroup][channelNumber]['physicalUnit']
        else:
            self.unit = ''
        if 'signalDescription' in info['CNBlock'][dataGroup][channelGroup][channelNumber]:
            self.desc = info['CNBlock'][dataGroup][channelGroup][channelNumber]['signalDescription']
        else:
            self.desc = ''
        self.conversion = info['CCBlock'][dataGroup][channelGroup][channelNumber]

    def __str__(self):
        output = [self.channelNumber, ' ', self.name, ' ', self.signalDataType, ' ',
                  self.channelType, ' ', self.RecordFormat, ' ', self.bitOffset, ' ',
                  self.byteOffset, ' ', 'unit ', self.unit, 'description ', self.desc]
        return ''.join(output)

    def changeChannelName(self, channelGroup):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channelGroup : int
            channelGroup bumber
        """
        self.name = '{0}_{1}'.format(self.name, channelGroup)
        self.recAttributeName = _convertName(self.name)
        self.RecordFormat = (('{}_title'.format(self.recAttributeName), self.recAttributeName), self.dataFormat)


def _datatypeformat3(signalDataType, numberOfBytes, ByteOrder):
    """ function returning C format string from channel data type and number of bits

    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    dataType : str
        C format used by fread to read channel raw data
    """
    if signalDataType in (0, 9, 13):  # unsigned
        if numberOfBytes == 1:
            dataType = 'B'
        elif numberOfBytes == 2:
            dataType = 'H'
        elif numberOfBytes <= 4:
            dataType = 'I'
        elif numberOfBytes <= 8:
            dataType = 'Q'
        else:
            warn('Unsupported number of bits for unsigned int {}'.format(signalDataType))

    elif signalDataType in (1, 10, 14):  # signed int
        if numberOfBytes == 1:
            dataType = 'b'
        elif numberOfBytes == 2:
            dataType = 'h'
        elif numberOfBytes <= 4:
            dataType = 'i'
        elif numberOfBytes <= 8:
            dataType = 'q'
        else:
            warn('Unsupported number of bits for signed int {}'.format(signalDataType))

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBytes == 4:
            dataType = 'f'
        elif numberOfBytes == 8:
            dataType = 'd'
        else:
            warn('Unsupported number of bit for floating point {}'.format(signalDataType))

    elif signalDataType == 7:  # string
        dataType = str(numberOfBytes) + 's'
    elif signalDataType == 8:  # array of bytes
        dataType = str(numberOfBytes) + 's'
    else:
        warn('Unsupported Signal Data Type {0} nBits {1}'.format(signalDataType, numberOfBytes))

    # deal with byte order
    if signalDataType in (0, 1, 2, 3):
        if ByteOrder:
            dataType = '>{}'.format(dataType)
        else:
            dataType = '<{}'.format(dataType)
    elif signalDataType in (13, 14, 15, 16):  # low endian
        dataType = '<{}'.format(dataType)
    elif signalDataType in (9, 10, 11, 12):  # big endian
        dataType = '>{}'.format(dataType)

    return dataType


def _arrayformat3(signalDataType, numberOfBytes, ByteOrder):
    """ function returning numpy style string from channel data type and number of bits
    Parameters
    ----------------
    signalDataType : int
        channel data type according to specification
    numberOfBytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    dataType : str
        numpy dtype format used by numpy.core.records to read channel raw data
    """
    # Formats used by numpy

    if signalDataType in (0, 9, 13):  # unsigned
        if numberOfBytes == 1:
            dataType = 'u1'
        elif numberOfBytes == 2:
            dataType = 'u2'
        elif numberOfBytes <= 4:
            dataType = 'u4'
        elif numberOfBytes <= 8:
            dataType = 'u8'
        else:
            warn('Unsupported number of bits for unsigned int {} nBits '.format(signalDataType, numberOfBytes))

    elif signalDataType in (1, 10, 14):  # signed int
        if numberOfBytes == 1:
            dataType = 'i1'
        elif numberOfBytes == 2:
            dataType = 'i2'
        elif numberOfBytes <= 4:
            dataType = 'i4'
        elif numberOfBytes <= 8:
            dataType = 'i8'
        else:
            warn('Unsupported number of bits for signed int {0} nBits {1}'.format(signalDataType, numberOfBytes))

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBytes == 4:
            dataType = 'f4'
        elif numberOfBytes == 8:
            dataType = 'f8'
        else:
            warn('Unsupported number of bit for floating point {0} nBits {1}'.format(signalDataType, numberOfBytes))

    elif signalDataType == 7:  # string
        dataType = 'S{}'.format(numberOfBytes)  # not directly processed
    elif signalDataType == 8:  # array of bytes
        dataType = 'V{}'.format(numberOfBytes)  # not directly processed
    else:
        warn('Unsupported Signal Data Type {0} nBits {1}'.format(signalDataType, numberOfBytes))

    nativeDataType = dataType
    # deal with byte order
    if signalDataType in (0, 1, 2, 3):
        if ByteOrder:
            dataType = '>{}'.format(dataType)
        else:
            dataType = '<{}'.format(dataType)
    elif signalDataType in (13, 14, 15, 16):  # low endian
        dataType = '<{}'.format(dataType)
    elif signalDataType in (9, 10, 11, 12):  # big endian
        dataType = '>{}'.format(dataType)

    return (dataType, nativeDataType)
