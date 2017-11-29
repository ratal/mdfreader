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
from struct import Struct
from sys import path, stderr
from os.path import dirname, abspath
_root = dirname(abspath(__file__))
path.append(_root)
from mdfinfo4 import ATBlock
from mdf import _bits_to_bytes, _convertName
from numpy import right_shift, bitwise_and


class channel4(object):
    __slots__ = ['channelNumber', 'channelGroup', 'dataGroup',
                 'type', 'name', 'VLSD_CG_Flag']
    """ channel class gathers all about channel structure in a record

    Attributes
    --------------
    name : str
        Name of channel
    type : str
        channel type. Can be 'std', 'NestCA',
        'CAN' or 'Inv'
    channelNumber : int
        channel number corresponding to mdfinfo4.info4 class
    channelGroup : int
        channel group number corresponding to mdfinfo4.info4 class
    dataGroup : int
        data group number corresponding to mdfinfo4.info4 class
    VLSD_CG_Flag : bool
        flag when Channel Group VLSD is used

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
    nBytes : int
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
    """

    def __init__(self):
        """ channel class constructor

        Attributes
        --------------

        name : str
            Name of channel
        type : str, default 'std'
            channel type. Can be 'std', 'NestCA',
            'CAN' or 'Inv'
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
        self.type = 'std'
        self.channelNumber = 0
        self.dataGroup = 0
        self.channelGroup = 0
        self.VLSD_CG_Flag = False

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
        if not self.type == 'Inv':
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
        if self.type in ('std', 'CA', 'NestCA'):
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_bit_count']
        elif self.type == 'CAN':
            if self.name == 'ms':
                if self.signalDataType(info) == 13:
                    return 16
                else:
                    return 32
            elif self.name == 'days':
                return 16
            else:
                return 8
        elif self.type == 'Inv':
            return self.nBytes(info) * 8
        else:
            print('Not found channel type')

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
            return 0  # in case of invaldi bytes channel

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
            print('No CA Block in this channel')

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
            return 0  # in case of invaldi bytes channel

    def nBytes(self, info):
        """ calculates channel bytes number

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        number of bytes integer
        """
        if not self.type == 'Inv':
            nBytes = _bits_to_bytes(info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_bit_count'])
            if self.type in ('CA', 'NestCA'):
                nBytes *= self.CABlock(info)['PNd']
                Block = self.CABlock(info)
                while 'CABlock' in Block:  # nested array
                    Block = Block['CABlock']
                    nBytes *= Block['PNd']
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
        if self.type == 'Inv':
            dataformat = '{}V'.format(self.nBytes(info))
        elif info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_composition'] and \
                'CABlock' in info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]:  # channel array
            CABlock = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            dataformat = arrayformat4(info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'],
                                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count'])
            # calculates total array size in bytes
            array_desc = CABlock['ca_dim_size']
            Block = CABlock
            while 'CABlock' in Block:  # nested array
                Block = Block['CABlock']
                Block['ca_dim_size']
                if isinstance(array_desc, list):
                    array_desc.append(Block['ca_dim_size'])
                else:
                    array_desc = [array_desc, Block['ca_dim_size']]
            if isinstance(array_desc, list):
                array_desc = str(tuple(array_desc))
            else:
                array_desc = str(array_desc)
            dataformat = array_desc + dataformat
        elif self.type == 'CAN':
            if self.name == 'ms':
                if self.signalDataType(info) == 13:
                    dataformat = '<u2'
                else:
                    dataformat = '<u4'
            elif self.name == 'days':
                dataformat = '<u2'
            else:
                dataformat = '<u1'
        else:  # not channel array
            dataformat = arrayformat4(info['CN'][self.dataGroup][self.channelGroup]\
                            [self.channelNumber]['cn_data_type'], info['CN'][self.dataGroup]\
                            [self.channelGroup][self.channelNumber]['cn_bit_count'])
        return dataformat

    def nativedataFormat(self, info):
        return self.dataFormat(info).lstrip('<').lstrip('>')

    def RecordFormat(self, info, recAttributeName=None):
        """ recarray dtype string

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        recarray dtype
        """
        if recAttributeName is None:
            recAttributeName = self.recAttributeName(info)
        return ((u'{}_title'.format(recAttributeName),
                 recAttributeName), self.dataFormat(info))

    def nativeRecordFormat(self, info):
        """ channel numpy.core.records data format without endian (< or >)

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data format
        """
        recAttributeName = self.recAttributeName(info)
        return ((u'{}_title'.format(recAttributeName), recAttributeName),
                self.dataFormat(info).lstrip('<').lstrip('>'))

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
        if self.type in ('std', 'CA', 'NestCA'):
            signalDataType = self.signalDataType(info)
            if signalDataType not in (13, 14):
                if not self.channelType(info) == 1:  # if not VSLD
                    return datatypeformat4(signalDataType, self.bitCount(info))
                else:  # VLSD
                    return datatypeformat4(0, self.bitCount(info))
        elif self.type == 'CAN':
            if self.name == 'ms':
                if self.signalDataType == 13:
                    return 'H'
                else:
                    return 'I'
            elif self.name == 'days':
                return 'H'
            else:
                return 'B'
        elif self.type == 'Inv':
            return '{}s'.format(self.nBytes(info))
        else:
            print('Not found channel type')

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

    def CANOpenOffset(self, info):
        """ CANopen channel bytes offset

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes offset
        """
        if self.name == 'ms':
            return 0
        elif self.name == 'days':
            return 4
        else:
            if self.name == 'minute':
                return 2
            elif self.name == 'hour':
                return 3
            elif self.name == 'day':
                return 4
            elif self.name == 'month':
                return 5
            elif self.name == 'year':
                return 6
            else:
                print('CANopen type not understood')

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
        if self.type in ('std', 'CA', 'NestCA'):
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_bit_offset']
        elif self.type == 'CAN':
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_bit_offset'] + self.CANOpenOffset(info)*8
        elif self.type == 'Inv':
            return 0
        else:
            print('Not found channel type')

    def byteOffset(self, info):
        """ channel data bytes offset in record (without record id)

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes offset
        """
        if self.type in ('std', 'CA', 'NestCA'):
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_byte_offset']
        elif self.type == 'CAN':
            return info['CN'][self.dataGroup][self.channelGroup]\
                        [self.channelNumber]['cn_byte_offset'] + self.CANOpenOffset(info)
        elif self.type == 'Inv':
            return info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        else:
            print('Not found channel type')

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
        return self.recordIDsize(info) + self.byteOffset(info)

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
        return self.posByteBeg(info) + self.nBytes(info)

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
        if not self.type == 'CAN':
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

    def invalid_bit(self, info):
        """ extrzcts from info4 the channels valid bits positions

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        dict of channels valid bits positions
        """
        invalid_bit = {}
        for channelNumber in info['CN'][self.dataGroup][self.channelGroup]:
            name = info['CN'][self.dataGroup][self.channelGroup][channelNumber]['name']
            invalid_bit[name] = info['CN'][self.dataGroup][self.channelGroup]\
                                    [channelNumber]['cn_invalid_bit_pos']
        return invalid_bit

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
        recordIDsize : int
            size of record ID in Bytes
        """
        self.name = info['CN'][dataGroup][channelGroup][channelNumber]['name']
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.type = 'std'
        if info['CN'][dataGroup][channelGroup][channelNumber]['cn_composition'] and \
                'CABlock' in info['CN'][dataGroup][channelGroup][channelNumber]:
            # channel array
            self.type = 'CA'
            Block = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            if 'CABlock' in Block:  # nested array
                self.type = 'NestCA'

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
        recordIDsize : int
            size of record ID in Bytes
        name : str
            name of channel. Should be in ('ms', 'day', 'days', 'hour',
            'month', 'minute', 'year')
        """
        self.type = 'CAN'
        self.name = name
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup

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
        recordIDsize : int
            size of record ID in Bytes
        byte_aligned : Bool
            Flag for byte alignement
        """
        self.type = 'Inv'
        self.name = 'invalid_bytes{}'.format(dataGroup)
        self.channelNumber = channelNumber
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup

    def validity_channel(self, info, invalid_bytes):
        """ extract channel validity bits

        Parameters
        ----------

        info : mdfinfo4.info4 class
        invalid_bytes : bytes
            bytes from where to extract validity bit array

        Return
        -------
        Numpy vector of validity

        """
        if self.type == 'Inv':
            return bitwise_and(right_shift(invalid_bytes, self.invalid_bit(info)), 1)
        else:
            print('asking for invalid byte array but channel is not invalid byte type')
    
    def changeChannelName(self, channelGroup):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channelGroup : int
            channelGroup bumber
        """
        self.name = '{0}_{1}'.format(self.name, channelGroup)


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

    if signalDataType == 0:  # unsigned, low endian
        if numberOfBits <= 8:
            dataType = '<u1'
        elif numberOfBits <= 16:
            dataType = '<u2'
        elif numberOfBits <= 32:
            dataType = '<u4'
        elif numberOfBits <= 64:
            dataType = '<u8'
        else:
            dataType = '<{}V'.format(_bits_to_bytes(numberOfBits) // 8)

    elif signalDataType == 2:  # signed int, low endian
        if numberOfBits <= 8:
            dataType = '<i1'
        elif numberOfBits <= 16:
            dataType = '<i2'
        elif numberOfBits <= 32:
            dataType = '<i4'
        elif numberOfBits <= 64:
            dataType = '<i8'
        else:
            print('Unsupported number of bits for signed int {}'.format(numberOfBits))

    elif signalDataType == 4:  # floating point, low endian
        if numberOfBits == 32:
            dataType = '<f4'
        elif numberOfBits == 64:
            dataType = '<f8'
        else:
            print('Unsupported number of bit for floating point ' + str(numberOfBits))

    elif signalDataType == 1:  # unsigned, big endian
        if numberOfBits <= 8:
            dataType = '>u1'
        elif numberOfBits <= 16:
            dataType = '>u2'
        elif numberOfBits <= 32:
            dataType = '>u4'
        elif numberOfBits <= 64:
            dataType = '>u8'
        else:
            dataType = '>{}V'.format(_bits_to_bytes(numberOfBits) // 8)

    elif signalDataType == 3:  # signed int, big endian
        if numberOfBits <= 8:
            dataType = '>i1'
        elif numberOfBits <= 16:
            dataType = '>i2'
        elif numberOfBits <= 32:
            dataType = '>i4'
        elif numberOfBits <= 64:
            dataType = '>i8'
        else:
            print('Unsupported number of bits for signed int {}'.format(numberOfBits))

    elif signalDataType == 5:  # floating point, big endian
        if numberOfBits == 32:
            dataType = '>f4'
        elif numberOfBits == 64:
            dataType = '>f8'
        else:
            print('Unsupported number of bit for floating point ' + str(numberOfBits))

    elif signalDataType == 6:  # string ISO-8859-1 Latin
        dataType = 'S{}'.format(numberOfBits // 8)
    elif signalDataType == 7:  # UTF-8
        dataType = 'S{}'.format(numberOfBits // 8)
    elif signalDataType == 8:  # UTF-16 low endian
        dataType = '<S{}'.format(numberOfBits // 8)
    elif signalDataType == 9:  # UTF-16 big endian
        dataType = '>S{}'.format(numberOfBits // 8)
    elif signalDataType == 10:  # bytes array
        dataType = 'V{}'.format(numberOfBits // 8)
    elif signalDataType in (11, 12):  # MIME sample or MIME stream
        dataType = 'V{}'.format(int(numberOfBits / 8))
    elif signalDataType in (13, 14):  # CANOpen date or time
        dataType = None
    else:
        print('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits)

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
            dataType = '{}s'.format(_bits_to_bytes(numberOfBits) // 8)

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
        dataType = '{}s'.format(numberOfBits // 8)
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits))
    # deal with byte order
    if signalDataType in (0, 2, 4, 8):  # low endian
        dataType = '<{}'.format(dataType)
    elif signalDataType in (1, 3, 5, 9):  # big endian
        dataType = '>{}'.format(dataType)

    return dataType


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
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfBits']
        ByteOrder = info['IDBlock']['ByteOrder']
        (self.dataFormat, self.nativedataFormat) = \
            _arrayformat3(self.signalDataType, self.bitCount, ByteOrder)
        self.CFormat = Struct(_datatypeformat3(self.signalDataType, self.bitCount, ByteOrder))
        self.nBytes = _bits_to_bytes(self.bitCount)
        self.posBitBeg = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfTheFirstBits']
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.byteOffset = self.posBitBeg // 8  # + info['CNBlock'][dataGroup][channelGroup][channelNumber]['ByteOffset']
        self.bitOffset = self.posBitBeg % 8
        self.embedding_channel_bitOffset = self.bitOffset  # for channel containing other channels
        self.posByteBeg = recordIDnumber + self.byteOffset
        self.posByteEnd = recordIDnumber + self.byteOffset + self.nBytes
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


def _datatypeformat3(signalDataType, numberOfBits, ByteOrder):
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
    if signalDataType in (0, 9, 13):  # unsigned
        if numberOfBits <= 8:
            dataType = 'B'
        elif numberOfBits <= 16:
            dataType = 'H'
        elif numberOfBits <= 32:
            dataType = 'I'
        elif numberOfBits <= 64:
            dataType = 'Q'
        else:
            print('Unsupported number of bits for unsigned int {}'.format(signalDataType))

    elif signalDataType in (1, 10, 14):  # signed int
        if numberOfBits <= 8:
            dataType = 'b'
        elif numberOfBits <= 16:
            dataType = 'h'
        elif numberOfBits <= 32:
            dataType = 'i'
        elif numberOfBits <= 64:
            dataType = 'q'
        else:
            print('Unsupported number of bits for signed int {}'.format(signalDataType))

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBits == 32:
            dataType = 'f'
        elif numberOfBits == 64:
            dataType = 'd'
        else:
            print('Unsupported number of bit for floating point {}'.format(signalDataType))

    elif signalDataType == 7:  # string
        dataType = str(numberOfBits // 8) + 's'
    elif signalDataType == 8:  # array of bytes
        dataType = str(numberOfBits // 8) + 's'
    else:
        print('Unsupported Signal Data Type {0} nBits {1}'.format(signalDataType, numberOfBits))

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


def _arrayformat3(signalDataType, numberOfBits, ByteOrder):
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
    # Formats used by numpy

    if signalDataType in (0, 9, 13):  # unsigned
        if numberOfBits <= 8:
            dataType = 'u1'
        elif numberOfBits <= 16:
            dataType = 'u2'
        elif numberOfBits <= 32:
            dataType = 'u4'
        elif numberOfBits <= 64:
            dataType = 'u8'
        else:
            print('Unsupported number of bits for unsigned int {} nBits '.format(signalDataType, numberOfBits))

    elif signalDataType in (1, 10, 14):  # signed int
        if numberOfBits <= 8:
            dataType = 'i1'
        elif numberOfBits <= 16:
            dataType = 'i2'
        elif numberOfBits <= 32:
            dataType = 'i4'
        elif numberOfBits <= 64:
            dataType = 'i8'
        else:
            print('Unsupported number of bits for signed int {0} nBits {1}'.format(signalDataType, numberOfBits))

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBits == 32:
            dataType = 'f4'
        elif numberOfBits == 64:
            dataType = 'f8'
        else:
            print('Unsupported number of bit for floating point {0} nBits {1}'.format(signalDataType, numberOfBits))

    elif signalDataType == 7:  # string
        dataType = 'S{}'.format(numberOfBits // 8)  # not directly processed
    elif signalDataType == 8:  # array of bytes
        dataType = 'V{}'.format(numberOfBits // 8)  # not directly processed
    else:
        print('Unsupported Signal Data Type {0} nBits {1}'.format(signalDataType, numberOfBits))

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
