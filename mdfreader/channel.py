# -*- coding: utf-8 -*-
""" Measured Data Format file reader module.

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Wed Oct 04 21:13:28 2017

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.4+

channel
--------------------------

"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from struct import Struct
from warnings import warn
from .mdfinfo4 import ATBlock
from .mdf import _bits_to_bytes_aligned, _bits_to_bytes_not_aligned, _convert_name

CAN_open_offset = {'ms': 0, 'days': 4, 'minute': 2, 'hour': 3, 'day': 4, 'month': 5, 'year': 6}


class Channel4(object):
    __slots__ = ['channelNumber', 'channelGroup', 'dataGroup',
                 'type', 'name', 'VLSD_CG_Flag', 'nBytes_aligned', 'byteOffset']
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
    nBytes_aligned : int
        number of bytes taken by channel at each sampling
    byteOffset : int
        byte offset in record
    bit_masking_needed : boolean


    Methods
    ------------
    __init__()
        constructor
    __str__()
        to print class attributes
    attachment(fid, info)
        in case of sync channel attached
    set(info, data_group, channel_group, channel_number, record_id_size)
        standard channel initialisation
    set_CANOpen(info, data_group, channel_group, channel_number,
                record_id_size, name)
        CANOpen channel initialisation
    set_invalid_bytes(info, data_group, channel_group, record_id_size, byte_aligned)
        Invalid Bytes channel initialisation
    rec_attribute_name : str
        Name of channel compatible with python attribute name conventions
    unit : str, default empty string
        channel unit
    desc : str
        channel description
    conversion : info class
        conversion dictionnary
    cn_block : info class
        Channel Block info class
    signal_data_type : int
        signal type according to specification
    bit_count : int
        number of bits used to store channel record
    calc_bytes : int
        number of bytes (1 byte = 8 bits) taken by channel record
    little_endian : Bool
        flag to inform of channel data endian
    data_format : str
        numpy dtype as string
    c_format :
        C format understood by fread
    c_format_structure : struct class instance
        struct instance to convert from C Format
    byte_offset : int
        position of channel record in complete record in bytes
    bit_offset : int
        bit position of channel value inside byte in case of channel
        having bit count below 8
    record_format : nested tuple of str
        dtype format used for numpy.core.records functions
        ((name_title,name),str_stype)
    native_record_format : nested tuple of str
        same as RecordFormat but using recAttributeName instead of name
    channel_type : int
        channel type ; 0 fixed length data, 1 VLSD, 2 master, 3 virtual master,
        4 sync, 5 MLSD, 6 virtual data
    channel_sync_type : int
        channel synchronisation type ; 0 None, 1 Time, 2 Angle,
        3 Distance, 4 Index
    pos_byte_beg : int
        start position in number of byte of channel record in complete record
    pos_byte_end : int
        end position in number of byte of channel record in complete record
    pos_bit_beg : int
        start position in number of bit of channel record in complete record
    pos_bit_end : int
        end position in number of bit of channel record in complete record
    ca_block : CABlock class
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
        nBytes : int, default 0
            number of bytes for this channel
        byteOffset : int, default
            position of channel in record, not including record ID
        """
        self.name = ''
        self.type = 0
        self.channelNumber = 0
        self.dataGroup = 0
        self.channelGroup = 0
        self.VLSD_CG_Flag = False
        self.nBytes_aligned = 0
        self.byteOffset = 0
        # self.bit_masking_needed = True

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
            return ATBlock(fid, info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data'])
        except KeyError:
            print('No Attachment block for this channel')

    def data(self, info):
        """ returns data block pointer for VLSD, MLD or sync channels"""
        try:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data']
        except KeyError:
            return None

    def cn_block(self, info):
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

    def signal_data_type(self, info, byte_aligned=True):
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
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type']
        else:
            if byte_aligned:
                return 10  # byte array
            else:
                return 0  # uint LE

    def bit_count(self, info):
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
                if self.signal_data_type(info) == 13:
                    return 16
                else:
                    return 32
            elif self.name == 'days':
                return 16
            else:
                return 8
        elif self.type == 4:  # Invalid bit channel
            return self.nBytes_aligned * 8
        else:
            warn('Not found channel type')

    def channel_sync_type(self, info):
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
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_sync_type']
        except KeyError:
            return 0  # in case of invalid bytes channel

    def ca_block(self, info):
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

    def is_ca_block(self, info):
        try:
            info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            return True
        except KeyError:
            return False

    def record_id_size(self, info):
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

    def channel_type(self, info):
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

    def calc_bytes(self, info, aligned=True):
        """ calculates channel aligned bytes number

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks
        aligned : boolean
            with or without aligned bytes

        Returns
        -----------
        number of bytes integer
        """
        if not self.type == 4:  # not channel containing invalid
            if aligned:
                n_bytes = _bits_to_bytes_aligned(
                    info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count']
                    + info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_offset'],
                    self.isnumeric(info))
            else:
                n_bytes = _bits_to_bytes_not_aligned(
                    info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count']
                    + info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_offset'])
            if self.type in (1, 2):  # array channel
                n_bytes *= self.ca_block(info)['PNd']
                block = self.ca_block(info)
                while 'CABlock' in block:  # nested array
                    block = block['CABlock']
                    n_bytes *= block['PNd']
            if self.type == 3:  # CAN channel
                if self.name == 'ms':
                    if info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'] == 13:
                        n_bytes = 2
                    else:
                        n_bytes = 4
                elif self.name == 'days':
                    n_bytes = 2
                else:
                    n_bytes = 1
            return n_bytes
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
        if info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'] in (1, 3, 5, 9):
            return False
        else:
            return True

    def record_attribute_name(self):
        """ clean up channel name from unauthorised characters

        Returns
        -----------
        channel name compliant to python attributes names (for recarray)
        """
        return _convert_name(self.name)

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
            data_format = '{}V'.format(self.nBytes_aligned)
        elif info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_composition'] and \
                'CABlock' in info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]:  # channel array
            ca_block = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            endian, data_format = array_format4(
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'],
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count'] // 8)
            # calculates total array size in bytes
            array_desc = ca_block['ca_dim_size']
            Block = ca_block
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
            data_format = array_desc + data_format
        elif self.type == 3:  # CAN channel
            if self.name == 'ms':
                if self.signal_data_type(info) == 13:
                    data_format = 'u2'
                else:
                    data_format = 'u4'
            elif self.name == 'days':
                data_format = 'u2'
            else:
                data_format = 'u1'
            endian = '<'
        else:  # not channel array
            endian, data_format = array_format4(
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_data_type'], self.nBytes_aligned)
        return endian, data_format

    def data_format(self, info):
        """ channel numpy.core.records data format

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data format
        """
        endian, _data_format = self.numpy_format(info)
        return ''.join([endian, _data_format])

    def native_data_format(self, info):
        endian, _native_data_format = self.numpy_format(info)
        return _native_data_format

    def c_format(self, info):
        """ channel data C format

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data C format
        """
        signal_data_type = self.signal_data_type(info)
        if self.type == 0:  # standard channel
            if signal_data_type not in (13, 14):
                if not self.channel_type(info) == 1:  # if not VSLD
                    endian, data_type = data_type_format4(signal_data_type, self.nBytes_aligned)
                else:  # VLSD
                    endian, data_type = data_type_format4(0, self.nBytes_aligned)
                return '{}{}'.format(endian, data_type)
        elif self.type in (1, 2):  # array channel
            ca = self.ca_block(info)
            n_bytes = _bits_to_bytes_aligned(
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_count'] +
                info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_bit_offset'],
                self.isnumeric(info))
            endian, data_type = data_type_format4(signal_data_type, n_bytes)
            return '{}{}{}'.format(endian, ca['PNd'], data_type)
        elif self.type == 3:  # CAN channel
            if self.name == 'ms':
                if signal_data_type == 13:
                    return 'H'
                else:
                    return 'I'
            elif self.name == 'days':
                return 'H'
            else:
                return 'B'
        elif self.type == 4:  # Invalid bit channel
            return '{}s'.format(self.nBytes_aligned)
        else:
            warn('Not found channel type')

    def c_format_structure(self, info):
        """ channel data C format struct object

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        string data C format struct object
        """
        return Struct(self.c_format(info))

    def CANOpen_offset(self):
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

    def bit_offset(self, info):
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
                   + self.CANOpen_offset() * 8
        elif self.type == 4:  # Invalid bit channel
            return 0
        else:
            warn('Not found channel type')

    def calc_byte_offset(self, info):
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
                   + self.CANOpen_offset()
        elif self.type == 4:  # Invalid bit channel
            return info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        else:
            warn('Not found channel type')

    def pos_byte_beg(self, info):
        """ channel data bytes starting position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes starting position
        """
        return self.record_id_size(info) + self.byteOffset

    def pos_byte_end(self, info):
        """ channel data bytes ending position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bytes ending position
        """
        return self.pos_byte_beg(info) + self.nBytes_aligned

    def pos_bit_beg(self, info):
        """ channel data bit starting position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bit starting position
        """
        return self.pos_byte_beg(info) * 8 + self.bit_offset(info)

    def pos_bit_end(self, info):
        """ channel data bit ending position in record

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        integer, channel bit ending position
        """
        return self.pos_bit_beg(info) + self.bit_count(info)

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
                    desc = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['Comment']
                    if (desc is not None) and isinstance(desc, dict):
                        if 'description' in desc:
                            desc = desc['description']
                        elif 'name' in desc:
                            desc = desc['name']
                else:
                    desc = ''
                return desc
            else:
                return 'Invalid Bytes DG group {0} CG group {1}'.format(self.dataGroup,
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

    def set(self, info, data_group, channel_group, channel_number):
        """ channel initialisation

        Parameters
        ------------

        info : mdfinfo4.info4 class
        data_group : int
            data group number in mdfinfo4.info4 class
        channel_group : int
            channel group number in mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class
        """
        self.name = info['CN'][data_group][channel_group][channel_number]['name']
        self.channelNumber = channel_number
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.type = 0
        if info['CN'][data_group][channel_group][channel_number]['cn_composition'] and \
                'CABlock' in info['CN'][data_group][channel_group][channel_number]:
            # channel array
            self.type = 1
            block = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['CABlock']
            if 'CABlock' in block:  # nested array
                self.type = 2
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)

    def set_CANOpen(self, info, data_group, channel_group, channel_number, name):
        """ CANOpen channel intialisation

        Parameters
        ------------

        info : mdfinfo4.info4 class
        data_group : int
            data group number in mdfinfo4.info4 class
        channel_group : int
            channel group number in mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class
        name : str
            name of channel. Should be in ('ms', 'day', 'days', 'hour',
            'month', 'minute', 'year')
        """
        self.type = 3
        self.name = name
        self.channelNumber = channel_number
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)

    def set_invalid_bytes(self, info, data_group, channel_group, channel_number):
        """ invalid_bytes channel initialisation

        Parameters
        ----------

        info : mdfinfo4.info4 class
        data_group : int
            data group number in mdfinfo4.info4 class
        channel_group : int
            channel group number in mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class
        """
        self.type = 4
        self.name = 'invalid_bytes{}'.format(data_group)
        self.channelNumber = channel_number
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)

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

    def change_channel_name(self, channel_group):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channel_group : int
            channelGroup number
        """
        self.name = '{0}_{1}'.format(self.name, channel_group)

    def bit_masking_need(self, info):
        """ Valid if bit masking need

        Parameters
        ----------------

        info : mdfinfo4.info4 class
            info4 class containing all MDF Blocks

        Returns
        -----------
        boolean True if channel needs bit masking, otherwise False

        """
        if not self.nBytes_aligned == self.bit_count(info) / 8:
            self.bit_masking_needed = True
        else:
            self.bit_masking_needed = False


def array_format4(signal_data_type, number_of_bytes):
    """ function returning numpy style string from channel data type and number of bits

    Parameters
    ----------------
    signal_data_type : int
        channel data type according to specification
    number_of_bytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    endian, data_type : str
        numpy dtype format used by numpy.core.records to read channel raw data
    """

    if signal_data_type == 0:  # unsigned, low endian
        if number_of_bytes == 1:
            data_type = 'u1'
        elif number_of_bytes == 2:
            data_type = 'u2'
        elif number_of_bytes <= 4:
            data_type = 'u4'
        elif number_of_bytes <= 8:
            data_type = 'u8'
        else:
            data_type = '{}V'.format(number_of_bytes)
        endian = '<'

    elif signal_data_type == 2:  # signed int, low endian
        if number_of_bytes == 1:
            data_type = 'i1'
        elif number_of_bytes == 2:
            data_type = 'i2'
        elif number_of_bytes <= 4:
            data_type = 'i4'
        elif number_of_bytes <= 8:
            data_type = 'i8'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(number_of_bytes))
        endian = '<'

    elif signal_data_type == 4:  # floating point, low endian
        if number_of_bytes == 4:
            data_type = 'f4'
        elif number_of_bytes == 8:
            data_type = 'f8'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(number_of_bytes))
        endian = '<'

    elif signal_data_type == 1:  # unsigned, big endian
        if number_of_bytes == 1:
            data_type = 'u1'
        elif number_of_bytes == 2:
            data_type = 'u2'
        elif number_of_bytes <= 4:
            data_type = 'u4'
        elif number_of_bytes <= 8:
            data_type = 'u8'
        else:
            data_type = '{}V'.format(number_of_bytes)
        endian = '>'

    elif signal_data_type == 3:  # signed int, big endian
        if number_of_bytes == 1:
            data_type = 'i1'
        elif number_of_bytes == 2:
            data_type = 'i2'
        elif number_of_bytes <= 4:
            data_type = 'i4'
        elif number_of_bytes <= 8:
            data_type = 'i8'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(number_of_bytes))
        endian = '>'

    elif signal_data_type == 5:  # floating point, big endian
        if number_of_bytes == 4:
            data_type = 'f4'
        elif number_of_bytes == 8:
            data_type = 'f8'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(number_of_bytes))
        endian = '>'

    elif signal_data_type == 6:  # string ISO-8859-1 Latin
        data_type = 'S{}'.format(number_of_bytes)
        endian = ''
    elif signal_data_type == 7:  # UTF-8
        data_type = 'S{}'.format(number_of_bytes)
        endian = ''
    elif signal_data_type == 8:  # UTF-16 low endian
        data_type = 'S{}'.format(number_of_bytes)
        endian = '<'
    elif signal_data_type == 9:  # UTF-16 big endian
        data_type = 'S{}'.format(number_of_bytes)
        endian = '>'
    elif signal_data_type == 10:  # bytes array
        data_type = 'V{}'.format(number_of_bytes)
        endian = ''
    elif signal_data_type in (11, 12):  # MIME sample or MIME stream
        data_type = 'V{}'.format(number_of_bytes)
        endian = ''
    elif signal_data_type in (13, 14):  # CANOpen date or time
        data_type = ''
        endian = ''
    else:
        warn('Unsupported Signal Data Type {} {}'.format(signal_data_type, number_of_bytes))

    return endian, data_type


def data_type_format4(signal_data_type, number_of_bytes):
    """ function returning C format string from channel data type and number of bits

    Parameters
    ----------------
    signal_data_type : int
        channel data type according to specification
    number_of_bytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    data_type : str
        C format used by fread to read channel raw data
    """

    if signal_data_type == 0:  # unsigned int
        if number_of_bytes == 1:
            data_type = 'B'
        elif number_of_bytes == 2:
            data_type = 'H'
        elif number_of_bytes <= 4:
            data_type = 'I'
        elif number_of_bytes <= 8:
            data_type = 'Q'
        else:
            data_type = '{}s'.format(number_of_bytes)
        endian = '<'

    elif signal_data_type == 1:  # unsigned int
        if number_of_bytes == 1:
            data_type = 'B'
        elif number_of_bytes == 2:
            data_type = 'H'
        elif number_of_bytes <= 4:
            data_type = 'I'
        elif number_of_bytes <= 8:
            data_type = 'Q'
        else:
            data_type = '{}s'.format(number_of_bytes)
        endian = '>'

    elif signal_data_type == 2:  # signed int
        if number_of_bytes == 1:
            data_type = 'b'
        elif number_of_bytes == 2:
            data_type = 'h'
        elif number_of_bytes <= 4:
            data_type = 'i'
        elif number_of_bytes <= 8:
            data_type = 'q'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(signal_data_type))
        endian = '<'

    elif signal_data_type == 3:  # signed int
        if number_of_bytes == 1:
            data_type = 'b'
        elif number_of_bytes == 2:
            data_type = 'h'
        elif number_of_bytes <= 4:
            data_type = 'i'
        elif number_of_bytes <= 8:
            data_type = 'q'
        else:
            warn('Unsupported number of bytes for signed int {}'.format(signal_data_type))
        endian = '>'

    elif signal_data_type == 4:  # floating point
        if number_of_bytes == 4:
            data_type = 'f'
        elif number_of_bytes == 8:
            data_type = 'd'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(signal_data_type))
        endian = '<'

    elif signal_data_type == 5:  # floating point
        if number_of_bytes == 4:
            data_type = 'f'
        elif number_of_bytes == 8:
            data_type = 'd'
        else:
            warn('Unsupported number of bytes for floating point {}'.format(signal_data_type))
        endian = '>'

    elif signal_data_type in (6, 7, 10, 11, 12):  # string/bytes
        data_type = '{}s'.format(number_of_bytes)
        endian = ''

    elif signal_data_type == 8:  # UTF16 string/bytes
        data_type = '{}s'.format(number_of_bytes)
        endian = '<'

    elif signal_data_type == 9:  # UTF16 string/bytes
        data_type = '{}s'.format(number_of_bytes)
        endian = '>'

    else:
        warn('Unsupported Signal Data Type {} {}'.format(signal_data_type, number_of_bytes))

    return endian, data_type


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
    nBytes_aligned : int
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
    change_channel_name(channel_group)
        rename duplicated channel name within unsorted channel groups
    """

    def __init__(self, info, data_group, channel_group, channel_number, record_id_number):
        """ Channel class constructor

        Parameters
        ------------
        info : mdfinfo3.info3 class
        data_group : int
            data group number in mdfinfo3.info3 class
        channel_group : int
            channel group number in mdfinfo3.info3 class
        channel_number : int
            channel number in mdfinfo3.info3 class
        record_id_number : int
            Number of record ID, each one Byte
        """
        self.name = info['CNBlock'][data_group][channel_group][channel_number]['signalName']
        self.channelNumber = channel_number
        self.signalDataType = info['CNBlock'][data_group][channel_group][channel_number]['signalDataType']
        if self.signalDataType not in (7, 8):
            numeric = True
        else:
            numeric = False
        self.bitCount = info['CNBlock'][data_group][channel_group][channel_number]['numberOfBits']
        byte_order = info['IDBlock']['ByteOrder']
        self.posBitBeg = info['CNBlock'][data_group][channel_group][channel_number]['numberOfTheFirstBits']
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.byteOffset = self.posBitBeg // 8
        self.bitOffset = self.posBitBeg % 8
        self.nBytes_aligned = _bits_to_bytes_aligned(self.bitCount + self.bitOffset, numeric)
        self.nBytes_not_aligned = _bits_to_bytes_not_aligned(self.bitCount + self.bitOffset)
        (self.dataFormat, self.nativedataFormat) = \
            _array_format3(self.signalDataType, self.nBytes_aligned, byte_order)
        self.CFormat = Struct(_data_type_format3(self.signalDataType, self.nBytes_aligned, byte_order))
        self.embedding_channel_bitOffset = self.bitOffset  # for channel containing other channels
        self.posByteBeg = record_id_number + self.byteOffset
        self.posByteEnd = record_id_number + self.byteOffset + self.nBytes_aligned
        if not self.nBytes_aligned == self.bitCount / 8:
            self.bit_masking_needed = True
        else:
            self.bit_masking_needed = False
        self.channelType = info['CNBlock'][data_group][channel_group][channel_number]['channelType']
        if 'physicalUnit' in info['CCBlock'][data_group][channel_group][channel_number]:
            self.unit = info['CCBlock'][data_group][channel_group][channel_number]['physicalUnit']
        else:
            self.unit = ''
        if 'signalDescription' in info['CNBlock'][data_group][channel_group][channel_number]:
            self.desc = info['CNBlock'][data_group][channel_group][channel_number]['signalDescription']
        else:
            self.desc = ''
        self.conversion = info['CCBlock'][data_group][channel_group][channel_number]

    def __str__(self):
        output = [self.channelNumber, ' ', self.name, ' ', self.signalDataType, ' ',
                  self.channelType, ' ', self.RecordFormat, ' ', self.bitOffset, ' ',
                  self.byteOffset, ' ', 'unit ', self.unit, 'description ', self.desc]
        return ''.join(output)

    def change_channel_name(self, channel_group):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channel_group : int
            channelGroup number
        """
        self.name = '{0}_{1}'.format(self.name, channel_group)
        self.recAttributeName = _convert_name(self.name)
        self.RecordFormat = (('{}_title'.format(self.recAttributeName), self.recAttributeName), self.dataFormat)


def _data_type_format3(signal_data_type, number_of_bytes, byte_order):
    """ function returning C format string from channel data type and number of bits

    Parameters
    ----------------
    signal_data_type : int
        channel data type according to specification
    number_of_bytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    dataType : str
        C format used by fread to read channel raw data
    """
    if signal_data_type in (0, 9, 13):  # unsigned
        if number_of_bytes == 1:
            data_type = 'B'
        elif number_of_bytes == 2:
            data_type = 'H'
        elif number_of_bytes <= 4:
            data_type = 'I'
        elif number_of_bytes <= 8:
            data_type = 'Q'
        else:
            warn('Unsupported number of bits for unsigned int {}'.format(signal_data_type))

    elif signal_data_type in (1, 10, 14):  # signed int
        if number_of_bytes == 1:
            data_type = 'b'
        elif number_of_bytes == 2:
            data_type = 'h'
        elif number_of_bytes <= 4:
            data_type = 'i'
        elif number_of_bytes <= 8:
            data_type = 'q'
        else:
            warn('Unsupported number of bits for signed int {}'.format(signal_data_type))

    elif signal_data_type in (2, 3, 11, 12, 15, 16):  # floating point
        if number_of_bytes == 4:
            data_type = 'f'
        elif number_of_bytes == 8:
            data_type = 'd'
        else:
            warn('Unsupported number of bit for floating point {}'.format(signal_data_type))

    elif signal_data_type == 7:  # string
        data_type = str(number_of_bytes) + 's'
    elif signal_data_type == 8:  # array of bytes
        data_type = str(number_of_bytes) + 's'
    else:
        warn('Unsupported Signal Data Type {0} nBits {1}'.format(signal_data_type, number_of_bytes))

    # deal with byte order
    if signal_data_type in (0, 1, 2, 3):
        if byte_order:
            data_type = '>{}'.format(data_type)
        else:
            data_type = '<{}'.format(data_type)
    elif signal_data_type in (13, 14, 15, 16):  # low endian
        data_type = '<{}'.format(data_type)
    elif signal_data_type in (9, 10, 11, 12):  # big endian
        data_type = '>{}'.format(data_type)

    return data_type


def _array_format3(signal_data_type, number_of_bytes, byte_order):
    """ function returning numpy style string from channel data type and number of bits
    Parameters
    ----------------
    signal_data_type : int
        channel data type according to specification
    number_of_bytes : int
        number of bytes taken by channel data in a record

    Returns
    -----------
    dataType : str
        numpy dtype format used by numpy.core.records to read channel raw data
    """
    # Formats used by numpy

    if signal_data_type in (0, 9, 13):  # unsigned
        if number_of_bytes == 1:
            data_type = 'u1'
        elif number_of_bytes == 2:
            data_type = 'u2'
        elif number_of_bytes <= 4:
            data_type = 'u4'
        elif number_of_bytes <= 8:
            data_type = 'u8'
        else:
            warn('Unsupported number of bits for unsigned int {} nBits '.format(signal_data_type, number_of_bytes))

    elif signal_data_type in (1, 10, 14):  # signed int
        if number_of_bytes == 1:
            data_type = 'i1'
        elif number_of_bytes == 2:
            data_type = 'i2'
        elif number_of_bytes <= 4:
            data_type = 'i4'
        elif number_of_bytes <= 8:
            data_type = 'i8'
        else:
            warn('Unsupported number of bits for signed int {0} nBits {1}'.format(signal_data_type, number_of_bytes))

    elif signal_data_type in (2, 3, 11, 12, 15, 16):  # floating point
        if number_of_bytes == 4:
            data_type = 'f4'
        elif number_of_bytes == 8:
            data_type = 'f8'
        else:
            warn('Unsupported number of bit for floating point {0} nBits {1}'.format(signal_data_type, number_of_bytes))

    elif signal_data_type == 7:  # string
        data_type = 'S{}'.format(number_of_bytes)  # not directly processed
    elif signal_data_type == 8:  # array of bytes
        data_type = 'V{}'.format(number_of_bytes)  # not directly processed
    else:
        warn('Unsupported Signal Data Type {0} nBits {1}'.format(signal_data_type, number_of_bytes))

    native_data_type = data_type
    # deal with byte order
    if signal_data_type in (0, 1, 2, 3):
        if byte_order:
            data_type = '>{}'.format(data_type)
        else:
            data_type = '<{}'.format(data_type)
    elif signal_data_type in (13, 14, 15, 16):  # low endian
        data_type = '<{}'.format(data_type)
    elif signal_data_type in (9, 10, 11, 12):  # big endian
        data_type = '>{}'.format(data_type)

    return (data_type, native_data_type)
