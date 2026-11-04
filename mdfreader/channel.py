# -*- coding: utf-8 -*-
""" Measured Data Format file reader module.

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Wed Oct 04 21:13:28 2017

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 3.4+

channel
--------------------------

"""
from struct import Struct
from warnings import warn
from .mdfinfo4 import ATBlock
from .mdf import _bits_to_bytes_aligned, _bits_to_bytes_not_aligned, _convert_name

CAN_open_offset = {'ms': 0, 'days': 4, 'minute': 2, 'hour': 3, 'day': 4, 'month': 5, 'year': 6}


class Channel4(object):
    __slots__ = ['channelNumber', 'channelGroup', 'dataGroup',
                 'type', 'name', 'VLSD_CG_Flag', 'nBytes_aligned', 'byteOffset', 'pos_bit_beg',
                 '_cn_info', '_cn_type', '_cn_data_type', '_cn_bit_count',
                 '_cn_bit_offset', '_cn_byte_offset', '_cn_flags', '_cn_composition']
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

    def __init__(self, data_group, channel_group, channel_number):
        """ channel class constructor

        Parameters
        ------------

        data_group : int
            data group number in mdfinfo4.info4 class
        channel_group : int
            channel group number in mdfinfo4.info4 class
        channel_number : int
            channel number in mdfinfo4.info4 class

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
        pos_bit_beg : int, default 0
            bit position in record
        """
        self.name = ''
        self.type = 0
        self.channelNumber = channel_number
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.VLSD_CG_Flag = False
        self.nBytes_aligned = 0
        self.byteOffset = 0
        self.pos_bit_beg = 0

    def __str__(self):
        """ channel object attributes print

        """
        return f'{self.name} {self.type} {self.dataGroup} {self.channelGroup} {self.channelNumber} {self.VLSD_CG_Flag}'

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
            warn('No Attachment block for this channel')
            return None

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
            return self._cn_info
        except AttributeError:
            return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]

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
        15 LE Complex
        16 BE Complex
        """
        if not self.type == 4:
            return self._cn_data_type
        else:
            if byte_aligned:
                return 10
            else:
                return 0

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
        if self.type in (0, 1, 2):
            return self._cn_bit_count
        elif self.type == 3:
            if self.name == 'ms':
                if self.signal_data_type(info) == 13:
                    return 16
                else:
                    return 32
            elif self.name == 'days':
                return 16
            else:
                return 8
        elif self.type == 4:
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
            return self._cn_info['cn_sync_type']
        except (AttributeError, KeyError):
            return 0

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
            return self._cn_info['CABlock']
        except (AttributeError, KeyError):
            return None

    def is_ca_block(self, info):
        try:
            self._cn_info['CABlock']
            return True
        except (AttributeError, KeyError, TypeError):
            return False

    def record_id_size(self, info):
        return info['DG'][self.dataGroup]['dg_rec_id_size']

    def channel_type(self, info):
        try:
            return self._cn_type
        except AttributeError:
            try:
                return info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]['cn_type']
            except KeyError:
                return 0

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
        if self._cn_data_type < 6:
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
        if not self.type == 4:
            bit_count = self._cn_bit_count
            bit_offset = self._cn_bit_offset
            if aligned:
                n_bytes = _bits_to_bytes_aligned(bit_count + bit_offset, self.isnumeric(info))
            else:
                n_bytes = _bits_to_bytes_not_aligned(bit_count + bit_offset)
            if self.type in (1, 2):
                n_bytes *= self.ca_block(info)['PNd']
                block = self.ca_block(info)
                while 'CABlock' in block:
                    block = block['CABlock']
                    n_bytes *= block['PNd']
            if self.type == 3:
                if self.name == 'ms':
                    if self._cn_data_type == 13:
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
        if self._cn_data_type in (1, 3, 5, 9):
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
        if self.type == 4:
            data_format = f'{self.nBytes_aligned}V'
        elif self._cn_composition and 'CABlock' in self._cn_info:
            ca_block = self._cn_info['CABlock']
            endian, data_format = array_format4(self._cn_data_type, self._cn_bit_count // 8)
            array_desc = ca_block['ca_dim_size']
            Block = ca_block
            while 'CABlock' in Block:
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
        elif self.type == 3:
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
        else:
            if self.channel_type(info) in (1, 7):
                endian, data_format = array_format4(0, self.nBytes_aligned)
            else:
                endian, data_format = array_format4(self._cn_data_type, self.nBytes_aligned)
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
        if self.type == 0:
            if signal_data_type not in (13, 14):
                if self.channel_type(info) in (1, 7):
                    endian, data_type = data_type_format4(0, self.nBytes_aligned)
                else:
                    endian, data_type = data_type_format4(signal_data_type, self.nBytes_aligned)
                return f'{endian}{data_type}'
        elif self.type in (1, 2):
            ca = self.ca_block(info)
            n_bytes = _bits_to_bytes_aligned(self._cn_bit_count + self._cn_bit_offset, self.isnumeric(info))
            endian, data_type = data_type_format4(signal_data_type, n_bytes)
            return f"{endian}{ca['PNd']}{data_type}"
        elif self.type == 3:
            if self.name == 'ms':
                if signal_data_type == 13:
                    return 'H'
                else:
                    return 'I'
            elif self.name == 'days':
                return 'H'
            else:
                return 'B'
        elif self.type == 4:
            return f'{self.nBytes_aligned}s'
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
        if self.type in (0, 1, 2):
            return self._cn_bit_offset
        elif self.type == 3:
            return self._cn_bit_offset
        elif self.type == 4:
            return 0
        else:
            warn('Not found channel type')

    def calc_byte_offset(self, info):
        if self.type in (0, 1, 2):
            return self._cn_byte_offset
        elif self.type == 3:
            return self._cn_byte_offset + self.CANOpen_offset()
        elif self.type == 4:
            return info['CG'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        else:
            warn('Not found channel type')

    def pos_byte_beg(self, info):
        return self.record_id_size(info) + self.byteOffset

    def pos_byte_end(self, info):
        return self.pos_byte_beg(info) + self.nBytes_aligned

    def pos_bit_begin(self, info):
        return self.pos_byte_beg(info) * 8 + self.bit_offset(info)

    def pos_bit_end(self, info):
        return self.pos_bit_beg + self.bit_count(info)

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
        elif 'unit' in self._cn_info:
            unit = self._cn_info['unit']
        else:
            unit = ''
        if 'Comment' in unit:
            unit = unit['Comment']
        return unit

    def desc(self, info):
        if not self.type == 3:
            if 'Comment' in self._cn_info:
                desc = self._cn_info['Comment']
                if (desc is not None) and isinstance(desc, dict):
                    if 'description' in desc:
                        desc = desc['description']
                    elif 'name' in desc:
                        desc = desc['name']
                else:
                    desc = ''
                return desc
            else:
                return f'Invalid Bytes DG group {self.dataGroup} CG group {self.channelGroup}'
        else:
            return self.name

    def conversion(self, info):
        try:
            return info['CC'][self.dataGroup][self.channelGroup][self.channelNumber]
        except KeyError:
            return None

    def set(self, info):
        cn_info = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]
        self._cn_info = cn_info
        self._cn_type = cn_info['cn_type']
        self._cn_data_type = cn_info['cn_data_type']
        self._cn_bit_count = cn_info['cn_bit_count']
        self._cn_bit_offset = cn_info['cn_bit_offset']
        self._cn_byte_offset = cn_info['cn_byte_offset']
        self._cn_flags = cn_info['cn_flags']
        self._cn_composition = cn_info['cn_composition']
        self.name = cn_info['name']
        self.type = 0
        if cn_info['cn_composition'] and 'CABlock' in cn_info:
            self.type = 1
            block = cn_info['CABlock']
            if 'CABlock' in block:
                self.type = 2
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)
        self.pos_bit_beg = self.pos_bit_begin(info)

    def set_CANOpen(self, info, name):
        self.type = 3
        self.name = name
        cn_info = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]
        self._cn_info = cn_info
        self._cn_type = cn_info['cn_type']
        self._cn_data_type = cn_info['cn_data_type']
        self._cn_bit_count = cn_info['cn_bit_count']
        self._cn_bit_offset = cn_info['cn_bit_offset']
        self._cn_byte_offset = cn_info['cn_byte_offset']
        self._cn_flags = cn_info.get('cn_flags', 0)
        self._cn_composition = cn_info['cn_composition']
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)
        self.pos_bit_beg = self.pos_bit_begin(info)

    def set_invalid_bytes(self, info):
        self.type = 4
        self.name = f'invalid_bytes{self.dataGroup}'
        try:
            cn_info = info['CN'][self.dataGroup][self.channelGroup][self.channelNumber]
            self._cn_info = cn_info
            self._cn_type = cn_info.get('cn_type', 0)
            self._cn_data_type = cn_info.get('cn_data_type', 0)
            self._cn_bit_count = cn_info.get('cn_bit_count', 0)
            self._cn_bit_offset = cn_info.get('cn_bit_offset', 0)
            self._cn_byte_offset = cn_info.get('cn_byte_offset', 0)
            self._cn_flags = cn_info.get('cn_flags', 0)
            self._cn_composition = cn_info.get('cn_composition', False)
        except KeyError:
            self._cn_info = None
            self._cn_type = 4
            self._cn_data_type = 10
            self._cn_bit_count = 0
            self._cn_bit_offset = 0
            self._cn_byte_offset = 0
            self._cn_flags = 0
            self._cn_composition = False
        self.nBytes_aligned = self.calc_bytes(info)
        self.byteOffset = self.calc_byte_offset(info)
        self.pos_bit_beg = self.pos_bit_begin(info)

    def invalid_bit(self, info):
        try:
            return self._cn_info['cn_invalid_bit_pos']
        except (AttributeError, KeyError):
            return 0

    def has_invalid_bit(self, info):
        try:
            flags = self._cn_flags
        except AttributeError:
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
        self.name = f'{self.name}_{channel_group}'

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


def _int_size_code(number_of_bytes, codes):
    """ picks the smallest of the (1, 2, 4, 8) byte integer codes able to hold number_of_bytes

    Parameters
    ----------------
    number_of_bytes : int
        number of bytes taken by channel data in a record
    codes : tuple of str
        codes for respectively 1, 2, 4 and 8 bytes

    Returns
    -----------
    code : str or None
        None if number_of_bytes is larger than 8 bytes
    """
    if number_of_bytes == 1:
        return codes[0]
    if number_of_bytes == 2:
        return codes[1]
    if number_of_bytes <= 4:  # includes the 3 bytes case, stored on 4 bytes
        return codes[2]
    if number_of_bytes <= 8:
        return codes[3]
    return None


# numpy dtype codes for 1, 2, 4 and 8 bytes integers
_NP_UINT = ('u1', 'u2', 'u4', 'u8')
_NP_INT = ('i1', 'i2', 'i4', 'i8')
# struct codes for 1, 2, 4 and 8 bytes integers
_STRUCT_UINT = ('B', 'H', 'I', 'Q')
_STRUCT_INT = ('b', 'h', 'i', 'q')
# data type -> (endian, kind) for numeric types, kind being 'uint', 'int', 'float' or 'complex'
_NUMERIC_TYPES4 = {
    0: ('<', 'uint'), 1: ('>', 'uint'),
    2: ('<', 'int'), 3: ('>', 'int'),
    4: ('<', 'float'), 5: ('>', 'float'),
    15: ('<', 'complex'), 16: ('>', 'complex'),
}
# floating point / complex sizes (bytes) -> numpy dtype
_NP_FLOAT = {2: 'f2', 4: 'f4', 8: 'f8'}
_NP_COMPLEX = {4: '(2,)f2', 8: 'c8', 16: 'c16'}  # 4 bytes: half-precision [real_f16, imag_f16]
# floating point / complex sizes (bytes) -> struct format
_STRUCT_FLOAT = {2: 'e', 4: 'f', 8: 'd'}
_STRUCT_COMPLEX = {4: '2e', 8: '2f', 16: '2d'}
# string / bytes types -> endian
_STRING_ENDIAN4 = {6: '', 7: '', 8: '<', 9: '>', 10: '', 11: '', 12: ''}


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
    raw = f'V{number_of_bytes}'
    if signal_data_type in _NUMERIC_TYPES4:
        endian, kind = _NUMERIC_TYPES4[signal_data_type]
        if kind == 'uint':
            data_type = _int_size_code(number_of_bytes, _NP_UINT) or f'{number_of_bytes}V'
        elif kind == 'int':
            data_type = _int_size_code(number_of_bytes, _NP_INT)
            if data_type is None:
                warn(f'Unsupported number of bytes for signed int {number_of_bytes}')
                data_type = raw
        elif kind == 'float':
            data_type = _NP_FLOAT.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for floating point {number_of_bytes}')
                data_type = raw
        else:  # complex
            data_type = _NP_COMPLEX.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for complex {number_of_bytes}')
                data_type = raw
    elif signal_data_type in (6, 7, 8, 9):  # strings
        endian, data_type = _STRING_ENDIAN4[signal_data_type], f'S{number_of_bytes}'
    elif signal_data_type in (10, 11, 12):  # bytes array, MIME sample or MIME stream
        endian, data_type = '', raw
    elif signal_data_type in (13, 14):  # CANOpen date or time
        endian, data_type = '', ''
    else:
        warn(f'Unsupported Signal Data Type {signal_data_type} {number_of_bytes}')
        endian, data_type = '', raw
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
    raw = f'{number_of_bytes}s'
    if signal_data_type in _NUMERIC_TYPES4:
        endian, kind = _NUMERIC_TYPES4[signal_data_type]
        if kind == 'uint':
            data_type = _int_size_code(number_of_bytes, _STRUCT_UINT) or raw
        elif kind == 'int':
            data_type = _int_size_code(number_of_bytes, _STRUCT_INT)
            if data_type is None:
                warn(f'Unsupported number of bytes for signed int {number_of_bytes}')
                data_type = raw
        elif kind == 'float':
            data_type = _STRUCT_FLOAT.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for floating point {number_of_bytes}')
                data_type = raw
        else:  # complex
            data_type = _STRUCT_COMPLEX.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for complex {number_of_bytes}')
                data_type = raw
    elif signal_data_type in _STRING_ENDIAN4:  # string/bytes
        endian, data_type = _STRING_ENDIAN4[signal_data_type], raw
    else:
        warn(f'Unsupported Signal Data Type {signal_data_type} {number_of_bytes}')
        endian, data_type = '', raw
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
    __slots__ = [
        'name', 'channelNumber', 'signalDataType', 'bitCount', 'posBitBeg',
        'posBitEnd', 'byteOffset', 'bitOffset', 'nBytes_aligned',
        'nBytes_not_aligned', 'dataFormat', 'nativedataFormat', 'CFormat',
        'embedding_channel_bitOffset', 'posByteBeg', 'posByteEnd',
        'bit_masking_needed', 'channelType', 'unit', 'desc', 'conversion',
        'recAttributeName', 'RecordFormat'
    ]

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
        cn = info['CNBlock'][data_group][channel_group][channel_number]
        cc = info['CCBlock'][data_group][channel_group][channel_number]
        self.name = cn['signalName']
        self.channelNumber = channel_number
        sdt = cn['signalDataType']
        self.signalDataType = sdt
        numeric = sdt not in (7, 8)
        bc = cn['numberOfBits']
        self.bitCount = bc
        byte_order = info['IDBlock']['ByteOrder']
        pbb = cn['numberOfTheFirstBits']
        self.posBitBeg = pbb
        self.posBitEnd = pbb + bc
        self.byteOffset = pbb // 8
        bo = pbb % 8
        self.bitOffset = bo
        self.nBytes_aligned = _bits_to_bytes_aligned(bc + bo, numeric)
        self.nBytes_not_aligned = _bits_to_bytes_not_aligned(bc + bo)
        (self.dataFormat, self.nativedataFormat) = \
            _array_format3(sdt, self.nBytes_aligned, byte_order)
        self.CFormat = Struct(_data_type_format3(sdt, self.nBytes_aligned, byte_order))
        self.embedding_channel_bitOffset = bo  # for channel containing other channels
        self.posByteBeg = record_id_number + self.byteOffset
        self.posByteEnd = record_id_number + self.byteOffset + self.nBytes_aligned
        self.bit_masking_needed = (self.nBytes_aligned != bc / 8)
        self.channelType = cn['channelType']
        self.unit = cc.get('physicalUnit', '')
        self.desc = cn.get('signalDescription', '')
        self.conversion = cc

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
        self.name = f'{self.name}_{channel_group}'
        self.recAttributeName = _convert_name(self.name)
        self.RecordFormat = (((f'{self.recAttributeName}_title', self.recAttributeName)), self.dataFormat)


# MDF3 numeric type -> kind ('uint', 'int', 'float')
_NUMERIC_TYPES3 = {
    0: 'uint', 9: 'uint', 13: 'uint',
    1: 'int', 10: 'int', 14: 'int',
    2: 'float', 3: 'float', 11: 'float', 12: 'float', 15: 'float', 16: 'float',
}
_STRUCT_FLOAT3 = {4: 'f', 8: 'd'}
_NP_FLOAT3 = {4: 'f4', 8: 'f8'}


def _mdf3_endian(signal_data_type, byte_order):
    if signal_data_type in (0, 1, 2, 3):
        return '>' if byte_order else '<'
    elif signal_data_type in (13, 14, 15, 16):
        return '<'
    elif signal_data_type in (9, 10, 11, 12):
        return '>'
    return ''


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
    raw = f'{number_of_bytes}s'
    if signal_data_type in _NUMERIC_TYPES3:
        kind = _NUMERIC_TYPES3[signal_data_type]
        if kind == 'uint':
            data_type = _int_size_code(number_of_bytes, _STRUCT_UINT)
            if data_type is None:
                warn(f'Unsupported number of bytes for unsigned int {number_of_bytes}')
                data_type = raw
        elif kind == 'int':
            data_type = _int_size_code(number_of_bytes, _STRUCT_INT)
            if data_type is None:
                warn(f'Unsupported number of bytes for signed int {number_of_bytes}')
                data_type = raw
        else:
            data_type = _STRUCT_FLOAT3.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for floating point {number_of_bytes}')
                data_type = raw
    elif signal_data_type in (7, 8):
        data_type = raw
    else:
        warn(f'Unsupported Signal Data Type {signal_data_type} nBytes {number_of_bytes}')
        data_type = raw

    endian = _mdf3_endian(signal_data_type, byte_order)
    return endian + data_type


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
    raw = f'V{number_of_bytes}'
    if signal_data_type in _NUMERIC_TYPES3:
        kind = _NUMERIC_TYPES3[signal_data_type]
        if kind == 'uint':
            data_type = _int_size_code(number_of_bytes, _NP_UINT)
            if data_type is None:
                warn(f'Unsupported number of bytes for unsigned int {number_of_bytes}')
                data_type = raw
        elif kind == 'int':
            data_type = _int_size_code(number_of_bytes, _NP_INT)
            if data_type is None:
                warn(f'Unsupported number of bytes for signed int {number_of_bytes}')
                data_type = raw
        else:
            data_type = _NP_FLOAT3.get(number_of_bytes)
            if data_type is None:
                warn(f'Unsupported number of bytes for floating point {number_of_bytes}')
                data_type = raw
    elif signal_data_type == 7:
        data_type = f'S{number_of_bytes}'
    elif signal_data_type == 8:
        data_type = raw
    else:
        warn(f'Unsupported Signal Data Type {signal_data_type} nBytes {number_of_bytes}')
        data_type = raw

    native_data_type = data_type
    endian = _mdf3_endian(signal_data_type, byte_order)
    return (endian + data_type, native_data_type)
