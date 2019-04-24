import numpy as np
cimport numpy as np
from sys import byteorder
#cimport cython

from cpython.bytes cimport PyBytes_AsString
from libc.string cimport memcpy

#@cython.boundscheck(False)
#@cython.wraparound(False)
def data_read(bytes tmp, unsigned short bit_count,
        unsigned short signal_data_type, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned char bit_offset,
        unsigned long pos_byte_beg, unsigned long n_bytes, array):
    """dataRead function to read in cython a channel from a byte stream

    Parameters
    ------------
    tmp : bytes
        byte stream
    bit_count : unsigned short
        number of bit taken by the channel in the record
    signal_data_type : unsigned short
        int to describe data type
    record_format : string
        basic numpy dtype description of data type, used to create
        returned numpy ndarray
    number_of_records : unsigned long long
        number of records in byte stream
    record_byte_size : unsigned long
        number of bytes taken by one record repeated in byte stream
    bit_offset : unsigned char
        bit offset of data in C aligned bytes
    pos_byte_beg : unsigned long
        beginning byte position of channel in record
    n_bytes : unsigned long
        bytes length of channel in record
    array : boolean
        reads an array, not a vector

    Returns
    -------
    ndarray of type record_format with number_of_records records.
    Byte order is swapped if necessary to match machine byte order before bits offset and masking
    """
    cdef char* bit_stream = PyBytes_AsString(tmp)
    if not array:
        if 'V' in record_format or 'S' in record_format or record_format is None:
            return read_byte(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset)
        elif signal_data_type in (4, 5) and n_bytes == 4:  # float
            if (byteorder == 'little' and signal_data_type == 4) or \
                    (byteorder == 'big' and signal_data_type == 5):
                return read_float(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, 0)
            else: #  swap bytes
                return read_float(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, 1)
        elif signal_data_type in (4, 5) and n_bytes == 8:  # double
            if (byteorder == 'little' and signal_data_type == 4) or \
                    (byteorder == 'big' and signal_data_type == 5):
                return read_double(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 0)
            else: #  swap bytes
                return read_double(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, 1)
        elif signal_data_type in (0, 1, 13) and n_bytes == 1:  # unsigned char
            return read_unsigned_char(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, bit_count, bit_offset)
        elif signal_data_type in (2, 3) and n_bytes == 1:  # signed char
            return read_signed_char(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, bit_count, bit_offset)
        elif signal_data_type in (0, 1, 13, 14) and n_bytes <= 2:  # unsigned short
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_short(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, bit_count, bit_offset, 0)
            else: #  swap bytes
                return read_unsigned_short(bit_stream, record_format, number_of_records,
                                      record_byte_size, pos_byte_beg, bit_count, bit_offset, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 2:  # signed short
            if (byteorder == 'little' and signal_data_type == 2) or \
                    (byteorder == 'big' and signal_data_type == 3):
                return read_signed_short(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, bit_count, bit_offset, 0)
            else: #  swap bytes
                return read_signed_short(bit_stream, record_format, number_of_records,
                                     record_byte_size, pos_byte_beg, bit_count, bit_offset, 1)
        elif signal_data_type in (0, 1, 14) and n_bytes <= 4:  # unsigned int
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_int(bit_stream, record_format, number_of_records,
                                    record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_unsigned_int(bit_stream, record_format, number_of_records,
                                    record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 4:  # signed int
            if (byteorder == 'little' and signal_data_type == 2) or \
                    (byteorder == 'big' and signal_data_type == 3):
                return read_signed_int(bit_stream, record_format, number_of_records,
                                   record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_signed_int(bit_stream, record_format, number_of_records,
                                   record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (0, 1) and n_bytes <= 8:  # unsigned long long
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_unsigned_longlong(bit_stream, record_format, number_of_records,
                                         record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_unsigned_longlong(bit_stream, record_format, number_of_records,
                                         record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        elif signal_data_type in (2, 3) and n_bytes <= 8:  # signed long long
            if (byteorder == 'little' and signal_data_type == 0) or \
                    (byteorder == 'big' and signal_data_type == 1):
                return read_signed_longlong(bit_stream, record_format, number_of_records,
                                        record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 0)
            else: #  swap bytes
                return read_signed_longlong(bit_stream, record_format, number_of_records,
                                        record_byte_size, pos_byte_beg, bit_count, bit_offset, n_bytes, 1)
        else:
            return read_byte(bit_stream, record_format, number_of_records,
                                record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset)
    else: # array
        if (byteorder == 'little' and signal_data_type in (0, 2, 4)) or \
                    (byteorder == 'big' and signal_data_type in (1, 3, 5)):
            return read_array(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset, 0)
        else: #  swap bytes
            return read_array(bit_stream, record_format, number_of_records,
                                 record_byte_size, pos_byte_beg, n_bytes, bit_count, bit_offset, 1)


cdef inline read_float(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.float32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef float temp_float = 0
    cdef char temp[4]
    for i in range(number_of_records):
        memcpy(&temp_float, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
        buf[i] = temp_float
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
    
cdef inline read_double(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned char swap):
    cdef np.ndarray[np.float64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef double temp_double = 0
    cdef char temp[8]
    for i in range(number_of_records):
        memcpy(&temp_double, &bit_stream[pos_byte_beg + record_byte_size * i], 8)
        buf[i] = temp_double
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
    
cdef inline read_unsigned_char(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray[np.uint8_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned char mask = ((1 << bit_count) - 1)
    cdef unsigned char temp1byte = 0
    if bit_count == 8:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            # right shift 
            if bit_offset > 0:
                temp1byte = temp1byte >> bit_offset
            # mask left part
            temp1byte &= mask
            buf[i] = temp1byte
    return buf

cdef inline read_signed_char(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray[np.int8_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef char mask = ((1 << bit_count) - 1)
    cdef char temp1byte = 0
    cdef char sign_bit = 0
    cdef char sign_bit_mask = (1 << (bit_count-1))
    cdef char sign_extend = ((1 << (8 - bit_count)) - 1) << bit_count
    if bit_count == 8:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(number_of_records):
            memcpy(&temp1byte, &bit_stream[pos_byte_beg + record_byte_size * i], 1)
            # right shift 
            if bit_offset > 0:
                temp1byte = temp1byte >> bit_offset
            # mask left part
            temp1byte &= mask
            sign_bit = temp1byte & sign_bit_mask
            if sign_bit: #  negative value, sign extend
                temp1byte |= sign_extend
            buf[i] = temp1byte
    return buf

cdef inline read_unsigned_short(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray[np.uint16_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bit_count) - 1)
    cdef unsigned short temp2byte = 0
    cdef char temp[2]
    if bit_count == 16:
        for i in range(number_of_records):
            memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                if bit_count < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        else:
            for i in range(number_of_records):
                memcpy(&temp, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                temp2byte = swap16(temp)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                if bit_count < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        return buf

cdef inline read_signed_short(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray[np.int16_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef short mask = ((1 << bit_count) - 1)
    cdef short temp2byte = 0
    cdef short sign_bit = 0
    cdef short sign_bit_mask = (1 << (bit_count-1))
    cdef short sign_extend = ((1 << (16 - bit_count)) - 1) << bit_count
    cdef char temp[2]
    if bit_count == 16:
        for i in range(number_of_records):
            memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp2byte, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                temp2byte &= mask
                sign_bit = temp2byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp2byte |= sign_extend
                buf[i] = temp2byte
        else:
            for i in range(number_of_records):
                memcpy(&temp, &bit_stream[pos_byte_beg + record_byte_size * i], 2)
                temp2byte = swap16(temp)
                # right shift
                if bit_offset > 0:
                    temp2byte = temp2byte >> bit_offset
                # mask left part
                temp2byte &= mask
                sign_bit = temp2byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp2byte |= sign_extend
                buf[i] = temp2byte
        return buf

cdef inline swap(char* temp, const char byte1, const char byte2):
    cdef char t
    t = temp[byte1]
    temp[byte1] = temp[byte2]
    temp[byte2] = t
    return temp

cdef inline swap16(char* temp):
    return swap(temp, 0, 1)

cdef inline read_unsigned_int(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.uint32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bit_count) - 1)
    cdef unsigned int temp4byte = 0
    cdef char temp4[4]
    cdef char temp3[3]
    if bit_count == 32:
        for i in range(number_of_records):
            memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 4:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp4, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = swap32(temp4) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp3, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = swap24(temp3)  #  swap bytes
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    
cdef inline swap24(char* temp):
    return swap(temp, 0, 2)

cdef inline swap32(char* temp):
    cdef char t
    t = temp[0]
    temp[0] = temp[3]
    temp[3] = t
    t = temp[1]
    temp[1] = temp[2]
    temp[2] = t
    return temp

cdef inline read_signed_int(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.int32_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef int mask = ((1 << bit_count) - 1)
    cdef int temp4byte = 0
    cdef int sign_bit = 0
    cdef int sign_bit_mask = (1 << (bit_count-1))
    cdef int sign_extend = ((1 << (32 - bit_count)) - 1) << bit_count
    cdef char temp4[4]
    cdef char temp3[3]
    if bit_count == 32:
        for i in range(number_of_records):
            memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], 4)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 4:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp4, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = swap32(temp4)  #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp4byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        else:
            for i in range(number_of_records):
                memcpy(&temp3, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp4byte = swap24(temp3)  #  swap bytes
                # right shift 
                if bit_offset > 0:
                    temp4byte = temp4byte >> bit_offset
                # mask left part
                if bit_count < 24:
                    temp4byte &= mask
                sign_bit = temp4byte & sign_bit_mask # assumes return in little endian, to be reviewed
                if sign_bit: #  negative value, sign extend
                    temp4byte |= sign_extend
                buf[i] = temp4byte
        return buf

cdef inline read_unsigned_longlong(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.uint64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bit_count) - 1)
    cdef unsigned long long temp8byte = 0
    cdef char temp8[8]
    cdef char temp7[7]
    cdef char temp6[6]
    cdef char temp5[5]
    if bit_count == 64:
        for i in range(number_of_records):
            memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 8:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp8, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap64(temp8) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 7:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp7, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap56(temp7) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 6:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp6, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap48(temp6) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif n_bytes == 5:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp5, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap40(temp5) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
    return buf

cdef inline swap40(char* temp):
    cdef char t
    t = temp[0]
    temp[0] = temp[4]
    temp[4] = t
    t = temp[1]
    temp[1] = temp[3]
    temp[3] = t
    return temp

cdef inline swap48(char* temp):
    cdef char t
    t = temp[0]
    temp[0] = temp[5]
    temp[5] = t
    t = temp[1]
    temp[1] = temp[4]
    temp[4] = t
    t = temp[2]
    temp[2] = temp[3]
    temp[3] = t
    return temp

cdef inline swap56(char* temp):
    cdef char t
    t = temp[0]
    temp[0] = temp[6]
    temp[6] = t
    t = temp[1]
    temp[1] = temp[5]
    temp[5] = t
    t = temp[2]
    temp[2] = temp[4]
    temp[4] = t
    return temp

cdef inline swap64(char* temp):
    cdef char t
    t = temp[0]
    temp[0] = temp[7]
    temp[7] = t
    t = temp[1]
    temp[1] = temp[6]
    temp[6] = t
    t = temp[2]
    temp[2] = temp[5]
    temp[5] = t
    t = temp[3]
    temp[3] = temp[4]
    temp[4] = t
    return temp

cdef inline read_signed_longlong(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg,
        unsigned long bit_count, unsigned char bit_offset, unsigned long n_bytes, unsigned char swap):
    cdef np.ndarray[np.int64_t] buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef long long mask = ((1 << bit_count) - 1)
    cdef long long temp8byte = 0
    cdef long sign_bit = 0
    cdef long long sign_bit_mask = (1 << (bit_count-1))
    cdef long long sign_extend = ((1 << (64 - bit_count)) - 1) << bit_count
    cdef char temp8[8]
    cdef char temp7[7]
    cdef char temp6[6]
    cdef char temp5[5]
    if bit_count == 64:
        for i in range(number_of_records):
            memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif n_bytes == 8:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp8, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap64(temp8) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 64:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 7:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp7, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap56(temp7) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 56:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 6:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp6, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap48(temp6) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 48:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    elif n_bytes == 5:
        if swap == 0:
            for i in range(number_of_records):
                memcpy(&temp8byte, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 40:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
        else:
            for i in range(number_of_records):
                memcpy(&temp5, &bit_stream[pos_byte_beg + record_byte_size * i], n_bytes)
                temp8byte = swap40(temp5) #  swap bytes
                # right shift
                if bit_offset > 0:
                    temp8byte = temp8byte >> bit_offset
                # mask left part
                if bit_count < 40:
                    temp8byte &= mask
                sign_bit = temp8byte & sign_bit_mask
                if sign_bit: #  negative value, sign extend
                    temp8byte |= sign_extend
                buf[i] = temp8byte
    return buf

cdef inline read_byte(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned long n_bytes,
        unsigned long bit_count, unsigned char bit_offset):
    cdef np.ndarray buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long pos_byte_end = pos_byte_beg + n_bytes
    for i in range(number_of_records):
            buf[i] = bytes(bit_stream[pos_byte_beg + record_byte_size * i:\
                    pos_byte_end + record_byte_size * i])
    return buf

cdef inline read_array(const char* bit_stream, str record_format, unsigned long long number_of_records,
        unsigned long record_byte_size, unsigned long pos_byte_beg, unsigned long n_bytes,
        unsigned long bit_count, unsigned char bit_offset, unsigned char swap):
    cdef np.ndarray buf = np.empty(number_of_records, dtype=record_format)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long pos_byte_end = pos_byte_beg + n_bytes
    for i in range(number_of_records):
        buf[i] = np.fromstring(bit_stream[pos_byte_beg + record_byte_size * i:\
            pos_byte_end + record_byte_size * i], dtype=record_format)
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
