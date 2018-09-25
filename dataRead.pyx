import numpy as np
cimport numpy as np
from sys import byteorder
#cimport cython

from cpython.bytes cimport PyBytes_AsString
from libc.string cimport memcpy
    
#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(bytes tmp, unsigned short bitCount,
        unsigned short signalDataType, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned char bitOffset,
        unsigned long posByteBeg, unsigned long nBytes, array):
    """dataRead function to read in cython a channel from a byte stream

    Parameters
    ------------
    tmp : bytes
        byte stream
    bitCount : unsigned short
        number of bit taken by the channel in the record
    signalDataType : unsigned short
        int to describe data type
    RecordFormat : string
        basic numpy dtype description of data type, used to create
	    returned numpy ndarray
    numberOfRecords : unsigned long long
        number of records in byte stream
    record_byte_size : unsigned long
        number of bytes taken by one record repeated in byte stream
    bitOffset : unsigned char
        bit offset of data in C aligned bytes
    posByteBeg : unsigned long
        beginning byte position of channel in record
    nBytes : unsigned long
        bytes length of channel in record
    array : boolean
        reads an array, not a vector

    Return
    -------
    ndarray of type RecordFormat with numberOfRecords records.
    Byte order is swapped if necessary to match machine byte order before bits offset and masking
    """
    cdef char* bita = PyBytes_AsString(tmp)
    if not array:
        if 'V' in RecordFormat or 'S' in RecordFormat or RecordFormat is None:
            return dataReadByte(bita, RecordFormat, numberOfRecords,
                record_byte_size, posByteBeg, nBytes, bitCount, bitOffset)
        elif signalDataType in (4, 5) and nBytes == 4:  # float
            if (byteorder == 'little' and signalDataType == 4) or \
                    (byteorder == 'big' and signalDataType == 5):
                return dataReadFloat(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, 0)
            else: #  swap bytes
                return dataReadFloat(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, 1)
        elif signalDataType in (4, 5) and nBytes == 8:  # double
            if (byteorder == 'little' and signalDataType == 4) or \
                    (byteorder == 'big' and signalDataType == 5):
                return dataReadDouble(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, 0)
            else: #  swap bytes
                return dataReadDouble(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, 1)
        elif signalDataType in (0, 1, 13) and nBytes == 1:  # unsigned char
            return dataReadUChar(bita, RecordFormat, numberOfRecords,
                record_byte_size, posByteBeg, bitCount, bitOffset)
        elif signalDataType in (2, 3) and nBytes == 1:  # signed char
            return dataReadChar(bita, RecordFormat, numberOfRecords,
                record_byte_size, posByteBeg, bitCount, bitOffset)
        elif signalDataType in (0, 1, 13, 14) and nBytes <= 2:  # unsigned short
            if (byteorder == 'little' and signalDataType == 0) or \
                    (byteorder == 'big' and signalDataType == 1):
                return dataReadUShort(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, 0)
            else: #  swap bytes
                return dataReadUShort(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, 1)
        elif signalDataType in (2, 3) and nBytes <= 2:  # signed short
            if (byteorder == 'little' and signalDataType == 2) or \
                    (byteorder == 'big' and signalDataType == 3):
                return dataReadShort(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, 0)
            else: #  swap bytes
                return dataReadShort(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, 1)
        elif signalDataType in (0, 1, 14) and nBytes <= 4:  # unsigned int
            if (byteorder == 'little' and signalDataType == 0) or \
                    (byteorder == 'big' and signalDataType == 1):
                return dataReadUInt(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 0)
            else: #  swap bytes
                return dataReadUInt(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 1)
        elif signalDataType in (2, 3) and nBytes <= 4:  # signed int
            if (byteorder == 'little' and signalDataType == 2) or \
                    (byteorder == 'big' and signalDataType == 3):
                return dataReadInt(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 0)
            else: #  swap bytes
                return dataReadInt(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 1)
        elif signalDataType in (0, 1) and nBytes <= 8:  # unsigned long long
            if (byteorder == 'little' and signalDataType == 0) or \
                    (byteorder == 'big' and signalDataType == 1):
                return dataReadULongLong(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 0)
            else: #  swap bytes
                return dataReadULongLong(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 1)
        elif signalDataType in (2, 3) and nBytes <= 8:  # signed long long
            if (byteorder == 'little' and signalDataType == 0) or \
                    (byteorder == 'big' and signalDataType == 1):
                return dataReadLongLong(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 0)
            else: #  swap bytes
                return dataReadLongLong(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, bitCount, bitOffset, nBytes, 1)
        else:
            return dataReadByte(bita, RecordFormat, numberOfRecords,
                record_byte_size, posByteBeg, nBytes, bitCount, bitOffset)
    else: # array
        if (byteorder == 'little' and signalDataType in (0, 2, 4)) or \
                    (byteorder == 'big' and signalDataType in (1, 3, 5)):
            return dataReadArray(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, nBytes, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadArray(bita, RecordFormat, numberOfRecords,
                    record_byte_size, posByteBeg, nBytes, bitCount, bitOffset, 1)


cdef inline dataReadFloat(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned char swap):
    cdef np.ndarray[np.float32_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef float tempfloat = 0
    cdef char temp[4]
    for i in range(numberOfRecords):
        memcpy(&tempfloat, &bita[posByteBeg + record_byte_size * i], 4)
        buf[i] = tempfloat
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
    
cdef inline dataReadDouble(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned char swap):
    cdef np.ndarray[np.float64_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef double tempDouble = 0
    cdef char temp[8]
    for i in range(numberOfRecords):
        memcpy(&tempDouble, &bita[posByteBeg + record_byte_size * i], 8)
        buf[i] = tempDouble
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
    
cdef inline dataReadUChar(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray[np.uint8_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned char mask = ((1 << bitCount) - 1)
    cdef unsigned char temp1byte = 0
    if bitCount == 8:
        for i in range(numberOfRecords):
            memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
            # right shift 
            if bitOffset > 0:
                temp1byte = temp1byte >> bitOffset
            # mask left part
            temp1byte &= mask
            buf[i] = temp1byte
    return buf

cdef inline dataReadChar(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray[np.int8_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef char mask = ((1 << bitCount) - 1)
    cdef char temp1byte = 0
    cdef char signBit = 0
    cdef char signBitMask = (1 << (bitCount-1))
    cdef char signExtend = ((1 << (8 - bitCount)) - 1) << bitCount
    if bitCount == 8:
        for i in range(numberOfRecords):
            memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
            buf[i] = temp1byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
            # right shift 
            if bitOffset > 0:
                temp1byte = temp1byte >> bitOffset
            # mask left part
            temp1byte &= mask
            signBit = temp1byte & signBitMask
            if signBit: #  negative value, sign extend
                temp1byte |= signExtend
            buf[i] = temp1byte
    return buf

cdef inline dataReadUShort(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.uint16_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bitCount) - 1)
    cdef unsigned short temp2byte = 0
    cdef char temp[2]
    if bitCount == 16:
        for i in range(numberOfRecords):
            memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
                # right shift
                if bitOffset > 0:
                    temp2byte = temp2byte >> bitOffset
                # mask left part
                if bitCount < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp, &bita[posByteBeg + record_byte_size * i], 2)
                temp2byte = temp[0]<<8 | temp[1]  #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp2byte = temp2byte >> bitOffset
                # mask left part
                if bitCount < 16:
                    temp2byte &= mask
                buf[i] = temp2byte
        return buf

cdef inline dataReadShort(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.int16_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef short mask = ((1 << bitCount) - 1)
    cdef short temp2byte = 0
    cdef short signBit = 0
    cdef short signBitMask = (1 << (bitCount-1))
    cdef short signExtend = ((1 << (16 - bitCount)) - 1) << bitCount
    cdef char temp[2]
    if bitCount == 16:
        for i in range(numberOfRecords):
            memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
            buf[i] = temp2byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
                # right shift
                if bitOffset > 0:
                    temp2byte = temp2byte >> bitOffset
                # mask left part
                temp2byte &= mask
                signBit = temp2byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp2byte |= signExtend
                buf[i] = temp2byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp, &bita[posByteBeg + record_byte_size * i], 2)
                temp2byte = temp[0]<<8 | temp[1]  #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp2byte = temp2byte >> bitOffset
                # mask left part
                temp2byte &= mask
                signBit = temp2byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp2byte |= signExtend
                buf[i] = temp2byte
        return buf

cdef inline dataReadUInt(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned long nBytes, unsigned char swap):
    cdef np.ndarray[np.uint32_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bitCount) - 1)
    cdef unsigned int temp4byte = 0
    cdef char temp4[4]
    cdef char temp3[3]
    if bitCount == 32:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], 2)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif nBytes == 4:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp4, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp4[0]<<24 | temp4[1]<<16 | temp4[2]<<8 | temp4[3]  #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift 
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp3, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp3[0]<<16 | temp3[1]<<8 | temp3[2]  #  swap bytes
                # right shift 
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 24:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    

cdef inline dataReadInt(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned long nBytes, unsigned char swap):
    cdef np.ndarray[np.int32_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef int mask = ((1 << bitCount) - 1)
    cdef int temp4byte = 0
    cdef int signBit = 0
    cdef int signBitMask = (1 << (bitCount-1))
    cdef int signExtend = ((1 << (32 - bitCount)) - 1) << bitCount
    cdef char temp4[4]
    cdef char temp3[3]
    if bitCount == 32:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], 2)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif nBytes == 4:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp4byte &= mask
                signBit = temp4byte & signBitMask # assumes return in little endian, to be reviewed
                if signBit: #  negative value, sign extend
                    temp4byte |= signExtend
                buf[i] = temp4byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp4, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp4[0]<<24 | temp4[1]<<16 | temp4[2]<<8 | temp4[3]  #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp4byte &= mask
                signBit = temp4byte & signBitMask # assumes return in little endian, to be reviewed
                if signBit: #  negative value, sign extend
                    temp4byte |= signExtend
                buf[i] = temp4byte
        return buf
    else:  # on 3 bytes
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift 
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 24:
                    temp4byte &= mask
                signBit = temp4byte & signBitMask # assumes return in little endian, to be reviewed
                if signBit: #  negative value, sign extend
                    temp4byte |= signExtend
                buf[i] = temp4byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp3, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp3[0]<<16 | temp3[1]<<8 | temp3[2]  #  swap bytes
                # right shift 
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 24:
                    temp4byte &= mask
                signBit = temp4byte & signBitMask # assumes return in little endian, to be reviewed
                if signBit: #  negative value, sign extend
                    temp4byte |= signExtend
                buf[i] = temp4byte
        return buf

cdef inline dataReadULongLong(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned long nBytes, unsigned char swap):
    cdef np.ndarray[np.uint64_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bitCount) - 1)
    cdef unsigned long long temp8byte = 0
    cdef char temp8[8]
    cdef char temp7[7]
    cdef char temp6[6]
    cdef char temp5[5]
    if bitCount == 64:
        for i in range(numberOfRecords):
            memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif nBytes == 8:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp8, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp8[0]<<56 | temp8[1]<<48 | temp8[2]<<40 | temp8[3]<<32 | \
                            temp8[4]<<24 | temp8[5]<<16 | temp8[6]<<8 | temp8[7] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 64:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif nBytes == 7:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp7, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp7[0]<<48 | temp7[1]<<40 | temp7[2]<<32 | \
                            temp7[3]<<24 | temp7[4]<<16 | temp7[5]<<8 | temp7[6] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 56:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif nBytes == 6:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp6, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp6[0]<<40 | temp6[1]<<32 | temp6[2]<<24 | \
                            temp6[3]<<16 | temp6[4]<<8 | temp6[5] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 48:
                    temp8byte &= mask
                buf[i] = temp8byte
    elif nBytes == 5:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp5, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp5[0]<<32 | temp5[1]<<24 | \
                            temp5[2]<<16 | temp5[3]<<8 | temp5[4] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp8byte &= mask
                buf[i] = temp8byte
    return buf

cdef inline dataReadLongLong(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg,
        unsigned long bitCount, unsigned char bitOffset, unsigned long nBytes, unsigned char swap):
    cdef np.ndarray[np.int64_t] buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef long long mask = ((1 << bitCount) - 1)
    cdef long long temp8byte = 0
    cdef long signBit = 0
    cdef long long signBitMask = (1 << (bitCount-1))
    cdef long long signExtend = ((1 << (64 - bitCount)) - 1) << bitCount
    cdef char temp8[8]
    cdef char temp7[7]
    cdef char temp6[6]
    cdef char temp5[5]
    if bitCount == 64:
        for i in range(numberOfRecords):
            memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp8byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    elif nBytes == 8:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 64:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp8, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp8[0]<<56 | temp8[1]<<48 | temp8[2]<<40 | temp8[3]<<32 | \
                            temp8[4]<<24 | temp8[5]<<16 | temp8[6]<<8 | temp8[7] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 64:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
    elif nBytes == 7:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 56:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp7, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp7[0]<<48 | temp7[1]<<40 | temp7[2]<<32 | \
                            temp7[3]<<24 | temp7[4]<<16 | temp7[5]<<8 | temp7[6] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 56:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
    elif nBytes == 6:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 48:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp6, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp6[0]<<40 | temp6[1]<<32 | temp6[2]<<24 | \
                            temp6[3]<<16 | temp6[4]<<8 | temp6[5] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 48:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
    elif nBytes == 5:
        if swap == 0:
            for i in range(numberOfRecords):
                memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 40:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
        else:
            for i in range(numberOfRecords):
                memcpy(&temp5, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp8byte = temp5[0]<<32 | temp5[1]<<24 | \
                            temp5[2]<<16 | temp5[3]<<8 | temp5[4] #  swap bytes
                # right shift
                if bitOffset > 0:
                    temp8byte = temp8byte >> bitOffset
                # mask left part
                if bitCount < 40:
                    temp8byte &= mask
                signBit = temp8byte & signBitMask
                if signBit: #  negative value, sign extend
                    temp8byte |= signExtend
                buf[i] = temp8byte
    return buf

cdef inline dataReadByte(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long nBytes,
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long posByteEnd = posByteBeg + nBytes
    for i in range(numberOfRecords):
            buf[i] = bytes(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])
    return buf

cdef inline dataReadArray(const char* bita, str RecordFormat, unsigned long long numberOfRecords,
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long nBytes,
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray buf = np.empty(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long posByteEnd = posByteBeg + nBytes
    for i in range(numberOfRecords):
        buf[i] = np.fromstring(bita[posByteBeg + record_byte_size * i:\
            posByteEnd + record_byte_size * i], dtype=RecordFormat)
    if swap == 0:
        return buf
    else:
        return buf.byteswap()
