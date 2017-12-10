import numpy as np
cimport numpy as np
from sys import byteorder
#cimport cython

from cpython.bytes cimport PyBytes_AsString
from libc.string cimport memcpy
    
#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(bytes tmp, unsigned short bitCount, \
        unsigned short signalDataType, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned char bitOffset, \
        unsigned long posByteBeg, unsigned long posByteEnd):
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
    posByteEnd : unsigned long
        ending byte position of channel in record

    Return
    -------
    ndarray of type RecordFormat with numberOfRecords records.
    Byte order is swapped if necessary to match machine byte order before bits offset and masking
    """
    cdef char* bita = PyBytes_AsString(tmp)
    if 'V' in RecordFormat or 'S' in RecordFormat or RecordFormat is None:
        return dataReadByte(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, bitCount, bitOffset)
    elif signalDataType in (4, 5) and bitCount == 32:  # float
        if (byteorder == 'little' and signalDataType == 4) or \
                (byteorder == 'big' and signalDataType == 5):
            return dataReadFloat(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, 0)
        else: #  swap bytes
            return dataReadFloat(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, 1)
    elif signalDataType in (4, 5) and bitCount == 64:  # double
        if (byteorder == 'little' and signalDataType == 4) or \
                (byteorder == 'big' and signalDataType == 5):
            return dataReadDouble(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, 0)
        else: #  swap bytes
            return dataReadDouble(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, 1)
    elif signalDataType in (0, 1, 13) and bitCount <= 8:  # unsigned char
        return dataReadUChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 8:  # signed char
        return dataReadChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (0, 1, 13, 14) and bitCount <=16:  # unsigned short
        if (byteorder == 'little' and signalDataType == 0) or \
                (byteorder == 'big' and signalDataType == 1):
            return dataReadUShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadUShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    elif signalDataType in (2, 3) and bitCount <= 16:  # signed short
        if (byteorder == 'little' and signalDataType == 2) or \
                (byteorder == 'big' and signalDataType == 3):
            return dataReadShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    elif signalDataType in (0, 1, 14) and bitCount <=32:  # unsigned int
        if (byteorder == 'little' and signalDataType == 0) or \
                (byteorder == 'big' and signalDataType == 1):
            return dataReadUInt(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadUInt(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    elif signalDataType in (2, 3) and bitCount <= 32:  # signed int
        if (byteorder == 'little' and signalDataType == 2) or \
                (byteorder == 'big' and signalDataType == 3):
            return dataReadInt(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadInt(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    elif signalDataType in (0, 1) and bitCount <=64:  # unsigned long long
        if (byteorder == 'little' and signalDataType == 0) or \
                (byteorder == 'big' and signalDataType == 1):
            return dataReadULongLong(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadULongLong(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    elif signalDataType in (2, 3) and bitCount <= 64:  # signed long long
        if (byteorder == 'little' and signalDataType == 0) or \
                (byteorder == 'big' and signalDataType == 1):
            return dataReadLongLong(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 0)
        else: #  swap bytes
            return dataReadLongLong(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset, 1)
    else:
        return dataReadByte(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, bitCount, bitOffset)
        
cdef inline dataReadFloat(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned char swap):
    cdef np.ndarray[np.float32_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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
    
cdef inline dataReadDouble(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned char swap):
    cdef np.ndarray[np.float64_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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
    
cdef inline dataReadUChar(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray[np.uint8_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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

cdef inline dataReadChar(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray[np.int8_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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

cdef inline dataReadUShort(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.uint16_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bitCount) - 1)
    cdef unsigned short temp2byte = 0
    cdef char temp[2]
    if bitCount == 16:
        for i in range(numberOfRecords):
            memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
            buf[i] = temp2byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
            # right shift 
            if bitOffset > 0:
                temp2byte = temp2byte >> bitOffset
            # mask left part
            if bitCount < 16:
                temp2byte &= mask
            buf[i] = temp2byte
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline dataReadShort(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.int16_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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
    else:
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
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline dataReadUInt(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.uint32_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bitCount) - 1)
    cdef unsigned int temp4byte = 0
    cdef unsigned char nBytes = 0
    cdef char temp[3]
    if bitCount > 24:
        nBytes = 4
    else:
        nBytes = 3
    if nBytes == 4:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:  # on 3 bytes
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
                memcpy(&temp, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp[0]<<16 | temp[1]<<8 | temp[2]  #  swap bytes
                # right shift 
                if bitOffset > 0:
                    temp4byte = temp4byte >> bitOffset
                # mask left part
                if bitCount < 32:
                    temp4byte &= mask
                buf[i] = temp4byte
        return buf
    

cdef inline dataReadInt(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.int32_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef int mask = ((1 << bitCount) - 1)
    cdef int temp4byte = 0
    cdef unsigned char nBytes = 0
    cdef int signBit = 0
    cdef int signBitMask = (1 << (bitCount-1))
    cdef int signExtend = ((1 << (32 - bitCount)) - 1) << bitCount
    cdef char temp[3]
    if bitCount > 24:
        nBytes = 4
    else:
        nBytes = 3
    if nBytes == 4:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp4byte
        if swap == 0:
            return buf
        else:
            return buf.byteswap()
    else:  # on 3 bytes
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
                memcpy(&temp, &bita[posByteBeg + record_byte_size * i], nBytes)
                temp4byte = temp[0]<<16 | temp[1]<<8 | temp[2]  #  swap bytes
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

cdef inline dataReadULongLong(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.uint64_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bitCount) - 1)
    cdef unsigned long long temp8byte = 0
    cdef unsigned char nBytes = bitCount // 8
    cdef char temp[8]
    if bitCount % 8 > 0:
        nBytes += 1
    if bitCount == 64:
        for i in range(numberOfRecords):
            memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp8byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            # right shift 
            if bitOffset > 0:
                temp8byte = temp8byte >> bitOffset
            # mask left part
            if bitCount < 64:
                temp8byte &= mask
            buf[i] = temp8byte
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline dataReadLongLong(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset, unsigned char swap):
    cdef np.ndarray[np.int64_t] buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef long long mask = ((1 << bitCount) - 1)
    cdef long long temp8byte = 0
    cdef unsigned char nBytes = bitCount // 8
    cdef long signBit = 0
    cdef long long signBitMask = (1 << (bitCount-1))
    cdef long long signExtend = ((1 << (64 - bitCount)) - 1) << bitCount
    cdef char temp[8]
    if bitCount % 8 > 0:
        nBytes += 1
    if bitCount == 64:
        for i in range(numberOfRecords):
            memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp8byte
    else:
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
    if swap == 0:
        return buf
    else:
        return buf.byteswap()

cdef inline dataReadByte(const char* bita, str RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    for i in range(numberOfRecords):
            buf[i] = bytes(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])
    return buf
