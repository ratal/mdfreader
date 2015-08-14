import numpy as np
cimport numpy as np
import sys
#cimport cython

cdef extern from *:
    ctypedef void const_void "const void"
cdef extern from "string.h" nogil:
    void *memcpy  (void *TO, const_void *FROM, size_t SIZE)

#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(const char[:] bita, unsigned short bitCount, \
        unsigned short signalDataType, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned char bitOffset, \
        unsigned long posByteBeg, unsigned long posByteEnd):
    order = sys.byteorder
    if signalDataType in (4, 5) and bitCount == 32:  # float
        if (order == 'little' and signalDataType == 4) or \
                (order == 'big' and signalDataType == 5):
            return dataReadFloat(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg)
        else: #  swap bytes
            return dataReadFloatSwapped(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg)
    elif signalDataType in (4, 5) and bitCount == 64:  # double
        if (order == 'little' and signalDataType == 4) or \
                (order == 'big' and signalDataType == 5):
            return dataReadDouble(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg)
        else:
            return dataReadDoubleSwapped(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg)
    elif signalDataType in (0, 1) and bitCount <= 8:  # unsigned char
        return dataReadUChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 8:  # signed char
        return dataReadChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=16:  # unsigned short
        if (order == 'little' and signalDataType == 0) or \
                (order == 'big' and signalDataType == 1):
            return dataReadUShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset)
        else: #  swap bytes
            return dataReadUShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 16:  # signed short
        if (order == 'little' and signalDataType == 2) or \
                (order == 'big' and signalDataType == 3):
            return dataReadShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset)
        else:
            return dataReadShort(bita, RecordFormat, numberOfRecords, \
                record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=32:  # unsigned int
        return dataReadUInt(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 32:  # signed int
        return dataReadInt(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=64:  # unsigned long long
        return dataReadULongLong(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 64:  # signed long long
        return dataReadLongLong(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, bitCount, bitOffset)
    else:
        return dataReadByte(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, bitCount, bitOffset)
        
cdef inline dataReadFloat(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef float tempfloat = 0
    for i in range(numberOfRecords):
        memcpy(&tempfloat, &bita[posByteBeg + record_byte_size * i], 4)
        buf[i] = tempfloat
    return buf

cdef inline dataReadFloatSwapped(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef float tempfloat = 0
    cdef char temp[4]
    for i in range(numberOfRecords):
        memcpy(&temp, &bita[posByteBeg + record_byte_size * i], 4)
        tempfloat = temp[0]<<24 | temp[1]<<16 | temp[2]<<8 | temp[3] #  swap bytes
        buf[i] = tempfloat
    return buf
    
cdef inline dataReadDouble(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef double tempDouble = 0
    for i in range(numberOfRecords):
        memcpy(&tempDouble, &bita[posByteBeg + record_byte_size * i], 8)
        buf[i] = tempDouble
    return buf

cdef inline dataReadDoubleSwapped(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef double tempDouble = 0
    cdef char temp[8]
    for i in range(numberOfRecords):
        memcpy(&temp, &bita[posByteBeg + record_byte_size * i], 8)
        tempDouble = temp[0]<<56 | temp[1]<<48 | temp[2]<<40 | temp[3]<<32 | \
                temp[4]<<24 | temp[5]<<16 | temp[6]<<8 | temp[7] #  swap bytes
        buf[i] = tempDouble
    return buf
    
cdef inline dataReadUChar(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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

cdef inline dataReadChar(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
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

cdef inline dataReadUShort(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned short mask = ((1 << bitCount) - 1)
    cdef unsigned short temp2byte = 0
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
    return buf

cdef inline dataReadShort(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef short mask = ((1 << bitCount) - 1)
    cdef short temp2byte = 0
    cdef short signBit = 0
    cdef short signBitMask = (1 << (bitCount-1))
    cdef short signExtend = ((1 << (16 - bitCount)) - 1) << bitCount
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
    return buf

cdef inline dataReadUInt(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned int mask = ((1 << bitCount) - 1)
    cdef unsigned int temp4byte = 0
    cdef unsigned char nBytes = 0
    if bitCount > 24:
        nBytes = 4
    else:
        nBytes = 3
    if bitCount == 32:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp4byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            # right shift 
            if bitOffset > 0:
                temp4byte = temp4byte >> bitOffset
            # mask left part
            if bitCount < 32:
                temp4byte &= mask
            buf[i] = temp4byte
    return buf

cdef inline dataReadInt(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef int mask = ((1 << bitCount) - 1)
    cdef int temp4byte = 0
    cdef unsigned char nBytes = 0
    cdef int signBit = 0
    cdef int signBitMask = (1 << (bitCount-1))
    cdef int signExtend = ((1 << (32 - bitCount)) - 1) << bitCount
    if bitCount > 24:
        nBytes = 4
    else:
        nBytes = 3
    if bitCount == 32:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            buf[i] = temp4byte
    else:
        for i in range(numberOfRecords):
            memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], nBytes)
            # right shift 
            if bitOffset > 0:
                temp4byte = temp4byte >> bitOffset
            # mask left part
            if bitCount < 32:
                temp4byte &= mask
            signBit = temp4byte & signBitMask
            if signBit: #  negative value, sign extend
                temp4byte |= signExtend
            buf[i] = temp4byte
    return buf

cdef inline dataReadULongLong(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long long mask = ((1 << bitCount) - 1)
    cdef unsigned long long temp8byte = 0
    cdef unsigned char nBytes = bitCount // 8
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
    return buf

cdef inline dataReadLongLong(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef long long mask = ((1 << bitCount) - 1)
    cdef long long temp8byte = 0
    cdef unsigned char nBytes = bitCount // 8
    cdef long signBit = 0
    cdef long long signBitMask = (1 << (bitCount-1))
    cdef long long signExtend = ((1 << (64 - bitCount)) - 1) << bitCount
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
    return buf

cdef inline dataReadByte(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long bitCount, unsigned char bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    for i in range(numberOfRecords):
            buf[i] = bytes(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])
    return buf
