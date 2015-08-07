import numpy as np
cimport numpy as np
cimport cython
from struct import unpack, Struct
#cdef extern from "Python.h":
#    char* PyByteArray_AsString(object bytearray) except NULL
#cdef extern from "stdint.h":
#    ctypedef char int8_t
#    ctypedef unsigned char uint8_t
#    ctypedef short int16_t
#    ctypedef unsigned short uint16_t
#    ctypedef long int32_t
#    ctypedef unsigned long uint32_t
#    ctypedef long long int64_t
#    ctypedef unsigned long long uint64_t

#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(unsigned char[:] bita, unsigned short bitCount, \
        unsigned short signalDataType, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, \
        unsigned long bitOffset, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd):
    # initialise variables
    cdef unsigned long trailBits = posByteEnd * 8 - posBitEnd
    #cdef char* temp = PyByteArray_AsString(bita)
    # slice stream in array
    if signalDataType in (4, 5) and bitCount == 32:  # float
        return dataReadFloat(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (4, 5) and bitCount == 64:  # double
        return dataReadDouble(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (0, 1) and bitCount <= 8:  # unsigned char
        return dataReadUChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 8:  # signed char
        return dataReadChar(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=16:  # unsigned short
        return dataReadUShort(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 16:  # signed short
        return dataReadShort(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=32:  # unsigned int
        return dataReadUInt(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 32:  # signed int
        return dataReadInt(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (0, 1) and bitCount <=64:  # unsigned long long
        return dataReadULongLong(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    elif signalDataType in (2, 3) and bitCount <= 64:  # signed long long
        return dataReadLongLong(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
    else:
        return dataReadByte(bita, RecordFormat, numberOfRecords, \
            record_byte_size, posByteBeg, posByteEnd, posBitEnd, trailBits, bitOffset)
        
cdef dataReadFloat(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef long mask = ((1 << maskBit) - 1)
    cdef long temp4bytes = 0
    if trailBits == 0 and bitOffset==0:
        CFormat = Struct('f')
        for i in range(numberOfRecords):
            buf[i] = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
    else:
        CFormat = Struct('i')
        for i in range(numberOfRecords):
            temp4bytes = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
            if trailBits > 0:
                temp4bytes = temp4bytes & mask
            # left shift 
            if bitOffset > 0:
                temp4bytes = temp4bytes << bitOffset
            buf[i] = <float>temp4bytes
    return buf

cdef dataReadDouble(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef long long mask = ((1 << maskBit) - 1)
    CFormat = Struct('d')
    for i in range(numberOfRecords):
        buf[i] = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                posByteEnd + record_byte_size * i])[0]
    return buf

cdef dataReadUChar(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned char mask = ((1 << maskBit) - 1)
    cdef unsigned char temp1byte = 0
    CFormat = Struct('B')
    for i in range(numberOfRecords):
        temp1byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp1byte = temp1byte & mask
        # left shift 
        if bitOffset > 0:
            temp1byte = temp1byte << bitOffset
        buf[i] = temp1byte
    return buf

cdef dataReadChar(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef char mask = ((1 << maskBit) - 1)
    cdef char temp1byte = 0
    CFormat = Struct('b')
    for i in range(numberOfRecords):
        temp1byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp1byte = temp1byte & mask
        # left shift 
        if bitOffset > 0:
            temp1byte = temp1byte << bitOffset
        buf[i] = temp1byte
    return buf

cdef dataReadUShort(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned short mask = ((1 << maskBit) - 1)
    cdef unsigned short temp2byte = 0
    CFormat = Struct('H')
    for i in range(numberOfRecords):
        temp2byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp2byte = temp2byte & mask
        # left shift 
        if bitOffset > 0:
            temp2byte = temp2byte << bitOffset
        buf[i] = temp2byte
    return buf

cdef dataReadShort(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef short mask = ((1 << maskBit) - 1)
    cdef short temp2byte = 0
    CFormat = Struct('h')
    for i in range(numberOfRecords):
        temp2byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp2byte = temp2byte & mask
        # left shift 
        if bitOffset > 0:
            temp2byte = temp2byte << bitOffset
        buf[i] = temp2byte
    return buf

cdef dataReadUInt(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned int mask = ((1 << maskBit) - 1)
    cdef unsigned int temp4byte = 0
    CFormat = Struct('I')
    for i in range(numberOfRecords):
        temp4byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp4byte = temp4byte & mask
        # left shift 
        if bitOffset > 0:
            temp4byte = temp4byte << bitOffset
        buf[i] = temp4byte
    return buf

cdef dataReadInt(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef int mask = ((1 << maskBit) - 1)
    cdef int temp4byte = 0
    CFormat = Struct('i')
    for i in range(numberOfRecords):
        temp4byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp4byte = temp4byte & mask
        # left shift 
        if bitOffset > 0:
            temp4byte = temp4byte << bitOffset
        buf[i] = temp4byte
    return buf

cdef dataReadULongLong(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned long long mask = ((1 << maskBit) - 1)
    cdef unsigned long long temp8byte = 0
    CFormat = Struct('Q')
    for i in range(numberOfRecords):
        temp8byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp8byte = temp8byte & mask
        # left shift 
        if bitOffset > 0:
            temp8byte = temp8byte << bitOffset
        buf[i] = temp8byte
    return buf

cdef dataReadLongLong(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef long long mask = ((1 << maskBit) - 1)
    cdef long long temp4byte = 0
    CFormat = Struct('q')
    for i in range(numberOfRecords):
        temp8byte = CFormat.unpack(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])[0]
        # mask right part
        if trailBits > 0:
            temp8byte = temp8byte & mask
        # left shift 
        if bitOffset > 0:
            temp8byte = temp8byte << bitOffset
        buf[i] = temp8byte
    return buf

cdef dataReadByte(unsigned char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    for i in range(numberOfRecords):
            buf[i] = bytes(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])
    return buf
