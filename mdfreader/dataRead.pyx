import numpy as np
cimport numpy as np
#cimport cython

#cdef extern from "Python.h":
#    char* PyByteArray_AsString(object bytearray) except NULL
#    char* PyBytes_AsString(object bytearray) except NULL
cdef extern from *:
    ctypedef void const_void "const void"
cdef extern from "string.h" nogil:
    void *memcpy  (void *TO, const_void *FROM, size_t SIZE)

#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(const char[:] bita, unsigned short bitCount, \
        unsigned short signalDataType, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, \
        unsigned long bitOffset, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd):
    # initialise variables
    cdef unsigned long trailBits = posByteEnd * 8 - posBitEnd
    #cdef const char* bita = PyByteArray_AsString(temp)
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
        
cdef inline dataReadFloat(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef long mask = ((1 << maskBit) - 1)
    cdef long temp4bytes = 0
    cdef float tempfloat = 0
    if trailBits == 0 and bitOffset==0:
        for i in range(numberOfRecords):
            memcpy(&temp4bytes, &bita[posByteBeg + record_byte_size * i], 4)
            buf[i] = temp4bytes
    else:
        for i in range(numberOfRecords):
            memcpy(&temp4bytes, &bita[posByteBeg + record_byte_size * i], 4)
            if trailBits > 0:
                temp4bytes = temp4bytes & mask
            # left shift 
            if bitOffset > 0:
                temp4bytes = temp4bytes << bitOffset
            memcpy(&tempfloat, &temp4bytes, 4)
            buf[i] = tempfloat
    return buf

cdef inline dataReadDouble(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef long temp8bytes = 0
    for i in range(numberOfRecords):
        memcpy(&temp8bytes, &bita[posByteBeg + record_byte_size * i], 8)
        buf[i] = temp8bytes
    return buf

cdef inline dataReadUChar(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned char mask = ((1 << maskBit) - 1)
    cdef unsigned char temp1byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
        # mask right part
        if trailBits > 0:
            temp1byte = temp1byte & mask
        # left shift 
        if bitOffset > 0:
            temp1byte = temp1byte << bitOffset
        buf[i] = temp1byte
    return buf

cdef inline dataReadChar(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef char mask = ((1 << maskBit) - 1)
    cdef char temp1byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp1byte, &bita[posByteBeg + record_byte_size * i], 1)
        # mask right part
        if trailBits > 0:
            temp1byte = temp1byte & mask
        # left shift 
        if bitOffset > 0:
            temp1byte = temp1byte << bitOffset
        buf[i] = temp1byte
    return buf

cdef inline dataReadUShort(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned short mask = ((1 << maskBit) - 1)
    cdef unsigned short temp2byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
        # mask right part
        if trailBits > 0:
            temp2byte = temp2byte & mask
        # left shift 
        if bitOffset > 0:
            temp2byte = temp2byte << bitOffset
        buf[i] = temp2byte
    return buf

cdef inline dataReadShort(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef short mask = ((1 << maskBit) - 1)
    cdef short temp2byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp2byte, &bita[posByteBeg + record_byte_size * i], 2)
        # mask right part
        if trailBits > 0:
            temp2byte = temp2byte & mask
        # left shift 
        if bitOffset > 0:
            temp2byte = temp2byte << bitOffset
        buf[i] = temp2byte
    return buf

cdef inline dataReadUInt(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned int mask = ((1 << maskBit) - 1)
    cdef unsigned int temp4byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], 4)
        # mask right part
        if trailBits > 0:
            temp4byte = temp4byte & mask
        # left shift 
        if bitOffset > 0:
            temp4byte = temp4byte << bitOffset
        buf[i] = temp4byte
    return buf

cdef inline dataReadInt(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef int mask = ((1 << maskBit) - 1)
    cdef int temp4byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp4byte, &bita[posByteBeg + record_byte_size * i], 4)
        # mask right part
        if trailBits > 0:
            temp4byte = temp4byte & mask
        # left shift 
        if bitOffset > 0:
            temp4byte = temp4byte << bitOffset
        buf[i] = temp4byte
    return buf

cdef inline dataReadULongLong(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned long long mask = ((1 << maskBit) - 1)
    cdef unsigned long long temp8byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], 8)
        # mask right part
        if trailBits > 0:
            temp8byte = temp8byte & mask
        # left shift 
        if bitOffset > 0:
            temp8byte = temp8byte << bitOffset
        buf[i] = temp8byte
    return buf

cdef inline dataReadLongLong(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef long long mask = ((1 << maskBit) - 1)
    cdef long long temp8byte = 0
    for i in range(numberOfRecords):
        memcpy(&temp8byte, &bita[posByteBeg + record_byte_size * i], 8)
        # mask right part
        if trailBits > 0:
            temp8byte = temp8byte & mask
        # left shift 
        if bitOffset > 0:
            temp8byte = temp8byte << bitOffset
        buf[i] = temp8byte
    return buf

cdef inline dataReadByte(const char[:] bita, RecordFormat, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd, unsigned long trailBits, unsigned long bitOffset):
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array
    cdef unsigned long long i
    for i in range(numberOfRecords):
            buf[i] = bytes(bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i])
    return buf
