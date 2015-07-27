from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free

#@cython.boundscheck(False)
#@cython.wraparound(False)
def dataRead(char* bita, np.dtype RecordFormat, \
        char Format, unsigned long long numberOfRecords, \
        unsigned long record_byte_size, \
        unsigned long bitOffset, unsigned long posByteBeg, unsigned long posByteEnd, \
        unsigned long posBitEnd):
    """ reads stream of record bytes needed for not byte aligned data
    
    Parameters
    ------------
    bita : stream
        stream of bytes
    channelList : List of str, optional
        list of channel to read
    
    Returns
    --------
    rec : numpy ndarray
        contains a matrix of raw data in a ndarray (attributes corresponding to channel name)
    """
    # initialise variables
    cdef np.ndarray buf = np.zeros(numberOfRecords, dtype=RecordFormat)  # return numpy array 
    cdef unsigned long trailBits = posByteEnd * 8 - posBitEnd
    cdef unsigned int nbytes = posByteEnd - posByteBeg  # size of each values
    cdef char *mask=<char *>malloc(nbytes)
    cdef unsigned long maskBit =  posBitEnd - posByteBeg * 8
    cdef unsigned long long i
    cdef signed char temp1byte
    cdef unsigned char utemp1byte
    cdef short temp2bytes
    cdef unsigned short utemp2bytes
    cdef long temp4bytes
    cdef unsigned long utemp4bytes
    cdef long long temp8bytes
    cdef unsigned long long utemp8bytes
    # slice stream in array
    if Format == 'f':  # float
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp4bytes = <long>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <long>mask
            # left shift 
            if bitOffset > 0:
                temp4bytes = temp4bytes << bitOffset
            buf[i] = <float>temp4bytes
    elif Format == 'd':  # double
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp8bytes = <long long>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <long long>mask
            # left shift 
            if bitOffset > 0:
                temp8bytes = temp8bytes << bitOffset
            buf[i] = <double>temp8bytes
    elif Format == 'H':  # unsigned short
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                utemp2bytes = <unsigned short>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <unsigned short>mask
            # left shift 
            if bitOffset > 0:
                utemp2bytes = utemp2bytes << bitOffset
            buf[i] = utemp2bytes
    elif Format == 'h':  # signed short
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp2bytes = <short>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <short>mask
            # left shift 
            if bitOffset > 0:
                temp2bytes = temp2bytes << bitOffset
            buf[i] = temp2bytes
    elif Format == 'B':  # unsigned char
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp1byte = <unsigned char>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <unsigned char>mask
            # left shift 
            if bitOffset > 0:
                temp1byte = temp1byte << bitOffset
            buf[i] = temp1byte
    elif Format == 'b':  # signed char
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                utemp1byte = <signed char>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <signed char>mask
            # left shift 
            if bitOffset > 0:
                utemp1byte = utemp1byte << bitOffset
            buf[i] = utemp1byte
    elif Format == 'I':  # unsigned int
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                utemp4bytes = <unsigned int>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <unsigned int>mask
            # left shift 
            if bitOffset > 0:
                utemp4bytes = utemp4bytes << bitOffset
            buf[i] = utemp4bytes
    elif Format == 'i':  # signed int
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp4bytes = <int>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <int>mask
            # left shift 
            if bitOffset > 0:
                temp4bytes = temp4bytes << bitOffset
            buf[i] = temp4bytes
    elif Format == 'Q':  # unsigned long long
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                utemp8bytes = <unsigned long long>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <unsigned long long>mask
            # left shift 
            if bitOffset > 0:
                utemp8bytes = utemp8bytes << bitOffset
            buf[i] = utemp8bytes
    elif Format == 'q':  # signed long long
        for i in range(numberOfRecords):
            # mask right part
            if trailBits > 0:
                mask = <char *>((1 << maskBit) - 1)
                temp8bytes = <long long>bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i] & <long long>mask
            # left shift 
            if bitOffset > 0:
                temp8bytes = temp8bytes << bitOffset
            buf[i] = temp8bytes
    else:
        for i in range(numberOfRecords):
            buf[i] = bita[posByteBeg + record_byte_size * i:\
                    posByteEnd + record_byte_size * i]
    free(mask)
    return buf
    
