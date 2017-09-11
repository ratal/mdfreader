# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 3.x

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Sun Oct 10 12:57:28 2010

Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>
- Sympy to convert channels with formula

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.2+

mdf3reader module
--------------------------
"""
from __future__ import print_function

from numpy import average, right_shift, bitwise_and, diff, max, min, interp
from numpy import asarray, zeros, recarray, array, searchsorted
from numpy.core.records import fromfile, fromarrays
from numpy.core.defchararray import encode as ncode
from mdf import mdf_skeleton, _open_MDF, _bits_to_bytes, _convertName,\
    dataField, conversionField, compressed_data
from mdfinfo3 import info3
from math import log, exp
from time import strftime, time, gmtime
from struct import pack, Struct
from io import open  # for python 3 and 2 consistency
from sys import platform, exc_info, version_info, stderr, path
from os.path import dirname, abspath
import os
_root = dirname(abspath(__file__))
path.append(_root)

PythonVersion = version_info
PythonVersion = PythonVersion[0]


def linearConv(data, conv):  # 0 Parametric, Linear: Physical =Integer*P2 + P1
    """ apply linear conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conv['P2'] == 1.0 and conv['P1'] in (0.0, -0.0):
        return data  # keeps dtype probably more compact than float64
    else:
        return data * conv['P2'] + conv['P1']


def tabInterpConv(data, conv):  # 1 Tabular with interpolation
    """ apply Tabular interpolation conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    tmp = array([(key, val['int'], val['phys'])
                for (key, val) in conv.items()])
    return interp(data, tmp[:, 1], tmp[:, 2])


def tabConv(data, conv):  # 2 Tabular
    """ apply Tabular conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    tmp = array([(key, val['int'], val['phys'])
                for (key, val) in conv.items()])
    indexes = searchsorted(tmp[:, 1], data)
    return tmp[indexes, 2]


def polyConv(data, conv):  # 6 Polynomial
    """ apply polynomial conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    return (conv['P2'] - conv['P4'] * (data - conv['P5'] - conv['P6'])) \
        / (conv['P3'] * (data - conv['P5'] - conv['P6']) - conv['P1'])


def expConv(data, conv):  # 7 Exponential
    """ apply exponential conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conv['P4'] == 0 and conv['P1'] != 0 and conv['P2'] != 0:
        return exp(((data - conv['P7']) * conv['P6'] - conv['P3'])
                   / conv['P1']) / conv['P2']
    elif conv['P1'] == 0 and conv['P4'] != 0 and conv['P5'] != 0:
        return exp((conv['P3'] / (data - conv['P7']) - conv['P6'])
                   / conv['P4']) / conv['P5']
    else:
        print('Non possible conversion parameters for channel ', file=stderr)


def logConv(data, conv):  # 8 Logarithmic
    """ apply logarithmic conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conv['P4'] == 0 and conv['P1'] != 0 and conv['P2'] != 0:
        return log(((data - conv['P7']) * conv['P6'] - conv['P3'])
                   / conv['P1']) / conv['P2']
    elif conv['P1'] == 0 and conv['P4'] != 0 and conv['P5'] != 0:
        return log((conv['P3'] / (data - conv['P7']) - conv['P6'])
                   / conv['P4']) / conv['P5']
    else:
        print('Non possible conversion parameters for channel ', file=stderr)


def rationalConv(data, conv):  # 9 rational
    """ apply rational conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    return (conv['P1'] * data * data + conv['P2'] * data + conv['P3'])\
        / (conv['P4'] * data * data + conv['P5'] * data + conv['P6'])


def formulaConv(data, conv):  # 10 Text Formula
    """ apply formula conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value

    Notes
    --------
    Requires sympy module
    """
    try:
        from sympy import lambdify, symbols
        X = symbols('X')  # variable is X
        formula = conv['textFormula']
        # remove trailing text after 0
        formula = formula[:formula.find('\x00')]
        # adapt ASAM-MCD2 syntax to sympy
        formula = formula.replace('pow(', 'power(')
        # formula to function for evaluation
        expr = lambdify(X, formula, modules='numpy', dummify=False)
        return expr(data)
    except:
        print('Please install sympy to convert channel ', file=stderr)
        print('Failed to convert formulae ' + conv['textFormula'], file=stderr)


def textRangeTableConv(data, conv):  # 12 Text range table
    """ apply text range table conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conv : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    try:
        npair = len(conv)
        lower = [conv[pair]['lowerRange'] for pair in range(npair)]
        upper = [conv[pair]['upperRange'] for pair in range(npair)]
        text = [conv[pair]['Textrange'] for pair in range(npair)]
        temp = []
        for Lindex in range(len(data)):
            value = text[0]  # default value
            for pair in range(1, npair):
                if lower[pair] <= data[Lindex] <= upper[pair]:
                    value = text[pair]
                    break
            temp.append(value)
        try:
            temp = asarray(temp)  # try to convert to numpy
        except:
            pass
        return temp
    except:
        print('Failed to convert text to range table', file=stderr)


class Channel():

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
        (self.dataFormat, nativeDataType) = \
            _arrayformat3(self.signalDataType, self.bitCount, ByteOrder[0])
        self.CFormat = Struct(_datatypeformat3(self.signalDataType, self.bitCount, ByteOrder[0]))
        self.nBytes = _bits_to_bytes(self.bitCount)
        self.posBitBeg = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfTheFirstBits']
        self.posBitEnd = self.posBitBeg + self.bitCount
        self.byteOffset = self.posBitBeg // 8  # + info['CNBlock'][dataGroup][channelGroup][channelNumber]['ByteOffset']
        self.bitOffset = self.posBitBeg % 8
        self.embedding_channel_bitOffset = self.bitOffset  # for channel containing other channels
        self.posByteBeg = recordIDnumber + self.byteOffset
        self.posByteEnd = recordIDnumber + self.byteOffset + self.nBytes
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['channelType']
        if self.channelType:  # master time channel
            self.name +=  '_' + str(dataGroup)
        self.recAttributeName = _convertName(self.name)
        self.RecordFormat = ((self.recAttributeName + '_title', self.recAttributeName), self.dataFormat)
        self.nativeRecordFormat = ((self.recAttributeName + '_title', self.recAttributeName), nativeDataType)
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
        output = str(self.channelNumber) + ' '
        output += self.name + ' '
        output += str(self.signalDataType) + ' '
        output += str(self.channelType) + ' '
        output += str(self.RecordFormat) + ' '
        output += str(self.bitOffset) + ' '
        output += str(self.byteOffset)
        output += 'unit ' + self.unit
        output += 'description ' + self.desc
        return output

    def changeChannelName(self, channelGroup):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channelGroup : int
            channelGroup bumber
        """
        self.name += '_' + str(channelGroup)
        self.recAttributeName = _convertName(self.name)
        self.RecordFormat = ((self.recAttributeName + '_title', self.recAttributeName), self.dataFormat)


class record(list):

    """ record class lists Channel classes,
        it is representing a channel group

    Attributes
    --------------
    CGrecordLength : int
        length of record from channel group block information in Byte
    recordLength : int
        length of record from channels information in Byte
    numberOfRecords : int
        number of records in data block
    recordID : int
        recordID corresponding to channel group
    recordIDnumber : int
        size of recordID
    dataGroup : int:
        data group number
    channelGroup : int
        channel group number
    numpyDataRecordFormat : list
        list of numpy (dtype) for each channel
    dataRecordName : list
        list of channel names used for recarray attribute definition
    master : dict
        define name and number of master channel
    recordToChannelMatching : dict
        helps to identify nested bits in byte
    channelNames : set
        channel names to be stored, useful for low memory consumption but slow
    hiddenBytes : Bool, False by default
        flag in case of non declared channels in record
    byte_aligned : Bool, True by default
        flag for byte aligned record

    Methods
    ------------
    addChannel(info, channelNumber)
    loadInfo(info)
    readSortedRecord(fid, pointer, channelSet=None)
    readRecordBuf(buf, channelSet=None)
    readRecordBits(bita, channelSet=None)
    changeChannelName(channelName)
    """

    def __init__(self, dataGroup, channelGroup):
        self.CGrecordLength = 0
        self.recordLength = 0
        self.dataBlockLength = 0
        self.numberOfRecords = 0
        self.recordID = 0
        self.recordIDnumber = 0
        self.dataGroup = dataGroup
        self.channelGroup = channelGroup
        self.numpyDataRecordFormat = []
        self.dataRecordName = []
        self.master = {}
        self.master['name'] = 'master_' + str(dataGroup)
        self.master['number'] = None
        self.recordToChannelMatching = {}
        self.channelNames = set()
        self.hiddenBytes = False
        self.byte_aligned = True

    def __repr__(self):
        output = 'Channels :\n'
        for chan in self.channelNames:
            output += chan + '\n'
        output += 'Datagroup number : ' + str(self.dataGroup) + '\n'
        if self.master['name'] is not None:
            output += 'Master channel : ' + self.master['name'] + '\n'
        output += 'Numpy records format : \n'
        for record in self.numpyDataRecordFormat:
            output += str(record) + '\n'
        return output

    def addChannel(self, info, channelNumber):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo3.info3 class
        channelNumber : int
            channel number in mdfinfo3.info3 class

        """
        self.append(Channel(info, self.dataGroup, self.channelGroup,
                    channelNumber, self.recordIDnumber))
        self.channelNames.add(self[-1].name)

    def loadInfo(self, info):
        """ gathers records related from info class

        Parameters
        ----------------
        info : mdfinfo3.info3 class

        """
        self.recordIDnumber = info['DGBlock'][self.dataGroup]['numberOfRecordIDs']
        self.recordID = info['CGBlock'][self.dataGroup][self.channelGroup]['recordID']
        self.CGrecordLength = info['CGBlock'][self.dataGroup][self.channelGroup]['dataRecordSize']
        self.numberOfRecords = info['CGBlock'][self.dataGroup][self.channelGroup]['numberOfRecords']
        self.dataBlockLength = self.CGrecordLength * self.numberOfRecords
        if self.recordIDnumber > 0:  # record ID existing at beginning of record
            self.dataRecordName.append('RecordID' + str(self.channelGroup))
            format = (self.dataRecordName[-1] + '_title', self.dataRecordName[-1])
            self.numpyDataRecordFormat.append((format, 'uint8'))
            self.dataBlockLength = (self.CGrecordLength + 1) * self.numberOfRecords
        embedding_channel = None
        for channelNumber in range(info['CGBlock'][self.dataGroup][self.channelGroup]['numberOfChannels']):
            channel = Channel(info, self.dataGroup, self.channelGroup,
                              channelNumber, self.recordIDnumber)
            if self.master['number'] is None or channel.channelType == 1:  # master channel found
                self.master['name'] = channel.name
                self.master['number'] = channelNumber
            self.append(channel)  # adds channel in record list
            self.channelNames.add(channel.name)
            # Checking if several channels are embedded in bytes
            if len(self) > 1:
                # all channels are already ordered in record based on byte_offset and bit_offset
                # so just comparing with previous channel
                prev_chan = self[-2]
                prev_chan_includes_curr_chan = channel.posBitBeg >= 8 * prev_chan.byteOffset \
                        and channel.posBitEnd <= 8 * (prev_chan.byteOffset + prev_chan.nBytes)
                if embedding_channel is not None:
                    embedding_channel_includes_curr_chan = channel.posBitEnd <= embedding_channel.posBitEnd
                else:
                    embedding_channel_includes_curr_chan = False
                if channel.byteOffset >= prev_chan.byteOffset and \
                        channel.posBitBeg < 8 * (prev_chan.byteOffset + prev_chan.nBytes) and \
                        channel.posBitEnd > 8 * (prev_chan.byteOffset + prev_chan.nBytes):
                    # not byte aligned
                    self.byte_aligned = False
                if embedding_channel is not None and \
                        channel.posBitEnd > embedding_channel.posBitEnd:
                    embedding_channel = None
                if prev_chan_includes_curr_chan or \
                        embedding_channel_includes_curr_chan:  # bit(s) in byte(s)
                    if embedding_channel is None and prev_chan_includes_curr_chan:
                        embedding_channel = prev_chan  # new embedding channel detected
                    if self.recordToChannelMatching:  # not first channel
                        self.recordToChannelMatching[channel.recAttributeName] = \
                            self.recordToChannelMatching[prev_chan.recAttributeName]
                        channel.embedding_channel_bitOffset = channel.posBitBeg - embedding_channel.posBitBeg
                    else:  # first channels
                        self.recordToChannelMatching[channel.recAttributeName] = \
                            channel.recAttributeName
                        self.numpyDataRecordFormat.append(channel.RecordFormat)
                        self.dataRecordName.append(channel.recAttributeName)
                        self.recordLength += channel.nBytes
            if embedding_channel is None:  # adding bytes
                self.recordToChannelMatching[channel.recAttributeName] = \
                    channel.recAttributeName
                self.numpyDataRecordFormat.append(channel.RecordFormat)
                self.dataRecordName.append(channel.recAttributeName)
                self.recordLength += channel.nBytes
        if self.recordIDnumber == 2:  # second record ID at end of record
            self.dataRecordName.append('RecordID' + str(self.channelGroup) + '_2')
            format = (self.dataRecordName[-1] + '_title', self.dataRecordName[-1])
            self.numpyDataRecordFormat.append((format, 'uint8'))
            self.dataBlockLength = (self.CGrecordLength + 2) * self.numberOfRecords
        # check for hidden bytes
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True
        # check record length consitency
        elif self.CGrecordLength < self.recordLength:
            # forces to use dataRead instead of numpy records.
            self.byte_aligned = False

    def readSortedRecord(self, fid, pointer, channelSet=None):
        """ reads record, only one channel group per datagroup

        Parameters
        ----------------
        fid : float
            file identifier
        pointer
            position in file of data block beginning
        channelSet : Set of str, optional
            list of channel to read

        Returns
        -----------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes
            corresponding to channel name)

        Notes
        --------
        If channelSet is None, read data using numpy.core.records.fromfile
        that is rather quick. However, in case of large file, you can use
        channelSet to load only interesting channels or only one channel
        on demand, but be aware it might be much slower.

        """
        fid.seek(pointer)
        if channelSet is None and not self.hiddenBytes and self.byte_aligned:
            # reads all, quickest but memory consuming
            return fromfile(fid, dtype=self.numpyDataRecordFormat,
                            shape=self.numberOfRecords,
                            names=self.dataRecordName)
        else:  # reads only some channels from a sorted data block
            if channelSet is None:
                channelSet = self.channelNames
            # memory efficient but takes time
            # are channelSet in this dataGroup
            if len(channelSet & self.channelNames) > 0:
                # check if master channel is in the list
                if not self.master['name'] in channelSet:
                    channelSet.add(self.master['name'])  # adds master channel
                try:  # use rather cython compiled code for performance
                    from dataRead import dataRead
                    # converts data type from mdf 3.x to 4.x
                    convertDataType3to4 = {0: 0, 1: 2, 2: 4, 3: 4,
                                           7: 6, 8: 10,
                                           9: 1, 10: 3, 11: 5, 12: 5,
                                           13: 0, 14: 2, 15: 4, 16: 4}
                    bita = fid.read(self.dataBlockLength)
                    format = []
                    for channel in self:
                        if channel.recAttributeName in channelSet:
                            format.append(channel.nativeRecordFormat)
                    buf = recarray(self.numberOfRecords, format)
                    for chan in range(len(self)):
                        if self[chan].recAttributeName in channelSet:
                            buf[self[chan].recAttributeName] = \
                                dataRead(bytes(bita),
                                         self[chan].bitCount,
                                         convertDataType3to4[self[chan].signalDataType],
                                         self[chan].nativeRecordFormat[1],
                                         self.numberOfRecords,
                                         self.CGrecordLength,
                                         self[chan].bitOffset,
                                         self[chan].posByteBeg,
                                         self[chan].posByteEnd)
                    return buf
                except:
                    print('Unexpected error:', exc_info(), file=stderr)
                    print('dataRead crashed, back to python data reading',
                          file=stderr)
                    rec = {}
                    recChan = []
                    numpyDataRecordFormat = []
                    for channel in channelSet:  # initialise data structure
                        rec[channel] = 0
                    for channel in self:  # list of Channels from channelSet
                        if channel.recAttributeName in channelSet:
                            recChan.append(channel)
                            numpyDataRecordFormat.append(channel.RecordFormat)
                    rec = zeros((self.numberOfRecords, ), dtype=numpyDataRecordFormat)
                    recordLength = self.recordIDnumber + self.CGrecordLength
                    for r in range(self.numberOfRecords):  # for each record,
                        buf = fid.read(recordLength)
                        for channel in recChan:
                            rec[channel.recAttributeName][r] = \
                                channel.CFormat.unpack(buf[channel.posByteBeg:
                                                           channel.posByteEnd])[0]
                    return rec.view(recarray)

    def readRecordBuf(self, buf, channelSet=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        channelSet : Set of str, optional
            list of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        temp = {}
        if channelSet is None:
            channelSet = self.channelNames
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channelSet:
                temp[self.recordToChannelMatching[Channel.recAttributeName]] = \
                    Channel.CFormat.unpack(buf[Channel.posByteBeg:
                                               Channel.posByteEnd])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def readRecordBits(self, bita, channelSet=None):
        """ read stream of record bits by bits in case of not aligned or hidden bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        channelSet : Set of str, optional
            list of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        from bitarray import bitarray
        B = bitarray(endian="little")  # little endian by default
        B.frombytes(bytes(bita))

        def signedInt(temp, extension):
            """ extend bits of signed data managing two's complement
            """
            extension.setall(False)
            extensionInv = bitarray(extension, endian='little')
            extensionInv.setall(True)
            signBit = temp[-1]
            if not signBit:  # positive value, extend with 0
                temp.extend(extension)
            else:  # negative value, extend with 1
                signBit = temp.pop(-1)
                temp.extend(extensionInv)
                temp.append(signBit)
            return temp
        # read data
        temp = {}
        if channelSet is None:
            channelSet = self.channelNames
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channelSet:
                temp[Channel.recAttributeName] = B[Channel.posBitBeg: Channel.posBitEnd]
                nbytes = len(temp[Channel.recAttributeName].tobytes())
                if not nbytes == Channel.nBytes:
                    byte = bitarray(8 * (Channel.nBytes - nbytes), endian='little')
                    byte.setall(False)
                    if Channel.signalDataType not in (1, 10, 14):  # not signed integer
                        temp[Channel.recAttributeName].extend(byte)
                    else:  # signed integer (two's complement), keep sign bit and extend with bytes
                        temp[Channel.recAttributeName] = signedInt(temp[Channel.recAttributeName], byte)
                nTrailBits = Channel.nBytes*8 - Channel.bitCount
                if Channel.signalDataType in (1, 10, 14) and \
                        nbytes == Channel.nBytes and \
                        nTrailBits > 0:  # Ctype byte length but signed integer
                    trailBits = bitarray(nTrailBits, endian='little')
                    temp[Channel.recAttributeName] = signedInt(temp[Channel.recAttributeName], trailBits)
                if 's' not in Channel.dataFormat:
                    temp[Channel.recAttributeName] = Channel.CFormat.unpack(temp[Channel.recAttributeName].tobytes())[0]
                else:
                    temp[Channel.recAttributeName] = temp[Channel.recAttributeName].tobytes()
        return temp  # returns dictionary of channel with its corresponding values

    def changeChannelName(self, channelName):
        """ In case of duplicate channel names within several channel groups
        for unsorted data, rename channel name

        Parameters
        ------------
        channelName : str
            channelName to be modified to avoid duplicates in unsorted
            channel groups
        """
        ind = self.dataRecordName.index(channelName) - 1
        OldrecAttributeName = self[ind].recAttributeName
        OldName = self[ind].name
        self[ind].changeChannelName(self.channelGroup)
        recAttributeName = self[ind].recAttributeName
        if self.recordToChannelMatching[OldrecAttributeName] == OldrecAttributeName:
            self.recordToChannelMatching.pop(OldrecAttributeName)
            self.recordToChannelMatching[recAttributeName] = recAttributeName
        if self.master['name'] == OldName:  # channel to change is master
            self.master['name'] = self[ind].name
        self.dataRecordName[ind] = recAttributeName
        self.channelNames.remove(OldName)
        self.channelNames.add(self[ind].name)


class DATA(dict):

    """ DATA class is organizing record classes itself made of channel.
    This class inherits from dict. Keys are corresponding to channel
    group recordID. A DATA class corresponds to a data block, a dict
    of record classes (one per channel group). Each record class contains
    a list of channel class representing the structure of channel record.

    Attributes
    --------------
    fid : io.open
        file identifier
    pointerToData : int
        position of Data block in mdf file
    BlockLength : int
        total size of data block

    Methods
    ------------
    addRecord(record)
        Adds a new record in DATA class dict
    read(channelSet)
        Reads data block
    loadSorted(record, nameList=None)
        Reads sorted data block from record definition
    loadUnSorted(nameList=None)
        Reads unsorted data block, not yet implemented
    """

    def __init__(self, fid, pointer):
        self.fid = fid
        self.pointerToData = pointer
        self.BlockLength = 0

    def addRecord(self, record):
        """Adds a new record in DATA class dict

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        """
        self[record.recordID] = {}
        self[record.recordID]['record'] = record
        self.BlockLength += record.dataBlockLength

    def read(self, channelSet):
        """Reads data block

        Parameters
        ----------------
        channelSet : set of str, optional
            list of channel names
        """
        if len(self) == 1:  # sorted dataGroup
            recordID = list(self.keys())[0]
            self[recordID]['data'] = \
                self.loadSorted(self[recordID]['record'],
                                nameList=channelSet)
        elif len(self) >= 2:  # unsorted DataGroup
            data = self.loadUnSorted(nameList=channelSet)
            for recordID in list(self.keys()):
                self[recordID]['data'] = {}
                for channel in self[recordID]['record']:
                    self[recordID]['data'][channel.recAttributeName] = \
                        data[self[recordID]['record'].
                             recordToChannelMatching[channel.recAttributeName]]
        else:  # empty data group
            pass

    def loadSorted(self, record, nameList=None):  # reads sorted data
        """Reads sorted data block from record definition

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        channelSet : set of str, optional
            list of channel names

        Returns
        -----------
        numpy recarray of data
        """
        return record.readSortedRecord(self.fid, self.pointerToData, nameList)

    def loadUnSorted(self, nameList=None):
        """Reads unsorted data block from record definition

        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        channelSet : set of str, optional
            list of channel names

        Returns
        -----------
        numpy recarray of data
        """
        self.fid.seek(self.pointerToData)
        stream = self.fid.read(self.BlockLength)
        # reads only the channels using offset functions, channel by channel.
        buf = {}
        position = 0
        recordIdCFormat = Struct('B')
        # several channel groups per data block, check for duplicate channel names
        temp = set()
        for recordID in self:
            intersec = set(self[recordID]['record'].dataRecordName) & temp
            for name in intersec: # duplicate channels found
                self[recordID]['record'].changeChannelName(name)
            temp.update(self[recordID]['record'].dataRecordName)
        # initialise data structure
        for recordID in self:
            for channelName in self[recordID]['record'].dataRecordName:
                buf[channelName] = []
        # read data
        while position < len(stream):
            recordID = recordIdCFormat.unpack(stream[position:position + 1])[0]
            if not self[recordID]['record'].hiddenBytes and self[recordID]['record'].byte_aligned:
                temp = self[recordID]['record'].readRecordBuf(stream[position:position + self[recordID]['record'].CGrecordLength + 1], nameList)
            else:  # do read bytes but bits in record
                temp = self[recordID]['record'].readRecordBits(stream[position:position + self[recordID]['record'].CGrecordLength + 1], nameList)
            # recordId is only unit8
            position += self[recordID]['record'].CGrecordLength + 1
            for channelName in temp:
                buf[channelName].append(temp[channelName])  # to remove append
        # convert list to array
        for chan in buf:
            buf[chan] = array(buf[chan])
        return buf


class mdf3(mdf_skeleton):

    """ mdf file version 3.0 to 3.3 class

    Attributes
    --------------
    fileName : str
        file name
    MDFVersionNumber : int
        mdf file version number
    masterChannelList : dict
        Represents data structure: a key per master channel with corresponding value containing
        a list of channels
        One key or master channel represents then a data group having same sampling interval.
    multiProc : bool
        Flag to request channel conversion multi processed for performance improvement.
        One thread per data group.
    convertAfterRead : bool
        flag to convert raw data to physical just after read
    filterChannelNames : bool
        flag to filter long channel names from its module names separated by '.'
    file_metadata : dict
        file metadata with minimum keys: author, organisation, project, subject, comment, time, date

    Methods
    ------------
    read3( fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True)
        Reads mdf 3.x file data and stores it in dict
    _getChannelData3(channelName)
        Returns channel numpy array
    _convertChannel3(channelName)
        converts specific channel from raw to physical data according to CCBlock information
    _convertAllChannel3()
        Converts all channels from raw data to converted data according to CCBlock information
    write3(fileName=None)
        Writes simple mdf 3.3 file
    """

    def read3(self, fileName=None, info=None, multiProc=False, channelList=None, \
        convertAfterRead=True, filterChannelNames=False, compression=False):
        """ Reads mdf 3.x file data and stores it in dict

        Parameters
        ----------------
        fileName : str, optional
            file name

        info : mdfinfo3.info3 class
            info3 class containing all MDF Blocks

        multiProc : bool
            flag to activate multiprocessing of channel data conversion

        channelList : list of str, optional
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory.
            Can be used to read big files

        convertAfterRead : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false, all data from channels
            will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .getChannelData()

        compression : bool, optional
            falg to activate data compression with blosc
        """
        self.multiProc = multiProc
        if platform == 'win32':
            self.multiProc = False  # no multiprocessing for windows platform

        if self.fileName is None and info is not None:
            self.fileName = info.fileName
            self._info = info
        elif fileName is not None:
            self.fileName = fileName

        if channelList is None:
            channelSetFile = None
        else:
            channelSetFile = set(channelList)

        # Read information block from file
        if self._info is None:
            if info is None:
                self._info = info3(self.fileName, None)
            else:
                self._info = info

        if self._info.fid.closed:
            self._info.fid = open(self.fileName, 'rb')

        # reads metadata
        try:
            comment = self._info['HDBlock']['TXBlock']['Text']
        except:
            comment = ''
        # converts date to be compatible with ISO8601
        day, month, year = self._info['HDBlock']['Date'].split(':')
        ddate = year + '-' + month + '-' + day
        self.add_metadata(author=self._info['HDBlock']['Author'], \
                organisation=self._info['HDBlock']['Organization'], \
                project=self._info['HDBlock']['ProjectName'], \
                subject=self._info['HDBlock']['Subject'], comment=comment, \
                    date=ddate, time=self._info['HDBlock']['Time'])

        try:
            fid = open(self.fileName, 'rb')
        except IOError:
            raise Exception('Can not find file ' + self.fileName)


        # Read data from file
        for dataGroup in self._info['DGBlock']:
            channelSet = channelSetFile
            if self._info['DGBlock'][dataGroup]['numberOfChannelGroups'] > 0:  # data exists
                # Pointer to data block
                pointerToData = self._info['DGBlock'][dataGroup]['pointerToDataRecords']
                buf = DATA(fid, pointerToData)

                for channelGroup in range(self._info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                    temp = record(dataGroup, channelGroup)  # create record class
                    temp.loadInfo(self._info)  # load all info related to record

                    if temp.numberOfRecords != 0:  # continue if there are at least some records
                        buf.addRecord(temp)
                        if self._noDataLoading and len(channelSet & \
                                buf[temp.recordID]['record'].channelNames) > 0:
                            channelSet = None  # will load complete datagroup

                buf.read(channelSet) # reads datablock potentially containing several channel groups

                for recordID in buf:
                    if 'record' in buf[recordID]:
                        master_channel = buf[recordID]['record'].master['name']
                        if master_channel in self and self[master_channel][dataField] is not None:
                            master_channel += '_' + str(dataGroup)

                        channels = (c for c in buf[recordID]['record']
                                    if channelSet is None or c.name in channelSet)

                        for chan in channels: # for each recordchannel
                            # in case record is used for several channels
                            recordName = buf[recordID]['record'].\
                                    recordToChannelMatching[chan.recAttributeName]
                            temp = buf[recordID]['data'][recordName]

                            if len(temp) != 0:
                                # Process concatenated bits inside uint8
                                if not chan.bitCount // 8.0 == chan.bitCount / 8.0 \
                                        and not buf[recordID]['record'].hiddenBytes \
                                        and buf[recordID]['record'].byte_aligned:
                                    # if channel data do not use complete bytes
                                    if chan.signalDataType in (0, 1, 9, 10, 13, 14):  # integers
                                        temp = right_shift(temp, chan.embedding_channel_bitOffset)
                                        mask = int(pow(2, chan.bitCount) - 1)  # masks isBitUint8
                                        temp = bitwise_and(temp, mask)
                                    else:  # should not happen
                                        print('bit count and offset not applied \
                                              to correct data type', file=stderr)
                                self.add_channel(dataGroup, chan.name, temp,
                                                 master_channel,
                                                 master_type=1,
                                                 unit=chan.unit,
                                                 description=chan.desc,
                                                 conversion=chan.conversion,
                                                 info=None,
                                                 compression=compression)
                del buf
        fid.close()  # close file
        if convertAfterRead and not compression:
            self._noDataLoading = False
            self._convertAllChannel3()

    def _getChannelData3(self, channelName):
        """Returns channel numpy array

        Parameters
        ----------------
        channelName : str
            channel name

        Returns:
        -----------
        numpy array
            converted, if not already done, data corresponding to channel name

        Notes
        ------
        This method is the safest to get channel data as numpy array from 'data' dict key might contain raw data
        """
        if channelName in self:
            vect = self.getChannel(channelName)[dataField]
            if vect is None: # noDataLoading reading argument flag activated
                if self._info.fid is None or (self._info.fid is not None and self._info.fid.closed):
                    (self._info.fid, self._info.fileName, self._info.zipfile) = _open_MDF(self.fileName)
                self.read3(fileName=None, info=self._info, channelList=[channelName], convertAfterRead=False)
            return self._convert3(channelName)
        else:
            raise KeyError('Channel ' + channelName + ' not in mdf dictionary')
            return channelName

    def _convert3(self, channelName):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channelName : str
            Name of channel

        Returns
        -----------
        numpy array
            returns numpy array converted to physical values according to conversion type
        """
        
        if self[channelName][dataField] is None:
            vect = self[channelName][dataField]
        else:
            if isinstance(self[channelName][dataField], compressed_data):
                vect = self[channelName][dataField].decompression()  # uncompress blosc
            else:
                vect = self[channelName][dataField][:]  # to have bcolz uncompressed data
        if conversionField in self[channelName]:  # there is conversion property
            conversion = self[channelName][conversionField]
            if conversion['type'] == 0:
                return linearConv(vect, conversion['parameters'])
            elif conversion['type'] == 1:
                return tabInterpConv(vect, conversion['parameters'])
            elif conversion['type'] == 2:
                return tabConv(vect, conversion['parameters'])
            elif conversion['type'] == 6:
                return polyConv(vect, conversion['parameters'])
            elif conversion['type'] == 7:
                return expConv(vect, conversion['parameters'])
            elif conversion['type'] == 8:
                return logConv(vect, conversion['parameters'])
            elif conversion['type'] == 9:
                return rationalConv(vect, conversion['parameters'])
            elif conversion['type'] == 10:
                return formulaConv(vect, conversion['parameters'])
            elif conversion['type'] == 12:
                return textRangeTableConv(vect, conversion['parameters'])
            else:
                return vect
        else:
            return vect

    def _convertChannel3(self, channelName):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channelName : str
            Name of channel
        """
        self.setChannelData(channelName, self._convert3(channelName))
        self.remove_channel_conversion(channelName)

    def _convertAllChannel3(self):
        """Converts all channels from raw data to converted data according to CCBlock information
        Converted data will take more memory.
        """
        if self._noDataLoading:  # no data loaded, load everything
            self.read3(self.fileName, convertAfterRead=True)
        else:
            for channel in self:
                self._convertChannel3(channel)

    def write3(self, fileName=None):
        """Writes simple mdf 3.3 file

        Parameters
        ----------------
        fileName : str, optional
            Name of file
            If file name is not input, written file name will be the one read with
            appended '_new' string before extension

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        # put master channel in first position for each datagroup if not already the case
        for master in self.masterChannelList:
            masterList = sorted(self.masterChannelList[master])
            masterPosition = masterList.index(master)
            masterList.pop(masterPosition)  # remove  master channel
            masterList.insert(0, master)  # insert at first position master channel
            self.masterChannelList[master] = masterList

        pointers = {}  # records pointers of blocks when writing

        # writes characters
        def writeChar(f, value, size=None):
            if size is None:
                temp = value
            else:
                if len(value) > size:
                    temp = value[:size]
                else:
                    temp = value + '\0' * (size - len(value))
                temp += '\0'
            if self.MDFVersionNumber < 400:
                if PythonVersion >= 3:
                    temp = temp.encode('latin1', 'replace')
                f.write(pack('<' + str(len(temp)) + 's', temp))
            else:
                temp = temp.encode('latin1', 'replace')
                f.write(pack('<' + str(len(temp)) + 's', temp))

        # write pointer of block and come back to current stream position
        def writePointer(f, pointer, value):
            currentPosition = f.tell()
            f.seek(pointer)
            f.write(pack('I', value))
            f.seek(currentPosition)

        fid = open(fileName, 'wb')  # buffering should automatically be set
        # Starts first to write ID and header
        head = (b'MDF     ', b'3.30    ', b'MDFreadr', 0, 0, 330,
                28591, b'\0' * 32)
        # code page ISO2859-1 latin 1 western europe
        fid.write(pack('<8s8s8s4H32s', *head))

        # Header Block
        pointers['HD'] = {}
        pointers['HD']['DG'] = 68
        pointers['HD']['TX'] = 72
        pointers['HD']['PR'] = 76
        ndataGroup = len(self.masterChannelList)
        if self.file_metadata['author'] is not None:  # Author
            author = '{:\x00<32}'.format(self.file_metadata['author']).encode('latin-1')
        elif PythonVersion >= 3:
            author = '{:\x00<32}'.format(os.getlogin()).encode('latin-1')
        else:
            author = b'\x00' * 32
        if self.file_metadata['organisation'] is not None:  # Organization
            organization = '{:\x00<32}'.format(self.file_metadata['organisation']).encode('latin-1')
        else:
            organization = b'\x00' * 32
        if self.file_metadata['project'] is not None:  # Project
            project = '{:\x00<32}'.format(self.file_metadata['project']).encode('latin-1')
        else:
            project = b'\x00' * 32
        if self.file_metadata['subject'] is not None:  # Subject
            subject = '{:\x00<32}'.format(self.file_metadata['subject']).encode('latin-1')
        else:
            subject = b'\x00' * 32
        # Header Block, block size
        # first Data block pointer, pointer to TX Block file comment, pointer to PR Block
        # number of data groups, date, time
        # Time Stamp, UTC time offset, Time quality, Timer identification
        timeOffset = gmtime()
        head = (b'HD', 208, 272, 0, 0, ndataGroup,
                '{:\x00<10}'.format(strftime("%d:%m:%Y", timeOffset)).encode('latin-1'),
                '{:\x00<8}'.format(strftime("%H:%M:%S", timeOffset)).encode('latin-1'),
                author, organization, project, subject,
                int(time() * 1000000000), 1, 0,
                b'Local PC Reference Time         ')
        fid.write(pack('<2sH3IH10s8s32s32s32s32sQhH32s', *head))

        # write DG block
        pointers['DG'] = {}
        pointers['CG'] = {}
        pointers['CN'] = {}
        dataGroup = 0

        for masterChannel in self.masterChannelList:
            # writes dataGroup Block
            pointers['DG'][dataGroup] = {}
            position = fid.tell()
            if 0 < dataGroup:  # not possible for first DG
                # previous datagroup pointer to this new datagroup
                writePointer(fid, pointers['DG'][dataGroup - 1]['nextDG'], position)
            else:
                # first datagroup pointer in header block
                writePointer(fid, pointers['HD']['DG'], position)
            pointers['DG'][dataGroup]['nextDG'] = position + 4
            pointers['DG'][dataGroup]['CG'] = position + 8
            pointers['DG'][dataGroup]['data'] = position + 16
            # DG block size, pointer to next DataGroup,
            # pointer to channel group, pointer to trigger block, pointer to data block
            # number of channel group, number of record IDs, reserved
            head = (b'DG', 28, 0, 0, 0, 0, 1, 0, b'\x00' * 32)
            fid.write(pack('<2sH4I2H4s', *head))

            # sorted data so only one channel group
            pointers['CG'][dataGroup] = {}
            # write first CG pointer in datagroup
            position = fid.tell()
            writePointer(fid, pointers['DG'][dataGroup]['CG'], position)
            pointers['CG'][dataGroup]['firstCN'] = position + 8
            pointers['CG'][dataGroup]['TX'] = position + 12
            pointers['CG'][dataGroup]['dataRecordSize'] = position + 20
            numChannels = len(self.masterChannelList[masterChannel])
            masterData = self.getChannelData(masterChannel)
            nRecords = len(masterData)
            # CG block size, pointer to next Channel Group
            # pointer to first channel block, pointer to TX block
            # No record ID no need for sorted data, Number of channels
            # Size of data record, Number of records, pointer to sample reduction block
            head = (b'CG', 30, 0, 0, 0, 0, numChannels, 0, nRecords, 0)
            fid.write(pack('<2sH3I3H2I', *head))

            # Channel blocks writing
            pointers['CN'][dataGroup] = {}
            dataList = ()
            dataTypeList = ''
            recordNumberOfBits = 0
            preceedingChannel = None
            bitOffset = 0
            byteOffset = 0
            # first channel bock pointer from CG
            position = fid.tell()
            writePointer(fid, pointers['CG'][dataGroup]['firstCN'], position)
            for channel in self.masterChannelList[masterChannel]:
                position = fid.tell()
                pointers['CN'][dataGroup][channel] = {}
                pointers['CN'][dataGroup][channel]['beginCN'] = position
                pointers['CN'][dataGroup][channel]['nextCN'] = position + 4
                pointers['CN'][dataGroup][channel]['CC'] = position + 8
                pointers['CN'][dataGroup][channel]['TX'] = position + 16
                if preceedingChannel is not None:  # not possible for first CN
                    # pointer in previous CN
                    writePointer(fid, pointers['CN'][dataGroup][preceedingChannel]['nextCN'],
                                 pointers['CN'][dataGroup][channel]['beginCN'])
                preceedingChannel = channel
                if channel not in set(self.masterChannelList):
                    masterFlag = 0  # data channel
                else:
                    masterFlag = 1  # master channel
                desc = self.getChannelDesc(channel)
                #if bitOffset exceeds two byte limit, we start using the byte offset field
                if bitOffset > 0xFFFF:
                    bitOffset -= 0x10000
                    byteOffset += 8192
                data = self.getChannelData(channel)  # channel data
                temp = data
                if PythonVersion >= 3 and data.dtype.kind in ['S', 'U']:
                    temp = ncode(temp, encoding='latin1', errors='replace')
                dataList = dataList + (temp, )
                if data.dtype in ('float64', 'int64', 'uint64'):
                    numberOfBits = 64
                elif data.dtype in ('float32', 'int32', 'uint32'):
                    numberOfBits = 32
                elif data.dtype in ('uint16', 'int16'):
                    numberOfBits = 16
                elif data.dtype in ('uint8', 'int8', 'bool'):
                    numberOfBits = 8
                else:
                    numberOfBits = 8  # if string, considered later
                if data.dtype == 'float64':
                    dataType = 3
                elif data.dtype in ('uint8', 'uint16', 'uint32', 'uint64', 'bool'):
                    dataType = 0
                elif data.dtype in ('int8', 'int16', 'int32', 'int64'):
                    dataType = 1
                elif data.dtype == 'float32':
                    dataType = 2
                elif data.dtype.kind in ['S', 'U']:
                    dataType = 7
                else:
                    raise Exception('Not recognized dtype')
                    return data.dtype
                if data.dtype.kind not in ['S', 'U']:
                    dataTypeList += data.dtype.char
                else:
                    dataTypeList += str(data.dtype.itemsize) + 's'
                    numberOfBits = 8 * data.dtype.itemsize
                recordNumberOfBits += numberOfBits
                if data.dtype.kind not in ['S', 'U']:
                    valueRangeValid = 1
                    if len(data) > 0:
                        maximum = max(data)
                        minimum = min(data)
                    else:
                        maximum = 0
                        minimum = 0
                else:
                    valueRangeValid = 0  # No value range valid
                    minimum = 0  # Min value
                    maximum = 0  # Max value
                pointers['CN'][dataGroup][channel]['longChannelName'] = position + 218
                # CN,  block size
                # pointer to next channel block, pointer to conversion block
                # pointer to source depending block, pointer to dependency block
                # pointer to comment TX, check if master channel
                # channel name, channel description,
                # bit position, Number of bits, Signal data type
                # value range valid, minimum, maximum, Sampling rate
                # pointer to long channel name
                # pointer to channel display name
                # additional byte offset
                try:
                    description = '{:\x00<128}'.format(desc).encode('latin-1')
                except:
                    description = b'\x00' * 128
                head = (b'CN', 228, 0, 0, 0, 0, 0, masterFlag,
                        '{:\x00<32}'.format(channel).encode('latin-1'),
                        description,
                        bitOffset, numberOfBits, dataType, valueRangeValid,
                        minimum, maximum, 0, 0, 0, byteOffset)
                fid.write(pack('<2sH5IH32s128s4H3d2IH', *head))
                bitOffset += numberOfBits

                # TXblock for long channel name
                writePointer(fid, pointers['CN'][dataGroup][channel]['longChannelName'], fid.tell())
                head = (b'TX', len(channel) + 4 + 1, channel.encode('latin-1') + b'\x00')
                fid.write(pack('<2sH{}s'.format(len(channel)+1), *head))

                # Conversion blocks writing
                writePointer(fid, pointers['CN'][dataGroup][channel]['CC'], fid.tell())
                # channel description
                # conversion already done during reading
                # additional size information, not necessary for 65535 conversion type ?
                try:
                    unit = '{:\x00<20}'.format(self.getChannelUnit(channel)\
                            .encode('latin-1', 'replace'))
                except:
                    unit = b'\x00' * 20
                head = (b'CC', 46, valueRangeValid, minimum, maximum,
                        unit, 65535, 0)
                fid.write(pack('<2shH2d20s2H', *head))

            # number of channels in CG
            currentPosition = fid.tell()
            fid.seek(pointers['CG'][dataGroup]['dataRecordSize'])
            fid.write(pack('H', int(recordNumberOfBits / 8)))  # Size of data record
            fid.seek(currentPosition)

            # data writing
            # write data pointer in datagroup
            writePointer(fid, pointers['DG'][dataGroup]['data'], fid.tell())
            # dumps data vector from numpy
            fid.write(fromarrays(dataList).tobytes(order='F'))

            dataGroup += 1

        # print(pointers, file=stderr)
        fid.close()


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
            print(('Unsupported number of bits for unsigned int ' + str(signalDataType)), file=stderr)

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
            print(('Unsupported number of bits for signed int ' + str(signalDataType)), file=stderr)

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBits == 32:
            dataType = 'f'
        elif numberOfBits == 64:
            dataType = 'd'
        else:
            print(('Unsupported number of bit for floating point ' + str(signalDataType)), file=stderr)

    elif signalDataType == 7:  # string
        dataType = str(numberOfBits // 8) + 's'
    elif signalDataType == 8:  # array of bytes
        dataType = str(numberOfBits // 8) + 's'
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits), file=stderr)

    # deal with byte order
    if signalDataType in (0, 1, 2, 3):
        if ByteOrder:
            dataType = '>' + dataType
        else:
            dataType = '<' + dataType
    elif signalDataType in (13, 14, 15, 16):  # low endian
        dataType = '<' + dataType
    elif signalDataType in (9, 10, 11, 12):  # big endian
        dataType = '>' + dataType

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
            print('Unsupported number of bits for unsigned int ' + str(signalDataType) + ' nBits ', numberOfBits, file=stderr)

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
            print('Unsupported number of bits for signed int ' + str(signalDataType) + ' nBits ', numberOfBits, file=stderr)

    elif signalDataType in (2, 3, 11, 12, 15, 16):  # floating point
        if numberOfBits == 32:
            dataType = 'f4'
        elif numberOfBits == 64:
            dataType = 'f8'
        else:
            print('Unsupported number of bit for floating point ' + str(signalDataType) + ' nBits ', numberOfBits, file=stderr)

    elif signalDataType == 7:  # string
        dataType = 'S' + str(numberOfBits // 8)  # not directly processed
    elif signalDataType == 8:  # array of bytes
        dataType = 'V' + str(numberOfBits // 8)  # not directly processed
    else:
        print('Unsupported Signal Data Type ' + str(signalDataType) + ' nBits ', numberOfBits, file=stderr)

    nativeDataType = dataType
    # deal with byte order
    if signalDataType in (0, 1, 2, 3):
        if ByteOrder:
            dataType = '>' + dataType
        else:
            dataType = '<' + dataType
    elif signalDataType in (13, 14, 15, 16):  # low endian
        dataType = '<' + dataType
    elif signalDataType in (9, 10, 11, 12):  # big endian
        dataType = '>' + dataType

    return (dataType, nativeDataType)
