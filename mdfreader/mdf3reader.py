# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 3.x

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Created on Sun Oct 10 12:57:28 2010


Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both
    python 2.6+ and 3.2+


mdf3reader
--------------------------
"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from __future__ import print_function
from numpy import right_shift, bitwise_and, interp
from numpy import max as npmax, min as npmin
from numpy import asarray, recarray, array, searchsorted, vectorize
from numpy import issubdtype, number as numpy_number
from numpy.core.records import fromstring, fromarrays
from collections import defaultdict
from math import log, exp
from time import strftime, time, gmtime
from struct import pack, Struct
from io import open  # for python 3 and 2 consistency
from sys import platform, exc_info, version_info
from warnings import warn
import os
from warnings import simplefilter
from .mdf import MdfSkeleton, _open_mdf, \
    dataField, conversionField, idField, CompressedData
from .mdfinfo3 import Info3
from .channel import Channel3
if os.name == 'posix':
    from os import getlogin


PythonVersion = version_info
PythonVersion = PythonVersion[0]

chunk_size_reading = 100000000  # reads by chunk of 100Mb, can be tuned for best performance


def _linear_conversion(data, conversion):  # 0 Parametric, Linear: Physical =Integer*P2 + P1
    """ apply linear conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conversion['P2'] == 1.0 and conversion['P1'] in (0.0, -0.0):
        return data  # keeps dtype probably more compact than float64
    else:
        return data * conversion['P2'] + conversion['P1']


def _tab_interp_conversion(data, conversion):  # 1 Tabular with interpolation
    """ apply Tabular interpolation conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    tmp = array([(key, val['int'], val['phys'])
                 for (key, val) in conversion.items()])
    return interp(data, tmp[:, 1], tmp[:, 2])


def _tab_conversion(data, conversion):  # 2 Tabular
    """ apply Tabular conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    tmp = array([(key, val['int'], val['phys'])
                 for (key, val) in conversion.items()])
    indexes = searchsorted(tmp[:, 1], data)
    return tmp[indexes, 2]


def _polynomial_conversion(data, conversion):  # 6 Polynomial
    """ apply polynomial conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    return (conversion['P2'] - conversion['P4'] * (data - conversion['P5'] - conversion['P6'])) \
        / (conversion['P3'] * (data - conversion['P5'] - conversion['P6']) - conversion['P1'])


def _exponential_conversion(data, conversion):  # 7 Exponential
    """ apply exponential conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conversion['P4'] == 0 and conversion['P1'] != 0 and conversion['P2'] != 0:
        return exp(((data - conversion['P7']) * conversion['P6'] - conversion['P3'])
                   / conversion['P1']) / conversion['P2']
    elif conversion['P1'] == 0 and conversion['P4'] != 0 and conversion['P5'] != 0:
        return exp((conversion['P3'] / (data - conversion['P7']) - conversion['P6'])
                   / conversion['P4']) / conversion['P5']
    else:
        warn('Non possible exponential parameters for channel')


def _log_conversion(data, conversion):  # 8 Logarithmic
    """ apply logarithmic conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    if conversion['P4'] == 0 and conversion['P1'] != 0 and conversion['P2'] != 0:
        return log(((data - conversion['P7']) * conversion['P6'] - conversion['P3'])
                   / conversion['P1']) / conversion['P2']
    elif conversion['P1'] == 0 and conversion['P4'] != 0 and conversion['P5'] != 0:
        return log((conversion['P3'] / (data - conversion['P7']) - conversion['P6'])
                   / conversion['P4']) / conversion['P5']
    else:
        warn('Non possible logarithmic conversion parameters for channel')


def _rational_conversion(data, conversion):  # 9 rational
    """ apply rational conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    return (conversion['P1'] * data * data + conversion['P2'] * data + conversion['P3'])\
        / (conversion['P4'] * data * data + conversion['P5'] * data + conversion['P6'])


def _formula_conversion(data, conversion):  # 10 Text Formula
    """ apply formula conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

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
        formula = conversion['textFormula']
        # remove trailing text after 0
        formula = formula[:formula.find('\x00')]
        # adapt ASAM-MCD2 syntax to sympy
        formula = formula.replace('pow(', 'power(')
        # formula to function for evaluation
        expr = lambdify(X, formula, modules='numpy', dummify=False)
        return expr(data)
    except:
        warn('Failed to convert formulae ' + conversion['textFormula'] +
             ' Sympy is correctly installed ?')


def _text_table_conversion(data, conversion):  # 11 Text table
    """ apply text table conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    conversion_table = dict()
    for pair in conversion:
        conversion_table[conversion[pair]['int']] = conversion[pair]['text']
    return vectorize(conversion_table.__getitem__)(data)


def _text_range_table_conversion(data, conversion):  # 12 Text range table
    """ apply text range table conversion to data

    Parameters
    ----------------
    data : numpy 1D array
        raw data to be converted to physical value
    conversion : mdfinfo3.info3 conversion block ('CCBlock') dict

    Returns
    -----------
    converted data to physical value
    """
    try:
        n_pair = len(conversion)
        lower = [conversion[pair]['lowerRange'] for pair in range(n_pair)]
        upper = [conversion[pair]['upperRange'] for pair in range(n_pair)]
        text = {}
        for pair in range(n_pair):
            text[pair] = conversion[pair]['Textrange']
            if text[pair] is None:
                continue
            if 'LINEAR_CONV' in text[pair]:  # linear conversion from CANape
                from sympy import lambdify, symbols
                X = symbols('X')  # variable is X
                left = text[pair].find('"')
                right = text[pair].rfind('"')
                text[pair] = text[pair][left + 1: right].replace('{', '').replace('}', '')
                text[pair] = lambdify(X, text[pair], modules='numpy', dummify=False)
        temp = []
        for l_index in range(len(data)):
            value = text[0]  # default value
            for pair in range(1, n_pair):
                if lower[pair] <= data[l_index] <= upper[pair]:
                    value = text[pair]
                    break
            if callable(value):
                value = value(data[l_index])
            temp.append(value)
        try:
            temp = asarray(temp)  # try to convert to numpy
        except:
            pass
        return temp
    except:
        warn('Failed to convert text to range table')


class Record(list):

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
    """

    def __init__(self, data_group, channel_group):
        self.CGrecordLength = 0
        self.recordLength = 0
        self.dataBlockLength = 0
        self.numberOfRecords = 0
        self.recordID = 0
        self.recordIDnumber = 0
        self.dataGroup = data_group
        self.channelGroup = channel_group
        self.numpyDataRecordFormat = []
        self.dataRecordName = []
        self.master = {}
        self.master['name'] = 'master_{}'.format(data_group)
        self.master['number'] = None
        self.recordToChannelMatching = {}
        self.channelNames = set()
        self.hiddenBytes = False
        self.byte_aligned = True

    def __repr__(self):
        output = list()
        output.append('Channels :\n')
        for chan in self.channelNames:
            output.append(''.join([chan, '\n']))
        output.append('Data group number : {}\n'.format(self.dataGroup))
        if self.master['name'] is not None:
            output.append(''.join(['Master channel : ', self.master['name'], '\n']))
        output.append('Numpy records format : \n')
        for record in self.numpyDataRecordFormat:
            output.append('{}\n'.format(record))
        return ''.join(output)

    def add_channel(self, info, channel_number):
        """ add a channel in class

        Parameters
        ----------------
        info : mdfinfo3.info3 class
        channel_number : int
            channel number in mdfinfo3.info3 class

        """
        self.append(Channel3(info, self.dataGroup, self.channelGroup, channel_number, self.recordIDnumber))
        self.channelNames.add(self[-1].name)

    def load_info(self, info):
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
            self.dataRecordName.append('RecordID{}'.format(self.channelGroup))
            self.numpyDataRecordFormat.append('uint8')
            self.dataBlockLength = (self.CGrecordLength + 1) * self.numberOfRecords
        embedding_channel = None
        for channelNumber in range(info['CGBlock'][self.dataGroup][self.channelGroup]['numberOfChannels']):
            channel = Channel3(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDnumber)
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
                    and channel.posBitEnd <= 8 * (prev_chan.byteOffset + prev_chan.nBytes_aligned)
                if embedding_channel is not None:
                    embedding_channel_includes_curr_chan = \
                        channel.posBitEnd <= embedding_channel.posByteEnd * 8
                else:
                    embedding_channel_includes_curr_chan = False
                if channel.byteOffset >= prev_chan.byteOffset and \
                        channel.posBitBeg < 8 * (prev_chan.byteOffset + prev_chan.nBytes_aligned) and \
                        channel.posBitEnd > 8 * (prev_chan.byteOffset + prev_chan.nBytes_aligned):
                    # not byte aligned
                    self.byte_aligned = False
                if embedding_channel is not None and \
                        channel.posBitEnd > embedding_channel.posByteEnd * 8:
                    embedding_channel = None
                if prev_chan_includes_curr_chan or \
                        embedding_channel_includes_curr_chan:  # bit(s) in byte(s)
                    if embedding_channel is None and prev_chan_includes_curr_chan:
                        embedding_channel = prev_chan  # new embedding channel detected
                    if self.recordToChannelMatching:  # not first channel
                        self.recordToChannelMatching[channel.name] = \
                            self.recordToChannelMatching[prev_chan.name]
                        channel.embedding_channel_bitOffset = \
                            channel.posBitBeg - embedding_channel.posBitBeg
                    else:  # first channels
                        self.recordToChannelMatching[channel.name] = \
                            channel.name
                        self.numpyDataRecordFormat.append(channel.dataFormat)
                        self.dataRecordName.append(channel.name)
                        self.recordLength += channel.nBytes_aligned
            if embedding_channel is None:  # adding bytes
                self.recordToChannelMatching[channel.name] = \
                    channel.name
                self.numpyDataRecordFormat.append(channel.dataFormat)
                self.dataRecordName.append(channel.name)
                self.recordLength += channel.nBytes_aligned
        if self.recordIDnumber == 2:  # second record ID at end of record
            self.dataRecordName.append('RecordID{}_2'.format(self.channelGroup))
            self.numpyDataRecordFormat.append('uint8')
            self.dataBlockLength = (self.CGrecordLength + 2) * self.numberOfRecords
        # check for hidden bytes
        if self.CGrecordLength > self.recordLength:
            self.hiddenBytes = True
        # check record length consistency
        elif self.CGrecordLength < self.recordLength:
            # forces to use dataRead instead of numpy records.
            self.byte_aligned = False

    def read_sorted_record(self, fid, pointer, channel_set=None):
        """ reads record, only one channel group per data group

        Parameters
        ----------------
        fid : float
            file identifier
        pointer
            position in file of data block beginning
        channel_set : Set of str, optional
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
        n_chunks = self.dataBlockLength // chunk_size_reading + 1
        chunk_length = self.dataBlockLength // n_chunks
        n_record_chunk = chunk_length // self.CGrecordLength
        chunks = [(n_record_chunk, self.CGrecordLength * n_record_chunk)] * n_chunks
        n_record_chunk = self.numberOfRecords - n_record_chunk * n_chunks
        if n_record_chunk > 0:
            chunks.append((n_record_chunk, self.CGrecordLength * n_record_chunk))
        previous_index = 0
        if channel_set is None and not self.hiddenBytes and self.byte_aligned:
            # reads all, quickest but memory consuming
            buf = recarray(self.numberOfRecords, dtype={'names': self.dataRecordName,
                                                        'formats': self.numpyDataRecordFormat})  # initialise array
            simplefilter('ignore', FutureWarning)
            for n_record_chunk, chunk_size in chunks:
                buf[previous_index: previous_index + n_record_chunk] = \
                    fromstring(fid.read(chunk_size),
                               dtype={'names': self.dataRecordName,
                                      'formats': self.numpyDataRecordFormat},
                               shape=n_record_chunk)
                previous_index += n_record_chunk
            return buf
        else:  # reads only some channels from a sorted data block
            if channel_set is None:
                channel_set = self.channelNames
            # memory efficient but takes time
            # are channelSet in this dataGroup
            if len(channel_set & self.channelNames) > 0:
                # check if master channel is in the list
                if not self.master['name'] in channel_set:
                    channel_set.add(self.master['name'])  # adds master channel
                rec_chan = []
                numpy_data_record_format = []
                data_record_name = []
                for channel in self:  # list of Channels from channelSet
                    if channel.name in channel_set:
                        rec_chan.append(channel)
                        data_record_name.append(channel.name)
                        numpy_data_record_format.append(channel.nativedataFormat)
                rec = recarray(self.numberOfRecords, dtype={'names': data_record_name,
                                                            'formats': numpy_data_record_format})
                try:  # use rather cython compiled code for performance
                    from dataRead import data_read
                    # converts data type from mdf 3.x to 4.x
                    convertDataType3to4 = {0: 0, 1: 2, 2: 4, 3: 4,
                                           7: 6, 8: 10,
                                           9: 1, 10: 3, 11: 5, 12: 5,
                                           13: 0, 14: 2, 15: 4, 16: 4}
                    for n_record_chunk, chunk_size in chunks:
                        bit_stream = fid.read(chunk_size)
                        for id, chan in enumerate(rec_chan):
                            rec[chan.name][previous_index: previous_index + n_record_chunk] = \
                                data_read(bytes(bit_stream),
                                          chan.bitCount,
                                          convertDataType3to4[chan.signalDataType],
                                          chan.nativedataFormat,
                                          n_record_chunk,
                                          self.CGrecordLength,
                                          chan.bitOffset,
                                          chan.posByteBeg,
                                          chan.nBytes_not_aligned, 0)
                            # masking already considered in dataRead
                            self[rec_chan[id].channelNumber].bit_masking_needed = False
                        previous_index += n_record_chunk
                    return rec
                except:
                    warn('Unexpected error: {}'.format(exc_info()))
                    warn('dataRead crashed, back to python data reading')
                    record_length = self.recordIDnumber + self.CGrecordLength
                    for r in range(self.numberOfRecords):  # for each record,
                        buf = fid.read(record_length)
                        for channel in rec_chan:
                            rec[channel.name][r] = \
                                channel.CFormat.unpack(buf[channel.posByteBeg:
                                                           channel.posByteEnd])[0]
                    return rec.view(recarray)

    def read_record_buf(self, buf, channel_set=None):
        """ read stream of record bytes

        Parameters
        ----------------
        buf : stream
            stream of bytes read in file
        channel_set : Set of str, optional
            list of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        temp = {}
        if channel_set is None:
            channel_set = self.channelNames
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channel_set:
                temp[self.recordToChannelMatching[Channel.name]] = \
                    Channel.CFormat.unpack(buf[Channel.posByteBeg:
                                               Channel.posByteEnd])[0]
        return temp  # returns dictionary of channel with its corresponding values

    def read_record_bits(self, bit_stream, channel_set=None):
        """ read stream of record bits by bits in case of not aligned or hidden bytes

        Parameters
        ----------------
        bit_stream : stream
            stream of bytes read in file
        channel_set : Set of str, optional
            list of channel to read

        Returns
        -----------
        rec : dict
            returns dictionary of channel with its corresponding values

        """
        from bitarray import bitarray
        bit_array = bitarray(endian="little")  # little endian by default
        bit_array.frombytes(bytes(bit_stream))

        def signed_int(temp, extension):
            """ extend bits of signed data managing two's complement
            """
            extension.setall(False)
            extension_inv = bitarray(extension, endian='little')
            extension_inv.setall(True)
            sign_bit = temp[-1]
            if not sign_bit:  # positive value, extend with 0
                temp.extend(extension)
            else:  # negative value, extend with 1
                sign_bit = temp.pop(-1)
                temp.extend(extension_inv)
                temp.append(sign_bit)
            return temp
        # read data
        temp = {}
        if channel_set is None:
            channel_set = self.channelNames
        for Channel in self:  # list of channel classes from channelSet
            if Channel.name in channel_set:
                temp[Channel.name] = bit_array[Channel.posBitBeg: Channel.posBitEnd]
                n_bytes = len(temp[Channel.name].tobytes())
                if not n_bytes == Channel.nBytes_aligned:
                    delta_byte = bitarray(8 * (Channel.nBytes_aligned - n_bytes), endian='little')
                    delta_byte.setall(False)
                    if Channel.signalDataType not in (1, 10, 14):  # not signed integer
                        temp[Channel.name].extend(delta_byte)
                    else:  # signed integer (two's complement), keep sign bit and extend with bytes
                        temp[Channel.name] = signed_int(temp[Channel.name], delta_byte)
                n_trail_bits = Channel.nBytes_aligned*8 - Channel.bitCount
                if Channel.signalDataType in (1, 10, 14) and \
                        n_bytes == Channel.nBytes_aligned and \
                        n_trail_bits > 0:  # Ctype byte length but signed integer
                    trail_bits = bitarray(n_trail_bits, endian='little')
                    temp[Channel.name] = signed_int(temp[Channel.name], trail_bits)
                if 's' not in Channel.dataFormat:
                    temp[Channel.name] = Channel.CFormat.unpack(temp[Channel.name].tobytes())[0]
                else:
                    temp[Channel.name] = temp[Channel.name].tobytes()
        return temp  # returns dictionary of channel with its corresponding values


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
    add_record(record)
        Adds a new record in DATA class dict
    read(channelSet)
        Reads data block
    load_sorted(record, nameList=None)
        Reads sorted data block from record definition
    load_unsorted(nameList=None)
        Reads unsorted data block, not yet implemented
    """

    def __init__(self, fid, pointer):
        self.fid = fid
        self.pointerToData = pointer
        self.BlockLength = 0

    def add_record(self, record):
        """Adds a new record in DATA class dict

        Parameters
        ----------------
        record: class
            channel group definition listing record channel classes
        """
        self[record.recordID] = {}
        self[record.recordID]['record'] = record
        self.BlockLength += record.dataBlockLength

    def read(self, channel_set, file_name):
        """Reads data block

        Parameters
        ----------------
        channel_set : set of str, optional
            list of channel names
        file_name : str
            name of file
        """
        # checks if file is closed
        if self.fid is None or self.fid.closed:
            self.fid = open(file_name, 'rb')
        if len(self) == 1:  # sorted dataGroup
            record_id = list(self.keys())[0]
            self[record_id]['data'] = \
                self.load_sorted(self[record_id]['record'], name_list=channel_set)
        elif len(self) >= 2:  # unsorted DataGroup
            data = self.load_unsorted(name_list=channel_set)
            for record_id in list(self.keys()):
                self[record_id]['data'] = {}
                for channel in self[record_id]['record']:
                    self[record_id]['data'][channel.name] = \
                        data[self[record_id]['record'].
                             recordToChannelMatching[channel.name]]

    def load_sorted(self, record, name_list=None):  # reads sorted data
        """Reads sorted data block from record definition

        Parameters
        ----------------
        record: class
            channel group definition listing record channel classes
        name_list : set of str, optional
            list of channel names

        Returns
        -----------
        numpy recarray of data
        """
        return record.read_sorted_record(self.fid, self.pointerToData, name_list)

    def load_unsorted(self, name_list=None):
        """Reads unsorted data block from record definition

        Parameters
        ----------------
        name_list : set of str, optional
            list of channel names

        Returns
        -----------
        numpy recarray of data
        """
        self.fid.seek(self.pointerToData)
        stream = self.fid.read(self.BlockLength)
        # reads only the channels using offset functions, channel by channel.
        buf = defaultdict(list)
        position = 0
        record_id_c_format = Struct('B')
        # initialise data structure
        for record_id in self:
            for channelName in self[record_id]['record'].dataRecordName:
                buf[channelName] = []
            if self[record_id]['record'].hiddenBytes or not self[record_id]['record'].byte_aligned:
                for ind, chan in enumerate(self[record_id]['record']):
                    # will already extract bits, no need of masking later
                    self[record_id]['record'][ind].bit_masking_needed = False
        # read data
        while position < len(stream):
            record_id = record_id_c_format.unpack(stream[position:position + 1])[0]
            if not self[record_id]['record'].hiddenBytes and self[record_id]['record'].byte_aligned:
                temp = self[record_id]['record'].read_record_buf(
                    stream[position:position + self[record_id]['record'].CGrecordLength + 1], name_list)
            else:  # do not read bytes but bits in record
                temp = self[record_id]['record'].read_record_bits(
                    stream[position:position + self[record_id]['record'].CGrecordLength + 1], name_list)
            # recordId is only unit8
            position += self[record_id]['record'].CGrecordLength + 1
            for channelName in temp:
                buf[channelName].append(temp[channelName])  # to remove append
        # convert list to array
        for chan in buf:
            buf[chan] = array(buf[chan])
        return buf


class Mdf3(MdfSkeleton):

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
    fileMetadata : dict
        file metadata with minimum keys: author, organisation, project, subject, comment, time, date

    Methods
    ------------
    read3( fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True)
        Reads mdf 3.x file data and stores it in dict
    _get_channel_data3(channelName)
        Returns channel numpy array
    _convert_channel3(channelName)
        converts specific channel from raw to physical data according to CCBlock information
    _convert_all_channel3()
        Converts all channels from raw data to converted data according to CCBlock information
    write3(fileName=None)
        Writes simple mdf 3.3 file
    """

    def read3(self, file_name=None, info=None, multi_processed=False, channel_list=None, convert_after_read=True,
              filter_channel_names=False, compression=False, metadata=2):
        """ Reads mdf 3.x file data and stores it in dict

        Parameters
        ----------------
        file_name : str, optional
            file name

        info : mdfinfo3.info3 class
            info3 class containing all MDF Blocks

        multi_processed : bool
            flag to activate multiprocessing of channel data conversion

        channel_list : list of str, optional
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory.
            Can be used to read big files

        convert_after_read : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false, all data from channels
            will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .get_channel_data()

        filter_channel_names : bool, optional
            flag to filter long channel names from its module names separated by '.'

        compression : bool, optional
            flag to activate data compression with blosc

        metadata: int, optional, default = 2
            Reading metadata has impact on performance, especially for mdf 4.x using xml.
            2: minimal metadata reading (mostly channel blocks)
            1: used for noDataLoading
            0: all metadata reading

        """
        self.multiProc = multi_processed
        if platform == 'win32':
            self.multiProc = False  # no multiprocessing for windows platform

        if self.fileName is None and info is not None:
            self.fileName = info.fileName
        elif file_name is not None and self.fileName is None:
            self.fileName = file_name

        minimal = metadata  # always reads minimum info by default

        if channel_list is None:
            channel_set_file = None
        else:
            channel_set_file = set(channel_list)
            minimal = 1  # reads at least CN to populate ChannelNamesByDG

        # Read information block from file
        if info is None:
            if self.info is None:
                info = Info3(self.fileName, fid=None, filter_channel_names=False, minimal=minimal)
            else:
                info = self.info

        if info.fid is None or info.fid.closed:
            try:
                info.fid = open(self.fileName, 'rb')
            except IOError:
                raise Exception('Can not find file ' + self.fileName)

        # reads metadata
        if not self._noDataLoading:
            try:
                comment = info['HDBlock']['TXBlock']
            except:
                comment = ''
            # converts date to be compatible with ISO8601
            day, month, year = info['HDBlock']['Date'].split(':')
            ddate = '-'.join([year, month, day])
            self.add_metadata(author=info['HDBlock']['Author'],
                              organisation=info['HDBlock']['Organization'],
                              project=info['HDBlock']['ProjectName'],
                              subject=info['HDBlock']['Subject'], comment=comment,
                              date=ddate, time=info['HDBlock']['Time'])

        data_groups = info['DGBlock']  # parse all data groups
        if self._noDataLoading and channel_list is not None:
            data_groups = [self[channel][idField][0] for channel in channel_list]

        # Read data from file
        for dataGroup in data_groups:
            channel_set = channel_set_file
            if info['DGBlock'][dataGroup]['numberOfChannelGroups'] > 0 and \
                    (channel_set is None or
                     len(channel_set & info['ChannelNamesByDG'][dataGroup]) > 0):  # data exists
                if minimal > 1 and not self._noDataLoading:  # load CG, CN and CC block info
                    info.read_cg_block(info.fid, dataGroup, minimal=minimal)
                # Pointer to data block
                pointer_to_data = info['DGBlock'][dataGroup]['pointerToDataRecords']

                if 'dataClass' not in info['DGBlock'][dataGroup]:
                    buf = DATA(info.fid, pointer_to_data)
                    for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                        temp = Record(dataGroup, channelGroup)  # create record class
                        temp.load_info(info)  # load all info related to record

                        if temp.numberOfRecords != 0:  # continue if there are at least some records
                            buf.add_record(temp)
                    if self._noDataLoading:
                        self.info['DGBlock'][dataGroup]['dataClass'] = buf
                else:
                    buf = self.info['DGBlock'][dataGroup]['dataClass']

                buf.read(channel_set, self.fileName)  # reads datablock potentially containing several channel groups

                channel_groups = buf
                if self._noDataLoading and channel_list is not None:
                    channel_groups = [info['CGBlock'][dataGroup][self[channel][idField][1]]['recordID']
                                      for channel in channel_list]

                for recordID in channel_groups:
                    if recordID in buf and 'record' in buf[recordID]:
                        master_channel = buf[recordID]['record'].master['name']

                        if channel_list is None or not self._noDataLoading:
                            channels = (c for c in buf[recordID]['record']
                                        if channel_set is None or c.name in channel_set)
                        else:
                            channels = buf[recordID]['record']
                            channels = [channels[self[channel][idField][2]] for channel in channel_list]

                        for chan in channels:  # for each channel
                            # in case record is used for several channels
                            if channel_set is None and not buf[recordID]['record'].hiddenBytes \
                                    and buf[recordID]['record'].byte_aligned:
                                record_name = buf[recordID]['record'].\
                                    recordToChannelMatching[chan.name]
                            else:
                                record_name = chan.name
                            temp = buf[recordID]['data'][record_name]

                            if len(temp) != 0:
                                # Process concatenated bits inside uint8
                                if chan.bit_masking_needed:
                                    # if channel data do not use complete bytes
                                    if chan.signalDataType in (0, 1, 9, 10, 13, 14):  # integers
                                        if chan.embedding_channel_bitOffset > 0:
                                            temp = right_shift(temp, chan.embedding_channel_bitOffset)
                                        mask = int(pow(2, chan.bitCount) - 1)  # masks isBitUint8
                                        temp = bitwise_and(temp, mask)
                                    else:  # should not happen
                                        warn('bit count and offset not applied to correct data type')
                                self.add_channel(chan.name, temp, master_channel, master_type=1, unit=chan.unit,
                                                 description=chan.desc, conversion=chan.conversion, info=None,
                                                 compression=compression)
                        buf[recordID].pop('data', None)
                del buf
                if minimal > 1:
                    # clean CN, CC and CG info to free memory
                    info.clean_dg_info(dataGroup)
        info.fid.close()  # close file
        if convert_after_read and not compression:
            self._noDataLoading = False
            self._convert_all_channel3()

    def _get_channel_data3(self, channel_name, raw_data=False):
        """Returns channel numpy array

        Parameters
        ----------------
        channel_name : str
            channel name
        raw_data: bool
            flag to return non converted data

        Returns
        -----------
        numpy array
            converted, if not already done, data corresponding to channel name

        Notes
        ------
        This method is the safest to get channel data as numpy array from 'data' dict key might contain raw data
        """
        if channel_name in self:
            vector = self.get_channel(channel_name)[dataField]
            if vector is None:  # noDataLoading reading argument flag activated
                if self.info.fid is None or (self.info.fid is not None and self.info.fid.closed):
                    (self.info.fid, self.info.fileName, zipfile) = _open_mdf(self.fileName)
                self.read3(file_name=None, info=self.info, channel_list=[channel_name], convert_after_read=False)
            if not raw_data:
                return self._convert3(channel_name, self.convertTables)
            else:
                return self.get_channel(channel_name)[dataField]
        else:
            return None

    def _convert3(self, channel_name, convert_tables=False):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channel_name : str
            Name of channel
        convert_tables : bool
            activates computation intensive loops for conversion with tables. Default is False

        Returns
        -----------
        numpy array
            returns numpy array converted to physical values according to conversion type
        """

        if self[channel_name][dataField] is None:
            vector = self[channel_name][dataField]
        else:
            if isinstance(self[channel_name][dataField], CompressedData):
                vector = self[channel_name][dataField].decompression()  # uncompress blosc
            else:
                vector = self[channel_name][dataField][:]  # to have bcolz uncompressed data
        if conversionField in self[channel_name]:  # there is conversion property
            conversion = self[channel_name][conversionField]
            if conversion['type'] == 0:
                return _linear_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 1:
                return _tab_interp_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 2:
                return _tab_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 6:
                return _polynomial_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 7:
                return _exponential_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 8:
                return _log_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 9:
                return _rational_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 10:
                return _formula_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 11 and convert_tables:
                return _text_table_conversion(vector, conversion['parameters'])
            elif conversion['type'] == 12 and convert_tables:
                return _text_range_table_conversion(vector, conversion['parameters'])
            else:
                return vector
        else:
            return vector

    def _convert_channel3(self, channel_name):
        """converts specific channel from raw to physical data according to CCBlock information

        Parameters
        ----------------
        channel_name : str
            Name of channel
        """
        self.set_channel_data(channel_name, self._convert3(channel_name, self.convertTables))
        self.remove_channel_conversion(channel_name)

    def _convert_all_channel3(self):
        """Converts all channels from raw data to converted data according to CCBlock information
        Converted data will take more memory.
        """
        if self._noDataLoading:  # no data loaded, load everything
            self.read3(self.fileName, convert_after_read=True)
        else:
            for channel in self:
                self._convert_channel3(channel)

    def write3(self, file_name=None):
        """Writes simple mdf 3.3 file

        Parameters
        ----------------
        file_name : str, optional
            Name of file
            If file name is not input, written file name will be the one read with
            appended '_new' string before extension

        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """

        # put master channel in first position for each datagroup if not already the case
        for master in self.masterChannelList:
            master_list = sorted(self.masterChannelList[master])
            master_position = master_list.index(master)
            master_list.pop(master_position)  # remove  master channel
            master_list.insert(0, master)  # insert at first position master channel
            self.masterChannelList[master] = master_list

        pointers = {}  # records pointers of blocks when writing

        # write pointer of block and come back to current stream position
        def write_pointer(f, pointer, value):
            current_position = f.tell()
            f.seek(pointer)
            f.write(pack('I', value))
            f.seek(current_position)

        fid = open(file_name, 'wb')  # buffering should automatically be set
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
        n_data_group = len(self.masterChannelList)
        if self.fileMetadata['author'] is not None:  # Author
            try:
                author = '{:\x00<32.31}'.format(self.fileMetadata['author']).encode('latin-1', 'ignore')
            except UnicodeEncodeError:
                author = b'\x00' * 32
        elif os.name == 'posix':
            author = '{:\x00<32.31}'.format(getlogin()).encode('latin-1')
        else:
            author = b'\x00' * 32
        if self.fileMetadata['organisation'] is not None:  # Organization
            try:
                organization = '{:\x00<32.31}'.format(self.fileMetadata['organisation']).encode('latin-1', 'ignore')
            except UnicodeEncodeError:
                organization = b'\x00' * 32
        else:
            organization = b'\x00' * 32
        if self.fileMetadata['project'] is not None:  # Project
            try:
                project = '{:\x00<32.31}'.format(self.fileMetadata['project']).encode('latin-1', 'ignore')
            except UnicodeEncodeError:
                project = b'\x00' * 32
        else:
            project = b'\x00' * 32
        if self.fileMetadata['subject'] is not None:  # Subject
            try:
                subject = '{:\x00<32.31}'.format(self.fileMetadata['subject']).encode('latin-1', 'ignore')
            except UnicodeEncodeError:
                subject = b'\x00' * 32
        else:
            subject = b'\x00' * 32
        # Header Block, block size
        # first Data block pointer, pointer to TX Block file comment, pointer to PR Block
        # number of data groups, date, time
        # Time Stamp, UTC time offset, Time quality, Timer identification
        time_offset = gmtime()
        head = (b'HD', 208, 272, 0, 0, n_data_group,
                '{:\x00<10}'.format(strftime("%d:%m:%Y", time_offset)).encode('latin-1'),
                '{:\x00<8}'.format(strftime("%H:%M:%S", time_offset)).encode('latin-1'),
                author, organization, project, subject,
                int(time() * 1000000000), 1, 0,
                b'Local PC Reference Time         ')
        fid.write(pack('<2sH3IH10s8s32s32s32s32sQhH32s', *head))

        # write DG block
        pointers['DG'] = {}
        pointers['CG'] = {}
        pointers['CN'] = {}
        data_group = 0

        for masterChannel in self.masterChannelList:
            # writes dataGroup Block
            pointers['DG'][data_group] = {}
            position = fid.tell()
            if 0 < data_group:  # not possible for first DG
                # previous datagroup pointer to this new datagroup
                write_pointer(fid, pointers['DG'][data_group - 1]['nextDG'], position)
            else:
                # first datagroup pointer in header block
                write_pointer(fid, pointers['HD']['DG'], position)
            pointers['DG'][data_group]['nextDG'] = position + 4
            pointers['DG'][data_group]['CG'] = position + 8
            pointers['DG'][data_group]['data'] = position + 16
            # DG block size, pointer to next DataGroup,
            # pointer to channel group, pointer to trigger block, pointer to data block
            # number of channel group, number of record IDs, reserved
            head = (b'DG', 28, 0, 0, 0, 0, 1, 0, b'\x00' * 32)
            fid.write(pack('<2sH4I2H4s', *head))

            # sorted data so only one channel group
            pointers['CG'][data_group] = {}
            # write first CG pointer in datagroup
            position = fid.tell()
            write_pointer(fid, pointers['DG'][data_group]['CG'], position)
            pointers['CG'][data_group]['firstCN'] = position + 8
            pointers['CG'][data_group]['TX'] = position + 12
            pointers['CG'][data_group]['dataRecordSize'] = position + 20
            num_channels = len(self.masterChannelList[masterChannel])
            master_data = self.get_channel_data(masterChannel)
            n_records = len(master_data)
            # CG block size, pointer to next Channel Group
            # pointer to first channel block, pointer to TX block
            # No record ID no need for sorted data, Number of channels
            # Size of data record, Number of records, pointer to sample reduction block
            head = (b'CG', 30, 0, 0, 0, 0, num_channels, 0, n_records, 0)
            fid.write(pack('<2sH3I3H2I', *head))

            # Channel blocks writing
            pointers['CN'][data_group] = {}
            data_list = ()
            data_type_list = ''
            record_number_of_bits = 0
            preceeding_channel = None
            bit_offset = 0
            byte_offset = 0
            # first channel bock pointer from CG
            position = fid.tell()
            write_pointer(fid, pointers['CG'][data_group]['firstCN'], position)
            for channel in self.masterChannelList[masterChannel]:
                position = fid.tell()
                pointers['CN'][data_group][channel] = {}
                pointers['CN'][data_group][channel]['beginCN'] = position
                pointers['CN'][data_group][channel]['nextCN'] = position + 4
                pointers['CN'][data_group][channel]['CC'] = position + 8
                pointers['CN'][data_group][channel]['TX'] = position + 16
                if preceeding_channel is not None:  # not possible for first CN
                    # pointer in previous CN
                    write_pointer(fid, pointers['CN'][data_group][preceeding_channel]['nextCN'],
                                  pointers['CN'][data_group][channel]['beginCN'])
                preceeding_channel = channel
                if channel not in set(self.masterChannelList):
                    master_flag = 0  # data channel
                else:
                    master_flag = 1  # master channel
                desc = self.get_channel_desc(channel)
                # if bitOffset exceeds two bytes limit, we start using the byte offset field
                if bit_offset > 0xFFFF:
                    bit_offset -= 0x10000
                    byte_offset += 8192
                data = self.get_channel_data(channel)  # channel data
                temp = data
                data_list = data_list + (temp, )
                cn_numpy_kind = data.dtype.kind
                cn_numpy_item_size = data.dtype.itemsize
                number_of_bits = cn_numpy_item_size * 8
                if cn_numpy_kind in ('u', 'b'):
                    data_type = 0
                elif cn_numpy_kind == 'i':
                    data_type = 1
                elif cn_numpy_kind == 'f':
                    if cn_numpy_item_size == 8:
                        data_type = 3
                    elif cn_numpy_item_size == 4:
                        data_type = 2
                    else:
                        raise Exception('Not recognized dtype')
                elif cn_numpy_kind in ('S', 'U', 'V'):
                    data_type = 7
                else:
                    raise Exception('Not recognized dtype')
                    return data.dtype
                if data.dtype.kind not in ['S', 'U']:
                    data_type_list = ''.join([data_type_list, data.dtype.char])
                else:
                    data_type_list = ''.join([data_type_list, '{}s'.format(data.dtype.itemsize)])
                    number_of_bits = 8 * data.dtype.itemsize
                record_number_of_bits += number_of_bits
                if data.dtype.kind not in ['S', 'U']:
                    value_range_valid = 1
                    if len(data) > 0 and issubdtype(data.dtype, numpy_number):
                        maximum = npmax(data)
                        minimum = npmin(data)
                    else:
                        maximum = 0
                        minimum = 0
                else:
                    value_range_valid = 0  # No value range valid
                    minimum = 0  # Min value
                    maximum = 0  # Max value
                pointers['CN'][data_group][channel]['longChannelName'] = position + 218
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
                    description = '{:\x00<128.127}'.format(desc).encode('latin-1')
                except:
                    description = b'\x00' * 128
                head = (b'CN', 228, 0, 0, 0, 0, 0, master_flag,
                        ('{:\x00<32.31}'.format(channel) + '\x00').encode('latin-1'),
                        description,
                        bit_offset, number_of_bits, data_type, value_range_valid,
                        minimum, maximum, 0, 0, 0, byte_offset)
                fid.write(pack('<2sH5IH32s128s4H3d2IH', *head))
                bit_offset += number_of_bits

                # TXblock for long channel name
                write_pointer(fid, pointers['CN'][data_group][channel]['longChannelName'], fid.tell())
                head = (b'TX', len(channel) + 4 + 1, channel.encode('latin-1') + b'\x00')
                fid.write(pack('<2sH{}s'.format(len(channel)+1), *head))

                # Conversion blocks writing
                write_pointer(fid, pointers['CN'][data_group][channel]['CC'], fid.tell())
                # channel description
                # conversion already done during reading
                # additional size information, not necessary for 65535 conversion type ?
                try:
                    unit = '{:\x00<20.19}'.format(self.get_channel_unit(channel)
                                                  .encode('latin-1', 'replace'))
                except:
                    unit = b'\x00' * 20
                head = (b'CC', 46, value_range_valid, minimum, maximum,
                        unit, 65535, 0)
                fid.write(pack('<2shH2d20s2H', *head))

            # number of channels in CG
            current_position = fid.tell()
            fid.seek(pointers['CG'][data_group]['dataRecordSize'])
            fid.write(pack('H', int(record_number_of_bits / 8)))  # Size of data record
            fid.seek(current_position)

            # data writing
            # write data pointer in datagroup
            write_pointer(fid, pointers['DG'][data_group]['data'], fid.tell())
            # dumps data vector from numpy
            fid.write(fromarrays(data_list).tobytes(order='F'))

            data_group += 1

        # print(pointers, file=stderr)
        fid.close()
