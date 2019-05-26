# -*- coding: utf-8 -*-
""" mdf_skeleton module describing basic mdf structure and methods

Created on Thu Sept 24 2015

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__


Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>


mdf
--------------------------
"""
from __future__ import absolute_import  # for consistency between python 2 and 3
from __future__ import print_function
from io import open
from zipfile import is_zipfile, ZipFile
from itertools import chain
from random import choice
from string import ascii_letters
from sys import version_info
from collections import OrderedDict, defaultdict
from warnings import warn
from numpy import array_repr, set_printoptions, recarray, fromstring
try:
    from pandas import set_option
except ImportError:
    warn('Pandas export will not be possible', ImportWarning)
set_printoptions(threshold=100, edgeitems=1)
_notAllowedChannelNames = set(dir(recarray))
try:
    CompressionPossible = True
    from blosc import compress, decompress
except ImportError:
    # Cannot compress data, please install bcolz and blosc
    CompressionPossible = False
PythonVersion = version_info
PythonVersion = PythonVersion[0]

descriptionField = 'description'
unitField = 'unit'
dataField = 'data'
masterField = 'master'
masterTypeField = 'masterType'
conversionField = 'conversion'
attachmentField = 'attachment'
idField = 'id'
invalidPosField = 'invalid_bit'
invalidChannel = 'invalid_channel'


class MdfSkeleton(dict):
    __slots__ = ['masterChannelList', 'fileName', 'MDFVersionNumber', 'multiProc',
                 'convertAfterRead', 'filterChannelNames', 'fileMetadata', 'convertTables',
                 '_pandasframe', 'info', '_compression_level', '_noDataLoading',
                 'fid', 'zipfile']
    """ MdfSkeleton class

    Attributes
    --------------
    fileName : str
        file name
    MDFVersionNumber : int
        mdf file version number
    masterChannelList : dict
        Represents data structure: a key per master channel with corresponding value containing a list of channels
        One key or master channel represents then a data group having same sampling interval.
    multiProc : bool
        Flag to request channel conversion multi processed for performance improvement.
        One thread per data group.
    convertAfterRead : bool
        flag to convert raw data to physical just after read
    filterChannelNames : bool
        flag to filter long channel names from its module names separated by '.'
    fileMetadata : dict
        file metadata with minimum keys : author, organisation, project, subject, comment, time, date

    Methods
    ------------
    add_channel(channel_name, data, master_channel, master_type=1, unit='', description='', conversion=None)
        adds channel to mdf dict
    remove_channel(channel_name)
        removes channel from mdf dict and returns its content
    rename_channel(channel_name, new_channel_name)
         renames a channel and returns its content
    copy()
        copy a mdf class
    add_metadata(author, organisation, project, subject, comment, date, time)
        adds basic metadata from file
    """

    def __init__(self, file_name=None, channel_list=None, convert_after_read=True,
                 filter_channel_names=False, no_data_loading=False,
                 compression=False, convert_tables=False, metadata=2):
        """ mdf_skeleton class constructor.

        Parameters
        ----------------
        file_name : str, optional
            file name

        channel_list : list of str, optional
            list of channel names to be read.
            If you use channelList, reading might be much slower but it will
            save you memory. Can be used to read big files.

        convert_after_read : bool, optional
            flag to convert channel after read, True by default.
            If you use convertAfterRead by setting it to false, all data from
            channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times
            memory footprint.
            To calculate value from channel, you can then use
            method .get_channel_data()

        filter_channel_names : bool, optional
            flag to filter long channel names from its module names
            separated by '.'

        compression : bool optional
            flag to compress data in memory.

        convert_tables : bool, optional, default False
            flag to convert or not only conversions with tables.
            These conversions types take generally long time and memory.
        """
        self.masterChannelList = OrderedDict()
        # flag to control multiprocessing, default deactivate,
        # giving priority to mdfconverter
        self.multiProc = False
        self.fileMetadata = dict()
        self.fileMetadata['author'] = ''
        self.fileMetadata['organisation'] = ''
        self.fileMetadata['project'] = ''
        self.fileMetadata['subject'] = ''
        self.fileMetadata['comment'] = ''
        self.fileMetadata['time'] = ''
        self.fileMetadata['date'] = ''
        self.MDFVersionNumber = 300
        self.filterChannelNames = filter_channel_names
        # by default, do not convert table conversion types, taking lot of time and memory
        self.convertTables = convert_tables
        self._pandasframe = False
        self.info = None
        self._compression_level = 9  # default compression level
        self._noDataLoading = False  # in case reading with this argument activated
        # clears class from previous reading and avoid to mess up
        self.clear()
        self.fileName = file_name
        if file_name is not None:
            self.read(file_name, channel_list=channel_list,
                      convert_after_read=convert_after_read,
                      filter_channel_names=filter_channel_names,
                      no_data_loading=no_data_loading,
                      compression=compression,
                      metadata=metadata)

    def add_channel(self, channel_name, data, master_channel, master_type=1, unit='', description='', conversion=None,
                    info=None, compression=False, identifier=None):
        """ adds channel to mdf dict.

        Parameters
        ----------------
        channel_name : str
            channel name
        data : numpy array
            numpy array of channel's data
        master_channel : str
            master channel name
        master_type : int, optional
            master channel type : 0=None, 1=Time, 2=Angle, 3=Distance, 4=index
        unit : str, optional
            unit description
        description : str, optional
            channel description
        conversion : info class, optional
            conversion description from info class
        info : info class for CNBlock, optional
            used for CABlock axis creation and channel conversion
        compression : bool
            flag to ask for channel data compression
        identifier : tuple
            tuple of int and str following below structure:
            (data group number, channel group number, channel number),
            (channel name, channel source, channel path),
            (group name, group source, group path)
        """
        if not self._noDataLoading:
            self[channel_name] = {}
            if master_channel not in self.masterChannelList:
                self.masterChannelList[master_channel] = []
            self.masterChannelList[master_channel].append(channel_name)
            self.set_channel_unit(channel_name, unit)
            self.set_channel_desc(channel_name, description)
            self.set_channel_master(channel_name, master_channel)
            if self.MDFVersionNumber < 400:  # mdf3
                self.set_channel_master_type(channel_name, 1)
            else:  # mdf4
                self.set_channel_master_type(channel_name, master_type)
        self.set_channel_data(channel_name, data, compression)
        if conversion is not None:
            self[channel_name]['conversion'] = {}
            self[channel_name]['conversion']['type'] = conversion['cc_type']
            self[channel_name]['conversion']['parameters'] = {}
            if self.MDFVersionNumber < 400:  # mdf3
                if 'conversion' in conversion:
                    self[channel_name]['conversion']['parameters'] = \
                        conversion['conversion']
                if conversion['cc_type'] == 0 and \
                        'P2' in self[channel_name]['conversion']['parameters'] and \
                        (self[channel_name]['conversion']['parameters']['P2'] == 1.0 and
                         self[channel_name]['conversion']['parameters']['P1'] in (0.0, -0.0)):
                    self[channel_name].pop('conversion')
            else:  # mdf4
                if 'cc_val' in conversion:
                    self[channel_name]['conversion']['parameters']['cc_val'] = \
                        conversion['cc_val']
                if 'cc_ref' in conversion:
                    self[channel_name]['conversion']['parameters']['cc_ref'] = \
                        conversion['cc_ref']
        if info is not None:  # axis from CABlock
            ca_block = info
            axis = []
            while 'CABlock' in ca_block:
                ca_block = ca_block['CABlock']
                if 'ca_axis_value' in ca_block:
                    if type(ca_block['ca_dim_size']) is list:
                        index = 0
                        for ndim in ca_block['ca_dim_size']:
                            axis.append(tuple(ca_block['ca_axis_value'][index:index+ndim]))
                            index += ndim
                    else:
                        axis = ca_block['ca_axis_value']
            self[channel_name]['axis'] = axis
        if identifier is not None:
            self[channel_name]['id'] = identifier

    def remove_channel(self, channel_name):
        """ removes channel from mdf dict.

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        value of mdf dict key=channel_name
        """
        self.masterChannelList[self.get_channel_master(channel_name)].remove(channel_name)
        return self.pop(channel_name)

    def rename_channel(self, channel_name, new_name):
        """Modifies name of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        new_name : str
            new channel name
        """
        if channel_name in self and new_name not in self:
            # add the new name to the same master
            self.masterChannelList[self.get_channel_master(channel_name)].append(new_name)
            # remove the old name
            self.masterChannelList[self.get_channel_master(channel_name)].remove(channel_name)
            self[new_name] = self.pop(channel_name)  # copy the data
            if channel_name in self.masterChannelList:  # it is a master channel
                self.masterChannelList[new_name] = self.pop(channel_name)
                for channel in self.masterChannelList[new_name]:
                    self.set_channel_master(channel, new_name)
            return self[new_name]
        else:
            return None

    def remove_channel_conversion(self, channel_name):
        """ removes conversion key from mdf channel dict.

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        removed value from dict
        """
        return self._remove_channel_field(channel_name, conversionField)

    def _remove_channel_field(self, channel_name, field):
        """general purpose function to remove key from channel dict in mdf

        Parameters
        ----------------
        channel_name : str
            channel name
        field : str
            channel dict key

        Returns
        -------
        removed value from dict
        """
        if field in self.get_channel(channel_name):
            return self[channel_name].pop(field)

    def get_channel_unit(self, channel_name):
        """Returns channel unit string
        Implemented for a future integration of pint

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -----------
        str
            unit string description
        """
        temp = self._get_channel_field(channel_name, field=unitField)
        if isinstance(temp, (dict, defaultdict)):
            try:
                return temp['Comment']
            except KeyError:
                return ''
        return temp

    def get_channel_desc(self, channel_name):
        """Extract channel description information from mdf structure

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        channel description string
        """
        return self._get_channel_field(channel_name, field=descriptionField)

    def get_channel_master(self, channel_name):
        """Extract channel master name from mdf structure

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        channel master name string
        """
        return self._get_channel_field(channel_name, field=masterField)

    def get_channel_master_type(self, channel_name):
        """Extract channel master type information from mdf structure

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        channel mater type integer : 0=None, 1=Time, 2=Angle, 3=Distance, 4=index
        """
        return self._get_channel_field(channel_name, field=masterTypeField)

    def get_channel_conversion(self, channel_name):
        """Extract channel conversion dict from mdf structure

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        channel conversion dict
        """
        return self._get_channel_field(channel_name, field=conversionField)

    def get_invalid_bit(self, channel_name):
        return self._get_channel_field(channel_name, field=invalidPosField)

    def get_invalid_channel(self, channel_name):
        return self._get_channel_field(channel_name, field=invalidChannel)

    def get_channel(self, channel_name):
        """Extract channel dict from mdf structure

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        channel dictionnary containing data, description, unit, etc.
        """
        try:
            return self[channel_name]
        except KeyError:
            return None

    def _get_channel_field(self, channel_name, field=None):
        """General purpose function to extract channel dict key value from mdf class

        Parameters
        ----------------
        channel_name : str
            channel name
        field : str
            channel dict key

        Returns
        -------
        channel description string
        """
        channel = self.get_channel(channel_name)
        if channel is not None:
            try:
                return channel[field]
            except KeyError:
                return ''
        else:
            return None

    def set_channel_unit(self, channel_name, unit):
        """Modifies unit of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        unit : str
            channel unit
        """
        self._set_channel(channel_name, unit, field=unitField)

    def set_channel_data(self, channel_name, data, compression=False):
        """Modifies data of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        data : numpy array
            channel data
        compression : bool or str
            trigger for data compression
        """
        if compression and CompressionPossible:
            temp = CompressedData()
            temp.compression(data)
            self._set_channel(channel_name, temp, field=dataField)
        else:
            self._set_channel(channel_name, data, field=dataField)

    def set_channel_desc(self, channel_name, desc):
        """Modifies description of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        desc : str
            channel description
        """
        self._set_channel(channel_name, desc, field=descriptionField)

    def set_channel_master(self, channel_name, master):
        """Modifies channel master name

        Parameters
        ----------------
        channel_name : str
            channel name
        master : str
            master channel name
        """
        self._set_channel(channel_name, master, field=masterField)

    def set_channel_master_type(self, channel_name, master_type):
        """Modifies master channel type

        Parameters
        ----------------
        channel_name : str
            channel name
        master_type : int
            master channel type
        """
        self._set_channel(channel_name, master_type, field=masterTypeField)

    def set_channel_conversion(self, channel_name, conversion):
        """Modifies conversion dict of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        conversion : dict
            conversion dictionary
        """
        self._set_channel(channel_name, conversion, field=conversionField)

    def set_channel_attachment(self, channel_name, attachment):
        """Modifies channel attachment

        Parameters
        ----------------
        channel_name : str
            channel name
        attachment
            channel attachment
        """
        self._set_channel(channel_name, attachment, field=attachmentField)

    def set_invalid_bit(self, channel_name, bit_position):
        """returns the invalid bit position of channel

        Parameters
        ----------------
        channel_name : str
            channel name
        bit_position
            invalid bit position of channel within invalid channel bytes
        Returns
        -------
        bit position
        """
        self[channel_name][invalidPosField] = bit_position

    def set_invalid_channel(self, channel_name, invalid_channel):
        self[channel_name][invalidChannel] = invalid_channel

    def _set_channel(self, channel_name, item, field=None):
        """General purpose method to modify channel values

        Parameters
        ----------------
        channel_name : str
            channel name
        item
            new replacing item
        field : str
            channel dict key of item
        """
        try:
            self[channel_name][field] = item
        except KeyError:
            warn('Channel {} not in dictionary'.format(channel_name))

    def _channel_in_mdf(self, channel_name):
        """Efficiently assess if channel is already in mdf

        Parameters
        ----------------
        channel_name : str
            channel name

        Returns
        -------
        bool
        """
        return channel_name in self.masterChannelList[masterField] or channel_name in self.masterChannelList

    def add_metadata(self, author='', organisation='', project='',
                     subject='', comment='', date='', time=''):
        """adds basic metadata to mdf class

        Parameters
        ----------------
        author : str
            author of file
        organisation : str
            organisation of author
        project : str
        subject : str
        comment : str
        date : str
        time : str

        Notes
        =====
        All fields are optional, default being empty string
        """
        self.fileMetadata['author'] = author
        self.fileMetadata['organisation'] = organisation
        self.fileMetadata['project'] = project
        self.fileMetadata['subject'] = subject
        self.fileMetadata['comment'] = comment
        self.fileMetadata['date'] = date
        self.fileMetadata['time'] = time

    def __str__(self):
        """representation a mdf_skeleton class data structure

        Returns
        ------------
        string of mdf class ordered as below
        "master_channel_name
            channel_name   description
            numpy_array    unit"
        """
        output = list()
        if self.fileName is not None:
            output.append('file name : {}\n'.format(self.fileName))
        else:
            output.append('')
        for m in self.fileMetadata.keys():
            if self.fileMetadata[m] is not None:
                output.append('{} : {}\n'.format(m, self.fileMetadata[m]))
        if not self._pandasframe:
            output.append('\nchannels listed by data groups:\n')
            for d in self.masterChannelList.keys():
                if d is not None:
                    output.append('{}\n'.format(d))
                for c in self.masterChannelList[d]:
                    output.append('  {} : '.format(c))
                    desc = self.get_channel_desc(c)
                    if desc is not None:
                        try:
                            output.append(str(desc))
                        except:
                            pass
                    output.append('\n    ')
                    data = self.get_channel_data(c)
                    # not byte, impossible to represent
                    if data.dtype.kind != 'V':
                        output.append(array_repr(data[:],
                                      precision=3, suppress_small=True))
                    unit = self.get_channel_unit(c)
                    if unit is not None:
                        output.append(' {}\n'.format(unit))
            return ''.join(output)
        else:
            set_option('max_rows', 3)
            set_option('expand_frame_repr', True)
            set_option('max_colwidth', 6)
            for master in self.masterGroups:
                output.append(master)
                output.append(str(self[master]))
            return ''.join(output)

    def copy(self):
        """copy a mdf class

        Returns
        ------------
        mdf_skeleton: class instance
            copy of a mdf_skeleton class
        """
        yop = MdfSkeleton()
        yop.multiProc = self.multiProc
        yop.fileName = self.fileName
        yop.masterChannelList = self.masterChannelList
        yop.fileMetadata = self.fileMetadata
        yop.MDFVersionNumber = self.MDFVersionNumber
        yop.filterChannelNames = self.filterChannelNames
        yop.convertTables = self.convert_tables
        for channel in self:
            yop[channel] = self[channel]
        return yop


def _open_mdf(file_name):
    """ Opens mdf, make a few checks and returns fid

    Parameters
    -----------
    file_name : str
        filename string

    Returns
    --------
    fid
        file identifier
    """

    try:
        fid = open(file_name, 'rb')
    except IOError:
        raise Exception('Can not find file {}'.format(file_name))
    zipfile = False
    # Check whether file is MDF file -- assumes that every MDF file starts
    # with the letters MDF
    if fid.read(3) not in ('MDF', b'MDF'):
        if is_zipfile(file_name):
            # this is .mfxz file, compressed zip file
            zipfile = True
            fid.close()
            zip_class = ZipFile(file_name, 'r')
            zip_name = zip_class.namelist()[0]  # there should be only one file
            zip_name = zip_class.extract(zip_name)  # locally extracts file
            fid = open(zip_name, 'rb')
            file_name = zip_name
        else:
            raise Exception('file {} is not an MDF file!'.format(file_name))
    return (fid, file_name, zipfile)


def _bits_to_bytes_aligned(n_bits, numeric=True):
    """ Converts number of bits into number of aligned bytes

    Parameters
    -------------
    n_bits : int
        number of bits
    numeric: bool
        flag to indicate channel is numeric

    Returns
    ----------
    number of equivalent aligned bytes
    """
    if numeric:
        if n_bits == 0:
            n_bytes = 0
        elif n_bits <= 8:
            n_bytes = 1
        elif n_bits <= 16:
            n_bytes = 2
        elif n_bits <= 32:
            n_bytes = 4
        elif n_bits <= 64:
            n_bytes = 8
        else:
            warn('error converting bits into bytes for a numeric channel, too many bits')
            n_bytes = _bits_to_bytes_not_aligned(n_bits)
    else:
        n_bytes = _bits_to_bytes_not_aligned(n_bits)
    return n_bytes


def _bits_to_bytes_not_aligned(n_bits):
    """ Converts number of bits into number of not aligned bytes

    Parameters
    -------------
    n_bits : int
        number of bits

    Returns
    ----------
    number of equivalent not aligned bytes
    """
    n_bytes = n_bits // 8
    if not n_bits % 8 == 0:
        n_bytes += 1
    return n_bytes


def _convert_name(channel_name):
    """ Check if channelName is valid python identifier
    to be removed with next function if no more need
    """

    if PythonVersion < 3:  # python 2
        channel_identifier = _sanitize_identifier(channel_name).encode('utf-8')
    else:  # python 3
        if channel_name.isidentifier():
            return channel_name
        else:
            channel_identifier = str(_sanitize_identifier(channel_name))
    # all characters of channel are not compliant to python
    if not channel_identifier:
        # generate random name for recarray
        channel_identifier = ''.join([choice(ascii_letters) for n in range(32)])
    if channel_identifier in _notAllowedChannelNames:
        # limitation from recarray object attribute
        channel_identifier = ''.join([channel_identifier, '_'])
    return channel_identifier


def _gen_valid_identifier(seq):
    # get an iterator
    itr = iter(seq)
    # pull characters until we get a legal one for first in identifier
    for ch in itr:
        if ch == '_' or ch.isalpha():
            yield ch
            break
        elif ch.isdigit():
            itr = chain(itr, ch)

    # pull remaining characters and yield legal ones for identifier
    for ch in itr:
        if ch == '_' or ch.isalpha() or ch.isdigit():
            yield ch
        else:
            yield '_'


def _sanitize_identifier(name):
    return ''.join(_gen_valid_identifier(name))


class CompressedData:
    __slots__ = ['data', 'dtype']
    """ class to represent compressed data by blosc
    """
    def __init__(self):
        """ data compression method

        Attributes
        -------------
        data : numpy array compressed
            compressed data
        dtype : numpy dtype object
            numpy array dtype
        """
        self.data = None
        self.dtype = None

    def compression(self, a):
        """ data compression method

        Parameters
        -------------
        a : numpy array
            data to be compresses
        """
        self.data = compress(a.tobytes())
        self.dtype = a.dtype

    def decompression(self):
        """ data decompression

        Returns
        -------------
        uncompressed numpy array
        """
        return fromstring(decompress(self.data), dtype=self.dtype)

    def __str__(self):
        """ prints compressed_data object content
        """
        return self.decompression()
