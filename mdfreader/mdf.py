# -*- coding: utf-8 -*-
""" mdf_skeleton module describing basic mdf structure and methods

Created on Thu Sept 24 2015

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

:Author: `Aymeric Rateau <https://github.com/ratal/mdfreader>`__

Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>

mdf module
--------------------------
"""
from numpy import array_repr, set_printoptions, recarray
import pandas as pd
from io import open
from zipfile import is_zipfile, ZipFile
from itertools import chain
from random import choice
from string import ascii_letters
set_printoptions(threshold=100, edgeitems=1)
from sys import version_info
PythonVersion = version_info
PythonVersion = PythonVersion[0]
_notAllowedChannelNames = set(dir(recarray))

descriptionField = 'description'
unitField = 'unit'
dataField = 'data'
masterField = 'master'
masterTypeField = 'masterType'
conversionField = 'conversion'
attachmentField = 'attachment'

class mdf_skeleton(dict):

    """ mdf_skeleton class

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
    file_metadata : dict
        file metadata with minimum keys : author, organisation, project, subject, comment, time, date

    Methods
    ------------
    add_channel(channel_name, data, master_channel, master_type=1, unit='', description='', conversion=None)
        adds channel to mdf dict
    remove_channel(channel_name)
        removes channel from mdf dict and returns its content
    copy()
        copy a mdf class
    add_metadata(author, organisation, project, subject, comment, date, time)
        adds basic metadata from file
    """

    def __init__(self, fileName=None, channelList=None, convertAfterRead=True, filterChannelNames=False, noDataLoading=False):
        """ mdf_skeleton class constructor.

        Parameters
        ----------------
        fileName : str, optional
            file name

        channelList : list of str, optional
            list of channel names to be read
            If you use channelList, reading might be much slower but it will save you memory. Can be used to read big files

        convertAfterRead : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false, all data from channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .getChannelData()

        filterChannelNames : bool, optional
            flag to filter long channel names from its module names separated by '.'
        """
        self.masterChannelList = {}
        self.multiProc = False  # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        self.file_metadata = {}
        self.file_metadata['author'] = ''
        self.file_metadata['organisation'] = ''
        self.file_metadata['project'] = ''
        self.file_metadata['subject'] = ''
        self.file_metadata['comment'] = ''
        self.file_metadata['time'] = ''
        self.file_metadata['date'] = ''
        self.MDFVersionNumber = 300
        self.filterChannelNames = filterChannelNames
        self.convert_tables = False
        self._pandasframe = False
        self._info = None
        # clears class from previous reading and avoid to mess up
        self.clear()
        self.fileName = fileName
        if fileName is not None:
            self.read(fileName, channelList=channelList, convertAfterRead=convertAfterRead, filterChannelNames=filterChannelNames, noDataLoading=noDataLoading)


    def add_channel(self, dataGroup, channel_name, data, master_channel, master_type=1, unit='', description='', conversion=None, info=None):
        """ adds channel to mdf dict.

        Parameters
        ----------------
        dataGroup : int
            dataGroup number. Is appended to master name for non unique channel names
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
        """
        if channel_name in self:
            if self[channel_name][dataField] is not None:
                channel_name += '_' + str(dataGroup)
                self[channel_name] = {}
            else:
                pass # using noDataLoading
        else:
            self[channel_name] = {}
        self.setChannelData(channel_name, data)
        self.setChannelUnit(channel_name, unit)
        self.setChannelDesc(channel_name, description)
        self.setChannelMaster(channel_name, master_channel)
        if master_channel not in self.masterChannelList:
            self.masterChannelList[master_channel] = []
        self.masterChannelList[master_channel].append(channel_name)
        if self.MDFVersionNumber < 400: #  mdf3
            self.setChannelMasterType(channel_name, 1)
        else: #  mdf4
            self.setChannelMasterType(channel_name, master_type)
        if conversion is not None: # 
            self[channel_name]['conversion'] = {}
            self[channel_name]['conversion']['type'] = conversion['cc_type']
            self[channel_name]['conversion']['parameters'] = {}
            if self.MDFVersionNumber < 400: #  mdf3
                if 'conversion' in conversion:
                    self[channel_name]['conversion']['parameters'] = conversion['conversion']
                if conversion['cc_type'] == 0 and \
                        (self[channel_name]['conversion']['parameters']['P2'] == 1.0 and \
                        self[channel_name]['conversion']['parameters']['P1'] in (0.0, -0.0)):
                    self[channel_name].pop('conversion')
            else: #  mdf4
                if 'cc_val' in conversion:
                    self[channel_name]['conversion']['parameters']['cc_val'] = conversion['cc_val']
                if 'cc_ref' in conversion:
                    self[channel_name]['conversion']['parameters']['cc_ref'] = conversion['cc_ref']
        if info is not None: # axis from CABlock
            CABlock = info
            axis=[]
            while 'CABlock' in CABlock:
                CABlock = CABlock['CABlock']
                if 'ca_axis_value' in CABlock:
                    if type(CABlock['ca_dim_size']) is list:
                        index = 0
                        for ndim in CABlock['ca_dim_size']:
                            axis.append(tuple(CABlock['ca_axis_value'][index:index+ndim]))
                            index += ndim
                    else:
                        axis = CABlock['ca_axis_value']
            self[channel_name]['axis'] = axis


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
        self.masterChannelList[self.getChannelMaster(channel_name)].remove(channel_name)
        return self.pop(channel_name)


    def remove_channel_conversion(self, channelName):
        """ removes conversion key from mdf channel dict.

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        removed value from dict
        """
        return self._remove_channel_field(channelName, conversionField)


    def _remove_channel_field(self, channelName, field):
        """general purpose function to remove key from channel dict in mdf

        Parameters
        ----------------
        channelName : str
            channel name
        field : str
            channel dict key
        Returns
        -------
        removed value from dict
        """
        if field in self.getChannel(channelName):
            return self[channelName].pop(field)


    def getChannelUnit(self, channelName):
        """Returns channel unit string
        Implemented for a future integration of pint

        Parameters
        ----------------
        channelName : str
            channel name

        Returns
        -----------
        str
            unit string description
        """
        return self._getChannelField(channelName, field = unitField)


    def getChannelDesc(self, channelName):
        """Extract channel description information from mdf structure

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        channel description string
        """
        return self._getChannelField(channelName, field = descriptionField)


    def getChannelMaster(self, channelName):
        """Extract channel master name from mdf structure

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        channel master name string
        """
        return self._getChannelField(channelName, field = masterField)


    def getChannelMasterType(self, channelName):
        """Extract channel master type information from mdf structure

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        channel mater type integer
        """
        return self._getChannelField(channelName, field = masterTypeField)


    def getChannelConversion(self, channelName):
        """Extract channel conversion dict from mdf structure

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        channel conversion dict
        """
        return self._getChannelField(channelName, field = conversionField)


    def getChannel(self, channelName):
        """Extract channel dict from mdf structure

        Parameters
        ----------------
        channelName : str
            channel name
        
        Returns
        -------
        channel dictionnary containing data, description, unit, etc.
        """
        if channelName in self:
            return self[channelName]
        else:
            return None


    def _getChannelField(self, channelName, field=None):
        """General purpose function to extract channel dict key value from mdf class

        Parameters
        ----------------
        channelName : str
            channel name
        field : str
            channel dict key
        Returns
        -------
        channel description string
        """
        channel = self.getChannel(channelName)
        if channelName in self:
            if field in channel:
                return channel[field]
            else:
                return ''
        else:
            return None


    def setChannelUnit(self, channelName, unit):
        """Modifies unit of channel

        Parameters
        ----------------
        channelName : str
            channel name
        unit : str
            channel unit
        """
        self._setChannel(channelName, unit, field = unitField)


    def setChannelData(self, channelName, data):
        """Modifies data of channel

        Parameters
        ----------------
        channelName : str
            channel name
        data : numpy array
            channel data
        """
        self._setChannel(channelName, data, field = dataField)


    def setChannelDesc(self, channelName, desc):
        """Modifies description of channel

        Parameters
        ----------------
        channelName : str
            channel name
        desc : str
            channel description
        """
        self._setChannel(channelName, desc, field = descriptionField)


    def setChannelMaster(self, channelName, master):
        """Modifies channel master name

        Parameters
        ----------------
        channelName : str
            channel name
        master : str
            master channel name
        """
        self._setChannel(channelName, master, field = masterField)


    def setChannelMasterType(self, channelName, masterType):
        """Modifies master channel type

        Parameters
        ----------------
        channelName : str
            channel name
        masterType : int
            master channel type
        """
        self._setChannel(channelName, masterType, field = masterTypeField)


    def setChannelConversion(self, channelName, conversion):
        """Modifies conversion dict of channel

        Parameters
        ----------------
        channelName : str
            channel name
        conversion : dict
            conversion dictionnary
        """
        self._setChannel(channelName, conversion, field = conversionField)


    def setChannelAttachment(self, channelName, attachment):
        """Modifies channel attachment

        Parameters
        ----------------
        channelName : str
            channel name
        attachment
            channel attachment
        """
        self._setChannel(channelName, attachment, field = attachmentField)


    def _setChannel(self, channelName, item, field=None):
        """General purpose method to modify channel values

        Parameters
        ----------------
        channelName : str
            channel name
        item
            new replacing item
        field : str
            channel dict key of item
        """
        if channelName in self:
            self[channelName][field] = item
        else:
            raise KeyError('Channel not in dictionary')

    def _channelInMDF(channelName):
        """Efficiently assess if channel is already in mdf

        Parameters
        ----------------
        channelName : str
            channel name
        
        Return
        -------
        bool
        """
        return channelName in self.masterChannelList[masterField] \
                or channelName in self.masterChannelList

    def add_metadata(self, author='', organisation='', project='', \
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

        Note
        =====
        All fields are optional, default being empty string
        """
        self.file_metadata['author'] = author
        self.file_metadata['organisation'] = organisation
        self.file_metadata['project'] = project
        self.file_metadata['subject'] = subject
        self.file_metadata['comment'] = comment
        self.file_metadata['date'] = date
        self.file_metadata['time'] = time


    def __str__(self):
        """representation a mdf_skeleton class data strucutre

        Returns:
        ------------
        string of mdf class ordered as below
        master_channel_name
            channel_name   description
            numpy_array    unit
        """
        if self.fileName is not None:
            output = 'file name : ' + self.fileName + '\n'
        else:
            output = ''
        for m in self.file_metadata.keys():
            if self.file_metadata[m] is not None:
                output += m + ' : ' + str(self.file_metadata[m]) + '\n'
        if not self._pandasframe:
            output += '\nchannels listed by data groups:\n'
            for d in self.masterChannelList.keys():
                if d is not None:
                    output += d + '\n'
                for c in self.masterChannelList[d]:
                    output += '  ' + c + ' : '
                    desc = self.getChannelDesc(c)
                    if desc is not None:
                        try:
                            output += str(desc)
                        except:
                            pass
                    output += '\n    '
                    data = self.getChannelData(c)
                    if data.dtype.kind != 'V': # not byte, impossible to represent
                        output += array_repr(data, \
                            precision=3, suppress_small=True)
                    unit = self.getChannelUnit(c)
                    if unit is not None:
                        output += ' ' + unit + '\n'
            return output
        else:
            pd.set_option('max_rows', 3)
            pd.set_option('expand_frame_repr', True)
            pd.set_option('max_colwidth', 6)
            for master in self.masterGroups:
                output += master
                output += str(self[master])
            return output


    def copy(self):
        """copy a mdf class

        Returns:
        ------------
        mdf_skeleton class instance
            copy of a mdf_skeleton class
        """
        yop = mdf_skeleton()
        yop.multiProc = self.multiProc
        yop.fileName = self.fileName
        yop.masterChannelList = self.masterChannelList
        yop.file_metadata = self.file_metadata
        yop.MDFVersionNumber = self.MDFVersionNumber
        yop.filterChannelNames = self.filterChannelNames
        yop.convert_tables = self.convert_tables
        for channel in list(self.keys()):
            yop[channel] = self[channel]
        return yop

def _open_MDF(fileName):
    """ Opens mdf, make a few checks and returns fid

    Parameters
    -----------
    filename : str
    filename string

    Returns:
    --------
    fid
    file identifier
    """

    try:
        fid = open(fileName, 'rb')
    except IOError:
        raise Exception('Can not find file ' + fileName)
    zipfile = False
    # Check whether file is MDF file -- assumes that every MDF file starts with the letters MDF
    if fid.read(3) not in ('MDF', b'MDF'):
        if is_zipfile(fileName):
            # this is .mfxz file, compressed zip file
            zipfile = True
            fid.close()
            zip_class = ZipFile(fileName, 'r')
            zip_name = zip_class.namelist()[0]  # there should be only one file
            zip_name = zip_class.extract(zip_name)  # locally extracts file
            fid = open(zip_name, 'rb')
            fileName = zip_name
        else:
            raise Exception('file ' + fileName + ' is not an MDF file!')
    return (fid, fileName, zipfile)

def _bits_to_bytes(nBits):
    """ Converts number of bits into number of aligned bytes
    
    Parameters
    -------------
    nBits : int
        number of bits
        
    Returns
    ----------
    number of equivalent bytes
    """
    if nBits <= 8:
        nBytes = 1
    elif nBits <= 16:
        nBytes = 2
    elif nBits <= 32:
        nBytes = 4
    elif nBits <= 64:
        nBytes = 8
    else:
        nBytes = nBits // 8
        if not nBits %8  == 0:
            nBytes += 1
    return nBytes

def _convertName(channelName):
    """ Check if channelName is valid python identifier"""


    if PythonVersion < 3: # python 2
        #channelIdentifier = channelName.encode('utf-8')
        channelIdentifier = _sanitize_identifier(channelName).encode('utf-8')
        #channelIdentifier = channelName.encode('utf-8')

    else: # python 3
        #channelIdentifier = str(channelName)
        channelIdentifier = str(_sanitize_identifier(channelName))
        #channelIdentifier = str(channelName)
    if not channelIdentifier: # all characters of channel ar not compliant to python
        channelIdentifier = ''.join([choice(ascii_letters) for n in range(32)])# generate random name for recarray
    if channelIdentifier in _notAllowedChannelNames:
        channelIdentifier += '_' #  limitation from recarray object attribute
    return channelIdentifier

def _gen_valid_identifier(seq):
    # get an iterator
    itr = iter(seq)
    # pull characters until we get a legal one for first in identifer
    for ch in itr:
        if ch == '_' or ch.isalpha():
            yield ch
            break
        elif ch.isdigit():
            itr = chain(itr,ch)

    # pull remaining characters and yield legal ones for identifier
    for ch in itr:
        if ch == '_' or ch.isalpha() or ch.isdigit():
            yield ch

def _sanitize_identifier(name):
    return ''.join(_gen_valid_identifier(name))