# -*- coding: utf-8 -*-
""" mdf module describing basic mdf structure and methods

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
from numpy import array_repr, set_printoptions
set_printoptions(threshold=100, edgeitems=1)

class mdf_skeleton(dict):

    """ mdf class

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
    """

    def __init__(self, fileName=None, channelList=None, convertAfterRead=True, filterChannelNames=False):
        """ mdf class constructor.

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
        # clears class from previous reading and avoid to mess up
        self.clear()
        self.fileName = fileName
        if fileName is not None:
            self.read(fileName, channelList=channelList, convertAfterRead=convertAfterRead, filterChannelNames=filterChannelNames)
    
    def add_channel(self, channel_name, data, master_channel, master_type=1, unit='', description='', conversion=None):
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
        """
        self[channel_name] = {}
        self[channel_name]['data'] = data
        self[channel_name]['unit'] = unit
        self[channel_name]['description'] = description
        self[channel_name]['master'] = master_channel
        if master_channel not in self.masterChannelList.keys():
            self.masterChannelList[master_channel] = []
        self.masterChannelList[master_channel].append(channel_name)
        if self.MDFVersionNumber < 400: #  mdf3
                self[channel_name]['masterType'] = 1
        else: #  mdf4
            self[channel_name]['masterType'] = master_type
        if conversion is not None:
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
        self.masterChannelList[self[channel_name]['master']].pop(channel_name)
        return self.pop(channel_name)
    
    def __repr__(self):
        output = 'file name : ' + self.fileName + '\n'
        for m in self.file_metadata.keys():
            output += m + ' : ' + str(self.file_metadata[m]) + '\n'
        output += '\nchannels listed by data groups:\n'
        for d in self.masterChannelList.keys():
            if d is not None:
                output += d + '\n'
            for c in self.masterChannelList[d]:
                output += '  ' + c + ' : ' + str(self[c]['description']) + '\n'
                output += '    ' + array_repr(self[c]['data'], precision=3, suppress_small=True) \
                    + ' ' + self[c]['unit'] + '\n'
        return output
