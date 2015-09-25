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

    def __init__(self, fileName=None, channelList=None, convertAfterRead=True, filterChannelNames=False):
        """ mdf class constructor.
        When mdf class is constructed, constructor can be called to directly reads file

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
    
    def add_channel(self, channel_name, data, master_channel, unit='', description='', conversion=None):
        self[channel_name] = {}
        self[channel_name]['data'] = data
        self[channel_name]['unit'] = unit
        self[channel_name]['description'] = description
        self[channel_name]['master'] = master_channel
        if master_channel not in self.masterChannelList.keys():
            self.masterChannelList[master_channel] = []
        self.masterChannelList[master_channel].append(channel_name)
        if conversion is not None and 'conversion' in conversion:
            self[channel_name]['conversion'] = {}
            self[channel_name]['conversion']['type'] = conversion['conversionFormulaIdentifier']
            self[channel_name]['conversion']['parameters'] = conversion['conversion']
            if conversion['conversionFormulaIdentifier'] == 0 and \
                    (self[channel_name]['conversion']['parameters']['P2'] == 1.0 and \
                    self[channel_name]['conversion']['parameters']['P1'] in (0.0, -0.0)):
                self[channel_name].pop('conversion')
    
    def remove_channel(self, channel_name):
        self.masterChannelList[self[channel_name]['master']].pop(channel_name)
        del self[channel_name]
