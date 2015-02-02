# -*- coding: utf-8 -*-
""" Measured Data Format file reader module for version 3.x

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

Created on Sun Oct 10 12:57:28 2010

Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>
- Sympy to convert channels with formula

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both python 2.6+ and 3.2+
    
mdf3reader module
--------------------------
"""
from numpy import average, right_shift, bitwise_and, all, diff, max, min, interp
from numpy import asarray, zeros, recarray, array, reshape
from numpy.core.records import fromfile
from math import log, exp
from sys import platform, version_info
from time import strftime, time
from struct import pack, Struct
from io import open # for python 3 and 2 consistency
from mdfinfo3 import info3
PythonVersion=version_info
PythonVersion=PythonVersion[0]

def processDataBlocks(Q, buf, info, dataGroup, channelList, multiProc ):
    """Put raw data from buf to a dict L and processes nested nBit channels
    
    Parameters
    ----------------
    Q : multiprocessing.Queue, optional
        Queue for multiprocessing
    buf : DATA class
        contains raw data
    info : info class
        contains infomation from MDF Blocks
    dataGroup : int
        data group number according to info class
    channelList : list of str, optional
        list of channel names to be processed
    multiProc : bool
        flag to return Queue or dict
        
    Returns
    -----------
    Q : multiprocessing.Queue
        updates Queue containing L dict
    L : dict
        dict of channels
    """
    L = {}

    if channelList is None:
        allChannel=True
    else:
        allChannel=False
    ## Processes Bits, metadata
    for recordID in buf.keys():
        for chan in buf[recordID]['record']:
            channelName=chan.name
            if (allChannel or channelName in channelList) and chan.signalDataType not in (132, 133):
                recordName=buf[recordID]['record'].recordToChannelMatching[channelName] # in case record is used for several channels
                temp = buf[recordID]['data'].__getattribute__( str(recordName)+'_title') # extract channel vector
                    
                if chan.channelType ==1: # master channel 
                    channelName = 'master' + str( dataGroup )

                # Process concatenated bits inside uint8
                if not chan.bitCount//8.0==chan.bitCount/8.0: # if channel data do not use complete bytes
                    mask = int(pow(2, chan.bitCount+1)-1) # masks isBitUnit8
                    if chan.signalDataType in (0,1, 9, 10, 13, 14): # integers
                        temp =  right_shift(temp,  chan.bitOffset)
                        temp =  bitwise_and(temp,  mask )
                        L[channelName] = temp
                    else: # should not happen
                        print('bit count and offset not applied to correct data type')
                        L[channelName] = temp
                else: #data using full bytes
                    L[channelName] = temp

    if multiProc:
        Q.put(L)
    else:
        return L
        
def linearConv(data, conv): # 0 Parametric, Linear: Physical =Integer*P2 + P1
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
        return data # keeps dtype probably more compact than float64
    else:
        return data * conv['P2'] + conv['P1']
def tabInterpConv(data, conv): # 1 Tabular with interpolation
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
    if all(diff(conv['int']) > 0):
        return interp(data, conv['int'], conv['phy'])
    else:
        print(( 'X values for interpolation of channel are not increasing'))
        return data
def tabConv(data, conv): # 2 Tabular
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
    if all(diff(conv['int']) > 0):
        return interp(data, conv['int'], conv['phy'])
    else:
        print(( 'X values for interpolation of channel are not increasing'))
        return data
def polyConv(data, conv): # 6 Polynomial
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
    return (conv['P2'] - conv['P4'] * (data - conv['P5'] - conv['P6'])) / (conv['P3'] * (data - conv['P5'] - conv['P6']) - conv['P1'])
def expConv(data, conv): # 7 Exponential
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
        return exp(((data - conv['P7']) * conv['P6'] - conv['P3']) / conv['P1']) / conv['P2']
    elif conv['P1'] == 0 and conv['P4'] != 0 and conv['P5'] != 0:
        return exp((conv['P3'] / (data - conv['P7']) - conv['P6']) / conv['P4']) / conv['P5']
    else:
        print('Non possible conversion parameters for channel ')
def logConv(data, conv): # 8 Logarithmic
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
    if conv['P4']  == 0 and conv['P1']  != 0 and conv['P2']  != 0:
        return log(((data - conv['P7'] ) * conv['P6']  - conv['P3'] ) / conv['P1'] ) / conv['P2'] 
    elif conv['P1']  == 0 and conv['P4']  != 0 and conv['P5']  != 0:
        return log((conv['P3']  / (data - conv['P7'] ) - conv['P6'] ) / conv['P4'] ) / conv['P5'] 
    else:
        print('Non possible conversion parameters for channel ')
def rationalConv(data, conv): # 9 rational
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
    return(conv['P1']*data * data + conv['P2']*data + conv['P3'])/(conv['P4']*data * data + conv['P5'] * data + conv['P6'])
def formulaConv(data, conv): # 10 Text Formula
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
        X = symbols('X') # variable is X
        formula = conv['textFormula']
        formula=formula[:formula.find('\x00')] # remove trailing text after 0
        formula = formula.replace('pow(', 'power(') # adapt ASAM-MCD2 syntax to sympy
        expr = lambdify(X,formula , 'numpy') # formula to function for evaluation
        return expr(data)
    except:
        print('Please install sympy to convert channel ')
        print('Failed to convert formulae '+conv['textFormula'])
def textRangeTableConv(data, conv): # 12 Text range table
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
        npair=len(conv)
        lower=[conv[pair]['lowerRange'] for pair in range(npair)]
        upper=[conv[pair]['upperRange'] for pair in range(npair)]
        text=[conv[pair]['Textrange'] for pair in range(npair)]
        temp=[]
        for Lindex in range(len(data)):
            value = text[0] # default value
            for pair in range(1, npair):
                if lower[pair] <= data[Lindex] <= upper[pair]:
                    value = text[pair]
                    break
            temp.append(value)
        try:
            temp=asarray(temp) # try to convert to numpy 
        except:
            pass
        return temp
    except:
        print('Failed to convert text to range table')

class recordChannel():
    """ recordChannel class gathers all about channel structure in a record
    
    Attributes
    --------------
    name : str
        Name of channel
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
        bit position of channel value inside byte in case of channel having bit count below 8
    RecordFormat : list of str
        dtype format used for numpy.core.records functions ((name,name_title),str_stype)
    channelType : int
        channel type
    posBeg : int
        start position in number of bit of channel record in complete record
    posEnd : int
        end position in number of bit of channel record in complete record
    
    Methods
    ------------
    __init__(info, dataGroup, channelGroup, channelNumber, recordIDsize)
        constructor
    __str__()
        to print class attributes
    """
    def __init__(self, info, dataGroup, channelGroup, channelNumber, recordIDsize):
        """ recordChannel class constructor
        
        Parameters
        ------------
        info : mdfinfo3.info3 class
        dataGroup : int
            data group number in mdfinfo3.info3 class
        channelGroup : int
            channel group number in mdfinfo3.info3 class
        channelNumber : int
            channel number in mdfinfo3.info3 class
        recordIDsize : int
            size of record ID in Bytes
        """
        self.name=info['CNBlock'][dataGroup][channelGroup][channelNumber]['signalName']
        self.channelNumber=channelNumber
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['signalDataType']
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfBits']
        self.dataFormat=arrayformat3( self.signalDataType, self.bitCount )
        self.CFormat=Struct(datatypeformat3( self.signalDataType, self.bitCount ))
        self.nBytes=self.bitCount // 8+1
        if self.bitCount%8==0:
            self.nBytes-= 1
        recordbitOffset=info['CNBlock'][dataGroup][channelGroup][channelNumber]['numberOfTheFirstBits']
        self.byteOffset=recordbitOffset // 8
        self.bitOffset=recordbitOffset % 8
        self.RecordFormat=((self.name, self.name+'_title'),  self.dataFormat)
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['channelType']
        self.posBeg=recordIDsize+self.byteOffset
        self.posEnd=recordIDsize+self.byteOffset+self.nBytes
    def __str__(self):
        output=str(self.channelNumber) + ' '
        output+=self.name+' '
        output+=str(self.signalDataType)+' '
        output+=str(self.channelType)+' '
        output+=str(self.RecordFormat)+' '
        output+=str(self.bitOffset)+' '
        output+=str(self.byteOffset)
        return output

class record(list):
    """ record class lists recordchannel classes, it is representing a channel group
    
    Attributes
    --------------
    recordLength : int
        length of record corresponding of channel group in Byte
    numberOfRecords : int
        number of records in data block
    recordID : int
        recordID corresponding to channel group
    recordIDsize : int
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
    channelNames : list
        channel names to be stored, useful for low memory consumption but slow
    
    Methods
    ------------
    addChannel(info, channelNumber)
    loadInfo(info)
    readSortedRecord(fid, pointer, channelList=None)
    readUnsortedRecord(buf, channelList=None)
    """
    def __init__(self, dataGroup, channelGroup):
        self.recordLength=0
        self.numberOfRecords=0
        self.recordID=0
        self.recordIDsize=0
        self.dataGroup=dataGroup
        self.channelGroup=channelGroup
        self.numpyDataRecordFormat=[]
        self.dataRecordName=[]
        self.master={}
        self.recordToChannelMatching={}
        self.channelNames=[]
    def addChannel(self, info, channelNumber):
        """ add a channel in class
        
        Parameters
        ----------------
        info : mdfinfo3.info3 class
        channelNumber : int
            channel number in mdfinfo3.info3 class
        """
        self.append(recordChannel(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize))
        self.channelNames.append(self[-1].name)
    def loadInfo(self, info):
        """ gathers records related from info class
        
        Parameters
        ----------------
        info : mdfinfo3.info3 class
        """
        self.recordIDsize=info['DGBlock'][self.dataGroup]['numberOfRecordIDs']
        self.recordID=info['CGBlock'][self.dataGroup][self.channelGroup]['recordID']
        if not self.recordIDsize==0: # record ID existing
            self.dataRecordName.append('RecordID'+str(self.channelGroup))
            format=(self.dataRecordName[-1], self.dataRecordName[-1]+'_title')
            self.numpyDataRecordFormat.append( ( format, 'uint8' ) )
        self.recordLength=info['CGBlock'][self.dataGroup][self.channelGroup]['dataRecordSize']
        self.numberOfRecords=info['CGBlock'][self.dataGroup][self.channelGroup]['numberOfRecords']
        for channelNumber in list(info['CNBlock'][self.dataGroup][self.channelGroup].keys()):
            channel=recordChannel(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize)
            if channel.channelType==1: # master channel found
                self.master['name']=channel.name
                self.master['number']=channelNumber
            self.append(channel)
            self.channelNames.append(channel.name)
            if len(self)>1 and channel.byteOffset==self[-2].byteOffset: # several channels in one byte, ubit1 or ubit2
                self.recordToChannelMatching[channel.name]=self.recordToChannelMatching[self[-2].name]
            else: # adding bytes
                self.recordToChannelMatching[channel.name]=channel.name
                self.numpyDataRecordFormat.append(channel.RecordFormat)
                self.dataRecordName.append(channel.name)
        if self.recordIDsize==2: # second record ID at end of record
            self.dataRecordName.append('RecordID'+str(self.channelGroup)+'_2')
            format=(self.dataRecordName[-1], self.dataRecordName[-1]+'_title')
            self.numpyDataRecordFormat.append( ( format, 'uint8' ) )

    def readSortedRecord(self, fid, pointer, channelList=None):
        """ reads record, only one channel group per datagroup
        Parameters
        ----------------
        fid : float
            file identifier
        pointer
            position in file of data block beginning
        channelList : list of str, optional
            list of channel to read
            
        Returns
        -----------
        rec : numpy recarray
            contains a matrix of raw data in a recarray (attributes corresponding to channel name)
            
        Notes
        --------
        If channelList is None, read data using numpy.core.records.fromfile that is rather quick.
        However, in case of large file, you can use channelList to load only interesting channels or 
        only one channel on demand, but be aware it might be much slower.
        """
        fid.seek(pointer)
        if channelList is None: # reads all, quickest but memory consuming
            return fromfile( fid, dtype = self.numpyDataRecordFormat, shape = self.numberOfRecords, names=self.dataRecordName)
        else: # reads only some channels from a sorted data block
            # memory efficient but takes time
            if len(list(set(channelList)&set(self.channelNames)))>0: # are channelList in this dataGroup
                # check if master channel is in the list
                if not self.master['name'] in channelList:
                    channelList.append(self.master['name']) # adds master channel
                rec={}
                recChan=[]
                numpyDataRecordFormat=[]
                for channel in channelList: # initialise data structure
                    rec[channel]=0
                for channel in self: # list of recordChannels from channelList
                    if channel.name in channelList:
                        recChan.append(channel)
                        numpyDataRecordFormat.append(channel.RecordFormat)
                rec=zeros((self.numberOfRecords, ), dtype=numpyDataRecordFormat)
                recordLength=self.recordIDsize+self.recordLength
                for r in range(self.numberOfRecords): # for each record,
                    buf=fid.read(recordLength)
                    for channel in recChan:
                        rec[channel.name][r]=channel.CFormat.unpack(buf[channel.posBeg:channel.posEnd])[0]
                return rec.view(recarray)
            
    def readUnsortedRecord(self, buf, channelList=None):
        """ Not implemented yet, no reference files available to test it
        """
        pass

class DATA(dict):
    """ DATA class is organizing record classes itself made of recordchannel.
    This class inherits from dict. Keys are corresponding to channel group recordID
    A DATA class corresponds to a data block, a dict of record classes (one per channel group) 
    Each record class contains a list of recordchannel class representing the structure of channel record.
    
    Attributes
    --------------
    fid : io.open
        file identifier
    pointerToData : int
        position of Data block in mdf file
    
    Methods
    ------------
    addRecord(record)
        Adds a new record in DATA class dict
    read(channelList, zip=None)
        Reads data block
    loadSorted(record, zip=None, nameList=None)
        Reads sorted data block from record definition
    load(nameList=None)
        Reads unsorted data block, not yet implemented
    """
    def __init__(self, fid, pointer):
        self.fid=fid
        self.pointerToData=pointer
    def addRecord(self, record):
        """Adds a new record in DATA class dict
        
        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        """
        self[record.recordID]={}
        self[record.recordID]['record']=record

    def read(self, channelList, zip=None):
        """Reads data block
        
        Parameters
        ----------------
        channelList : list of str, optional
            list of channel names
        zip : bool, optional
            flag to track if data block is compressed
        """
        if len(self)==1: #sorted dataGroup
            recordID=list(self.keys())[0]
            self[recordID]['data']=self.loadSorted( self[recordID]['record'], zip=None, nameList=channelList)
        else: # unsorted DataGroup
            self.load( nameList=channelList)
            
    def loadSorted(self, record, zip=None, nameList=None): # reads sorted data
        """Reads sorted data block from record definition
        
        Parameters
        ----------------
        record class
            channel group definition listing record channel classes
        zip : bool, optional
            flag to track if data block is compressed
        channelList : list of str, optional
            list of channel names
            
        Returns
        -----------
        numpy recarray of data
        """
        return record.readSortedRecord(self.fid, self.pointerToData, nameList)
    
    def load(self, nameList=None):
        """ not yet implemented
        """ 
        return None

class mdf3(dict):
    """ mdf file version 3.0 to 3.3 class
    
    Attributes
    --------------
    fileName : str
        file name
    VersionNumber : int
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
    author : str
    organisation : str
    project : str
    subject : str
    comment : str
    time : str
    date : str
    
    Methods
    ------------
    read3( fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True)
        Reads mdf 3.x file data and stores it in dict
    getChannelData3(channelName)
        Returns channel numpy array
    convertChannel3(channelName)
        converts specific channel from raw to physical data according to CCBlock information
    convertAllChannel3()
        Converts all channels from raw data to converted data according to CCBlock information
    write3(fileName=None)
        Writes simple mdf 3.3 file
    """
    
    def __init__( self, fileName=None, info=None, multiProc=False,  channelList=None, convertAfterRead=True, filterChannelNames=False):
        self.masterChannelList = {}
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        self.author=''
        self.organisation=''
        self.project=''
        self.subject=''
        self.comment=''
        self.time=''
        self.date=''
        self.VersionNumber=300
        self.filterChannelNames=False
        # clears class from previous reading and avoid to mess up
        self.clear()
        if fileName is None and info is not None:
            self.fileName = info.fileName
            self.read3(self.fileName, info, multiProc, channelList, convertAfterRead)
        elif fileName is not None:
            self.fileName = fileName
            self.read3(self.fileName, info, multiProc, channelList, convertAfterRead)

    def read3( self, fileName=None, info=None, multiProc=False, channelList=None, convertAfterRead=True):
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
            If you use channelList, reading might be much slower but it will save you memory. Can be used to read big files
        
        convertAfterRead : bool, optional
            flag to convert channel after read, True by default
            If you use convertAfterRead by setting it to false, all data from channels will be kept raw, no conversion applied.
            If many float are stored in file, you can gain from 3 to 4 times memory footprint
            To calculate value from channel, you can then use method .getChannelData()
        """
        self.multiProc = multiProc
        if platform == 'win32':
            self.multiProc = False # no multiprocessing for windows platform
        try:
            from multiprocessing import Queue, Process
        except:
            print('No multiprocessing module found')
            self.multiProc = False
        
        if self.fileName is None:
            self.fileName = info.fileName
        else:
            self.fileName = fileName
            
        # inttime = time.clock()
        ## Read information block from file
        if info is None:
            info = info3(self.fileName,  None,  self.filterChannelNames)
            
        # reads metadata
        self.author=info['HDBlock']['Author']
        self.organisation=info['HDBlock']['Organization']
        self.project=info['HDBlock']['ProjectName']
        self.subject=info['HDBlock']['Subject']
        try:
            self.comment=info['HDBlock']['TXBlock']['Text']
        except:
            pass
        self.time=info['HDBlock']['Time']
        self.date=info['HDBlock']['Date']
        # converts date to be compatible with ISO8601
        day, month, year=self.date.split(':')
        self.date=year+'-'+month+'-'+day

        try:
            fid = open(self.fileName, 'rb')
        except IOError:
            print('Can not find file'+self.fileName)
            raise

        # Look for the biggest group to process first, to reduce processing time when mutiprocessed
        dataGroupList = dict.fromkeys(list(range( info['HDBlock']['numberOfDataGroups'])))
        for dataGroup in list(dataGroupList.keys()):
            dataGroupList[dataGroup] = info['CGBlock'][dataGroup][0]['numberOfRecords']
        sortedDataGroup = sorted(dataGroupList, key=dataGroupList.__getitem__, reverse=True)
        
        if self.multiProc:
            # prepare multiprocessing of dataGroups
            proc = []
            Q = Queue()
        L = {}
        masterDataGroup={}  # datagroup name correspondence with its master channel
        ## Read data from file
        for dataGroup in sortedDataGroup:
            if info['DGBlock'][dataGroup]['numberOfChannelGroups']>0: # data exists
                #Pointer to data block
                pointerToData = info['DGBlock'][dataGroup]['pointerToDataRecords']
                buf=DATA(fid, pointerToData)

                for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                    temp=record(dataGroup, channelGroup) # create record class
                    temp.loadInfo(info) # load all info related to record

                    if temp.numberOfRecords != 0:  # continue if there are at least some records
                        buf.addRecord(temp)
                        for channel in range(info['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):
                            if info['CNBlock'][dataGroup][channelGroup][channel]['channelType'] ==1:
                                masterDataGroup[dataGroup]=info['CNBlock'][dataGroup][channelGroup][channel]['signalName']

                buf.read(channelList)

                if self.multiProc:
                    proc.append(Process(target=processDataBlocks,
                        args=(Q, buf, info, dataGroup, channelList, self.multiProc)))
                    proc[-1].start()
                else:  # for debugging purpose, can switch off multiprocessing
                    L.update(processDataBlocks( None, buf, info, dataGroup,  channelList, self.multiProc))
                del buf # free memory

        fid.close()  # close file

        if self.multiProc:
            for p in proc:
                L.update(Q.get())  # concatenate results of processes in dict
            for p in proc:
                p.join()
            del Q # free memory

        # After all processing of channels,
        # prepare final class data with all its keys
        for dataGroup in range(info['HDBlock']['numberOfDataGroups']):
            for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                for channel in range(info['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):
                    numberOfRecords = info['CGBlock'][dataGroup][channelGroup]['numberOfRecords']
                    if numberOfRecords != 0 :
                        channelName = info['CNBlock'][dataGroup][channelGroup][channel]['signalName']
                        if info['CNBlock'][dataGroup][channelGroup][channel]['channelType'] == 1:  # time channel
                            channelName = 'master' + str(dataGroup)
                        if channelName in L and len(L[channelName]) != 0:
                            if ('master' + str(dataGroup)) not in list(self.masterChannelList.keys()):
                                self.masterChannelList['master' + str(dataGroup)] = []
                            self.masterChannelList['master' + str(dataGroup)].append(channelName)
                            self[channelName] = {}
                            self[channelName]['master'] = 'master' + str(dataGroup) # master channel of channel
                            self[channelName]['unit'] = info['CCBlock'][dataGroup][channelGroup][channel]['physicalUnit']
                            self[channelName]['description'] = info['CNBlock'][dataGroup][channelGroup][channel]['signalDescription']
                            self[channelName]['data'] = L[channelName]
                            L.pop(channelName, None) # free memory
                            convType = info['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier']
                            if convType in (0, 1, 2, 6, 7, 8, 9, 10, 11, 12): # needs conversion
                                self[channelName]['conversion'] = {}
                                self[channelName]['conversion']['type'] = convType
                                self[channelName]['conversion']['parameters'] = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']
                            if convType == 0 and (self[channelName]['conversion']['parameters']['P2'] == 1.0 and self[channelName]['conversion']['parameters']['P1'] in (0.0, -0.0)):
                                self[channelName].pop('conversion')
        
        if convertAfterRead: 
            self.convertAllChannel3()
        #print( 'Finished in ' + str( time.clock() - inttime ) )
    
    def getChannelData3(self, channelName):
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
            return self.convert3(channelName)
        else:
            raise('Channel not in dictionary')
    
    def convert3(self, channelName):
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
        if 'conversion' in self[channelName]: # there is conversion property
            if self[channelName]['conversion']['type'] == 0:
                return linearConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 1:
                return tabInterpConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 2:
                return tabConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 6:
                return polyConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 7:
                return expConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 8:
                return logConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 9:
                return rationalConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 10:
                return formulaConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            elif self[channelName]['conversion']['type'] == 12:
                return textRangeTableConv(self[channelName]['data'], self[channelName]['conversion']['parameters'])
            else:
                return self[channelName]['data']
        else:
            return self[channelName]['data']
            
    def convertChannel3(self, channelName):
        """converts specific channel from raw to physical data according to CCBlock information
        
        Parameters
        ----------------
        channelName : str
            Name of channel
        """
        if 'conversion' in self[channelName]:
            self[channelName]['data'] = self.convert3(channelName)
            self[channelName].pop('conversion')
                
    def convertAllChannel3(self):
        """Converts all channels from raw data to converted data according to CCBlock information
        Converted data will take more memory.
        """
        for channel in self:
            self.convertChannel3(channel)

    def write3(self, fileName=None):
        """Writes simple mdf 3.3 file
        
        Parameters
        ----------------
        fileName : str, optional
            Name of file
            If file name is not input, written file name will be the one read with appended '_new' string before extension
        
        Notes
        --------
        All channels will be converted to physical data, so size might be bigger than original file
        """
        LINK = 'I'
        #CHAR = 'c'
        REAL = 'd'
        BOOL = 'h'
        #UINT8 = 'B'
        #BYTE =  'B'
        INT16 = 'h'
        UINT16 = 'H'
        UINT32 = 'I'
        #INT32 = 'i'
        UINT64 = 'Q'
        #INT64 = 'q'
        
        # put master channel in first position for each datagroup if not already the case
        for master in list(self.masterChannelList.keys()):
            masterList = self.masterChannelList[master]
            masterList.sort()  # alphabetically sort the channel names
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
                    temp = value+'\0'*(size-len(value))
                temp += '\0'
            if self.VersionNumber<400:
                if PythonVersion>=3:
                    temp=temp.encode('latin1', 'replace')
                f.write(pack('<'+str(len(temp))+'s', temp))
            else:
                temp=temp.encode('latin1', 'replace')
                f.write(pack('<'+str(len(temp))+'s', temp))

        # write pointer of block and come back to current stream position
        def writePointer(f, pointer, value):
            currentPosition = f.tell()
            f.seek(pointer)
            f.write(pack(LINK, value))
            f.seek(currentPosition)
        
        # Starts first to write ID and header
        fid = open(fileName, 'wb')  # buffering should automatically be set
        writeChar(fid, 'MDF     ')
        writeChar(fid, '3.30    ')
        writeChar(fid, 'MDFreadr')
        fid.write(pack(UINT16, 0))  # little endian
        fid.write(pack(UINT16, 0))  # floating format
        fid.write(pack(UINT16, 330))  # version 3.0
        fid.write(pack(UINT16, 28591))  # code page ISO2859-1 latin 1 western europe
        writeChar(fid, '\0'*32)  # reserved
        
        # Header Block
        writeChar(fid, 'HD')
        fid.write(pack(UINT16, 208))  # block size
        pointers['HD'] = {}
        pointers['HD']['DG'] = fid.tell()
        fid.write(pack(LINK, 272))  # first Data block pointer
        pointers['HD']['TX'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to TX Block file comment
        pointers['HD']['PR'] = fid.tell()
        fid.write(pack(LINK, 0))  # pointer to PR Block
        ndataGroup = len(self.masterChannelList)
        fid.write(pack(UINT16, ndataGroup))  # number of data groups
        writeChar(fid, strftime("%d:%m:%Y"))  # date
        writeChar(fid, strftime("%H:%M:%S"))  # time
        if self.author is not None:
            writeChar(fid, self.author,  size=31)  # Author
        else:
            writeChar(fid, ' ',  size=31)  # Author
        if self.organisation is not None:
            writeChar(fid, self.organisation,  size=31)  # Organization
        else:
            writeChar(fid, ' ',  size=31)  
        if self.project is not None:
            writeChar(fid, self.project,  size=31)  # Project
        else:
            writeChar(fid, ' ',  size=31)  
        if self.subject is not None:
            writeChar(fid, self.subject,  size=31)  # Subject
        else:
            writeChar(fid, ' ',  size=31)  
        fid.write(pack(UINT64, int(time()*1000000000)))  # Time Stamp
        fid.write(pack(INT16, 1))  # UTC time offset
        fid.write(pack(UINT16, 0))  # Time quality
        writeChar(fid, 'Local PC Reference Time         ')  # Timer identification
        
        # write DG block
        pointers['DG'] = {}
        pointers['CG'] = {}
        pointers['CN'] = {}
        
        for dataGroup in range(ndataGroup):
            # writes dataGroup Block
            pointers['DG'][dataGroup] = {}
            if 0 < dataGroup:  # not possible for first DG
                # previous datagroup pointer to this new datagroup
                writePointer(fid, pointers['DG'][dataGroup-1]['nextDG'], fid.tell())
            else:
                # first datagroup pointer in header block
                writePointer(fid, pointers['HD']['DG'], fid.tell())
            writeChar(fid, 'DG')
            fid.write(pack(UINT16, 28))  # DG block size
            pointers['DG'][dataGroup]['nextDG'] = fid.tell()
            # pointer to next DataGroup, 0 by default until it is known when creating new datagroup
            fid.write(pack(LINK, 0))
            pointers['DG'][dataGroup]['CG'] = fid.tell()
            fid.write(pack(LINK, 0))  # pointer to channel group, 0 until CG created
            fid.write(pack(LINK, 0))  # pointer to trigger block, not used
            pointers['DG'][dataGroup]['data'] = fid.tell()
            fid.write(pack(LINK, 0))  # pointer to data block
            fid.write(pack(UINT16, 1))  # number of channel group, 1 because sorted data
            fid.write(pack(UINT16, 0))  # number of record IDs
            writeChar(fid, '\0'*32)  # reserved
            
            # sorted data so only one channel group
            pointers['CG'][dataGroup] = {}
            # write first CG pointer in datagroup
            writePointer(fid, pointers['DG'][dataGroup]['CG'], fid.tell())
            writeChar(fid, 'CG')
            fid.write(pack(UINT16, 30))  # CG block size
            fid.write(pack(LINK, 0))  # pointer to next Channel Group but no other, one CG per DG
            pointers['CG'][dataGroup]['firstCN'] = fid.tell()
            fid.write(pack(LINK, 0))  # pointer to first channel block
            pointers['CG'][dataGroup]['TX'] = fid.tell()
            fid.write(pack(LINK, 0))  # pointer to TX block
            fid.write(pack(UINT16, 0))  # No record ID no need for sorted data
            masterChannel = list(self.masterChannelList.keys())[dataGroup]
            numChannels = len(self.masterChannelList[masterChannel])
            fid.write(pack(UINT16, numChannels))  # Number of channels
            pointers['CG'][dataGroup]['dataRecordSize'] = fid.tell()
            fid.write(pack(UINT16, 0))  # Size of data record
            nRecords = len(self[masterChannel]['data'])
            fid.write(pack(UINT32, nRecords))  # Number of records
            fid.write(pack(LINK, 0))  # pointer to sample reduction block, not used
            sampling = 0
            if self[masterChannel]['data'] is not None and len(self[masterChannel]['data'])>0 and self[masterChannel]['data'].dtype.kind not in ['S', 'U']:
                sampling = average(diff(self[masterChannel]['data']))            
            
            # Channel blocks writing
            pointers['CN'][dataGroup] = {}
            dataList = ()
            dataTypeList = ''
            recordNumberOfBits = 0
            preceedingChannel = None
            bitOffset=0
            writePointer(fid, pointers['CG'][dataGroup]['firstCN'], fid.tell())  # first channel bock pointer from CG
            for channel in self.masterChannelList[masterChannel]:
                pointers['CN'][dataGroup][channel] = {}
                pointers['CN'][dataGroup][channel]['beginCN'] = fid.tell()
                writeChar(fid, 'CN')
                fid.write(pack(UINT16, 228))  # CN block size
                pointers['CN'][dataGroup][channel]['nextCN'] = fid.tell()
                if preceedingChannel is not None:  # not possible for first CN
                    writePointer(fid, pointers['CN'][dataGroup][preceedingChannel]['nextCN'], pointers['CN'][dataGroup][channel]['beginCN'])  # pointer in previous cN
                preceedingChannel = channel
                fid.write(pack(LINK, 0))  # pointer to next channel block, 0 as not yet known
                pointers['CN'][dataGroup][channel]['CC'] = fid.tell()
                fid.write(pack(LINK, 0))  # pointer to conversion block
                fid.write(pack(LINK, 0))  # pointer to source depending block
                fid.write(pack(LINK, 0))  # pointer to dependency block
                pointers['CN'][dataGroup][channel]['TX'] = fid.tell()
                fid.write(pack(LINK, 0))  # pointer to comment TX, no comment
                # check if master channel
                if channel not in list(self.masterChannelList.keys()):
                    fid.write(pack(UINT16, 0))  # data channel
                else:
                    fid.write(pack(UINT16, 1))  # master channel
                # make channel name in 32 bytes
                writeChar(fid, channel, size=31)  # channel name
                # channel description
                desc = self[channel]['description']
                writeChar(fid, desc, size=127)  # channel description
                fid.write(pack(UINT16, bitOffset))  # bit position
                data = self[channel]['data']  # channel data
                temp=data
                if PythonVersion>=3 and data.dtype.kind in ['S', 'U']:
                    temp=temp.encode('latin1', 'replace')
                dataList = dataList+(temp, )

                if data.dtype in ('float64', 'int64', 'uint64'):
                    numberOfBits = 64
                elif data.dtype in ('float32', 'int32', 'uint32'):
                    numberOfBits = 32
                elif data.dtype in ('uint16', 'int16'):
                    numberOfBits = 16
                elif data.dtype in ('uint8', 'int8'):
                    numberOfBits = 8
                else:
                    numberOfBits = 8  # if string, not considered
                recordNumberOfBits += numberOfBits
                fid.write(pack(UINT16, numberOfBits))  # Number of bits
                bitOffset += numberOfBits
                if data.dtype == 'float64':
                    dataType = 3
                elif data.dtype in ('uint8', 'uint16', 'uint32', 'uint64'):
                    dataType = 0
                elif data.dtype in ('int8', 'int16', 'int32', 'int64'):
                    dataType = 1
                elif data.dtype == 'float32':
                    dataType = 2
                elif data.dtype.kind in ['S', 'U']: 
                    dataType = 7
                else:
                    print('Not recognized dtype')
                    raise
                if not data.dtype.kind in ['S', 'U']: 
                    dataTypeList += data.dtype.char
                else:
                    dataTypeList += str(data.dtype.itemsize)+'s'
                fid.write(pack(UINT16, dataType))  # Signal data type
                if not data.dtype.kind in ['S', 'U']:
                    fid.write(pack(BOOL, 1))  # Value range valid
                    if len(data)>0:
                        maximum=max(data)
                        minimum=min(data)
                    else:
                        maximum=0
                        minimum=0
                    fid.write(pack(REAL, minimum))  # Min value
                    fid.write(pack(REAL, maximum))  # Max value
                else:
                    fid.write(pack(BOOL, 0))  # No value range valid
                    fid.write(pack(REAL, 0))  # Min value
                    fid.write(pack(REAL, 0))  # Max value
                fid.write(pack(REAL, sampling))  # Sampling rate
                pointers['CN'][dataGroup][channel]['longChannelName'] = fid.tell()
                fid.write(pack(LINK, 0))  # pointer to long channel name
                fid.write(pack(LINK, 0))  # pointer to channel display name
                fid.write(pack(UINT16, 0))  # No Byte offset
                
                # TXblock for long channel name
                writePointer(fid, pointers['CN'][dataGroup][channel]['longChannelName'], fid.tell())
                writeChar(fid, 'TX')
                fid.write(pack(UINT16, len(channel)+4+1))  # TX block size
                writeChar(fid, channel+'\0')  # channel name that can be long, should ends by 0 (NULL)
                
                # Conversion blocks writing
                writePointer(fid, pointers['CN'][dataGroup][channel]['CC'], fid.tell())
                writeChar(fid, 'CC')
                fid.write(pack(UINT16, 46))  # CC block size
                if not data.dtype.kind in ['S', 'U']:
                    fid.write(pack(BOOL, 1))  # Value range valid
                    fid.write(pack(REAL, minimum))  # Min value
                    fid.write(pack(REAL, maximum))  # Max value
                else:
                    fid.write(pack(BOOL, 0))  # No value range valid
                    fid.write(pack(REAL, 0))  # Min value
                    fid.write(pack(REAL, 0))  # Max value
                writeChar(fid, self[channel]['unit'], size=19)  # channel description
                fid.write(pack(UINT16, 65535))  # conversion already done during reading
                fid.write(pack(UINT16, 0))  # additional size information, not necessary for 65535 conversion type ?
            # number of channels in CG
            currentPosition = fid.tell()
            fid.seek(pointers['CG'][dataGroup]['dataRecordSize'])
            fid.write(pack(UINT16, int(recordNumberOfBits/8)))  # Size of data record
            fid.seek(currentPosition)
            
            # data writing
            # write data pointer in datagroup
            writePointer(fid, pointers['DG'][dataGroup]['data'], fid.tell())
            records = array(dataList, object).T
            records = reshape(records, (1, len(self.masterChannelList[masterChannel])*nRecords), order='C')[0]  # flatten the matrix
            fid.write(pack('<'+dataTypeList*nRecords, *records))  # dumps data vector from numpy
                
        #print(pointers)
        fid.close()
        
    
def datatypeformat3(signalDataType, numberOfBits):
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
    if signalDataType in (0, 9, 10, 11):  # unsigned
        if numberOfBits <= 8:
            dataType = 'B'
        elif numberOfBits <= 16:
            dataType = 'H'
        elif numberOfBits <= 32:
            dataType = 'I'
        else:
            print(('Unsupported number of bits for unsigned int ' + str(signalDataType)))

    elif signalDataType == 1:  # signed int
        if numberOfBits <= 8:
            dataType = 'b'
        elif numberOfBits <= 16:
            dataType = 'h'
        elif numberOfBits <= 32:
            dataType = 'i'
        else:
            print(('Unsupported number of bits for signed int ' + str(signalDataType)))

    elif signalDataType in (2, 3):  # floating point
        if numberOfBits == 32:
            dataType = 'f'
        elif numberOfBits == 64:
            dataType = 'd'
        else:
            print(('Unsupported number of bit for floating point ' + str(signalDataType)))

    elif signalDataType == 7:  # string
        dataType = 's'
    elif signalDataType == 8:  # array of bytes
        dataType = 'B'  # take unit8 as default ?
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits))
    return dataType


def arrayformat3(signalDataType, numberOfBits):
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

    if signalDataType in (0, 9, 10, 11):  # unsigned
        if numberOfBits <= 8:
            dataType = 'uint8'
        elif numberOfBits == 16:
            dataType = 'uint16'
        elif numberOfBits == 32:
            dataType = 'uint32'
        elif numberOfBits == 1:
            dataType = 'uint8'  # not directly processed
        elif numberOfBits == 2:
            dataType = 'uint8'  # not directly processed
        else:
            print(('Unsupported number of bits for unsigned int ' + str(signalDataType)))

    elif signalDataType == 1:  # signed int
        if numberOfBits <= 8:
            dataType = 'int8'
        elif numberOfBits == 16:
            dataType = 'int16'
        elif numberOfBits == 32:
            dataType = 'int32'
        else:
            print(('Unsupported number of bits for signed int ' + str(signalDataType)))

    elif signalDataType in (2, 3):  # floating point
        if numberOfBits == 32:
            dataType = 'float32'
        elif numberOfBits == 64:
            dataType = 'float64'
        else:
            print(('Unsupported number of bit for floating point ' + str(signalDataType)))

    elif signalDataType == 7:  # string
        dataType = 'str'  # not directly processed
    elif signalDataType == 8:  # array of bytes
        dataType = 'buffer'  # not directly processed
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits))
    return dataType
