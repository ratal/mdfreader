# -*- coding: utf-8 -*-
""" Measured Data Format blocks parser for version 3.x

Created on Thu Dec 9 12:57:28 2014

Platform and python version
----------------------------------------
With Unix and Windows for python 2.6+ and 3.2+

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

Dependencies
-------------------
- Python >2.6, >3.2 <http://www.python.org>
- Numpy >1.6 <http://numpy.scipy.org>

Attributes
--------------
PythonVersion : float
    Python version currently running, needed for compatibility of both python 2.6+ and 3.2+
    
mdfinfo3 module
--------------------------
"""
from sys import version_info
PythonVersion=version_info
PythonVersion=PythonVersion[0]
from struct import calcsize, unpack

class info3(dict):
    """ mdf file info class version 3.x
    MDFINFO is a class information about an MDF (Measure Data Format) file
    Based on following specification http://powertrainnvh.com/nvh/MDFspecificationv03.pdf
    
    Attributes
    --------------
    filterChannelNames : bool, optional
        flag to filter long channel names including module names separated by a '.'
    fileName : str
        name of file
    
    Notes
    --------
    mdfinfo(FILENAME) contains a dict of structures, for
    each data group, containing key information about all channels in each
    group. FILENAME is a string that specifies the name of the MDF file. 
    General dictionary structure is the following
    
    - mdfinfo['HDBlock'] header block
    - mdfinfo['DGBlock'][dataGroup] Data Group block
    - mdfinfo['CGBlock'][dataGroup][channelGroup] Channel Group block
    - mdfinfo['CNBlock'][dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    - mdfinfo['CCBlock'][dataGroup][channelGroup][channel] Channel conversion information
    """
    def __init__(self, fileName=None, fid=None,  filterChannelNames = False):
        """ info3 class constructor
        
        Parameters
        ----------------
        fileName : str, optional
            name of file
        fid : float, optional
            file identifier
        filterChannelNames : bool, optional
            flag to filter long channel names including module names separated by a '.'

        Notes
        --------
        If fileName is given it will read file blocks directly by calling method readinfo3
        """
        self['IDBlock'] = {} # Identifier Block
        self['HDBlock'] = {} # Header Block
        self['DGBlock']= {} # Data Group Block
        self['CGBlock'] = {} # Channel Group Block
        self['CNBlock'] = {}# Channel Block
        self['CCBlock'] = {} # Conversion block
        self.filterChannelNames = filterChannelNames
        self.fileName = fileName
        if fileName != None and fid==None:
            try:
                fid = open( self.fileName, 'rb' )
            except IOError:
                print('Can not find file'+self.fileName)
                raise
            self.readinfo3( fid )
        elif fileName == None and fid!=None:
            self.readinfo3(fid)

    def readinfo3(self, fid):
        """ read all file blocks except data
        
        Parameters
        ----------------
        fid : float
            file identifier
        """
        # reads IDBlock
        fid.seek(24)
        self['IDBlock']['ByteOrder']=unpack( '<H', fid.read( 2 ) )
        self['IDBlock']['FloatingPointFormat']=unpack( '<H', fid.read( 2 ) )
        self['IDBlock']['id_ver']=unpack( '<H', fid.read( 2 ) )
        self['IDBlock']['id_ver']=self['IDBlock']['id_ver'][0]
        self['IDBlock']['CodePageNumber']=unpack( '<H', fid.read( 2 ) )
        self['IDBlock']['CodePageNumber']=self['IDBlock']['CodePageNumber'][0]
        
        ### Read header block (HDBlock) information
        # Read Header block info into structure, HD pointer at 64 as mentioned in specification
        self['HDBlock'] = self.mdfblockread3( self.blockformats3( 'HDFormat' ,  self['IDBlock']['id_ver']), fid, 64 )

        ### Read text block (TXBlock) information
        self['HDBlock']['TXBlock'] = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, self['HDBlock']['pointerToTXBlock'] )

        ### Read text block (PRBlock) information
        self['HDBlock']['PRBlock'] = self.mdfblockread3( self.blockformats3( 'PRFormat' ), fid, self['HDBlock']['pointerToPRBlock'] )

        ### Read Data Group blocks (DGBlock) information
        # Get pointer to first Data Group block
        DGpointer = self['HDBlock']['pointerToFirstDGBlock']
        for dataGroup in range( self['HDBlock']['numberOfDataGroups'] ):

            # Read data Data Group block info into structure
            self['DGBlock'][dataGroup] = self.mdfblockread3( self.blockformats3( 'DGFormat' ), fid, DGpointer )
            # Get pointer to next Data Group block
            DGpointer = self['DGBlock'][dataGroup]['pointerToNextDGBlock']

            ### Read Channel Group block (CGBlock) information - offset set already

            # Read data Channel Group block info into structure
            CGpointer = self['DGBlock'][dataGroup]['pointerToNextCGBlock']
            self['CGBlock'][dataGroup]={}
            self['CNBlock'][dataGroup]={}
            self['CCBlock'][dataGroup]={}
            for channelGroup in range( self['DGBlock'][dataGroup]['numberOfChannelGroups'] ):
                self['CNBlock'][dataGroup][channelGroup]={}
                self['CCBlock'][dataGroup][channelGroup]={}
                self['CGBlock'][dataGroup][channelGroup] = self.mdfblockread3( self.blockformats3( 'CGFormat' ), fid, CGpointer )
                CGpointer = self['CGBlock'][dataGroup][channelGroup]['pointerToNextCGBlock']

                CGTXpointer = self['CGBlock'][dataGroup][channelGroup]['pointerToChannelGroupCommentText']

                # Read data Text block info into structure
                self['CGBlock'][dataGroup][channelGroup]['TXBlock'] = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, CGTXpointer )

                # Get pointer to next first Channel block
                CNpointer = self['CGBlock'][dataGroup][channelGroup]['pointerToFirstCNBlock']

                # For each Channel
                for channel in range( self['CGBlock'][dataGroup][channelGroup]['numberOfChannels'] ):

                    ### Read Channel block (CNBlock) information
                    #self.numberOfChannels += 1
                    # Read data Channel block info into structure

                    self['CNBlock'][dataGroup][channelGroup][channel] = self.mdfblockread3( self.blockformats3( 'CNFormat' ), fid, CNpointer )
                    CNpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToNextCNBlock']

                    ### Read Channel text blocks (TXBlock)

                    # Clean signal name
                    shortSignalName = self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] # short name of signal
                    CNTXpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToASAMNameBlock']
                    if CNTXpointer > 0:
                        longSignalName = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, CNTXpointer ) # long name of signal
                        longSignalName = longSignalName['Text']
                        if len(longSignalName)>len(shortSignalName): # long name should be used
                            signalname = longSignalName
                        else:
                            signalname = shortSignalName
                    else:
                        signalname = shortSignalName
                    signalname = signalname.split( '\\' )
                    if len(signalname) > 1:
                        self['CNBlock'][dataGroup][channelGroup][channel]['deviceName'] = signalname[1]
                    signalname = signalname[0]
                    if self.filterChannelNames:
                        signalname = signalname.split('.')[-1] #filters channels modules
                    self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] = signalname
                    #self.channelNameList.append( signalname )

                    # Read channel description
                    CNTXpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToChannelCommentBlock']
                    self['CNBlock'][dataGroup][channelGroup][channel]['ChannelCommentBlock'] = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, CNTXpointer )

                    CNTXpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToSignalIdentifierBlock']
                    # Read data Text block info into structure
                    self['CNBlock'][dataGroup][channelGroup][channel]['SignalIdentifierBlock'] = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, CNTXpointer )

                    ### Read Channel Conversion block (CCBlock)

                    # Get pointer to Channel conversion block
                    CCpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToConversionFormula']

                    if CCpointer == 0: # If no conversion formula, set to 1:1
                        self['CCBlock'][dataGroup][channelGroup][channel] = {}
                        self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] = 65535
                    else: # Otherwise get conversion formula, parameters and physical units
                        # Read data Channel Conversion block info into structure
                        self['CCBlock'][dataGroup][channelGroup][channel] = self.mdfblockread3( self.blockformats3( 'CCFormat' ), fid, CCpointer )

                        # Extract Channel Conversion formula based on conversion
                        # type(conversionFormulaIdentifier)

                        # Get current file position
                        currentPosition = fid.tell()

                        if self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 0: # Parameteric, Linear: Physical =Integer*P2 + P1

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula0' ), fid, currentPosition )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 1: # Table look up with interpolation

                            # Get number of parameters sets
                            num = self['CCBlock'][dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula1' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 2: # table look up

                            # Get number of parameters sets
                            num = self['CCBlock'][dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula2' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 6: # Polynomial

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula6' ), fid, currentPosition )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 7: # Exponential

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula7' ), fid, currentPosition )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 8: # Logarithmic

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula8' ), fid, currentPosition )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 9: #  Rational

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula9' ), fid, currentPosition )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 10: #  Text Formula

                            # Read data Channel Conversion parameters info into structure
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula10' ), fid, currentPosition,  removeTrailing0=False )

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 11: #  Text table

                            # Get number of parameters sets
                            num = self['CCBlock'][dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula11' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 12: #  Text range table

                            # Get number of parameters sets
                            num = self['CCBlock'][dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self['CCBlock'][dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ): # first pair is default

                                # Read data Channel Conversion parameters info into structure
                                self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread3( self.blockformats3( 'CCFormatFormula12' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                            for pair in range( num ): # get text range
                                # Read corresponding text
                                temp = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair]['pointerToTXBlock'] )
                                try:
                                    self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair]['Textrange'] = temp['Text']
                                except:
                                    self['CCBlock'][dataGroup][channelGroup][channel]['conversion'][pair]['Textrange'] = ""

                        elif self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 65535: #  No conversion int=phys
                            pass
                        else:

                            # Give warning that conversion formula is not being
                            # made
                            print(( 'Conversion Formula type (conversionFormulaIdentifier=' + str( self['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] ) + ')not supported.' ))

        # CLose the file
        fid.close()
        
    def listChannels3( self, fileName = None ):
        """ reads data, channel group and channel blocks to list channel names
        
        Attributes
        --------------
        fileName : str
            file name
            
        Returns
        -----------
        list of channel names
        """
        # Read MDF file and extract its complete structure
        if not fileName == None:
            self.fileName = fileName
        # Open file
        fid = open( self.fileName, 'rb' )
        channelNameList=[]

        ### Read header block (HDBlock) information
        # Set file pointer to start of HDBlock
        HDpointer = 64
        # Read Header block info into structure
        self['HDBlock'] = self.mdfblockread3( self.blockformats3( 'HDFormat' ), fid, HDpointer )

        ### Read Data Group blocks (DGBlock) information
        # Get pointer to first Data Group block
        DGpointer = self['HDBlock']['pointerToFirstDGBlock']
        for dataGroup in range( self['HDBlock']['numberOfDataGroups'] ):

            # Read data Data Group block info into structure
            self['DGBlock'][dataGroup] = self.mdfblockread3( self.blockformats3( 'DGFormat' ), fid, DGpointer )
            # Get pointer to next Data Group block
            DGpointer = self['DGBlock'][dataGroup]['pointerToNextDGBlock']

            ### Read Channel Group block (CGBlock) information - offset set already

            # Read data Channel Group block info into structure
            CGpointer = self['DGBlock'][dataGroup]['pointerToNextCGBlock']
            self['CGBlock'][dataGroup]={}
            self['CNBlock'][dataGroup]={}
            self['CCBlock'][dataGroup]={}
            for channelGroup in range( self['DGBlock'][dataGroup]['numberOfChannelGroups'] ):
                self['CNBlock'][dataGroup][channelGroup]={}
                self['CCBlock'][dataGroup][channelGroup]={}
                self['CGBlock'][dataGroup][channelGroup] = self.mdfblockread3( self.blockformats3( 'CGFormat' ), fid, CGpointer )
                CGpointer = self['CGBlock'][dataGroup][channelGroup]['pointerToNextCGBlock']

                # Get pointer to next first Channel block
                CNpointer = self['CGBlock'][dataGroup][channelGroup]['pointerToFirstCNBlock']

                # For each Channel
                for channel in range( self['CGBlock'][dataGroup][channelGroup]['numberOfChannels'] ):

                    ### Read Channel block (CNBlock) information
                    #self.numberOfChannels += 1
                    # Read data Channel block info into structure
                    self['CNBlock'][dataGroup][channelGroup][channel] = self.mdfblockread3( self.blockformats3( 'CNFormat' ), fid, CNpointer )
                    CNpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToNextCNBlock']

                    ### Read Channel text blocks (TXBlock)

                    # Clean signal name
                    shortSignalName = self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] # short name of signal
                    CNTXpointer = self['CNBlock'][dataGroup][channelGroup][channel]['pointerToASAMNameBlock']
                    if CNTXpointer > 0:
                        longSignalName = self.mdfblockread3( self.blockformats3( 'TXFormat' ), fid, CNTXpointer ) # long name of signal
                        longSignalName = longSignalName['Text']
                        if len(longSignalName)>len(shortSignalName): # long name should be used
                            signalname = longSignalName
                        else:
                            signalname = shortSignalName
                    else:
                        signalname = shortSignalName
                    signalname = signalname.split( '\\' )
                    signalname = signalname[0]
                    if self.filterChannelNames:
                        signalname = signalname.split('.')[-1]
                    self['CNBlock'][dataGroup][channelGroup][channel]['signalName'] = signalname
                    channelNameList.append( signalname )

        # CLose the file
        fid.close()
        return channelNameList
        
    #######blockformats3####################
    @staticmethod
    def blockformats3( block , version=0):
        """ This function returns all the predefined formats for the different blocks in the MDF file
        
        Parameters
        ----------------
        block : str
            kind of block
        version : int
            mdf version
            
        Returns
        -----------
        nested list of str and int describing structure of block to be used by mdfblockread3 method
        """
        ## Data Type Definitions '<' for little-endian byteOrder
        LINK = '<I'
        CHAR = '<c'
        REAL = '<d'
        BOOL = '<h'
        #UINT8 = '<B'
        #BYTE =  '<B'
        INT16 = '<h'
        UINT16 = '<H'
        UINT32 = '<I'
        #INT32 = '<i'
        UINT64 = '<Q'
        #INT64 = '<q'
        
        formats = ()
        
        if block == 'TXFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( CHAR, 'Var', 'Text' ) )
        elif block == 'CNFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( LINK, 1, 'pointerToNextCNBlock' ),
                ( LINK, 1, 'pointerToConversionFormula' ),
                ( LINK, 1, 'pointerToCEBlock' ),
                ( LINK, 1, 'pointerToCDBlock' ),
                ( LINK, 1, 'pointerToChannelCommentBlock' ),
                ( UINT16, 1, 'channelType' ),
                ( CHAR , 32, 'signalName' ),
                ( CHAR, 128, 'signalDescription' ),
                ( UINT16, 1, 'numberOfTheFirstBits' ),
                ( UINT16, 1, 'numberOfBits' ),
                ( UINT16, 1, 'signalDataType' ),
                ( BOOL, 1, 'valueRangeKnown' ),
                ( REAL, 1, 'valueRangeMinimum' ),
                ( REAL, 1, 'valueRangeMaximum' ),
                ( REAL, 1, 'rateVariableSampled' ),
                ( LINK, 1, 'pointerToASAMNameBlock' ),
                ( LINK, 1, 'pointerToSignalIdentifierBlock' ),
                ( UINT16, 1, 'ByteOffset' ) )
        elif block == 'CCFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( BOOL, 1, 'valueRangeKnown' ),
                ( REAL, 1, 'valueRangeMinimum' ),
                ( REAL, 1, 'valueRangeMaximum' ),
                ( CHAR, 20, 'physicalUnit' ),
                ( UINT16, 1, 'conversionFormulaIdentifier' ),
                ( UINT16, 1, 'numberOfValuePairs' ) )
        elif block == 'CCFormatFormula0':
            formats = ( # parametric, linear
                ( REAL, 1, 'P1' ),
                ( REAL, 1, 'P2' ) )
        elif block in ( 'CCFormatFormula1', 'CCFormatFormula2' ):
            formats = ( # Tablular or Tabular with interp
                ( REAL, 1, 'int' ),
                ( REAL, 1, 'phys' ) )
        elif block == 'CCFormatFormula10':
            formats = ( # Text formula
                ( CHAR , 256, 'textFormula' ),  )
        elif block == 'CCFormatFormula11':
            formats = ( # ASAM-MCD2 text table
                ( REAL, 1, 'int' ),
                ( CHAR, 32, 'text' ) )
        elif block in ( 'CCFormatFormula6', 'CCFormatFormula9' ):
            formats = ( # polynomial
                ( REAL, 1, 'P1' ) ,
                ( REAL, 1, 'P2' ),
                ( REAL, 1, 'P3' ),
                ( REAL, 1, 'P4' ),
                ( REAL, 1, 'P5' ),
                ( REAL, 1, 'P6' ) )
        elif block in ( 'CCFormatFormula7', 'CCFormatFormula8' ):
            formats = ( # Exponential and logarithmic
                ( REAL, 1, 'P1' ) ,
                ( REAL, 1, 'P2' ),
                ( REAL, 1, 'P3' ),
                ( REAL, 1, 'P4' ),
                ( REAL, 1, 'P5' ),
                ( REAL, 1, 'P6' ),
                ( REAL, 1, 'P7' ) )
        elif block == 'CCFormatFormula12':
            formats = ( # ASAM-MCD2 text range table
                ( REAL, 1, 'lowerRange' ),
                ( REAL, 1, 'upperRange' ),
                ( LINK, 1, 'pointerToTXBlock' ) )
        elif block == 'PRFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( CHAR, 'Var', 'PRData' ) )
        elif block == 'CGFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( LINK, 1, 'pointerToNextCGBlock' ),
                ( LINK, 1, 'pointerToFirstCNBlock' ),
                ( LINK, 1, 'pointerToChannelGroupCommentText' ),
                ( UINT16, 1, 'recordID' ),
                ( UINT16, 1, 'numberOfChannels' ),
                ( UINT16, 1, 'dataRecordSize' ),
                ( UINT32, 1, 'numberOfRecords' ) )
                # last one missing
        elif block == 'DGFormat':
            formats = (
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ), #All DGBlocks Size
                ( LINK , 1, 'pointerToNextDGBlock' ),
                ( LINK , 1, 'pointerToNextCGBlock' ),
                ( LINK, 1, 'reserved' ),
                ( LINK, 1, 'pointerToDataRecords' ),
                ( UINT16, 1, 'numberOfChannelGroups' ),
                ( UINT16, 1, 'numberOfRecordIDs' ) )
        elif block == 'HDFormat':
            if version < 3.2:
                formats = (
                    ( CHAR , 2, 'BlockType' ),
                    ( UINT16, 1, 'BlockSize' ),
                    ( LINK, 1, 'pointerToFirstDGBlock' ),
                    ( LINK, 1, 'pointerToTXBlock' ),
                    ( LINK, 1, 'pointerToPRBlock' ),
                    ( UINT16, 1, 'numberOfDataGroups' ),
                    ( CHAR, 10, 'Date' ),
                    ( CHAR, 8, 'Time' ),
                    ( CHAR, 32, 'Author' ),
                    ( CHAR, 32, 'Organization' ),
                    ( CHAR, 32, 'ProjectName' ),
                    ( CHAR, 32, 'Subject' ) )
            else:
                formats = (
                    ( CHAR , 2, 'BlockType' ),
                    ( UINT16, 1, 'BlockSize' ),
                    ( LINK, 1, 'pointerToFirstDGBlock' ),
                    ( LINK, 1, 'pointerToTXBlock' ),
                    ( LINK, 1, 'pointerToPRBlock' ),
                    ( UINT16, 1, 'numberOfDataGroups' ),
                    ( CHAR, 10, 'Date' ),
                    ( CHAR, 8, 'Time' ),
                    ( CHAR, 32, 'Author' ),
                    ( CHAR, 32, 'Organization' ),
                    ( CHAR, 32, 'ProjectName' ),
                    ( CHAR, 32, 'Subject' ),
                    ( UINT64, 1, 'TimeStamp' ),
                    ( INT16, 1, 'UTCTimeOffset' ),
                    ( UINT16, 1, 'TimeQualityClass' ),
                    ( CHAR, 32, 'TimeIdentification' ) )
        elif block== 'IDFormat':
            formats = (
                ( CHAR, 8, 'FileID' ),
                ( CHAR, 8, 'FormatID' ),
                ( CHAR, 8, 'ProgramID' ),
                ( UINT16, 1, 'ByteOrder' ),
                ( UINT16, 1, 'FloatingPointFormat' ),
                ( UINT16, 1, 'VersionNumber' ),
                ( UINT16, 1, 'CodePageNumber' ))
        else:
            print( 'Block format name error ' )
        return formats

    #######mdfblockread3####################
    @staticmethod
    def mdfblockread3( blockFormat, fid, pointer,  removeTrailing0=True ):
        """ Extract block of data from MDF file in original data types. 
        Returns a dictionary with keys specified in data structure blockFormat
        
        Parameters
        ----------------
        blockFormat : nested list
            output of blockformats3 method
        fid : float
            file identifier
        pointer : int
            position of block in file
        removeTrailing0 : bool, optional
            removes or not the trailing 0 from strings
        
        Returns
        -----------
        Block content in a dict
        """
        Block = {}
        if pointer != 0 and not pointer == None:
            fid.seek( pointer )
            # Extract parameters
            for field in range( len( blockFormat ) ):
                fieldTuple = blockFormat[field]
                fieldFormat = fieldTuple[0]
                fieldNumber = fieldTuple[1]
                fieldName = fieldTuple[2]
                if fieldFormat == '<c': # Char
                    if fieldNumber != 'Var':
                        value = fid.read( fieldNumber )
                    elif Block['BlockSize'] and Block['BlockSize'] - 4 > 0:
                        value = fid.read( Block['BlockSize'] - 4 )
                    elif not Block['BlockSize'] :
                        value=''
                    if PythonVersion<3:
                        if removeTrailing0:
                            Block[fieldName] = value.replace('\x00', '')
                        else:
                            Block[fieldName] = value
                    else:
                        if removeTrailing0:
                            Block[fieldName] = value.decode('latin1', 'replace').replace('\x00', '') # decode bytes
                        else:
                            Block[fieldName] = value.decode('latin1', 'replace')
                else: # numeric
                    formatSize = calcsize( fieldFormat )
                    number=fid.read( formatSize * fieldNumber )
                    if number:
                        value = unpack( fieldFormat, number)
                        Block[fieldName] = value[0]
                    else:
                        Block[fieldName] = None
            return Block
