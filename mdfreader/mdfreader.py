# -*- coding: utf-8 -*-
""" Measured Data Format file reader
Created on Sun Oct 10 12:57:28 2010

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

This module contains 2 classes :
First class mdfinfo is meant to extract all blocks in file, giving like metadata of measurement channels
	Method read() will read file and add in the class all Blocks info as defined in MDF specification
	"Format Specification MDF Format Version 3.0", 14/11/2002

:Version: 2010.11.06

Second class mdf reads file and add dict type data
	Each channel is identified by its name the dict key.
	At a higher level this dict contains other keys :
		'data' : containing vector of data (numpy)
		'unit' : unit (string)
		'time' : name of the string key related to the channel (string)
		
	Method plot('channelName') will plot the channel versus corresponding time (not iterable for now)
	Method resample(samplingTime) will resample all channels by interpolation and make only one channel 'time'

Requirements
------------

* `Python 2.6 <http://www.python.org>`__
* `Numpy 1.4 <http://numpy.scipy.org>`__
* `Matplotlib 1.0 <http://matplotlib.sourceforge.net>`__

Examples
--------
>>> import mdfreader
>>> yop=mdfreader.mdf('NameOfFile')
>>> yop.keys() # list channels
>>> yop.plot('channelName')

"""
import io
import struct
import math
#import time
import multiprocessing
import numpy

def processDataBlocks( L, buf, info, channelToList, recordFormat, dataRecordFormat, numberOfRecords, dataGroup ):
	## Processes recorded data blocks
	# Outside of class to allow multiprocessing
	numberOfRecordIDs = info.DGBlock[dataGroup]['numberOfRecordIDs']

	if numberOfRecords != 0:
		buf = [struct.unpack( recordFormat, buf[i] ) for i in range( numberOfRecords )]

		# converts list of records into list of values (for each channel, except bits that are packed)
		buf = zip( *buf )

		np = [( numpy.array( buf[i], dtype = dataRecordFormat[i] ) ) for i in range( len( dataRecordFormat ) )]

		# Check if recordId are used for unsorted writing
		if numberOfRecordIDs >= 1:
			offset = 1
		else:
			offset = 0

		isBitInUnit8 = 0 # Initialise Bit counter
		## Processes Bits, metadata and conversion
		for channelGroup in range( info.DGBlock[dataGroup]['numberOfChannelGroups'] ):
			for channel in range( info.CGBlock[dataGroup][channelGroup]['numberOfChannels'] ):
				# Type of channel data, signed, unsigned, floating or string
				signalDataType = info.CNBlock[dataGroup][channelGroup][channel]['signalDataType']
				# Corresponding number of bits for this format
				numberOfBits = info.CNBlock[dataGroup][channelGroup][channel]['numberOfBits']

				channelName = info.CNBlock[dataGroup][channelGroup][channel]['signalName']
				if channelName == 'time': # assumes first channel is time
					channelName = channelName + str( dataGroup )
				# Process concatenated bits inside unint8
				if signalDataType not in ( 1, 2, 3, 8 ):
					if signalDataType == 0:
						if numberOfBits == 1: # One bit, considered Bool in numpy
							temp = buf[channelToList[channel + offset]] # used datatype is uint8
							mask = 1 << isBitInUnit8 # masks isBitUnit8
							temp = [( ( int( temp[record] ) & mask ) >> ( isBitInUnit8 ) ) for record in range( numberOfRecords )]
							L[channelName] = numpy.array( temp, dtype = 'uint8' )
							isBitInUnit8 += 1
						elif numberOfBits == 2:
							temp = buf[channelToList[channel + offset]] # used datatype is uint8
							mask = 1 << isBitInUnit8
							mask = mask | ( 1 << ( isBitInUnit8 + 1 ) ) # adds second bit in mask
							temp = [( ( int( temp[record] ) & mask ) >> ( isBitInUnit8 ) ) for record in range( numberOfRecords )]
							L[channelName] = numpy.array( temp, dtype = 'uint8' )
							isBitInUnit8 += 2
						else:
							L[channelName] = np[channelToList[channel + offset]]
							isBitInUnit8 = 0 # Initialise Bit counter
					else:
						L[channelName] = np[channelToList[channel + offset]]
						isBitInUnit8 = 0 # Initialise Bit counter
				else:
					L[channelName] = np[channelToList[channel + offset]]
					isBitInUnit8 = 0 # Initialise Bit counter

				## Process Conversion
				conversionFormulaIdentifier = info.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier']
				if conversionFormulaIdentifier == 0: # Parameteric, Linear: Physical =Integer*P2 + P1
					P1 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P1']
					P2 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P2']
					L[channelName] = L[channelName] * P2 + P1
				elif conversionFormulaIdentifier == 11: # Text Table
					pass # considers number representative enough, no need to convert into string, more complex
				elif conversionFormulaIdentifier == 65535: # No conversion, keep data
					pass
				elif conversionFormulaIdentifier == 1: # Tabular with interpolation
					intVal = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['int' ]
					physVal = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['phys' ]
					if numpy.all( numpy.diff( intVal ) > 0 ):
						L[channelName] = numpy.interp( L[channelName], intVal, physVal )
					else:
						print( 'X values for interpolation of channel ' + channelName + ' are not increasing' )
				elif conversionFormulaIdentifier == 2: # Tabular
					intVal = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['int' ]
					physVal = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['phys' ]
					if numpy.all( numpy.diff( intVal ) > 0 ):
						L[channelName] = numpy.interp( L[channelName], intVal, physVal )
					else:
						print( 'X values for interpolation of channel ' + channelName + ' are not increasing' )
				elif conversionFormulaIdentifier == 6: # Polynomial
					P1 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P1']
					P2 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P2']
					P3 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P3']
					P4 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P4']
					P5 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P5']
					P6 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P6']
					L[channelName] = ( P2 - P4 * ( L[channelName]['data'] - P5 - P6 ) ) / ( P3 * ( L[channelName] - P5 - P6 ) - P1 )
				elif conversionFormulaIdentifier == 7: # Exponential
					P1 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P1']
					P2 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P2']
					P3 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P3']
					P4 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P4']
					P5 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P5']
					P6 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P6']
					P7 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P7']
					if P4 == 0 and P1 != 0 and P2 != 0:
						L[channelName]['data'] = math.log( ( ( L[channelName] - P7 ) * P6 - P3 ) / P1 ) / P2
					elif P1 == 0 and P4 != 0 and P5 != 0:
						L[channelName] = math.log( ( P3 / ( L[channelName] - P7 ) - P6 ) / P4 ) / P5
					else:
						print( 'Non possible conversion parameters for channel ' + channelName )
				elif conversionFormulaIdentifier == 8: # Logarithmic
					P1 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P1']
					P2 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P2']
					P3 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P3']
					P4 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P4']
					P5 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P5']
					P6 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P6']
					P7 = info.CCBlock[dataGroup][channelGroup][channel]['conversion'] ['P7']
					if P4 == 0 and P1 != 0 and P2 != 0:
						L[channelName] = math.exp( ( ( L[channelName] - P7 ) * P6 - P3 ) / P1 ) / P2
					elif P1 == 0 and P4 != 0 and P5 != 0:
						L[channelName] = math.exp( ( P3 / ( L[channelName] - P7 ) - P6 ) / P4 ) / P5
					else:
						print( 'Non possible conversion parameters for channel ' + channelName )
				elif conversionFormulaIdentifier == 10: # Text Formula
					print( 'Conversion of formula : ' + info.CCBlock[dataGroup][channelGroup][channel]['conversion']['textFormula'] + 'not yet supported' )
				elif conversionFormulaIdentifier == 12: # Text Range Table
					pass # Not yet supported, practically not used format
		del np

class superDict( dict ):
    # dict with autovivification to allow easy dict nesting
    def __getitem__( self, item ):
        try:
            return dict.__getitem__( self, item )
        except KeyError:
            value = self[item] = type( self )()
            return value

class mdfinfo( superDict ):
    """ mdf file info class
    # MDFINFO is a class information about an MDF (Measure Data Format) file
    #	Based on following specification http://powertrainnvh.com/nvh/MDFspecificationv03.pdf
    #   mdfinfo(FILENAME) contains a dict of structures, for
    #   each data group, containing key information about all channels in each
    #   group. FILENAME is a string that specifies the name of the MDF file.
    #
    #   	mdfinfo.read(FIELNAME) will process Filename file
    #  General dictionary structure is the following :
    #  mdfinfo.HDBlock header block
    #  mdfinfo.DGBlock[dataGroup] Data Group block
    #  mdfinfo.CGBlock[dataGroup][channelGroup] Channel Group block
    #  mdfinfo.CNBlock[dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    #  mdfinfo.CCBlock[dataGroup][channelGroup][channel] Channel conversion information"""

    def __init__( self, fileName = None ):
        self.fileName = fileName
        self.byteOrder = 0 # default is little Endian
        self.version = 300 # default is version 3.00
        self.HDBlock = {} # Header Block
        self.DGBlock = {} # Data Group Block
        self.CGBlock = superDict() # Channel Group Block
        self.CNBlock = superDict() # Channel Block
        self.CCBlock = superDict() # Conversion block
        self.numberOfChannels = 0 # number of channels
        self.channelNameList = [] # List of channel names
        if fileName != None:
            self.readinfo( fileName )
    ## Reads block informations inside file
    def readinfo( self, fileName = None ):
        # Read MDF file and extract its complete structure
        if self.fileName == None:
            self.fileName = fileName
        # Open file
        fid = io.open( self.fileName, 'rb' )
		### Read some Identifier block (IDBlock) informations
        fid.seek( 24 ) # Look for byte order
        self.byteOrder=struct.unpack( '<H', fid.read( 2 ) )
        if self.byteOrder[0]==0:
        	self.byteOrder=0 #little Endian to be used
        else:
        	self.byteOrder=1 #Big Endian to be used
		fid.seek( 28 ) # Look for version number
        self.version=struct.unpack( '<H', fid.read( 2 ) )
        print self.version
        self.version=self.version[0]

        ### Read header block (HDBlock) information
        # Set file pointer to start of HDBlock
        HDpointer = 64
        # Read Header block info into structure
        self.HDBlock = self.mdfblockread( self.blockformats( 'HDFormat' ,self.version, self.byteOrder ), fid, HDpointer )

        ### Read text block (TXBlock) information
        self.HDBlock['TXBlock'] = self.mdfblockread( self.blockformats( 'TXFormat' ), fid, self.HDBlock['pointerToTXBlock'] )

        ### Read text block (PRBlock) information
        self.HDBlock['PRBlock'] = self.mdfblockread( self.blockformats( 'PRFormat' ), fid, self.HDBlock['pointerToPRBlock'] )

        ### Read Data Group blocks (DGBlock) information
        # Get pointer to first Data Group block
        DGpointer = self.HDBlock['pointerToFirstDGBlock']
        for dataGroup in range( self.HDBlock['numberOfDataGroups'] ):

            # Read data Data Group block info into structure
            self.DGBlock[dataGroup] = self.mdfblockread( self.blockformats( 'DGFormat' ), fid, DGpointer )
            # Get pointer to next Data Group block
            DGpointer = self.DGBlock[dataGroup]['pointerToNextDGBlock']

            ### Read Channel Group block (CGBlock) information - offset set already

            # Read data Channel Group block info into structure
            CGpointer = self.DGBlock[dataGroup]['pointerToNextCGBlock']
            for channelGroup in range( self.DGBlock[dataGroup]['numberOfChannelGroups'] ):
                self.CGBlock[dataGroup][channelGroup] = self.mdfblockread( self.blockformats( 'CGFormat' ), fid, CGpointer )
                CGpointer = self.CGBlock[dataGroup][channelGroup]['pointerToNextCGBlock']

                CGTXpointer = self.CGBlock[dataGroup][channelGroup]['pointerToChannelGroupCommentText']

                # Read data Text block info into structure
                self.CGBlock[dataGroup][channelGroup]['TXBlock'] = self.mdfblockread( self.blockformats( 'TXFormat' ), fid, CGTXpointer )

                # Get pointer to next first Channel block
                CNpointer = self.CGBlock[dataGroup][channelGroup]['pointerToFirstCNBlock']

                # For each Channel
                for channel in range( self.CGBlock[dataGroup][channelGroup]['numberOfChannels'] ):

                    ### Read Channel block (CNBlock) information
                    self.numberOfChannels += 1
                    # Read data Channel block info into structure
                    self.CNBlock[dataGroup][channelGroup][channel] = self.mdfblockread( self.blockformats( 'CNFormat' ), fid, CNpointer )
                    CNpointer = self.CNBlock[dataGroup][channelGroup][channel]['pointerToNextCNBlock']

                    ### Read Channel text blocks (TXBlock)

                    # Clean signal name
                    signalname = self.CNBlock[dataGroup][channelGroup][channel]['signalName']
                    slashlocation = signalname.find( '\\' )
                    if slashlocation > 0:
                        if slashlocation + 1 < len( signalname ):
                            self.CNBlock[dataGroup][channelGroup][channel]['deviceName'] = signalname[slashlocation + 1:len( signalname )]
                        # Clean signal name from device name
                        signalname = signalname[0:slashlocation]
                    self.CNBlock[dataGroup][channelGroup][channel]['signalName'] = signalname
                    self.channelNameList.append( signalname )

                    CNTXpointer = self.CNBlock[dataGroup][channelGroup][channel]['pointerToChannelCommentBlock']
                    # Read data Text block info into structure
                    self.CNBlock[dataGroup][channelGroup][channel]['ChannelCommentBlock'] = self.mdfblockread( self.blockformats( 'TXFormat' ), fid, CNTXpointer )

                    CNTXpointer = self.CNBlock[dataGroup][channelGroup][channel]['pointerToASAMNameBlock']
                    # Read data Text block info into structure
                    temp = self.mdfblockread( self.blockformats( 'TXFormat' ), fid, CNTXpointer )
                    signalname = temp['Text']
                    slashlocation = signalname.find( '\\' )
                    if slashlocation > 0:
                        if slashlocation + 1 < len( signalname ):
                            # Save device name
                            temp['DeviceName'] = signalname[slashlocation + 1:len( signalname )]
                        # Clean signal name from device name
                        signalname = signalname[0:slashlocation]
                    temp['Text'] = signalname
                    self.CNBlock[dataGroup][channelGroup][channel]['ASAMNameBlock'] = temp

                    CNTXpointer = self.CNBlock[dataGroup][channelGroup][channel]['pointerToSignalIdentifierBlock']
                    # Read data Text block info into structure
                    self.CNBlock[dataGroup][channelGroup][channel]['SignalIdentifierBlock'] = self.mdfblockread( self.blockformats( 'TXFormat' ), fid, CNTXpointer )

                    ### Read Channel Conversion block (CCBlock)

                    # Get pointer to Channel conversion block
                    CCpointer = self.CNBlock[dataGroup][channelGroup][channel]['pointerToConversionFormula']

                    if CCpointer == 0: # If no conversion formula, set to 1:1
                        self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] = 65535
                    else: # Otherwise get conversion formula, parameters and physical units
                        # Read data Channel Conversion block info into structure
                        self.CCBlock[dataGroup][channelGroup][channel] = self.mdfblockread( self.blockformats( 'CCFormat' ), fid, CCpointer )

                        # Extract Channel Conversion formula based on conversion
                        # type(conversionFormulaIdentifier)

                        # Get current file position
                        currentPosition = fid.tell()

                        if self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 0: # Parameteric, Linear: Physical =Integer*P2 + P1

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula0' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 1: # Table look up with interpolation

                            # Get number of parameters sets
                            num = self.CCBlock[dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self.CCBlock[dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread( self.blockformats( 'CCFormatFormula1' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 2: # table look up

                            # Get number of parameters sets
                            num = self.CCBlock[dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self.CCBlock[dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread( self.blockformats( 'CCFormatFormula2' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 6: # Polynomial

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula6' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 7: # Exponential

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula7' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 8: # Logarithmic

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula8' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 9: #  Rational

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula9' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 10: #  Text Formula

                            # Read data Channel Conversion parameters info into structure
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = self.mdfblockread( self.blockformats( 'CCFormatFormula10' ), fid, currentPosition )

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 11: #  Text table

                            # Get number of parameters sets
                            num = self.CCBlock[dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num ):

                                # Read data Channel Conversion parameters info into structure
                                self.CCBlock[dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread( self.blockformats( 'CCFormatFormula11' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 12: #  Text range table

                            # Get number of parameters sets
                            num = self.CCBlock[dataGroup][channelGroup][channel]['numberOfValuePairs']
                            self.CCBlock[dataGroup][channelGroup][channel]['conversion'] = {}

                            for pair in range( num + 1 ): # first 2 pairs are ignored

                                # Read data Channel Conversion parameters info into structure
                                self.CCBlock[dataGroup][channelGroup][channel]['conversion'][pair] = self.mdfblockread( self.blockformats( 'CCFormatFormula12' ), fid, currentPosition )
                                # Get current file position
                                currentPosition = fid.tell()

                            for pair in range( num ): # get text range    
                                # Read corresponding text
                                temp = self.mdfblockread( self.blockformats( 'TXBlock' ), fid, self.CCBlock[dataGroup][channelGroup][channel][pair]['pointerToTXBlock'] )
                                self.CCBlock[dataGroup][channelGroup][channel]['conversion'][pair]['Textrange'] = temp['Text']

                        elif self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] == 65535: #  No conversion int=phys
                            pass
                        else:

                            # Give warning that conversion formula is not being
                            # made
                            print( 'Conversion Formula type (conversionFormulaIdentifier=' + str( self.CCBlock[dataGroup][channelGroup][channel]['conversionFormulaIdentifier'] ) + ')not supported.' )

        # CLose the file
        fid.close()

    #######BLOCKFORMATS####################
    @staticmethod
    def blockformats( block, version=300, byteOrder=0):
        # This function returns all the predefined formats for the different blocks
        # in the MDF file as specified in "Format Specification MDF Format Version 3.0"
        # document version 2.0, 14/11/2002
        if byteOrder==0:
			endian='<'
        else:
			endian='>'
        ## Data Type Definitions '<' for little-endian byteOrder

        CHAR = endian+'c'
        REAL = endian+'d'
        BOOL = endian+'h'
        UINT8 = endian+'B'
        INT16 = endian+'h'
        UINT16 = endian+'H'
        UINT32 = endian+'I'
        UINT64 = endian+'Q'
        if version<320:
			LINK = endian+'i'
        else:
			LINK = endian+'I'
		
        if block == 'HDFormat':
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
                ( CHAR, 32, 'Vehicle' ) )
            if version>=320:
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
                ( CHAR, 32, 'Vehicle' ) ,
				( UINT64, 1, 'TimeStamp' ),
				( INT16, 1, 'UTCTiemOffset' ),
				( UINT16, 1, 'TimeQualityClass' ),
				( CHAR, 32, 'TimeIdentification' ))
        elif block == 'TXFormat':
            formats = ( 
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( CHAR, 'Var', 'Text' ) )
        elif block == 'PRFormat':
            formats = ( 
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ),
                ( CHAR, 'Var', 'PRData' ) )
        elif block == 'DGFormat':
            formats = ( 
                ( CHAR, 2, 'BlockType' ),
                ( UINT16, 1, 'BlockSize' ), #All DGBlocks Size
                ( LINK , 1, 'pointerToNextDGBlock' ),
                ( LINK , 1, 'pointerToNextCGBlock' ),
                ( LINK, 1, 'pointerToTRBlock' ),
                ( LINK, 1, 'pointerToDataRecords' ),
                ( UINT16, 1, 'numberOfChannelGroups' ),
                ( UINT16, 1, 'numberOfRecordIDs' ) )
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
                ( CHAR , 256, 'textFormula' ) )
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
        elif block == 'TRBlock':
			formats = ( # Trigger block
				( CHAR, 2, 'BlockType' ),
				( UINT16, 1, 'BlockSize' ),
				( LINK, 1, 'pointerToTXBlock' ),
				( UINT16, 1, 'numberOfTriggerEvents' ))
        else:
            print( 'Block format name error ' )
        return formats

    #######MDFBLOCKREAD####################
    @staticmethod
    def mdfblockread( blockFormat, fid, pointer ):
        # MDFBLOCKREAD Extract block of data from MDF file in original data types
        #   Block=MDFBLOCKREAD(BLOCKFORMAT, FID, OFFSET) returns a
        #   dictionary with keys specified in data structure BLOCKFORMAT, fid
        #   FID, at byte offset in the file OFFSET
        #
        # Example block format is:
        # formats=(
        #        (CHAR,2,'BlockType'),
        #        (UINT16,1,'BlockSize'),
        #        (CHAR,'Var',Text))
        #
        # Example function call is:
        # Block=mdfblockread(blockFormat, 1, 413)
        fid.seek( pointer )
        Block = {}
        if pointer != 0:
            # Extract parameters
            for field in range( len( blockFormat ) ):
                fieldTuple = blockFormat[field]
                fieldFormat = fieldTuple[0]
                fieldNumber = fieldTuple[1]
                fieldName = fieldTuple[2]
                if fieldFormat == '<c':
                    if fieldNumber != 'Var':
                        value = fid.read( fieldNumber )
                    elif Block['BlockSize'] - 4 > 0:
                        value = fid.read( Block['BlockSize'] - 4 )
                    value = value.strip( '\x00' )
                    Block[fieldName] = value
                else:
                    formatSize = struct.calcsize( fieldFormat )
                    value = struct.unpack( fieldFormat, fid.read( formatSize * fieldNumber ) )
                    Block[fieldName] = value[0]
            return Block


class mdf( mdfinfo, superDict ):
	" mdf file class"

	def __init__( self, fileName = None ):
		self.fileName = fileName
		self.timeChannelList = []
		if fileName != None:
			self.read( fileName )

	## reads mdf file
	def read( self, fileName ):
		# read mdf file
		if self.fileName == None:
		    self.fileName = fileName

		#inttime = time.clock()
		## Read information block from file
		info = mdfinfo( self.fileName )

		# Open file
		fid = io.open( fileName, 'rb' )
		# prepare multiprocessing of dataGroups
		proc = []
		manager = multiprocessing.Manager()
		L = manager.dict()

		## Defines record format
		for dataGroup in range( info.HDBlock['numberOfDataGroups'] ):
			# Number for records in data block, meaning number of samples
			numberOfRecordIDs = info.DGBlock[dataGroup]['numberOfRecordIDs']
			#Pointer to data block
			pointerToData = info.DGBlock[dataGroup]['pointerToDataRecords']
			fid.seek( pointerToData )
			# initialize counter of bits
			precedingNumberOfBits = 0
			dataRecordNumber = 0
			dataRecordFormat = []
			formatSize = 0
			channelToList = {}
			# Check if recordId are used for unsorted writing
			if numberOfRecordIDs >= 1:
				# start first with uint8
				recordFormat = '<B' # adds uint8
				formatSize += 1 # increases by 1byte recordFomat size
			else: # otherwise initialize recordFormat
				recordFormat = '<' # little-endian
				formatSize = 0
			# Prcess dataFormats
			for channelGroup in range( info.DGBlock[dataGroup]['numberOfChannelGroups'] ):
				# Reads size of each channel
				#dataRecordSize = info.CGBlock[dataGroup][channelGroup]['dataRecordSize']
				# Number of samples recorded
				numberOfRecords = info.CGBlock[dataGroup][channelGroup]['numberOfRecords']

				for channel in range( info.CGBlock[dataGroup][channelGroup]['numberOfChannels'] ):
					# Type of channel data, signed, unsigned, floating or string
					signalDataType = info.CNBlock[dataGroup][channelGroup][channel]['signalDataType']
					# Corresponding number of bits for this format
					numberOfBits = info.CNBlock[dataGroup][channelGroup][channel]['numberOfBits']
					# Defines data format
					datatype = self.datatypeformat( signalDataType, numberOfBits )

					if numberOfBits < 8: # adding bit format, ubit1 or ubit2
						if precedingNumberOfBits == 8: # 8 bit make a byte
						    recordFormat = recordFormat + 'B' # adds uint8
						    formatSize += 1 # increases by 1byte recordFomat size
						    channelToList[channel] = dataRecordNumber
						    dataRecordNumber += 1  # 1 additionnal item in future dataRecord
						    dataRecordFormat.append( self.arrayformat( signalDataType, numberOfBits ) )
						    precedingNumberOfBits = 0
						else:
						    precedingNumberOfBits += numberOfBits # counts successive bits
						    channelToList[channel] = dataRecordNumber
						    # Do not increase dataRecordNumber as bit channel inside byte uint8
					else: # adding bytes
						if precedingNumberOfBits != 0: # There was bits in previous channel
						    recordFormat = recordFormat + 'B' # adds uint8
						    formatSize += 1 # increases by 1byte recordFomat size
						    precedingNumberOfBits = 0

						if signalDataType == 7: # string
						    recordFormat = recordFormat + str( numberOfBits / 8 ) + datatype # append datatype to record format

						elif signalDataType == 8: # array of bytes
						    recordFormat = recordFormat + str( numberOfBits / 8 ) + datatype # append datatype to record format

						else: # append numerical byte datatype to record format
						    recordFormat = recordFormat + datatype

						formatSize = formatSize + numberOfBits / 8 # increases number of bytes
						channelToList[channel] = dataRecordNumber # channelToList converts index of channel into fread output list index
						dataRecordNumber += 1 # 1 additionnal item in future dataRecord
						dataRecordFormat.append( self.arrayformat( signalDataType, numberOfBits ) )

				if numberOfBits < 8: # last channel in record is bit inside byte unit8
					recordFormat = recordFormat + 'B' # adds uint8
					formatSize += 1 # increases by 1byte recordFomat size
					dataRecordNumber += 1  # 1 additionnal item in future dataRecord
					dataRecordFormat.append( self.arrayformat( signalDataType, numberOfBits ) )
					precedingNumberOfBits = 0

			# Offset of record id at the end of record
			if numberOfRecordIDs == 2:
				recordFormat = recordFormat + 'B' # adds unsigned int
				formatSize += 1 # increases by 1byte recordFomat size

			## reads recursively the records
			buf = [fid.read( formatSize ) for i in range( numberOfRecords ) ]

			proc.append( multiprocessing.Process( target = processDataBlocks,
				args = ( L, buf, info, channelToList, recordFormat, dataRecordFormat, numberOfRecords, dataGroup ) ) )
			proc[dataGroup].start()

		fid.close() # close file
		for i in range ( len( proc ) ): # Make sure all processes are finished
			try:
				proc[i].join()
			except:
				print ( 'process ' + str( i ) + ' can not be closed' )
				pass

		tru = []
		truc = ''
		for dataGroup in range( info.HDBlock['numberOfDataGroups'] ):
			for channelGroup in range( info.DGBlock[dataGroup]['numberOfChannelGroups'] ):
				for channel in range( info.CGBlock[dataGroup][channelGroup]['numberOfChannels'] ):
					numberOfRecords = info.CGBlock[dataGroup][channelGroup]['numberOfRecords']
					if numberOfRecords != 0 :
						channelName = info.CNBlock[dataGroup][channelGroup][channel]['signalName']
						if channelName == 'time': # assumes first channel is time
							channelName = channelName + str( dataGroup )
							self.timeChannelList.append( channelName )
						if channelName in L :
							# corresponding time channel of considered channel
							self[channelName]['time'] = 'time' + str( dataGroup )
							self[channelName]['unit'] = info.CCBlock[dataGroup][channelGroup][channel]['physicalUnit'] # Unit of channel
							self[channelName]['description'] = info.CNBlock[dataGroup][channelGroup][channel]['signalDescription']
							if not ( type( L[channelName] ) == type( tru ) or type( L[channelName] ) == type( truc ) ) :
								self[channelName]['data'] = L[channelName]

	@staticmethod
	def datatypeformat( signalDataType, numberOfBits ):
		# DATATYPEFORMAT Data type format precision to give to fread
		#   DATATYPEFORMAT(SIGNALDATATYPE,NUMBEROFBITS) is the precision string to
		#   give to fread for reading the data type specified by SIGNALDATATYPE and
		#   NUMBEROFBITS

		if signalDataType == 0: # unsigned
			if numberOfBits == 8:
				dataType = 'B'
			elif numberOfBits == 16:
				dataType = 'H'
			elif numberOfBits == 32:
				dataType = 'I'
			elif numberOfBits == 1:
				dataType = 'ubit1' # not directly processed
			elif numberOfBits == 2:
				dataType = 'ubit2' # not directly processed
			else:
				print( 'Unsupported number of bits for unsigned int ' + str( dataType ) )

		elif signalDataType == 1: # signed int
			if numberOfBits == 8:
				dataType = 'b'
			elif numberOfBits == 16:
				dataType = 'h'
			elif numberOfBits == 32:
				dataType = 'i'
			else:
				print( 'Unsupported number of bits for signed int ' + str( dataType ) )

		elif signalDataType in ( 2, 3 ): # floating point
			if numberOfBits == 32:
				dataType = 'f'
			elif numberOfBits == 64:
				dataType = 'd'
			else:
				print( 'Unsupported number of bit for floating point ' + str( dataType ) )

		elif signalDataType == 7: # string
		    dataType = 's'
		elif signalDataType == 8: # array of bytes
		    dataType = 'B' # take unit8 as default ?
		else:
		    print ( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits )
		return dataType

	@staticmethod
	def arrayformat( signalDataType, numberOfBits ):
		# Formats used by numpy

		if signalDataType == 0: # unsigned
			if numberOfBits == 8:
				dataType = 'uint8'
			elif numberOfBits == 16:
				dataType = 'uint16';
			elif numberOfBits == 32:
				dataType = 'uint32'
			elif numberOfBits == 1:
				dataType = 'uint8' # not directly processed
			elif numberOfBits == 2:
				dataType = 'uint8' # not directly processed
			else:
				print( 'Unsupported number of bits for unsigned int ' + str( dataType ) )

		elif signalDataType == 1: # signed int
			if numberOfBits == 8:
				dataType = 'int8'
			elif numberOfBits == 16:
				dataType = 'int16'
			elif numberOfBits == 32:
				dataType = 'int32'
			else:
				print( 'Unsupported number of bits for signed int ' + str( dataType ) )

		elif signalDataType in ( 2, 3 ): # floating point
			if numberOfBits == 32:
				dataType = 'float32'
			elif numberOfBits == 64:
				dataType = 'float64'
			else:
				print( 'Unsupported number of bit for floating point ' + str( dataType ) )

		elif signalDataType == 7: # string
		    dataType = 'c' # not directly processed
		elif signalDataType == 8: # array of bytes
		    dataType = 'B' # not directly processed
		else:
		    print( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits )
		return dataType


	def plot( self, channelName ):
		import matplotlib.pyplot as plt
		if self[channelName]['data'].dtype != '|S1': # if channel not a string
			# plot using matplotlib the channel versus time
			if self[channelName].has_key( 'time' ):
				timeName = self[channelName]['time']
				if timeName == '': # resampled channels, only one time channel called 'time'
					timeName = 'time'
			else:
				timeName = 'time'
			self.fig = plt.figure()
			plt.plot( self[timeName]['data'], self[channelName]['data'] )
			plt.title( self[channelName]['description'] )
			plt.xlabel( timeName + ' [' + self[timeName]['unit'] + ']' )
			plt.ylabel( channelName + ' [' + self[channelName]['unit'] + ']' )
			plt.grid( True )
			plt.show()

	def allPlot( self ):
		# plot all channels in the object, be careful for test purpose only,
		# can display many many many plots overloading your computer
		for Name in self.keys():
			try:
				self.plot( Name )
			except:
				print( Name )

	def resample( self, samplingTime ):
		# resample all channels to one sampling time vector
		channelNames = self.keys()
		try: # define new time
			newTime = numpy.arange( 0, self['time0']['data'][len( self['time0']['data'] ) - 1], samplingTime )
			self['time']['data'] = newTime
			self['time']['unit'] = self['time0']['unit']
		except:
			print( 'Cannot process time, already resampled or empty object' )
		for Name in channelNames:
			try:
				if Name not in self.timeChannelList:
					timevect = self[self[Name]['time']]['data']
					if self[Name]['data'].dtype != '|S1': # if channel not a string
						self[Name]['data'] = numpy.interp( newTime, timevect, self[Name]['data'] )
					del self[Name]['time']
			except:
				if len( timevect ) != len( self[Name]['data'] ):
					print( Name + ' and time channel ' + self[Name]['time'] + ' do not have same length' )
				elif numpy.all( numpy.diff( timevect ) > 0 ):
					print( Name + ' has non regularly increasing time channel ' + self[Name]['time'] )
		# remove time channels in timeChannelList
		for ind in range( len( self.timeChannelList ) ):
			del self[self.timeChannelList[ind]]
		self.timeChannelList = [] # empty list
		self.timeChannelList.append( 'time' )

if __name__ == '__main__': # to allow multiprocessing
	MDF = mdf()
	MDF.read()
