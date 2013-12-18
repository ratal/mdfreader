# -*- coding: utf-8 -*-
""" Measured Data Format file reader
Created on Sun Oct 10 12:57:28 2010

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

This module contains 2 classes :
First class mdfinfo is meant to extract all blocks in file, giving like metadata of measurement channels
    Method read() will read file and add in the class all Blocks info as defined in MDF specification
    "Format Specification MDF Format Version 3.0", 14/11/2002

:Version: 2013.12.16 r21

Second class mdf reads file and add dict type data
    Each channel is identified by its name the dict key.
    At a higher level this dict contains other keys :
        'data' : containing vector of data (numpy)
        'unit' : unit (string)
        'time' : name of the string key related to the channel (string)
        'description' : Description of channel

    Method plot('channelName') will plot the channel versus corresponding time
    Method resample(samplingTime) will resample all channels by interpolation and make only one channel 'time'

Requirements
------------

* `Python >2.6 <http://www.python.org>`__
* `Numpy >1.6 <http://numpy.scipy.org>`__
* `Matplotlib >1.0 <http://matplotlib.sourceforge.net>`__
* 'NetCDF
* 'h5py for the HDF5 export
* 'xlwt for the excel export (not existing for python3)
* 'openpyxl for the excel 2007 export
* 'scipy for the Matlab file conversion

Examples
--------
>>> import mdfreader
>>> yop=mdfreader.mdf('NameOfFile')
>>> yop.keys() # list channels names
>>> yop.plot('channelName')
>>> yop.resample(0.1)
>>> yop.exportoCSV(sampling=0.01)
>>> yop.exportNetCDF()
>>> yop.exporttoHDF5()
>>> yop.exporttoMatlab()
>>> yop.exporttoExcel()
>>> yop.exporttoXlsx()
>>> yop.keepChannels(ChannelListToKeep)

"""
from io import open
from struct import unpack
from math import ceil, log, exp
import numpy

from sys import platform, version_info
PythonVersion=version_info
PythonVersion=PythonVersion[0]
channelTime={'time','TIME'} # possible channel time writing, adjust if necessary

def processDataBlocks( Q, buf, info, numberOfRecords, dataGroup,  multiProc ):
    ## Processes recorded data blocks
    # Outside of class to allow multiprocessing
    #numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
    #procTime = time.clock()
    #print( 'Start process ' + str( dataGroup ) + ' Number of Channels ' + str( info['CGBlock'][dataGroup][0]['numberOfChannels'] ) + ' of length ' + str( len( buf ) ) + ' ' + str( time.clock() ) )
    L={}
    isBitInUnit8 = 0 # Initialize Bit counter
    previousChannelName=info['CNBlock'][0][0][0]['signalName'] # initialise channelName

    ## Processes Bits, metadata and conversion
    for channelGroup in range( info['DGBlock'][dataGroup]['numberOfChannelGroups'] ):
        for channel in range( info['CGBlock'][dataGroup][channelGroup]['numberOfChannels'] ):
            # Type of channel data, signed, unsigned, floating or string
            signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['signalDataType']
            # Corresponding number of bits for this format
            numberOfBits = info['CNBlock'][dataGroup][channelGroup][channel]['numberOfBits']

            cName = info['CNBlock'][dataGroup][channelGroup][channel]['signalName']
            channelName=cName
            try:
                temp = buf.__getattribute__( str(cName)) # extract channel vector
                previousChannelName=cName
            except: # if tempChannelName not in buf -> bits in unit8
                temp = buf.__getattribute__( str(previousChannelName) ) # extract channel vector

            if cName in channelTime: # assumes first channel is time
                channelName = 'time' + str( dataGroup )
            # Process concatenated bits inside unint8
            if signalDataType not in ( 1, 2, 3, 8 ):
                if signalDataType == 0:
                    if numberOfBits == 1: # One bit, considered Boolean in numpy
                        mask = 1 << isBitInUnit8 # masks isBitUnit8
                        temp = [( ( int( temp[record] ) & mask ) >> ( isBitInUnit8 ) ) for record in range( numberOfRecords )]
                        L[channelName] = numpy.array( temp, dtype = 'uint8' )
                        isBitInUnit8 += 1
                    elif numberOfBits == 2:
                        mask = 1 << isBitInUnit8
                        mask = mask | ( 1 << ( isBitInUnit8 + 1 ) ) # adds second bit in mask
                        temp = [( ( int( temp[record] ) & mask ) >> ( isBitInUnit8 ) ) for record in range( numberOfRecords )]
                        L[channelName] = numpy.array( temp, dtype = 'uint8' )
                        isBitInUnit8 += 2
                    else:
                        L[channelName] = temp
                        isBitInUnit8 = 0 # Initialize Bit counter
                else:
                    L[channelName] = temp
                    isBitInUnit8 = 0 # Initialize Bit counter
            else:
                L[channelName] = temp
                isBitInUnit8 = 0 # Initialize Bit counter

            ## Process Conversion
            conversionFormulaIdentifier = info['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier']
            if conversionFormulaIdentifier == 0: # Parametric, Linear: Physical =Integer*P2 + P1
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P2']
                L[channelName] = L[channelName] * P2 + P1
            elif conversionFormulaIdentifier == 1: # Tabular with interpolation
                intVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['int' ]
                physVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['phys' ]
                if numpy.all( numpy.diff( intVal ) > 0 ):
                    L[channelName] = numpy.interp( L[channelName], intVal, physVal )
                else:
                    print(( 'X values for interpolation of channel ' + channelName + ' are not increasing' ))
            elif conversionFormulaIdentifier == 2: # Tabular
                intVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['int' ]
                physVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['phys' ]
                if numpy.all( numpy.diff( intVal ) > 0 ):
                    L[channelName] = numpy.interp( L[channelName], intVal, physVal )
                else:
                    print(( 'X values for interpolation of channel ' + str(channelName) + ' are not increasing' ))
            elif conversionFormulaIdentifier == 6: # Polynomial
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P6']
                L[channelName] = ( P2 - P4 * ( L[channelName]['data'] - P5 - P6 ) ) / ( P3 * ( L[channelName] - P5 - P6 ) - P1 )
            elif conversionFormulaIdentifier == 7: # Exponential
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P6']
                P7 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P7']
                if P4 == 0 and P1 != 0 and P2 != 0:
                    L[channelName]['data'] = log( ( ( L[channelName] - P7 ) * P6 - P3 ) / P1 ) / P2
                elif P1 == 0 and P4 != 0 and P5 != 0:
                    L[channelName] = log( ( P3 / ( L[channelName] - P7 ) - P6 ) / P4 ) / P5
                else:
                    print(( 'Non possible conversion parameters for channel ' + str(channelName) ))
            elif conversionFormulaIdentifier == 8: # Logarithmic
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P6']
                P7 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P7']
                if P4 == 0 and P1 != 0 and P2 != 0:
                    L[channelName] = exp( ( ( L[channelName] - P7 ) * P6 - P3 ) / P1 ) / P2
                elif P1 == 0 and P4 != 0 and P5 != 0:
                    L[channelName] = exp( ( P3 / ( L[channelName] - P7 ) - P6 ) / P4 ) / P5
                else:
                    print(( 'Non possible conversion parameters for channel ' + str(channelName) ))
            elif conversionFormulaIdentifier == 9: # rationnal
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion'] ['P6']
                L[channelName] = (P1*L[channelName] *L[channelName] +P2*L[channelName] +P3)/(P4*L[channelName] *L[channelName] +P5*L[channelName] +P6)
            elif conversionFormulaIdentifier == 10: # Text Formula, not really supported yet
                print(( 'Conversion of formula : ' + str(info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['textFormula']) + 'not yet supported' ))
            elif conversionFormulaIdentifier == 11: # Text Table, not really supported yet
                pass # considers number representative enough, no need to convert into string, more complex
            elif conversionFormulaIdentifier == 12: # Text Range Table, not really supported yet
                print('Not supported text range tables')
                pass # Not yet supported, practically not used format
            elif conversionFormulaIdentifier == 132: # date, not really supported yet
                pass
            elif conversionFormulaIdentifier == 133: # time, not really supported yet
                pass
            elif conversionFormulaIdentifier == 65535: # No conversion, keep data
                pass
            else:
                print('Unrecognized conversion type')
    if multiProc:
        Q.put(L)
    else:
        return L
        #print( 'Process ' + str( dataGroup ) + ' Finished ' + str( time.clock() - procTime ) )

def convertMatlabName(channel):
# removes non allowed characters for a matlab variable name
    if PythonVersion<3:
        channel=channel.decode('utf-8')
    channelName=channel.replace('[', '_ls_')
    channelName=channelName.replace(']', '_rs_')
    channelName=channelName.replace('$', '')
    channelName=channelName.replace('.', 'p')
    channelName=channelName.replace('\\','_bs_')
    channelName=channelName.replace('/','_fs_')
    channelName=channelName.replace('(','_lp_')
    channelName=channelName.replace(',','_rp_')
    channelName=channelName.replace('@','_am_')
    channelName=channelName.replace(' ','_')
    channelName=channelName.replace(':','_co_')
    channelName=channelName.replace('-','_hy_')
    channelName=channelName.replace('-','_hy_')
    return channelName

class mdfinfo( dict):
    """ mdf file info class
    # MDFINFO is a class information about an MDF (Measure Data Format) file
    #   Based on following specification http://powertrainnvh.com/nvh/MDFspecificationv03.pdf
    #   mdfinfo(FILENAME) contains a dict of structures, for
    #   each data group, containing key information about all channels in each
    #   group. FILENAME is a string that specifies the name of the MDF file.
    #
    #       mdfinfo.readinfo(FILENAME) will process Filename file
    #  General dictionary structure is the following :
    #  mdfinfo['HDBlock'] header block
    #  mdfinfo['DGBlock'][dataGroup] Data Group block
    #  mdfinfo['CGBlock'][dataGroup][channelGroup] Channel Group block
    #  mdfinfo['CNBlock'][dataGroup][channelGroup][channel] Channel block including text blocks for comment and identifier
    #  mdfinfo['CCBlock'][dataGroup][channelGroup][channel] Channel conversion information"""

    def __init__( self, fileName = None ):
        self.fileName = fileName
        if fileName != None:
            self.readinfo( fileName )
    ## Reads block informations inside file
    def readinfo( self, fileName = None ):
        # Read MDF file and extract its complete structure
        if self.fileName == None:
            self.fileName = fileName
        # Open file
        try:
            fid = open( self.fileName, 'rb' )
        except IOError:
            print('Can not find file'+self.fileName)
            raise
        # read Identifier block
        fid.seek(28)
        VersionNumber=unpack( '<H', fid.read( 2 ) )
        VersionNumber=VersionNumber[0]
        if VersionNumber<400: # up to version 3.x not compatible with version 4.x
            from mdfinfo3 import info3
            self.update(info3(None, fid))
        else: #MDF version 4.x
            from mdfinfo4 import info4
            self.update(info4(None, fid))

    def listChannels( self, fileName = None ):
        # Read MDF file and extract its complete structure
        if self.fileName == None:
            self.fileName = fileName
        # Open file
        try:
            fid = open( self.fileName, 'rb' )
        except IOError:
            print('Can not find file'+self.fileName)
            raise
        # read Identifier block
        fid.seek(28)
        VersionNumber=unpack( '<H', fid.read( 2 ) )
        VersionNumber=VersionNumber[0]
        if VersionNumber<400: # up to version 3.x not compatible with version 4.x
            from mdfinfo3 import info3
            channelNameList=info3()
            channelNameList.listChannels(fileName)
        else:
            from mdfinfo4 import info4
            channelNameList=info4()
            channelNameList.listChannels(fileName)
        return channelNameList

class mdf( dict ):
    """ mdf file class
    It imports mdf files version 3.0 and partially above.
    To use : yop= mdfreader.mdf('FileName.dat')
    Some additional useful methods:
    Resample : yop.resample(SamplingRate_in_secs)
    plot a specific channel : yop.plot('ChannelName')
    export to csv file : yop.exportCSV() , specific filename can be input
    export to netcdf : yop.exportNetCDF() """

    def __init__( self, fileName = None ):
        self.fileName = None
        self.timeChannelList = []
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        # clears class from previous reading and avoid to mess up
        self.clear()
        if not fileName == None:
            self.read( fileName )

    ## reads mdf file
    def read( self, fileName = None, multiProc = False, channelList=None):
        # read mdf file
        if self.fileName == None:
            self.fileName = fileName
        self.multiProc = multiProc
        if platform=='win32':
            self.multiProc=False # no multiprocessing for windows platform
        try:
            from multiprocessing import Queue, Process
        except:
            print('No multiprocessing module found')
            self.multiProc = False

        #inttime = time.clock()
        ## Read information block from file
        info = mdfinfo( self.fileName )

        # Open file
        fid = open( self.fileName, 'rb', buffering = 65536 )
        print(self.fileName)

        # Look for the biggest group to process first, to reduce processing time when mutiprocessed
        dataGroupList = dict.fromkeys( list(range( info['HDBlock']['numberOfDataGroups'])) )
        for dataGroup in list(dataGroupList.keys()):
            dataGroupList[dataGroup] = info['CGBlock'][dataGroup][0]['numberOfRecords']
        sortedDataGroup = sorted( dataGroupList, key = dataGroupList.__getitem__, reverse = True )

        if self.multiProc:
            # prepare multiprocessing of dataGroups
            proc = []
            Q=Queue()
        L={}
        ## Defines record format
        for dataGroup in sortedDataGroup:
            # Number for records before and after the data records
            numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
            #Pointer to data block
            pointerToData = info['DGBlock'][dataGroup]['pointerToDataRecords']
            fid.seek( pointerToData )
            # initialize counter of bits
            precedingNumberOfBits = 0
            numpyDataRecordFormat = []
            dataRecordName = []
            if numberOfRecordIDs >= 1:
                # start first with uint8
                numpyDataRecordFormat.append( ( 'firstRecordID', 'uint8' ) )

            for channelGroup in range( info['DGBlock'][dataGroup]['numberOfChannelGroups'] ):
                # Reads size of each channel
                #dataRecordSize = info['CGBlock'][dataGroup][channelGroup]['dataRecordSize']
                # Number of samples recorded
                numberOfRecords = info['CGBlock'][dataGroup][channelGroup]['numberOfRecords']

                if numberOfRecords != 0: # continue if there are at least some records
                    for channel in range( info['CGBlock'][dataGroup][channelGroup]['numberOfChannels'] ):
                        # Type of channel data, signed, unsigned, floating or string
                        signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['signalDataType']
                        # Corresponding number of bits for this format
                        numberOfBits = info['CNBlock'][dataGroup][channelGroup][channel]['numberOfBits']
                        # Defines data format
                        #datatype = self.datatypeformat( signalDataType, numberOfBits )

                        if numberOfBits < 8: # adding bit format, ubit1 or ubit2
                            if precedingNumberOfBits == 8: # 8 bit make a byte
                                dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                                numpyDataRecordFormat.append( ( dataRecordName[-1], self.arrayformat( signalDataType, numberOfBits ) ) )
                                precedingNumberOfBits = 0
                            else:
                                precedingNumberOfBits += numberOfBits # counts successive bits

                        else: # adding bytes
                            if precedingNumberOfBits != 0: # There was bits in previous channel
                                precedingNumberOfBits = 0
                            dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                            numpyDataRecordFormat.append( ( dataRecordName[-1], self.arrayformat( signalDataType, numberOfBits ) ) )

                    if numberOfBits < 8: # last channel in record is bit inside byte unit8
                        dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                        numpyDataRecordFormat.append( ( dataRecordName[-1], self.arrayformat( signalDataType, numberOfBits ) ) )
                        precedingNumberOfBits = 0
                    # Offset of record id at the end of record
                    if numberOfRecordIDs == 2:
                        numpyDataRecordFormat.append( ( 'secondRecordID', 'uint8' ) )
                    #print( 'Group ' + str( dataGroup ) + ' ' + str( time.clock() - inttime ) )
                    ## reads the records
                    buf = numpy.core.records.fromfile( fid, dtype = numpyDataRecordFormat, shape = numberOfRecords , names=dataRecordName)

                    if self.multiProc:
                        proc.append( Process( target = processDataBlocks,
                            args = ( Q, buf, info, numberOfRecords, dataGroup , self.multiProc) ) )
                        proc[-1].start()
                    else: # for debugging purpose, can switch off multiprocessing
                        L.update(processDataBlocks( None, buf, info, numberOfRecords, dataGroup,  self.multiProc))

        fid.close() # close file

        if self.multiProc:
            for i in proc:
                L.update(Q.get()) # concatenate results of processes in dict

        # After all processing of channels,
        # prepare final class data with all its keys
        for dataGroup in range( info['HDBlock']['numberOfDataGroups'] ):
            for channelGroup in range( info['DGBlock'][dataGroup]['numberOfChannelGroups'] ):
                for channel in range( info['CGBlock'][dataGroup][channelGroup]['numberOfChannels'] ):
                    numberOfRecords = info['CGBlock'][dataGroup][channelGroup]['numberOfRecords']
                    if numberOfRecords != 0 :
                        channelName = info['CNBlock'][dataGroup][channelGroup][channel]['signalName']
                        if channelName in channelTime:# assumes first channel is time
                            channelName = 'time' + str( dataGroup )
                            if channelName in L and len( L[channelName] ) != 0:
                                self.timeChannelList.append( channelName )
                        if channelName in L and len( L[channelName] ) != 0:
                            self[channelName] = {}
                            self[channelName]['time'] = 'time' + str( dataGroup )
                            self[channelName]['unit'] = info['CCBlock'][dataGroup][channelGroup][channel]['physicalUnit']
                            self[channelName]['description'] = info['CNBlock'][dataGroup][channelGroup][channel]['signalDescription']
                            self[channelName]['data'] = L[channelName]
        if not channelList==None and len(channelList)>0:
            self.keepChannels(channelList)
        #print( 'Finished in ' + str( time.clock() - inttime ) )
    @staticmethod
    def datatypeformat( signalDataType, numberOfBits ):
        # DATATYPEFORMAT Data type format precision to give to fread
        #   DATATYPEFORMAT(SIGNALDATATYPE,NUMBEROFBITS) is the precision string to
        #   give to fread for reading the data type specified by SIGNALDATATYPE and
        #   NUMBEROFBITS

        if signalDataType in ( 0, 9, 10, 11 ): # unsigned
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
                print(( 'Unsupported number of bits for unsigned int ' + str( dataType ) ))

        elif signalDataType == 1: # signed int
            if numberOfBits == 8:
                dataType = 'b'
            elif numberOfBits == 16:
                dataType = 'h'
            elif numberOfBits == 32:
                dataType = 'i'
            else:
                print(( 'Unsupported number of bits for signed int ' + str( dataType ) ))

        elif signalDataType in ( 2, 3 ): # floating point
            if numberOfBits == 32:
                dataType = 'f'
            elif numberOfBits == 64:
                dataType = 'd'
            else:
                print(( 'Unsupported number of bit for floating point ' + str( dataType ) ))

        elif signalDataType == 7: # string
            dataType = 's'
        elif signalDataType == 8: # array of bytes
            dataType = 'B' # take unit8 as default ?
        else:
            print(( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits ))
        return dataType

    @staticmethod
    def arrayformat( signalDataType, numberOfBits ):

        # Formats used by numpy

        if signalDataType in ( 0, 9, 10, 11 ): # unsigned
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
                print(( 'Unsupported number of bits for unsigned int ' + str( dataType ) ))

        elif signalDataType == 1: # signed int
            if numberOfBits == 8:
                dataType = 'int8'
            elif numberOfBits == 16:
                dataType = 'int16'
            elif numberOfBits == 32:
                dataType = 'int32'
            else:
                print(( 'Unsupported number of bits for signed int ' + str( dataType ) ))

        elif signalDataType in ( 2, 3 ): # floating point
            if numberOfBits == 32:
                dataType = 'float32'
            elif numberOfBits == 64:
                dataType = 'float64'
            else:
                print(( 'Unsupported number of bit for floating point ' + str( dataType ) ))

        elif signalDataType == 7: # string
            dataType = 'str' # not directly processed
        elif signalDataType == 8: # array of bytes
            dataType = 'buffer' # not directly processed
        else:
            print(( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits ))
        return dataType


    def plot( self, channels ):
        try:
            import matplotlib.pyplot as plt
        except:
            print('matplotlib not found' )
            raise
        for channelName in channels:
            if channelName in self:
                if not self[channelName]['data'].dtype.kind in ['S', 'U']: # if channel not a string
                    self.fig = plt.figure()
                    # plot using matplotlib the channel versus time
                    if 'time' in self[channelName]: # Resampled signals
                        timeName = self[channelName]['time']
                        if timeName == '': # resampled channels, only one time channel called 'time'
                            timeName = 'time'
                        if timeName in self: # time channel properly defined
                            plt.plot( self[timeName]['data'], self[channelName]['data'] )
                            plt.xlabel( timeName + ' [' + self[timeName]['unit'] + ']' )
                        else: # no time channel found
                            plt.plot( self[channelName]['data'] )
                    else: # no time signal recognized,
                        if 'time' in self: # most probably resampled
                            plt.plot( self['time']['data'], self[channelName]['data'] )
                            plt.xlabel( 'time [' + self['time']['unit']+ ']' )
                        else:
                            plt.plot( self[channelName]['data'] )

                    plt.title( self[channelName]['description'])
                    if self[channelName]['unit'] == {}:
                        plt.ylabel( channelName )
                    else:
                        plt.ylabel( channelName + ' [' + self[channelName]['unit'] + ']' )
                    plt.grid( True )
                    plt.show()
            else:
                print(( 'Channel ' + channelName + ' not existing' ))

    def allPlot( self ):
        # plot all channels in the object, be careful for test purpose only,
        # can display many many many plots overloading your computer
        for Name in list(self.keys()):
            try:
                self.plot( Name )
            except:
                print( Name )

    def resample( self, samplingTime = 0.1 ):
        """ Resamples mdf channels inside object for one single time"""
        # resample all channels to one sampling time vector
        if 'time' not in self.timeChannelList: # Not yet resampled
            channelNames = list(self.keys())
            minTime = maxTime = []
            self['time'] = {}
            unit = ''
            for time in self.timeChannelList:
                if time in self and len( self[time]['data'] ) > 5: # consider groups having minimum size 
                    minTime.append( self[time]['data'][0] )
                    maxTime.append( self[time]['data'][len( self[time]['data'] ) - 1] )
                    if self[time]['unit'] != '':
                        unit = self[time]['unit']
            self['time']['data'] = numpy.arange( min( minTime ),max( maxTime ),samplingTime )
            self['time']['unit'] = unit
            self['time']['description'] = 'Unique time channel'

            # Interpolate channels
            timevect=[]
            for Name in channelNames:
                try:
                    if Name not in self.timeChannelList:
                        timevect = self[self[Name]['time']]['data']
                        if not self[Name]['data'].dtype.kind in ('S', 'U'): # if channel not array of string
                            self[Name]['data'] = numpy.interp( self['time']['data'], timevect, self[Name]['data'] )
                            if 'time' in self[Name]:
                                del self[Name]['time']
                except:
                    if len( timevect ) != len( self[Name]['data'] ):
                        print(( Name + ' and time channel ' + self[Name]['time'] + ' do not have same length' ))
                    elif not numpy.all( numpy.diff( timevect ) > 0 ):
                        print(( Name + ' has non regularly increasing time channel ' + self[Name]['time'] ))
            # remove time channels in timeChannelList
            for ind in self.timeChannelList:
                del self[ind]
            self.timeChannelList = [] # empty list
            self.timeChannelList.append( 'time' )
        else:
            pass

    def exportToCSV( self, filename = None, sampling = 0.1 ):
        # Export mdf file into csv
        # If no name defined, it will use original mdf name and path
        # By default, sampling is 0.1sec but can be changed
        import csv
        self.resample( sampling )
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.csv'
        if PythonVersion <3:
            f = open( filename, "wb")
        else:
            f = open( filename, "wt" , encoding='latin-1')
        writer = csv.writer( f, dialect = csv.excel )
        # writes header
        writer.writerow( [name for name in list(self.keys()) if self[name]['data'].dtype in ('float64','float32')] ) # writes channel names
        writer.writerow( [(self[name]['unit']) for name in list(self.keys()) if self[name]['data'].dtype in ('float64','float32')] ) # writes units
        # concatenate all channels
        buf = numpy.vstack( [self[name]['data'].transpose() for name in list(self.keys()) if self[name]['data'].dtype in ('float64','float32')] )
        buf = buf.transpose()
        # Write all rows
        r, c = buf.shape
        writer.writerows( [list( buf[i, :] ) for i in range( r )] )
        f.close()

    def exportToNetCDF( self, filename = None, sampling = None ):
        # Export mdf file into netcdf file
        try:
            from Scientific.IO import NetCDF
        except:
            print( 'Scientific.IO module not found' )
            raise
        def cleanName( name ):
            output = name.replace( '$', '' )
            ouput = output.replace( '[', '' )
            ouput = output.replace( ']', '' )
            return output
        if sampling != None:
            self.resample( sampling )
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.nc'
        f = NetCDF.NetCDFFile( filename, 'w' )
        info = mdfinfo( self.fileName ) # info for the metadata
        setattr( f, 'Author',( info['HDBlock']['Author']))
        setattr( f, 'Date', (info['HDBlock']['Date']))
        setattr( f, 'Time', (info['HDBlock']['Time']))
        setattr( f, 'Organization', (info['HDBlock']['Organization']))
        setattr( f, 'ProjectName', (info['HDBlock']['ProjectName']))
        setattr( f, 'Vehicle', (info['HDBlock']['Vehicle']))
        setattr( f, 'Comment', (info['HDBlock']['TXBlock']['Text']))
        # Create dimensions having name of all time channels
        for time in self.timeChannelList:
            f.createDimension( time, len( self[time]['data'] ) )
        # Create variables definition, dimension and attributes
        var = {}
        for name in list(self.keys()):
            if self[name]['data'].dtype == 'float64':
                type = 'd'
            elif self[name]['data'].dtype == 'float32':
                type = 'f'
            elif self[name]['data'].dtype in ['int8', 'int16', 'uint8', 'uint16']:
                type = 'i'
            elif self[name]['data'].dtype in ['int32', 'uint32']:
                type = 'l'
            elif self[name]['data'].dtype.kind in ['S', 'U'] :
                type = 'c'
            else:
                print(( 'Can not process numpy type ' + str(self[name]['data'].dtype) + ' of channel' ))
            # create variable
            CleanedName = cleanName( name )
            if len( self.timeChannelList ) == 1: # mdf resampled
                var[name] = f.createVariable( CleanedName, type, ( self.timeChannelList[0], ) )
            else: # not resampled
                var[name] = f.createVariable( CleanedName, type, ( self[name]['time'], ) )
            # Create attributes
            setattr( var[name], 'title', CleanedName )
            setattr( var[name], 'units', self[name]['unit'])
            setattr( var[name], 'Description', self[name]['description'])
            if name in self.timeChannelList:
                setattr( var[name], 'Type', 'Time Channel' )
                setattr( var[name], 'datatype', 'time' )
            else:
                setattr( var[name], 'Type', 'Data Channel' )
        # put data in variables
        for name in list(self.keys()):
            var[name] = self[name]['data']
        f.close()

    def exportToHDF5( self, filename = None, sampling = None ):
        # export class data structure into hdf5 file
        try:
            import h5py
            import os
        except:
            print( 'h5py not found' )
            raise
        if sampling != None:
            self.resample( sampling )
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.hdf'
        info = mdfinfo( self.fileName ) # info for the metadata
        f = h5py.File( filename, 'w' ) # create hdf5 file
        filegroup=f.create_group(os.path.basename(filename)) # create group in root associated to file
        attr = h5py.AttributeManager( filegroup ) # adds group/file metadata
        attr.create('Author',( info['HDBlock']['Author']))
        attr.create('Date', (info['HDBlock']['Date']))
        attr.create('Time', (info['HDBlock']['Time']))
        attr.create('Organization', (info['HDBlock']['Organization']))
        attr.create('ProjectName', (info['HDBlock']['ProjectName']))
        attr.create('Vehicle', (info['HDBlock']['Vehicle']))
        attr.create('Comment', (info['HDBlock']['TXBlock']['Text']))
        if len( self.timeChannelList ) > 1:
            # if several time groups of channels, not resampled
            groups = {}
            ngroups = 0
            grp = {}
            for channel in list(self.keys()):
                if self[channel]['time'] not in list(groups.keys()):
                    # create new time group
                    ngroups += 1
                    groups[self[channel]['time'] ] = ngroups
                    grp[ngroups] = filegroup.create_group( self[channel]['time'] )
                dset = grp[groups[self[channel]['time'] ]].create_dataset( channel, data = self[channel]['data'] )
                attr = h5py.AttributeManager( dset )
                attr.create( 'unit', (self[channel]['unit'] ))
                attr.create( 'description', (self[channel]['description']))
        else: # resampled or only one time for all channels : no groups
            for channel in list(self.keys()):
                channelName=convertMatlabName(channel)
                dset = filegroup.create_dataset( channelName, data = self[channel]['data'] )
                attr = h5py.AttributeManager( dset )
                attr.create( 'unit', (self[channel]['unit'] ))
                attr.create( 'description', (self[channel]['description'] ))
        f.close()

    def exportToMatlab( self, filename = None ):
        # export class data struture into .mat file
        try:
            from scipy.io import savemat
        except:
            print( 'scipy module not found' )
            raise
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.mat'
        # convert self into simple dict without and metadata
        temp = {}
        for channel in list(self.keys()):
            if not self[channel]['data'].dtype.kind in ('S', 'U'): # does not like special characters chains, skip
                channelName=convertMatlabName(channel)
                temp[channelName] = self[channel]['data']
        try: # depends of version used , compression can be used
            savemat( filename , temp, long_field_names = True,format='5',do_compression=True,oned_as='column' )
        except:
            savemat( filename , temp, long_field_names = True,format='5')

    def exportToExcel( self , filename = None ):
        # export to excel 95 to 2003
        # currently xlwt is not supporting python 3.x
        # finally long to process, to be multiprocessed
        try:
            if PythonVersion<3:
                import xlwt
            else:
                import xlwt3 as xlwt
        except:
            print( 'xlwt module missing' )
            raise
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.xls'
        styleText = xlwt.easyxf( 'font: name Times New Roman, color-index black, bold off' )
        coding='utf-8'
        wb = xlwt.Workbook(encoding=coding)
        channelList = list(self.keys())
        if PythonVersion<3:
            Units=[ self[channel]['unit'].decode(coding, 'replace') for channel in list(self.keys())]
        else:
            Units=[ self[channel]['unit'] for channel in list(self.keys())]
        # Excel 2003 limits
        maxCols = 255
        maxLines = 65535
        workbooknumber = int( ceil( len( channelList ) * 1.0 / ( maxCols * 1.0 ) ) )
        tooLongChannels = []
        # split colmuns in several worksheets if more than 256 cols
        for workbook in range( workbooknumber ):
            ws = wb.add_sheet( 'Sheet' + str( workbook ) ) #, cell_overwrite_ok = True )
            if workbook == workbooknumber - 1: # last sheet
                columnrange = list(range( workbook * maxCols, len( channelList )))
            elif workbook < workbooknumber - 1 and workbooknumber > 1: # first sheets
                columnrange = list(range( workbook * maxCols, ( workbook + 1 ) * maxCols))
            for col in columnrange:
                # write header
                ws.write( 0, col - workbook * maxCols, channelList[col] , styleText )
                ws.write( 1, col - workbook * maxCols, Units[col] , styleText )
                vect = self[channelList[col]]['data'] # data vector
                if not len( vect ) > maxLines :
                    if  vect.dtype.kind not in ['S', 'U']: # if not a string or unicode
                        [ws.row( row + 2 ).set_cell_number( col - workbook * maxCols, vect[row] ) for row in list(range( len( vect ) ))]
                    else: # it's a string, cannot write for the moment
                        if PythonVersion <3:
                            vect=vect.encode(coding)
                        [ws.row( row + 2 ).set_cell_text( col - workbook * maxCols, vect[row]) for row in list(range( len( vect ) ))]
                else: # channel too long, written until max Excel line limit
                    if vect.dtype.kind not in ['S', 'U']: # if not a string
                        [ws.row( row + 2 ).set_cell_number( col - workbook * maxCols, vect[row] ) for row in list(range( maxLines ))]
                    else: # it's a string, cannot write for the moment
                        if PythonVersion <3:
                            vect=vect.encode(coding)
                        [ws.row( row + 2 ).set_cell_text( col - workbook * maxCols, vect[row] ) for row in list(range( maxLines ))]
                    tooLongChannels.append( channelList[col] ) # to later warn user the channel is not completely written
        wb.save( filename ) # writes workbook on HDD
        if len( tooLongChannels ) > 0: # if not empty, some channels have been not processed
            print( 'Following channels were too long to be processed completely, maybe you should resample : ' )
            print( tooLongChannels )

    def exportToXlsx(self, filename=None):
        # export to excel 2007&2010
        # requires openpyxl
        # It is recommended to export resampled data for performances
        try:
            import openpyxl
        except:
            print('Module openpyxl missing')
            raise
        if filename == None:
            filename = self.fileName.replace( '.dat', '' )
            filename = filename.replace( '.DAT', '' )
            filename = filename + '.xlsx'
        channels=list(self.keys())
        maxRows=max([len(self[channel]['data']) for channel in list(self.keys())]) # find max column length
        maxCols=len(list(self.keys())) # number of columns
        print('Creating Excel sheet')
        if len( self.timeChannelList ) > 1: # not resampled data, can be long, writing cell by cell !
            wb=openpyxl.workbook.Workbook(encoding='utf-8')
            ws=wb.get_active_sheet()
            # write header
            if PythonVersion<3:
                for j in range(maxCols):
                    ws.cell(row=0, column=j).value=channels[j].decode('utf-8', 'ignore')
                    ws.cell(row=1, column=j).value=self[channels[j]]['unit'].decode('utf-8', 'ignore')
            else:
                for j in range(maxCols):
                    ws.cell(row=0, column=j).value=channels[j]
                    ws.cell(row=1, column=j).value=self[channels[j]]['unit']               
            for j in range(maxCols):
                print(channels[j])
                if self[channels[j]]['data'].dtype  in ['int8', 'int16', 'uint8', 'uint16']:
                     for r in range(len(self[channels[j]]['data'])):
                        ws.cell(row=r+2, column=j).value=numpy.float64(self[channels[j]]['data'][r])
                else:
                     for r in range(len(self[channels[j]]['data'])):
                        ws.cell(row=r+2, column=j).value=self[channels[j]]['data'][r]
        else: # resampled data
            wb=openpyxl.workbook.Workbook(optimized_write=True, encoding='utf-8')
            ws=wb.create_sheet()
            # write header
            ws.append(channels)
            ws.append([ self[channel]['unit'] for channel in list(self.keys())])
            # write data
            maxRows=max([len(self[channel]['data']) for channel in list(self.keys())]) # find max column length
            maxCols=len(list(self.keys())) # number of columns
            bigmat=numpy.zeros(maxRows) # create empty column
            buf=bigmat
            for col in range(maxCols):
                if not self[channels[col]]['data'].dtype.kind in ['S', 'U']:
                    chanlen=len(self[channels[col]]['data'])
                    if chanlen<maxRows:
                        buf[:]=None
                        buf[0:chanlen]=self[channels[col]]['data']
                        bigmat=numpy.vstack((bigmat, buf))
                    else:
                        bigmat=numpy.vstack((bigmat, self[channels[col]]['data']))
                else:
                    buf[:]=None
                    bigmat=numpy.vstack((bigmat, buf))
            bigmat=numpy.delete(bigmat, 0, 0)
            [ws.append(bigmat[:, row]) for row in range(maxRows)]
        print('Writing file, please wait')
        wb.save(filename)

    def keepChannels(self, channelList):
        # keep only list of channels and removes the rest
        channelList=[channel for channel in channelList]
        removeChannels=[]
        for channel in list(self.keys()):
            if channel not in channelList and not 'time'==channel[0:4] and channel not in self.timeChannelList :
                # avoid to remove time channels otherwise problems with resample
                removeChannels.append(channel)
        if not len(removeChannels)==0:
            [self.pop(channel) for channel in removeChannels]

    def copy(self):
        # copy a mdf class
        yop=mdf()
        yop.multiProc=self.multiProc
        yop.fileName=self.fileName
        yop.timeChannelList=self.timeChannelList
        for channel in list(self.keys()):
            yop[channel]=self[channel]
        return yop

    def mergeMdf(self, mdfClass):
        # merges data of 2 mdf classes
        # both must have been ressampled, otherwise, impossible to know time vector to match
        # create union of both lists
        unionedList=list(mdfClass.keys()) and list(self.keys())
        if 'time' not in self:
            print( 'Data not resampled')
            raise
        initialTimeSize=len(self['time']['data'])
        for channel in unionedList:
            if channel in mdfClass and channel in self: # channel exists in both class
                if not channel=='time':
                    self[channel]['data']=numpy.hstack((self[channel]['data'], mdfClass[channel]['data']))
                else:
                    offset=numpy.mean(numpy.diff(mdfClass[channel]['data'])) # sampling
                    offset=self[channel]['data'][-1]+offset # time offset
                    self[channel]['data']=numpy.hstack((self[channel]['data'], mdfClass[channel]['data']+offset))
            elif channel in mdfClass: # new channel for self from mdfClass
                self[channel]=mdfClass[channel] # initialise all fields, units, descriptions, etc.
                refill=numpy.empty(initialTimeSize)
                refill.fil(numpy.nan) # fill with NANs
                self[channel]['data']=numpy.hstack((refill,  mdfClass[channel]['data'])) # readjust against time
            else: #channel missing in mdfClass
                refill=numpy.empty(len(mdfClass['time']['data']))
                refill.fill(numpy.nan) # fill with NANs
                self[channel]['data']=numpy.hstack((self[channel]['data'], refill))

    def convertToPandas(self, merging=False, sampling=None):
        # convert data structure into pandas module
        try:
            import pandas as pd
        except:
            print('Module pandas missing')
            raise
        if sampling is not None:
            self.resample(sampling)
        info=mdfinfo(self.fileName)
        datetimeInfo=numpy.datetime64(info['HDBlock']['Date'].replace(':','-')+' '+info['HDBlock']['Time'])
        originalKeys=self.keys()   
        # first find groups of channels sharing same time vector
        # initialise frame description
        frame={}
        for group in self.timeChannelList:
            frame[group]=[]
        # create list of channel for each time group
        [frame[self[channel]['time']].append(channel) for channel in self.keys()]
        # create frames
        for group in frame.keys():
            temp={}
            # convert time channel into timedelta
            time=datetimeInfo+numpy.array(self[group]['data']*10E6, dtype='timedelta64[us]')
            for channel in frame[group]:
                temp[channel]=pd.Series(self[channel]['data'], index=time)
            self[group+'_group']=pd.DataFrame(temp)
            self[group+'_group'].pop(group) # delete time channel, no need anymore
        # clean rest of self from data and time channel information
        [self[channel].pop('data') for channel in originalKeys]
        [self[channel].pop('time') for channel in originalKeys]
        self.timeChannelList=[]
        self.timeGroups=[] # save time groups name in list
        [self.timeGroups.append(group+'_group') for group in frame.keys()]

if __name__ == "__main__":
    try:
        from multiprocessing import freeze_support
        freeze_support()
    except:
        None
    mdf()

