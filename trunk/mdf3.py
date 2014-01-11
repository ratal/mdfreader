# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Thu Dec 9 12:57:28 2014

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
import numpy
from math import log, exp
from sys import platform
from mdfinfo3 import info3
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

class mdf3(dict):
    """ mdf file class
    It imports mdf files version 3.0 to 3.3
    To use : yop= mdfreader.mdf('FileName.dat') """
    
    def __init__( self, fileName = None, info=None,multiProc = False,  channelList=None):
        self.timeChannelList = []
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        # clears class from previous reading and avoid to mess up
        self.clear()
        if fileName == None and info!=None:
            self.fileName=info.fileName
        else:
            self.fileName=fileName
        self.read(self.fileName, info, multiProc, channelList)

    ## reads mdf file
    def read( self, fileName=None, info = None, multiProc = False, channelList=None):
        # read mdf file
        self.multiProc = multiProc
        if platform=='win32':
            self.multiProc=False # no multiprocessing for windows platform
        try:
            from multiprocessing import Queue, Process
        except:
            print('No multiprocessing module found')
            self.multiProc = False
        
        if self.fileName == None:
            self.fileName=info.fileName
        else:
            self.fileName=fileName
            
        #inttime = time.clock()
        ## Read information block from file
        if info==None:
            info = info3( self.fileName,  None )

        try:
            fid = open( self.fileName, 'rb' )
        except IOError:
            print('Can not find file'+self.fileName)
            raise

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
