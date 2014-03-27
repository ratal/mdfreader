# -*- coding: utf-8 -*-
""" Measured Data Format blocks parser for version 4.x
Created on Thu Dec 9 12:57:28 2014

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
import numpy
from math import log, exp
from sys import platform
import time
from struct import pack
from mdfinfo3 import info3
channelTime = {'time', 'TIME'}  # possible channel time writing, adjust if necessary


def processDataBlocks(Q, buf, info, numberOfRecords, dataGroup,  multiProc ):
    ## Processes recorded data blocks
    # Outside of class to allow multiprocessing
    #numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
    #procTime = time.clock()
    #print( 'Start process ' + str( dataGroup ) )
    #print( ' Number of Channels ' + str( info['CGBlock'][dataGroup][0]['numberOfChannels'] ) )
    #print( ' of length ' + str( len( buf ) ) + ' ' + str( time.clock() ) )
    L = {}
    isBitInUnit8 = 0  # Initialize Bit counter
    previousChannelName = info['CNBlock'][0][0][0]['signalName']  # initialise channelName

    ## Processes Bits, metadata and conversion
    for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
        for channel in range(info['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):
            # Type of channel data, signed, unsigned, floating or string
            signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['signalDataType']
            # Corresponding number of bits for this format
            numberOfBits = info['CNBlock'][dataGroup][channelGroup][channel]['numberOfBits']

            cName = info['CNBlock'][dataGroup][channelGroup][channel]['signalName']
            channelName=cName
            try:
                temp = buf.__getattribute__( str(cName))  # extract channel vector
                previousChannelName = cName
            except:  # if tempChannelName not in buf -> bits in unit8
                temp = buf.__getattribute__(str(previousChannelName))  # extract channel vector

            if cName in channelTime:  # assumes first channel is time
                channelName = 'time' + str(dataGroup)
            # Process concatenated bits inside unint8
            if signalDataType not in (1, 2, 3, 8):
                if signalDataType == 0:
                    if numberOfBits == 1:  # One bit, considered Boolean in numpy
                        mask = 1 << isBitInUnit8  # masks isBitUnit8
                        temp = [((int(temp[record]) & mask) >> (isBitInUnit8)) for record in range(numberOfRecords)]
                        L[channelName] = numpy.array( temp, dtype='uint8')
                        isBitInUnit8 += 1
                    elif numberOfBits == 2:
                        mask = 1 << isBitInUnit8
                        mask = mask | ( 1 << (isBitInUnit8+1))  # adds second bit in mask
                        temp = [((int(temp[record]) & mask) >> ( isBitInUnit8 )) for record in range(numberOfRecords)]
                        L[channelName] = numpy.array(temp, dtype='uint8')
                        isBitInUnit8 += 2
                    else:
                        L[channelName] = temp
                        isBitInUnit8 = 0  # Initialize Bit counter
                else:
                    L[channelName] = temp
                    isBitInUnit8 = 0  # Initialize Bit counter
            else:
                L[channelName] = temp
                isBitInUnit8 = 0  # Initialize Bit counter

            ## Process Conversion
            conversionFormulaIdentifier = info['CCBlock'][dataGroup][channelGroup][channel]['conversionFormulaIdentifier']
            if conversionFormulaIdentifier == 0:  # Parametric, Linear: Physical =Integer*P2 + P1
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P2']
                L[channelName] = (L[channelName] * P2 + P1)  # .astype(L[channelName].dtype)
            elif conversionFormulaIdentifier == 1:  # Tabular with interpolation
                intVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['int']
                physVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['phys']
                if numpy.all(numpy.diff(intVal) > 0):
                    L[channelName] = numpy.interp(L[channelName], intVal, physVal)
                else:
                    print(( 'X values for interpolation of channel ' + channelName + ' are not increasing'))
            elif conversionFormulaIdentifier == 2:  # Tabular
                intVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['int']
                physVal = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['phys']
                if numpy.all(numpy.diff(intVal) > 0):
                    L[channelName] = numpy.interp(L[channelName], intVal, physVal)
                else:
                    print(('X values for interpolation of channel ' + str(channelName) + ' are not increasing'))
            elif conversionFormulaIdentifier == 6:  # Polynomial
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P6']
                L[channelName] = (P2 - P4 * (L[channelName]['data'] - P5 - P6)) / (P3 * (L[channelName] - P5 - P6) - P1)
            elif conversionFormulaIdentifier == 7:  # Exponential
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P6']
                P7 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P7']
                if P4 == 0 and P1 != 0 and P2 != 0:
                    L[channelName]['data'] = log(((L[channelName] - P7) * P6 - P3) / P1) / P2
                elif P1 == 0 and P4 != 0 and P5 != 0:
                    L[channelName] = log((P3 / (L[channelName] - P7) - P6) / P4) / P5
                else:
                    print(('Non possible conversion parameters for channel ' + str(channelName)))
            elif conversionFormulaIdentifier == 8:  # Logarithmic
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P6']
                P7 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P7']
                if P4 == 0 and P1 != 0 and P2 != 0:
                    L[channelName] = exp(((L[channelName] - P7) * P6 - P3) / P1) / P2
                elif P1 == 0 and P4 != 0 and P5 != 0:
                    L[channelName] = exp((P3 / (L[channelName] - P7) - P6) / P4) / P5
                else:
                    print(('Non possible conversion parameters for channel ' + str(channelName)))
            elif conversionFormulaIdentifier == 9:  # rationnal
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P1']
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P2']
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P3']
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P4']
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P5']
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['P6']
                L[channelName] = (P1*L[channelName] * L[channelName] + P2*L[channelName] + P3)/(P4*L[channelName] * L[channelName] + P5 * L[channelName] + P6)
            elif conversionFormulaIdentifier == 10:  # Text Formula, not really supported yet
                print(('Conversion of formula : ' + str(info['CCBlock'][dataGroup][channelGroup][channel]['conversion']['textFormula']) + 'not yet supported'))
            elif conversionFormulaIdentifier == 11:  # Text Table, not really supported yet
                pass # considers number representative enough, no need to convert into string, more complex
            elif conversionFormulaIdentifier == 12:  # Text Range Table, not really supported yet
                print('Not supported text range tables')
                pass # Not yet supported, practically not used format
            elif conversionFormulaIdentifier == 132:  # date, not really supported yet
                pass
            elif conversionFormulaIdentifier == 133:  # time, not really supported yet
                pass
            elif conversionFormulaIdentifier == 65535:  # No conversion, keep data
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
    
    def __init__( self, fileName=None, info=None,multiProc=False,  channelList=None):
        self.timeChannelList = []
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        self.author=None
        self.organisation=None
        self.project=None
        self.subject=None
        # clears class from previous reading and avoid to mess up
        self.clear()
        if fileName is None and info is not None:
            self.fileName = info.fileName
        else:
            self.fileName = fileName
        self.read(self.fileName, info, multiProc, channelList)

    ## reads mdf file
    def read( self, fileName=None, info=None, multiProc=False, channelList=None):
        # read mdf file
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
            info = info3(self.fileName,  None)
            
        # reads metadata
        self.author=info['HDBlock']['Author']
        self.organisation=info['HDBlock']['Organization']
        self.project=info['HDBlock']['ProjectName']
        self.subject=info['HDBlock']['Subject']

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
        ## Defines record format
        for dataGroup in sortedDataGroup:
            # Number for records before and after the data records
            numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
            #Pointer to data block
            pointerToData = info['DGBlock'][dataGroup]['pointerToDataRecords']
            fid.seek(pointerToData)
            # initialize counter of bits
            precedingNumberOfBits = 0
            numpyDataRecordFormat = []
            dataRecordName = []
            if numberOfRecordIDs >= 1:
                # start first with uint8
                numpyDataRecordFormat.append(('firstRecordID', 'uint8'))

            for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                # Reads size of each channel
                #dataRecordSize = info['CGBlock'][dataGroup][channelGroup]['dataRecordSize']
                # Number of samples recorded
                numberOfRecords = info['CGBlock'][dataGroup][channelGroup]['numberOfRecords']

                if numberOfRecords != 0:  # continue if there are at least some records
                    for channel in range(info['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):
                        # Type of channel data, signed, unsigned, floating or string
                        signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['signalDataType']
                        # Corresponding number of bits for this format
                        numberOfBits = info['CNBlock'][dataGroup][channelGroup][channel]['numberOfBits']
                        # Defines data format
                        #datatype = self.datatypeformat( signalDataType, numberOfBits )

                        if numberOfBits < 8:  # adding bit format, ubit1 or ubit2
                            if precedingNumberOfBits == 8:  # 8 bit make a byte
                                dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                                numpyDataRecordFormat.append((dataRecordName[-1], self.arrayformat(signalDataType, numberOfBits)))
                                precedingNumberOfBits = 0
                            else:
                                precedingNumberOfBits += numberOfBits  # counts successive bits

                        else: # adding bytes
                            if precedingNumberOfBits != 0:  # There was bits in previous channel
                                precedingNumberOfBits = 0
                            dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                            numpyDataRecordFormat.append((dataRecordName[-1], self.arrayformat(signalDataType, numberOfBits )))

                    if numberOfBits < 8: # last channel in record is bit inside byte unit8
                        dataRecordName.append(str(info['CNBlock'][dataGroup][channelGroup][channel]['signalName']))
                        numpyDataRecordFormat.append((dataRecordName[-1], self.arrayformat(signalDataType, numberOfBits)))
                        precedingNumberOfBits = 0
                    # Offset of record id at the end of record
                    if numberOfRecordIDs == 2:
                        numpyDataRecordFormat.append(('secondRecordID', 'uint8'))
                    #print( 'Group ' + str( dataGroup ) + ' ' + str( time.clock() - inttime ) )
                    ## reads the records
                    buf = numpy.core.records.fromfile(fid, dtype=numpyDataRecordFormat, shape=numberOfRecords, names=dataRecordName)

                    if self.multiProc:
                        proc.append(Process(target=processDataBlocks,
                            args=(Q, buf, info, numberOfRecords, dataGroup, self.multiProc)))
                        proc[-1].start()
                    else:  # for debugging purpose, can switch off multiprocessing
                        L.update(processDataBlocks( None, buf, info, numberOfRecords, dataGroup,  self.multiProc))

        fid.close()  # close file

        if self.multiProc:
            for i in proc:
                L.update(Q.get())  # concatenate results of processes in dict

        # After all processing of channels,
        # prepare final class data with all its keys
        for dataGroup in range(info['HDBlock']['numberOfDataGroups']):
            for channelGroup in range(info['DGBlock'][dataGroup]['numberOfChannelGroups']):
                for channel in range(info['CGBlock'][dataGroup][channelGroup]['numberOfChannels']):
                    numberOfRecords = info['CGBlock'][dataGroup][channelGroup]['numberOfRecords']
                    if numberOfRecords != 0 :
                        channelName = info['CNBlock'][dataGroup][channelGroup][channel]['signalName']
                        if channelName in channelTime:  # assumes first channel is time
                            channelName = 'time' + str(dataGroup)
                            if channelName in L and len(L[channelName]) != 0:
                                self.timeChannelList.append(channelName)
                        if channelName in L and len(L[channelName]) != 0:
                            self[channelName] = {}
                            self[channelName]['time'] = 'time' + str(dataGroup)
                            self[channelName]['unit'] = info['CCBlock'][dataGroup][channelGroup][channel]['physicalUnit']
                            self[channelName]['description'] = info['CNBlock'][dataGroup][channelGroup][channel]['signalDescription']
                            self[channelName]['data'] = L[channelName]
        if channelList is not None and len(channelList) > 0:
            self.keepChannels(channelList)
        #print( 'Finished in ' + str( time.clock() - inttime ) )

    def write(self, fileName=None):
        # write data in sorted format: 1 datagroup for 1  channel group and many channel with same sampling
        LINK = 'I'
        CHAR = 'c'
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
        
        # rediscover mdf class structure
        rasters = {}
        for master in self.timeChannelList:
            rasters[master] = []
        for channel in self.keys():
            rasters[self[channel]['time']].append(channel)  # group channels by master channel = datagroup&channelgroup
        pointers = {}  # records pointers of blocks when writing
        
        # writes characters 
        def writeChar(f, value, size=None):
            if size is None:
                temp = value
            else:
                if len(value) > size:
                    temp = value[:size]
                else:
                    temp = value+' '*(size-len(value))
                temp += '\0'
            f.write(pack('<'+CHAR*len(temp), *temp))

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
        fid.write(pack(UINT16, len(self.timeChannelList)))  # number of data groups
        writeChar(fid, time.strftime("%d:%m:%Y"))  # date
        writeChar(fid, time.strftime("%H:%M:%S"))  # time
        if self.author is not None:
            writeChar(fid, self.author,  size=31)  # Author
        else:
            writeChar(fid, ' ',  size=31)  # Author
        if self.organisation is not None:
            writeChar(fid, self.organisation,  size=31)  # Organisation
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
        fid.write(pack(UINT64, int(time.time()*1000000000)))  # Time Stamp
        fid.write(pack(INT16, 1))  # UTC time offset
        fid.write(pack(UINT16, 0))  # Time quality
        writeChar(fid, 'Local PC Reference Time         ')  # Timer identification
        
        # write DG block
        pointers['DG'] = {}
        pointers['CG'] = {}
        pointers['CN'] = {}
        ndataGroup = len(self.timeChannelList)
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
            masterChannel = self.timeChannelList[dataGroup]
            numChannels = len(rasters[masterChannel])
            fid.write(pack(UINT16, numChannels))  # Number of channels
            pointers['CG'][dataGroup]['dataRecordSize'] = fid.tell()
            fid.write(pack(UINT16, 0))  # Size of data record
            nRecords = len(self[self.timeChannelList[dataGroup]]['data'])
            fid.write(pack(UINT32, nRecords))  # Number of records
            fid.write(pack(LINK, 0))  # pointer to sample reduction block, not used
            sampling = numpy.average(numpy.diff(self[self.timeChannelList[dataGroup]]['data']))            
            
            # Channel blocks writing
            pointers['CN'][dataGroup] = {}
            dataList = ()
            dataTypeList = ''
            recordNumberOfBits = 0
            preceedingChannel = None
            writePointer(fid, pointers['CG'][dataGroup]['firstCN'], fid.tell())  # first channel bock pointer from CG
            for channel in rasters[masterChannel]:
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
                if channel not in self.timeChannelList:
                    fid.write(pack(UINT16, 0))  # data channel
                else:
                    fid.write(pack(UINT16, 1))  # master channel
                # make channel name in 32 bytes
                writeChar(fid, channel, size=31)  # channel name
                # channel description
                desc = self[channel]['description']
                writeChar(fid, desc, size=127)  # channel description
                fid.write(pack(UINT16, 0))  # no offset
                data = self[channel]['data']  # channel data
                dataList = dataList+(data, )

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
                    dataTypeList += 's'
                fid.write(pack(UINT16, dataType))  # Signal data type
                if not data.dtype.kind in ['S', 'U']:
                    fid.write(pack(BOOL, 1))  # Value range valid
                    fid.write(pack(REAL, min(data)))  # Min value
                    fid.write(pack(REAL, max(data)))  # Max value
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
                    fid.write(pack(REAL, min(data)))  # Min value
                    fid.write(pack(REAL, max(data)))  # Max value
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
            fid.write(pack(UINT16, recordNumberOfBits/8))  # Size of data record
            fid.seek(currentPosition)
            
            # data writing
            # write data pointer in datagroup
            writePointer(fid, pointers['DG'][dataGroup]['data'], fid.tell())
            records = numpy.array(dataList, object).T
            records = numpy.reshape(records, (1, len(dataTypeList)*nRecords), order='C')[0]  # flatten the matrix
            fid.write(pack('<'+dataTypeList*nRecords, *records))  # dumps data vector from numpy
                
        #print(pointers)
        fid.close()
        
    @staticmethod
    def datatypeformat(signalDataType, numberOfBits):
        # DATATYPEFORMAT Data type format precision to give to fread
        #   DATATYPEFORMAT(SIGNALDATATYPE,NUMBEROFBITS) is the precision string to
        #   give to fread for reading the data type specified by SIGNALDATATYPE and
        #   NUMBEROFBITS

        if signalDataType in (0, 9, 10, 11):  # unsigned
            if numberOfBits == 8:
                dataType = 'B'
            elif numberOfBits == 16:
                dataType = 'H'
            elif numberOfBits == 32:
                dataType = 'I'
            elif numberOfBits == 1:
                dataType = 'ubit1'  # not directly processed
            elif numberOfBits == 2:
                dataType = 'ubit2'  # not directly processed
            else:
                print(('Unsupported number of bits for unsigned int ' + str(signalDataType)))

        elif signalDataType == 1:  # signed int
            if numberOfBits == 8:
                dataType = 'b'
            elif numberOfBits == 16:
                dataType = 'h'
            elif numberOfBits == 32:
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

    @staticmethod
    def arrayformat(signalDataType, numberOfBits):

        # Formats used by numpy

        if signalDataType in (0, 9, 10, 11):  # unsigned
            if numberOfBits == 8:
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
            if numberOfBits == 8:
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
