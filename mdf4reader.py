# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Thu Dec 10 12:57:28 2013

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
import numpy
from struct import unpack
from sys import platform
from mdfinfo4 import info4, MDFBlock
from collections import OrderedDict
from sys import version_info
from time import gmtime, strftime
PythonVersion=version_info
PythonVersion=PythonVersion[0]

LINK='<Q' 
REAL='<d'
BOOL='<h'
UINT8='<B'
BYTE='<c'
INT16='<h'
UINT16='<H'
UINT32='<I'
INT32='<i'
UINT64='<Q'
INT64='<q'
CHAR='<c'

class DATABlock(MDFBlock):
    def __init__(self, fid,  pointer, numpyDataRecordFormat, numberOfRecords , dataRecordName, zip_type=None, channelList=None):
        # block header
        self.loadHeader(fid, pointer)
        if self['id'] in ('##DT', '##RD', b'##DT', b'##RD'): # normal data block
#            if channelList==None: # reads all blocks
            self['data']=numpy.core.records.fromfile( fid, dtype = numpyDataRecordFormat, shape = numberOfRecords , names=dataRecordName)
#            else: # channelList defined
#                # reads only the channels using offset functions, channel by channel
#                buf={}
#                # order offsets and names based on offset value
#                channelList=OrderedDict(sorted(channelList.items(), key=lambda t: t[1]))
#                for record in range(numberOfRecords):
#                    for name in dataRecordName:
#                        if name in channelList.keys(): 
#                            fid.seek()
#                            buf.append(numpy.core.records.fromfile( fid, dtype = numpyDataRecordFormat, shape = 1 , names=name,  offset=None ))
#                self['data']=buf
                
        elif self['id'] in ('##SD', b'##SD'): # Signal Data Block
            unpack('uint32', fid.read(4)) # length of data
            
        elif self['id'] in ('##DZ', b'##DZ'): # zipped data block
            self['dz_org_block_type']=self.mdfblockreadCHAR(fid, 2)
            self['dz_zip_type']=self.mdfblockread(fid, UINT8, 1)
            if zip_type==None: # HLBlock used
                self['dz_zip_type']=zip_type
            self['dz_reserved']=self.mdfblockreadBYTE(fid, 1)
            self['dz_zip_parameter']=self.mdfblockread(fid, UINT32, 1)
            self['dz_org_data_length']=self.mdfblockread(fid, UINT64, 1)
            self['dz_data_length']=self.mdfblockread(fid, UINT64, 1)
            self['data']=self.mdfblockreadBYTE(fid, self['link_count']-24)
            # uncompress data
            try:
                from zlib import decompress
                if self['dz_zip_type']==0:
                    self['data']=decompress(self['data'])
                elif self['dz_zip_type']==1: # data transposed
                    pass # not yet implemented
            except:
                print('zlib module not found or error while uncompressing')
            if channelList==None: # reads all blocks
                self['data']=numpy.core.records.fromstring(self['data'] , dtype = numpyDataRecordFormat, shape = numberOfRecords , names=dataRecordName)
            else:
                # reads only the channels using offset functions, channel by channel
                buf={}
                # order offsets and names based on offset value
                channelList=OrderedDict(sorted(channelList.items(), key=lambda t: t[1]))
                for name in dataRecordName:
                    if name in channelList.keys():
                        buf.append(numpy.core.records.fromstring(self['data'], dtype = numpyDataRecordFormat, shape = 1 , names=name,  offset=None ))
                self['data']=buf
     
class DATA(MDFBlock):
    def __init__(self, fid,  pointer, numpyDataRecordFormat, numberOfRecords , dataRecordName, zip=None, nameList=None):
        # block header
        self.loadHeader(fid, pointer)
        if not self['id'] in ('##DL', '##HL', b'##DL', b'##HL'):
            self['data']=DATABlock(fid, pointer, numpyDataRecordFormat, numberOfRecords , dataRecordName, zip_type=zip, channelList=nameList)
        elif self['id'] in ('##DL', b'##DL'): # data list block
            # link section
            self['dl_dl_next']=self.mdfblockread(fid, LINK, 1)
            self['dl_data']=self.mdfblockread(fid, LINK, self['link_count']-1)
            # data section
            self['dl_flags']=self.mdfblockread(fid, UINT8, 1)
            self['dl_reserved']=self.mdfblockreadBYTE(fid, 3)
            self['dl_count']=self.mdfblockread(fid, UINT32, 1)
            self['dl_equal_length']=self.mdfblockread(fid, UINT64, 1)
            self['dl_offset']=self.mdfblockread(fid, UINT64, 1)
            # read data list
            if self['dl_count']==1:
                self['data']=DATABlock(fid, self['dl_data'], numpyDataRecordFormat, numberOfRecords , dataRecordName, channelList=nameList)
            else:
                self['data']=numpy.vstack([DATABlock(fid, self['dl_data'][dl], numpyDataRecordFormat, numberOfRecords , dataRecordName, channelList=nameList) for dl in range(self['dl_count'])])
            
        elif self['id'] in ('##HL', b'##HL'): # header list block, if DZBlock used
            # link section
            self['hl_dl_first']=self.mdfblockread(fid, LINK, 1)
            # data section
            self['hl_flags']=self.mdfblockread(fid, UINT16, 1)
            self['hl_zip_type']=self.mdfblockread(fid, UINT8, 1)
            self['hl_reserved']=self.mdfblockreadBYTE(fid, 5)
            self['data']=DATABlock(fid, self['hl_dl_first'], numpyDataRecordFormat, numberOfRecords , dataRecordName, zip_type=self['hl_zip_type'], channelList=nameList)

class mdf4(dict):
    """ mdf file class
    It imports mdf files version 4.1
    To use : yop= mdfreader.mdf('FileName.mf4')"""
    
    def __init__( self, fileName = None, info=None,multiProc = False,  channelList=None):
        self.masterChannelList = {}
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        self.author=''
        self.organisation=''
        self.project=''
        self.subject=''
        self.comment=''
        self.time=''
        self.date=''
        # clears class from previous reading and avoid to mess up
        self.clear()
        if fileName == None and info!=None:
            self.fileName=info.fileName
        else:
            self.fileName=fileName
        self.read4(self.fileName, info, multiProc, channelList)

    ## reads mdf file
    def read4( self, fileName=None, info = None, multiProc = False, channelList=None):
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
            info = info4( self.fileName,  None )

        # reads metadata
        fileDateTime=gmtime(info['HDBlock']['hd_start_time_ns']/1000000000)
        self.date=strftime('%Y-%m-%d', fileDateTime)
        self.time=strftime('%H:%M:%S', fileDateTime)
        if 'Comment' in list(info['HDBlock'].keys()):
            if 'author' in list(info['HDBlock']['Comment'].keys()):
                self.author=info['HDBlock']['Comment']['author']
            if 'department' in list(info['HDBlock']['Comment'].keys()):
                self.organisation=info['HDBlock']['Comment']['department']
            if 'project' in list(info['HDBlock']['Comment'].keys()):
                self.project=info['HDBlock']['Comment']['project']
            if 'subject' in list(info['HDBlock']['Comment'].keys()):
                self.subject=info['HDBlock']['Comment']['subject']
            if 'TX' in list(info['HDBlock']['Comment'].keys()):
                self.comment=info['HDBlock']['Comment']['TX']

        try:
            fid = open( self.fileName, 'rb' )
        except IOError:
            print('Can not find file'+self.fileName)
            raise

        if self.multiProc:
            # prepare multiprocessing of dataGroups
            proc = []
            Q=Queue()
        L={}
        
        masterDataGroup={}  # datagroup name correspondence with its master channel
        for dataGroup in info['DGBlock'].keys():
            if not info['DGBlock'][dataGroup]['dg_data']==0: # data exists
                precedingNumberOfBits = 0
                numpyDataRecordFormat = []
                dataRecordName = []
                #Pointer to data block
                pointerToData =info['DGBlock'][dataGroup]['dg_data']
                # defines record ID if any
                if info['DGBlock'][dataGroup]['dg_rec_id_size']==0: # no record ID
                    pass
                elif info['DGBlock'][dataGroup]['dg_rec_id_size']==1:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint8' ) )
                elif info['DGBlock'][dataGroup]['dg_rec_id_size']==2:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint16' ) )
                elif info['DGBlock'][dataGroup]['dg_rec_id_size']==3:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint16' ) )
                elif info['DGBlock'][dataGroup]['dg_rec_id_size']==4:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint64' ) )
                
                # defines data record for each channel group 
                for channelGroup in info['CGBlock'][dataGroup].keys():
                    # check if this is a VLSD ChannelGroup
                    if info['CGBlock'][dataGroup][channelGroup]['cg_cn_first']==0 or info['CGBlock'][dataGroup][channelGroup]['cg_flags']&1:
                        print('VLSD ChannelGroup, not really implemented yet')
                    numberOfRecords=info['CGBlock'][dataGroup][channelGroup]['cg_cycle_count']
                    #channelListByte=[]
                    #dataBytes=info['CGBlock'][dataGroup][channelGroup]['cg_dataBytes']
                    #invalidBytes=info['CGBlock'][dataGroup][channelGroup]['cg_invalid_bytes']
                    #invalidBit=info['CGBlock'][dataGroup][channelGroup]['cg_invalid_bit_pos']
                    for channel in info['CNBlock'][dataGroup][channelGroup].keys():
                        # Type of channel data, signed, unsigned, floating or string
                        signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['cn_data_type']
                        # Corresponding number of bits for this format
                        numberOfBits = info['CNBlock'][dataGroup][channelGroup][channel]['cn_bit_count']
                        #channelByteOffset= info['CNBlock'][dataGroup][channelGroup][channel]['cn_byte_offset']+info['DGBlock'][dataGroup]['dg_rec_id_size']
                        channelName=info['CNBlock'][dataGroup][channelGroup][channel]['name']
                        # is it a master channel for the group
                        if info['CNBlock'][dataGroup][channelGroup][channel]['cn_type'] in (2, 3, 4): # master channel
                            masterDataGroup[dataGroup]=channelName
                        #if channelList is not None:
                         #   channelListByte[channelName]=channelByteOffset
                            
                        if numberOfBits < 8: # adding bit format, ubit1 or ubit2
                            if precedingNumberOfBits == 8: # 8 bit make a byte
                                dataRecordName.append(info['CNBlock'][dataGroup][channelGroup][channel]['name'])
                                numpyDataRecordFormat.append( ( (dataRecordName[-1], convertName(dataRecordName[-1])), arrayformat4( signalDataType, numberOfBits ) ) )
                                precedingNumberOfBits = 0
                            else:
                                precedingNumberOfBits += numberOfBits # counts successive bits

                        else: # adding bytes
                            if precedingNumberOfBits != 0: # There was bits in previous channel
                                precedingNumberOfBits = 0
                            dataRecordName.append(info['CNBlock'][dataGroup][channelGroup][channel]['name'])
                            numpyDataRecordFormat.append( ( (dataRecordName[-1], convertName(dataRecordName[-1])), arrayformat4( signalDataType, numberOfBits ) ) )

                    if numberOfBits < 8: # last channel in record is bit inside byte unit8
                        dataRecordName.append(info['CNBlock'][dataGroup][channelGroup][channel]['name'])
                        numpyDataRecordFormat.append( ( (dataRecordName[-1], convertName(dataRecordName[-1])), arrayformat4( signalDataType, numberOfBits ) ) )
                        precedingNumberOfBits = 0
                
                # converts channel group records into channels
                buf = DATA(fid, pointerToData, numpyDataRecordFormat, numberOfRecords , dataRecordName, nameList=channelList)
                
                if self.multiProc:
                    proc.append( Process( target = processDataBlocks4, args = ( Q, buf, info, numberOfRecords, dataGroup , self.multiProc) ) )
                    proc[-1].start()
                else: # for debugging purpose, can switch off multiprocessing
                    L.update(processDataBlocks4( None, buf, info, numberOfRecords, dataGroup,  self.multiProc))

        fid.close() # close file

        if self.multiProc:
            for i in proc:
                L.update(Q.get()) # concatenate results of processes in dict

        # After all processing of channels,
        # prepare final class data with all its keys
        for dataGroup in list(info['DGBlock'].keys()):
            self.masterChannelList[masterDataGroup[dataGroup]]=[]
            for channelGroup in list(info['CGBlock'][dataGroup].keys()):
                for channel in list(info['CNBlock'][dataGroup][channelGroup].keys()):
                    channelName = info['CNBlock'][dataGroup][channelGroup][channel]['name']
                    if channelName in L and len( L[channelName] ) != 0:
                        self.masterChannelList[masterDataGroup[dataGroup]].append(channelName)
                        self[channelName] = {}
                        self[channelName]['master'] = masterDataGroup[dataGroup]
                        # sync_types : 0 data, 1 time (sec), 2 angle (radians), 4 index
                        self[channelName]['masterType'] = info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type']
                        if 'unit' in list(info['CCBlock'][dataGroup][channelGroup][channel].keys()):
                            self[channelName]['unit'] = info['CCBlock'][dataGroup][channelGroup][channel]['unit']
                        elif 'unit' in list(info['CNBlock'][dataGroup][channelGroup][channel].keys()):
                            self[channelName]['unit'] = info['CNBlock'][dataGroup][channelGroup][channel]['unit']
                        else:
                            self[channelName]['unit'] = ''
                        if 'Comment' in list(info['CNBlock'][dataGroup][channelGroup][channel].keys()):
                            self[channelName]['description'] = info['CNBlock'][dataGroup][channelGroup][channel]['Comment']
                        else:
                            self[channelName]['description'] = ''
                        self[channelName]['data'] = L[channelName]

        if not channelList==None and len(channelList)>0:
            self.keepChannels(channelList)
        #print( 'Finished in ' + str( time.clock() - inttime ) )

def convertName(channelName):
    if PythonVersion<3:
        channelIdentifier=channelName.encode('ASCII')+'_title'
    else:
        channelIdentifier=str(channelName)+'_title'
    return channelIdentifier

def arrayformat4( signalDataType, numberOfBits ):

    # Formats used by numpy dtype

    if signalDataType in (0, 1): # unsigned
        if numberOfBits <= 8:
            dataType = 'uint8'
        elif numberOfBits == 16:
            dataType = 'uint16';
        elif numberOfBits == 32:
            dataType = 'uint32'
        elif numberOfBits == 64:
            dataType = 'uint64'
        else:
            print(( 'Unsupported number of bits for unsigned int ' + str( numberOfBits) ))

    elif signalDataType in (2, 3): # signed int
        if numberOfBits <= 8:
            dataType = 'int8'
        elif numberOfBits == 16:
            dataType = 'int16'
        elif numberOfBits == 32:
            dataType = 'int32'
        elif numberOfBits == 64:
            dataType = 'int64'
        else:
            print(( 'Unsupported number of bits for signed int ' + str( numberOfBits ) ))

    elif signalDataType in ( 4, 5 ): # floating point
        if numberOfBits == 32:
            dataType = 'float32'
        elif numberOfBits == 64:
            dataType = 'float64'
        else:
            print(( 'Unsupported number of bit for floating point ' + str( numberOfBits ) ))
    
    elif signalDataType == 6: # string ISO-8859-1 Latin
        dataType = 'str' # not directly processed
    elif signalDataType == 7: # UTF-8
        dataType = 'unicode' # not directly processed
    elif signalDataType == (8, 9): # UTF-16
        dataType = 'unicode16' # not directly processed
    elif signalDataType == 10: # bytes array
        dataType = 'buffer' # not directly processed
    elif signalDataType in (11, 12): # MIME sample or MIME stream
        dataType = ('unint16', 'unint8', 'unint8', 'unint8', 'unint8', 'unint8') # 7 Byte Date data structure
    elif signalDataType in (13, 14): # CANopen date and time
        dataType = ('unint32',  'unint16') # 6 Byte time data structure
    else:
        print(( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits ))
        
    # deal with byte order
    if signalDataType in (1, 3, 5, 9):
        dataType='>'+dataType
    
    return dataType

def processDataBlocks4( Q, buf, info, numberOfRecords, dataGroup,  multiProc ):
    ## Processes recorded data blocks
    # Outside of class to allow multiprocessing
    #numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
    #procTime = time.clock()
    #print( 'Start process ' + str( dataGroup ) + ' Number of Channels ' + str( info['CGBlock'][dataGroup][0]['numberOfChannels'] ) + ' of length ' + str( len( buf ) ) + ' ' + str( time.clock() ) )
    L={}
    isBitInUnit8 = 0 # Initialize Bit counter
    previousChannelName=info['CNBlock'][0][0][0]['name'] # initialise channelName

    ## Processes Bits, metadata and conversion
    for channelGroup in info['CGBlock'][dataGroup].keys():
        for channel in info['CNBlock'][dataGroup][channelGroup].keys():
            # Type of channel data, signed, unsigned, floating or string
            signalDataType = info['CNBlock'][dataGroup][channelGroup][channel]['cn_data_type']
            # Corresponding number of bits for this format
            numberOfBits =  info['CNBlock'][dataGroup][channelGroup][channel]['cn_bit_count']

            cName = info['CNBlock'][dataGroup][channelGroup][channel]['name']
            channelName=cName
            try:
                temp = buf['data']['data'].__getattribute__( str(cName)) # extract channel vector
                previousChannelName=cName
            except: # if tempChannelName not in buf -> bits in unit8
                temp = buf['data']['data'].__getattribute__( str(previousChannelName) ) # extract channel vector

            if info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type'] in (2, 3, 4):# master channel 
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
            conversionFormulaIdentifier = info['CCBlock'][dataGroup][channelGroup][channel]['cc_type']
            if conversionFormulaIdentifier == 0: # 1:1 conversion
                pass
            elif conversionFormulaIdentifier == 1: # Parametric, Linear: Physical =Integer*P2 + P1
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [0]
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [1]
                L[channelName] = L[channelName] * P2 + P1
            elif conversionFormulaIdentifier == 2: # rationnal
                P1 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [0]
                P2 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [1]
                P3 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [2]
                P4 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [3]
                P5 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [4]
                P6 = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'] [5]
                L[channelName] = (P1*L[channelName] *L[channelName] +P2*L[channelName] +P3)/(P4*L[channelName] *L[channelName] +P5*L[channelName] +P6)
            elif conversionFormulaIdentifier == 3: #  Algebraic conversion, needs simpy
                pass # not yet implemented
                print(info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']) # formulae
            elif conversionFormulaIdentifier in (4, 5): # value to value table with or without interpolation
                val_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /2
                cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                intVal = cc_val[range(0, val_count, 2)]
                physVal = cc_val[range(1, val_count, 2)]
                if numpy.all( numpy.diff( intVal ) > 0 ):
                    if conversionFormulaIdentifier == 4:
                        L[channelName] = numpy.interp( L[channelName], intVal, physVal ) # with interpolation
                    else:
                        from scipy import interpolate
                        f = interpolate.interp1d( intVal, physVal , kind='nearest', bounds_error=False) # nearest
                        L[channelName] =  f(L[channelName]) # fill with Nan out of bounds while should be bounds
                else:
                    print(( 'X values for interpolation of channel ' + channelName + ' are not increasing' ))
            elif conversionFormulaIdentifier == 6: # Value range to value table
                val_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /3
                cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                key_min = cc_val[range(0, val_count, 3)]
                key_max = cc_val[range(1, val_count, 3)]
                value = cc_val[range(1, val_count, 3)]
                # look up in range keys
                for Lindex in range(len(L[channelName] )):
                    key_index=0 # default index if not found
                    for i in range(val_count):
                        if  key_min[i] < L[channelName] [Lindex] < key_max[i]:
                            key_index=i
                            break
                    L[channelName] [Lindex]=value[key_index]
            elif conversionFormulaIdentifier == 7: # Value to text / scale conversion table
                val_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] 
                cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                for Lindex in range(len(L[channelName] )):
                    key_index=0 # default index if not found
                    for i in range(val_count):
                        if  L[channelName] [Lindex] == cc_val[i]:
                            key_index=i
                            break
                    L[channelName] [Lindex]=value[key_index]
            elif conversionFormulaIdentifier == 8: # Value range to Text / Scale conversion Table
                val_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /2
                cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                key_min = cc_val[range(0, val_count, 2)]
                key_max = cc_val[range(1, val_count, 2)]
                # look up in range keys
                for Lindex in range(len(L[channelName] )):
                    key_index=0 # default index if not found
                    for i in range(val_count):
                        if  key_min[i] < L[channelName] [Lindex] < key_max[i]:
                            key_index=i
                            break
                    L[channelName] [Lindex]=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref'][key_index]
            elif conversionFormulaIdentifier == 9: # Text to value table
                val_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] 
                cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                cc_ref=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
                for Lindex in range(len(L[channelName] )):
                    key_index=0 # default index if not found
                    for i in range(val_count):
                        if  L[channelName] [Lindex] == cc_ref[i]:
                            key_index=i
                            break
                    L[channelName] [Lindex]=cc_val[key_index]
            elif conversionFormulaIdentifier == 10: # Text to text table
                ref_count=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /2
                cc_ref=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
                for Lindex in range(len(L[channelName] )):
                    key_index=0 # default index if not found
                    for i in range(0, ref_count, 2):
                        if  L[channelName] [Lindex] == cc_ref[i]:
                            key_index=i
                            break
                    L[channelName] [Lindex]=cc_ref[key_index+1]
            else:
                print('Unrecognized conversion type')
    if multiProc:
        Q.put(L)
    else:
        return L
        #print( 'Process ' + str( dataGroup ) + ' Finished ' + str( time.clock() - procTime ) )
