# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Thu Dec 10 12:57:28 2013

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
import numpy
from struct import unpack, Struct
from math import pow
from sys import platform, version_info
from mdfinfo4 import info4, MDFBlock#,  ATBlock, CNBlock
from collections import OrderedDict
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
    def __init__(self, fid,  pointer, record, zip_type=None, channelList=None):
        # block header
        self.loadHeader(fid, pointer)
        if self['id'] in ('##DT', '##RD', b'##DT', b'##RD'): # normal data block
            record.numberOfRecords = (self['length']-24) // record.recordLength
            self['data']=record.readSortedRecord(fid, pointer, channelList)
                
        elif self['id'] in ('##SD', b'##SD'): # Signal Data Block
            unpack('uint32', fid.read(4)) # length of data
            
        elif self['id'] in ('##DZ', b'##DZ'): # zipped data block
            self['dz_org_block_type']=self.mdfblockreadCHAR(fid, 2)
            self['dz_zip_type']=self.mdfblockread(fid, UINT8, 1)
            if zip_type is not None: # HLBlock used
                self['dz_zip_type']=zip_type
            self['dz_reserved']=self.mdfblockreadBYTE(fid, 1)
            self['dz_zip_parameter']=self.mdfblockread(fid, UINT32, 1)
            self['dz_org_data_length']=self.mdfblockread(fid, UINT64, 1)
            record.numberOfRecords = self['dz_org_data_length'] // record.recordLength
            self['dz_data_length']=self.mdfblockread(fid, UINT64, 1)
            self['data']=fid.read( self['dz_data_length'] )
            # uncompress data
            try:
                from zlib import decompress
            except:
                raise('zlib module not found or error while uncompressing')
            self['data']=decompress(self['data']) # decompress data
            if self['dz_zip_type']==1: # data bytes transposed
                N = self['dz_zip_parameter']
                M = self['dz_org_data_length']%N
                if PythonVersion<3:
                    temp=numpy.frombuffer(self['data'][:M*N])
                    temp.reshape(M, N).T.reshape(1, M*N)
                    self['data']=numpy.getbuffer(temp)+self['data'][M*N:]
                else:
                    temp=numpy.asarray(self['data'][:M*N])
                    temp.reshape(M, N).T.reshape(1, M*N)
                    self['data']=numpy.getbuffer(temp)+self['data'][M*N:]
            if channelList is None: # reads all blocks
                self['data']=numpy.core.records.fromstring(self['data'] , dtype = record.numpyDataRecordFormat, shape = record.numberOfRecords , names=record.dataRecordName)
            else:
                # reads only the channels using offset functions, channel by channel. Not yet ready
                buf={}
                # order offsets and names based on offset value
                channelList=OrderedDict(sorted(channelList.items(), key=lambda t: t[1]))
                for name in record.dataRecordName:
                    if name in list(channelList.keys()):
                        buf.append(numpy.core.records.fromstring(self['data'], dtype = record.numpyDataRecordFormat, shape = 1 , names=name,  offset=None ))
                self['data']=buf
     
class DATA(dict):
    def __init__(self, fid, pointer):
        self.fid=fid
        self.pointerTodata=pointer
    def addRecord(self, record):
        self[record.recordID]={}
        self[record.recordID]['record']=record
    def read(self, channelList, zip=None):
        if len(self)==1: #sorted dataGroup
            recordID=list(self.keys())[0]
            self[recordID]['data']=self.loadSorted( self[recordID]['record'], zip=None, nameList=channelList)
        else: # unsorted DataGroup
            self.load( zip=None, nameList=channelList)
            
    def loadSorted(self, record, zip=None, nameList=None): # reads sorted data
        temps=MDFBlock()
        # block header
        temps.loadHeader(self.fid, self.pointerTodata)
        if temps['id'] in ('##DL', b'##DL'): # data list block
            # link section
            temps['dl_dl_next']=temps.mdfblockread(self.fid, LINK, 1)
            temps['dl_data']=[temps.mdfblockread(self.fid, LINK, 1) for Link in range(temps['link_count']-1)]
            # data section
            temps['dl_flags']=temps.mdfblockread(self.fid, UINT8, 1)
            temps['dl_reserved']=temps.mdfblockreadBYTE(self.fid, 3)
            temps['dl_count']=temps.mdfblockread(self.fid, UINT32, 1)
            if temps['dl_count']:
                if temps['dl_flags']: # equal length datalist
                    temps['dl_equal_length']=temps.mdfblockread(self.fid, UINT64, 1)
                else: # datalist defined by byte offset
                    temps['dl_offset']=temps.mdfblockread(self.fid, UINT64, temps['dl_count'])
                # read data list
                if temps['dl_count']==1:
                    temps['data']=DATABlock(self.fid, temps['dl_data'][0], record, channelList=nameList)
                elif temps['dl_count']==0:
                    raise('empty datalist')
                else:
                    # define size of each block to be read
                    temps['data']={}
                    temps['data']['data']=[]
                    temps['data']['data']=[DATABlock(self.fid, temps['dl_data'][dl], record, zip, channelList=nameList)['data'] for dl in range(temps['dl_count'])]
                    if temps['dl_dl_next']:
                         self.pointerTodata=temps['dl_dl_next']
                         temps['data']['data'].append(self.loadSorted(record, zip, nameList)['data'])
                    #Concatenate all data
                    temps['data']['data']=numpy.hstack(temps['data']['data']) # concatenate data list
                    temps['data']['data']=temps['data']['data'].view(numpy.recarray) # vstack output ndarray instead of recarray
            else: # empty datalist
                temps['data']=None
        elif temps['id'] in ('##HL', b'##HL'): # header list block for DZBlock 
            # link section
            temps['hl_dl_first']=temps.mdfblockread(self.fid, LINK, 1)
            # data section
            temps['hl_flags']=temps.mdfblockread(self.fid, UINT16, 1)
            temps['hl_zip_type']=temps.mdfblockread(self.fid, UINT8, 1)
            temps['hl_reserved']=temps.mdfblockreadBYTE(self.fid, 5)
            self.pointerTodata= temps['hl_dl_first']
            temps['data']=self.loadSorted(record, zip=temps['hl_zip_type'], nameList=nameList)#, record, zip=temps['hl_zip_type'], nameList=nameList)
        else:
            temps['data']=DATABlock(self.fid, self.pointerTodata, record, zip_type=zip, channelList=nameList)
        return temps['data']
    
    def load(self, record, zip=None, nameList=None):
        # read unsorted records in Datablock
        # identify record channels in buffer
        # iterate in buffer
        # identify record ID
        # reads the channels from record
        temps=MDFBlock()
        # block header
        temps.loadHeader(self.fid, self.pointerTodata)
        if temps['id'] in ('##DL', b'##DL'): # data list block
            # link section
            temps['dl_dl_next']=temps.mdfblockread(self.fid, LINK, 1)
            temps['dl_data']=[temps.mdfblockread(self.fid, LINK, 1) for Link in range(temps['link_count']-1)]
            # data section
            temps['dl_flags']=temps.mdfblockread(self.fid, UINT8, 1)
            temps['dl_reserved']=temps.mdfblockreadBYTE(self.fid, 3)
            temps['dl_count']=temps.mdfblockread(self.fid, UINT32, 1)
            if temps['dl_count']:
                if temps['dl_flags']: # equal length datalist
                    temps['dl_equal_length']=temps.mdfblockread(self.fid, UINT64, 1)
                else: # datalist defined by byte offset
                    temps['dl_offset']=temps.mdfblockread(self.fid, UINT64, temps['dl_count'])
                # read data list
                if temps['dl_count']==1:
                    temps['data']=DATABlock(self.fid, temps['dl_data'][0], record, nameList)
                elif temps['dl_count']==0:
                    raise('empty datalist')
                else:
                    # define size of each block to be read
                    temps['data']={}
                    temps['data']['data']=[]
                    temps['data']['data']=[DATABlock(self.fid, temps['dl_data'][dl], record, zip, nameList)['data'] for dl in range(temps['dl_count'])]
                    if temps['dl_dl_next']:
                         self.pointerTodata=temps['dl_dl_next']
                         temps['data']['data'].append(self.load(record, zip, nameList)['data'])
                    #Concatenate all data
                    temps['data']['data']=numpy.hstack(temps['data']['data']) # concatenate data list
                    temps['data']['data']=temps['data']['data'].view(numpy.recarray) # vstack output ndarray instead of recarray
            else: # empty datalist
                temps['data']=None
        elif temps['id'] in ('##HL', b'##HL'): # header list block for DZBlock 
            # link section
            temps['hl_dl_first']=temps.mdfblockread(self.fid, LINK, 1)
            # data section
            temps['hl_flags']=temps.mdfblockread(self.fid, UINT16, 1)
            temps['hl_zip_type']=temps.mdfblockread(self.fid, UINT8, 1)
            temps['hl_reserved']=temps.mdfblockreadBYTE(self.fid, 5)
            self.pointerTodata= temps['hl_dl_first']
            temps['data']=self.load(record, zip=temps['hl_zip_type'], nameList=nameList)#, record, zip=temps['hl_zip_type'], nameList=nameList)
        else:
            temps['data']=DATABlock(self.fid, self.pointerTodata, record, zip_type=zip, channelList=nameList)
        return temps['data']
        
class recordChannel():
    def __init__(self, info, dataGroup, channelGroup, channelNumber, recordIDsize):
        self.name=info['CNBlock'][dataGroup][channelGroup][channelNumber]['name']
        self.channelNumber=channelNumber
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data_type']
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_count']
        self.dataFormat=arrayformat4( self.signalDataType, self.bitCount )
        if not self.signalDataType in (13, 14): # processed by several channels
            self.CFormat=Struct(datatypeformat4( self.signalDataType, self.bitCount ))
        if self.bitCount%8==0:
            self.nBytes=self.bitCount // 8
        else:
            self.nBytes=self.bitCount // 8 + 1
        self.bitOffset=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_offset']
        self.byteOffset=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_byte_offset']
        self.RecordFormat=((self.name, convertName(self.name)),  self.dataFormat)
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_type']
        if self.channelType in (1, 4, 5): # if VSLD or sync or max length channel
            self.data=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data']
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
        self.VLSDFlag=0
        self.recordToChannelMatching={}
        self.channelNames=[]
    def addChannel(self, info, channelNumber):
        self.append(recordChannel(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize))
        self.channelNames.append(self[-1].name)
    def loadInfo(self, info):
        # gathers records related from info class
        self.recordIDsize=info['DGBlock'][self.dataGroup]['dg_rec_id_size']
        if not self.recordIDsize==0: # no record ID
            self.dataRecordName.append('RecordID'+str(self.channelGroup))
            format=(self.dataRecordName[-1], convertName(self.dataRecordName[-1]))
            if self.recordIDsize==1:
                self.numpyDataRecordFormat.append( ( format, 'uint8' ) )
            elif self.recordIDsize==2:
                self.numpyDataRecordFormat.append( ( format, 'uint16' ) )
            elif self.recordIDsize==3:
                self.numpyDataRecordFormat.append( ( format, 'uint32' ) )
            elif self.recordIDsize==4:
                self.numpyDataRecordFormat.append( ( format, 'uint64' ) )
        self.recordID=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_record_id']
        self.recordLength=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        self.recordLength+= info['CGBlock'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
        self.numberOfRecords=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        self.VLSDFlag=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_flags']
        for channelNumber in list(info['CNBlock'][self.dataGroup][self.channelGroup].keys()):
            channel=recordChannel(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize)
            if channel.channelType in (2, 3): # master channel found
                self.master['name']=channel.name
                self.master['number']=channelNumber
            if channel.channelType in (0, 2): # not virtual channel
                self.append(channel)
                self.channelNames.append(channel.name)
                if len(self)>1 and channel.byteOffset==self[-2].byteOffset: # several channels in one byte, ubit1 or ubit2
                    self.recordToChannelMatching[channel.name]=self.recordToChannelMatching[self[-2].name]
                else: # adding bytes
                    self.recordToChannelMatching[channel.name]=channel.name
                    if channel.signalDataType not in (13, 14):
                        self.numpyDataRecordFormat.append(channel.RecordFormat)
                        self.dataRecordName.append(channel.name)
                    elif channel.signalDataType == 13:
                        self.dataRecordName.append( 'ms' )
                        self.numpyDataRecordFormat.append( ( ('ms', 'ms_title'), '<u2') ) 
                        self.dataRecordName.append( 'min' )
                        self.numpyDataRecordFormat.append( ( ('min', 'min_title'), '<u1') ) 
                        self.dataRecordName.append( 'hour' )
                        self.numpyDataRecordFormat.append( ( ('hour', 'hour_title'), '<u1') ) 
                        self.dataRecordName.append( 'day' )
                        self.numpyDataRecordFormat.append( ( ('day', 'day_title'), '<u1') ) 
                        self.dataRecordName.append( 'month' )
                        self.numpyDataRecordFormat.append( ( ('month', 'month_title'), '<u1') )
                        self.dataRecordName.append( 'year' )
                        self.numpyDataRecordFormat.append( ( ('year', 'year_title'), '<u1') )
                    elif channel.signalDataType == 14:
                        self.dataRecordName.append( 'ms' )
                        self.numpyDataRecordFormat.append( ( ('ms', 'ms_title'), '<u4') ) 
                        self.dataRecordName.append( 'days' )
                        self.numpyDataRecordFormat.append( ( ('days', 'days_title'), '<u2') )
            elif channel.channelType in (3, 6): # virtual channel
                pass # channel calculated based on record index later in conversion function

    def readSortedRecord(self, fid, pointer, channelList=None):
        # reads record, only one channel group per datagroup
        fid.seek(pointer)
        if channelList is None: # reads all, quickest but memory consuming
            return numpy.core.records.fromfile( fid, dtype = self.numpyDataRecordFormat, shape = self.numberOfRecords, names=self.dataRecordName)
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
                rec=numpy.zeros((self.numberOfRecords, ), dtype=numpyDataRecordFormat)
                recordLength=self.recordIDsize+self.recordLength
                for r in range(self.numberOfRecords): # for each record,
                    buf=fid.read(recordLength)
                    for channel in recChan:
                        rec[channel.name][r]=channel.CFormat.unpack(buf[channel.posBeg:channel.posEnd])[0]
                return rec.view(numpy.recarray)
            
    def readUnsortedRecord(self, buf, channelList=None):
        # read stream of record bytes
        pass
        
    def readVLSDCGBlock(self, fid, cg_data_bytes):
        VLSDLength=unpack('<u4', fid.read(4))
        return fid.read(VLSDLength)
            
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
        if platform in ('win32', 'win64'):
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
                #Pointer to data block
                pointerToData =info['DGBlock'][dataGroup]['dg_data']
                buf=DATA(fid, pointerToData)
                
                for channelGroup in info['CGBlock'][dataGroup].keys():
                    temp=record(dataGroup, channelGroup) # create record class
                    temp.loadInfo(info) # load all info related to record
                    buf.addRecord(temp)
                    for channel in info['CNBlock'][dataGroup][channelGroup].keys():
                        if info['CNBlock'][dataGroup][channelGroup][channel]['cn_type'] in (2, 3, 4):
                            masterDataGroup[dataGroup]=info['CNBlock'][dataGroup][channelGroup][channel]['name']

                buf.read(channelList)

                # Convert channels to physical values
                OkBuf=len(buf)>0 and 'data' in buf[list(buf.keys())[0]] and buf[list(buf.keys())[0]]['data'] is not None
                if self.multiProc and OkBuf: 
                    proc.append( Process( target = processDataBlocks4, args = ( Q, buf, info, dataGroup, channelList, self.multiProc) ) )
                    proc[-1].start()
                elif OkBuf: # for debugging purpose, can switch off multiprocessing
                    L.update(processDataBlocks4( None, buf, info, dataGroup, channelList, self.multiProc))
                else:
                    print('No data in dataGroup '+ str(dataGroup))

        fid.close() # close file

        if self.multiProc and 'data' in buf:
            for i in proc:
                L.update(Q.get()) # concatenate results of processes in dict

        # After all processing of channels,
        # prepare final class data with all its keys
        for dataGroup in list(info['DGBlock'].keys()):
            if masterDataGroup: #master channel exist
                self.masterChannelList[masterDataGroup[dataGroup]]=[]
            for channelGroup in list(info['CGBlock'][dataGroup].keys()):
                for channel in list(info['CNBlock'][dataGroup][channelGroup].keys()):
                    channelName = info['CNBlock'][dataGroup][channelGroup][channel]['name']
                    if channelName in L and len( L[channelName] ) != 0:
                        if masterDataGroup: #master channel exist
                            self.masterChannelList[masterDataGroup[dataGroup]].append(channelName)
                        self[channelName] = {}
                        if masterDataGroup: #master channel exist
                            self[channelName]['master'] = masterDataGroup[dataGroup]
                        else:
                            self[channelName]['master'] = ''
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
                    elif info['CNBlock'][dataGroup][channelGroup][channel]['cn_data_type']==13: #CANopen date
                        identifier=['ms', 'min', 'hour', 'day', 'month', 'year']
                        for name in identifier:
                            self[name]={}
                            self[name]['data']=L[name]
                            self[name]['unit']=name
                            self[name]['description']=''
                            self[name]['masterType']=info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type']
                            if masterDataGroup: #master channel exist
                                self.masterChannelList[masterDataGroup[dataGroup]].append(name)
                    elif info['CNBlock'][dataGroup][channelGroup][channel]['cn_data_type']==14: #CANopen time
                        identifier=['ms', 'days']
                        for name in identifier:
                            self[name]={}
                            self[name]['data']=L[name]
                            self[name]['unit']=name
                            self[name]['description']=''
                            self[name]['masterType']=info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type']
                            if masterDataGroup: #master channel exist
                                self.masterChannelList[masterDataGroup[dataGroup]].append(name)
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
            dataType = 'u1'
        elif numberOfBits <= 16:
            dataType = 'u2';
        elif numberOfBits <= 32:
            dataType = 'u4'
        elif numberOfBits <= 64:
            dataType = 'u8'
        else:
            print(( 'Unsupported number of bits for unsigned int ' + str( numberOfBits) ))

    elif signalDataType in (2, 3): # signed int
        if numberOfBits <= 8:
            dataType = 'i1'
        elif numberOfBits <= 16:
            dataType = 'i2'
        elif numberOfBits <= 32:
            dataType = 'i4'
        elif numberOfBits <= 64:
            dataType = 'i8'
        else:
            print(( 'Unsupported number of bits for signed int ' + str( numberOfBits ) ))

    elif signalDataType in ( 4, 5 ): # floating point
        if numberOfBits == 32:
            dataType = 'f4'
        elif numberOfBits == 64:
            dataType = 'f8'
        else:
            print(( 'Unsupported number of bit for floating point ' + str( numberOfBits ) ))
    
    elif signalDataType == 6: # string ISO-8859-1 Latin
        dataType = 'str' # not directly processed
    elif signalDataType == 7: # UTF-8
        dataType = 'U2' # not directly processed
    elif signalDataType in (8, 9): # UTF-16
        dataType = 'U4' # not directly processed
    elif signalDataType == 10: # bytes array
        dataType = 'V'+str(int(numberOfBits/8)) # not directly processed
    elif signalDataType in (11, 12): # MIME sample or MIME stream
        dataType = 'V'+str(int(numberOfBits/8)) # not directly processed
    elif signalDataType in (13, 14): # CANOpen date or time
        dataType=None
    else:
        print(( 'Unsupported Signal Data Type ' + str( signalDataType ) + ' ', numberOfBits ))
        
    # deal with byte order
    if signalDataType in (1, 3, 5, 9): # by default low endian, but here big endian
        dataType='>'+dataType
    
    return dataType
    
def datatypeformat4(signalDataType, numberOfBits):
    # DATATYPEFORMAT Data type format precision to give to fread
    #   datatypeformat4(SIGNALDATATYPE,NUMBEROFBITS) is the precision string to
    #   give to fread for reading the data type specified by SIGNALDATATYPE and
    #   NUMBEROFBITS

    if signalDataType in (0, 1):  # unsigned int
        if numberOfBits <= 8:
            dataType = 'B'
        elif numberOfBits <= 16:
            dataType = 'H'
        elif numberOfBits <= 32:
            dataType = 'I'
        elif numberOfBits <= 64:
            dataType = 'L'  
        else:
            print(('Unsupported number of bits for unsigned int ' + str(signalDataType)))

    elif signalDataType in (2, 3):  # signed int
        if numberOfBits <= 8:
            dataType = 'b'
        elif numberOfBits <= 16:
            dataType = 'h'
        elif numberOfBits <= 32:
            dataType = 'i'
        elif numberOfBits <= 64:
            dataType = 'l'
        else:
            print(('Unsupported number of bits for signed int ' + str(signalDataType)))

    elif signalDataType in (4, 5):  # floating point
        if numberOfBits == 32:
            dataType = 'f'
        elif numberOfBits == 64:
            dataType = 'd'
        else:
            print(('Unsupported number of bit for floating point ' + str(signalDataType)))

    elif signalDataType in (6, 7, 8, 9, 10, 11, 12):  # string/bytes
        dataType = 's'
    else:
        print(('Unsupported Signal Data Type ' + str(signalDataType) + ' ', numberOfBits))
    # deal with byte order
    if signalDataType in (1, 3, 5, 9): # by default low endian, but here big endian
        dataType='>'+dataType
    
    return dataType
    
def processDataBlocks4( Q, buf, info, dataGroup,  channelList, multiProc ):
    ## Processes recorded data blocks
    # Outside of class to allow multiprocessing
    #numberOfRecordIDs = info['DGBlock'][dataGroup]['numberOfRecordIDs']
    #procTime = time.clock()
    #print( 'Start process ' + str( dataGroup ) + ' Number of Channels ' + str( info['CGBlock'][dataGroup][0]['numberOfChannels'] ) + ' of length ' + str( len( buf ) ) + ' ' + str( time.clock() ) )
    L={}

    if channelList is None:
        allChannel=True
    else:
        allChannel=False
    ## Processes Bits, metadata and conversion
    for recordID in buf.keys():
        channelGroup=buf[recordID]['record'].channelGroup
        for chan in buf[recordID]['record']:
            channelName=chan.name
            if (allChannel or channelName in channelList) and chan.signalDataType not in (13, 14):
                channel=chan.channelNumber
                recordName=buf[recordID]['record'].recordToChannelMatching[channelName] # in case record is used for several channels
                if not chan.channelType in (3, 6): # not virtual channel
                    temp = buf[recordID]['data']['data'].__getattribute__( str(recordName)+'_title') # extract channel vector
                else :# virtual channel 
                    temp=numpy.arange(buf[recordID]['record'].numberOfRecords)
                    
                if info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type'] in (2, 3, 4): # master channel 
                    channelName = 'master' + str( dataGroup )

                # Process concatenated bits inside uint8
                if not chan.bitCount//8.0==chan.bitCount/8.0: # if channel data do not use complete bytes
                    mask = int(pow(2, chan.bitCount+1)-1) # masks isBitUnit8
                    if chan.signalDataType in (0,1, 2, 3): # integers
                        temp =  numpy.right_shift(temp,  chan.bitOffset)
                        temp =  numpy.bitwise_and(temp,  mask )
                        L[channelName] = temp
                    else: # should not happen
                        print('bit count and offset not applied to correct data type')
                        L[channelName] = temp
                else: #data using bull bytes
                    L[channelName] = temp

                ## Process Conversion
                conversionFormulaIdentifier = info['CCBlock'][dataGroup][channelGroup][channel]['cc_type']
                if conversionFormulaIdentifier == 0: # 1:1 conversion
                    pass
                elif conversionFormulaIdentifier == 1: # Parametric, Linear: Physical =Integer*P2 + P1
                    L[channelName] = linearConv(L[channelName], info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'])
                elif conversionFormulaIdentifier == 2: # rationnal
                    L[channelName] = rationalConv(L[channelName], info['CCBlock'][dataGroup][channelGroup][channel]['cc_val'])
                elif conversionFormulaIdentifier == 3: #  Algebraic conversion, needs simpy
                    try:
                        L[channelName] = formulaConv(L[channelName], info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']['Comment'])
                    except:
                        print('Please install sympy to convert channel '+channelName)
                        print('Had problem to convert '+channelName+' formulae '+info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref'])
                elif conversionFormulaIdentifier in (4, 5): # value to value table with or without interpolation
                    val_count=2*int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /2)
                    cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                    intVal = [cc_val[i][0] for i in range(0, val_count, 2)]
                    physVal = [cc_val[i][0] for i in range(1, val_count, 2)]
                    if numpy.all( numpy.diff( intVal ) > 0 ):
                        if conversionFormulaIdentifier == 4:
                            L[channelName] = numpy.interp( L[channelName], intVal, physVal ) # with interpolation
                        else:
                            try:
                                from scipy import interpolate
                            except:
                                raise('Please install scipy to convert channel'+channelName)
                            f = interpolate.interp1d( intVal, physVal , kind='nearest', bounds_error=False) # nearest
                            L[channelName] =  f(L[channelName]) # fill with Nan out of bounds while should be bounds
                    else:
                        print(( 'X values for interpolation of channel ' + channelName + ' are not increasing' ))
                elif conversionFormulaIdentifier == 6: # Value range to value table
                    val_count=int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /3)
                    cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                    key_min = [cc_val[i][0] for i in range(0, 3*val_count+1, 3)]
                    key_max = [cc_val[i][0] for i in range(1, 3*val_count+1, 3)]
                    value = [cc_val[i][0] for i in range(2, 3*val_count+1, 3)]
                    # look up in range keys
                    for Lindex in range(len(L[channelName] )):
                        key_index=0 # default index if not found
                        for i in range(val_count):
                            if  key_min[i] < L[channelName] [Lindex] < key_max[i]:
                                key_index=i
                                break
                        L[channelName][Lindex]=value[key_index]
                elif conversionFormulaIdentifier == 7: # Value to text / scale conversion table
                    val_count=int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] )
                    cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                    value=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
                    temp=[]
                    for Lindex in range(len(L[channelName] )):
                        key_index=val_count # default index if not found
                        for i in range(val_count):
                            if  L[channelName] [Lindex] == cc_val[i]:
                                key_index=i
                                break
                        temp.append(value[key_index])
                    L[channelName]=temp
                elif conversionFormulaIdentifier == 8: # Value range to Text / Scale conversion Table
                    val_count=int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] /2)
                    cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                    key_min = [cc_val[i][0] for i in range(0, 2*val_count, 2)]
                    key_max = [cc_val[i][0] for i in range(1, 2*val_count, 2)]
                    # look up in range keys
                    temp=[]
                    for Lindex in range(len(L[channelName] )):
                        key_index=val_count # default index if not found
                        for i in range(val_count):
                            if  key_min[i] < L[channelName] [Lindex] < key_max[i]:
                                key_index=i
                                break
                        temp.append(info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref'][key_index])
                    L[channelName]=temp
                elif conversionFormulaIdentifier == 9: # Text to value table
                    ref_count=int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref_count'] )
                    cc_val=info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                    cc_ref=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
                    temp=[]
                    for Lindex in range(len(L[channelName] )):
                        key_index=ref_count # default index if not found
                        for i in range(ref_count):
                            if  L[channelName] [Lindex] == cc_ref[i]:
                                key_index=i
                                break
                        temp.append(cc_val[key_index])
                    L[channelName]=temp
                elif conversionFormulaIdentifier == 10: # Text to text table
                    val_count=int(info['CCBlock'][dataGroup][channelGroup][channel]['cc_val_count'] )
                    cc_ref=info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
                    for Lindex in range(len(L[channelName] )):
                        key_index=val_count+1 # default index if not found
                        for i in range(0, val_count, 2):
                            if  L[channelName] [Lindex] == cc_ref[i]:
                                key_index=i
                                break
                        L[channelName] [Lindex]=cc_ref[key_index]
                else:
                    print('Unrecognized conversion type')
            elif chan.signalDataType==13:
                L['ms']=buf[recordID]['data']['data'].__getattribute__( 'ms_title') 
                L['min']=buf[recordID]['data']['data'].__getattribute__( 'min_title') 
                L['hour']=buf[recordID]['data']['data'].__getattribute__( 'hour_title') 
                L['day']=buf[recordID]['data']['data'].__getattribute__( 'day_title') 
                L['month']=buf[recordID]['data']['data'].__getattribute__( 'month_title') 
                L['year']=buf[recordID]['data']['data'].__getattribute__( 'year_title') 
            elif chan.signalDataType==14:
                L['ms']=buf[recordID]['data']['data'].__getattribute__( 'ms_title') 
                L['days']=buf[recordID]['data']['data'].__getattribute__( 'days_title') 
    if multiProc:
        Q.put(L)
    else:
        return L
        #print( 'Process ' + str( dataGroup ) + ' Finished ' + str( time.clock() - procTime ) )
def linearConv(vect, cc_val):
    P1 = cc_val[0]
    P2 = cc_val[1]
    return vect* P2 + P1
def rationalConv(vect,  cc_val):
    P1 = cc_val[0]
    P2 = cc_val[1]
    P3 = cc_val[2]
    P4 = cc_val[3]
    P5 = cc_val[4]
    P6 = cc_val[5]
    return (P1*vect *vect +P2*vect +P3)/(P4*vect *vect +P5*vect +P6)
def formulaConv(vect, formula):
    from sympy import lambdify, symbols
    X=symbols('X')
    expr=lambdify(X, formula, 'numpy')
    return expr(vect) 
