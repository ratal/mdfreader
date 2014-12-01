# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Thu Dec 10 12:57:28 2013

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
from numpy.core.records import fromrecords, fromstring, fromfile
from numpy import array, recarray, append, asarray, empty, zeros, dtype
from numpy import hstack, arange, right_shift, bitwise_and, all, diff, interp
from struct import unpack, Struct
from math import pow
from sys import platform, version_info
from io import open # for python 3 and 2 consistency
from mdfinfo4 import info4, MDFBlock,  ATBlock#, CNBlock
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
    def __init__(self, fid,  pointer, record, zip_type=None, channelList=None, sortedFlag=True):
        # block header
        self.loadHeader(fid, pointer)
        if self['id'] in ('##DT', '##RD', b'##DT', b'##RD'): # normal data block
            if sortedFlag:
                record.numberOfRecords = (self['length']-24) // record.recordLength
                self['data']=record.readSortedRecord(fid, pointer, channelList)
            else: # unsorted reading
                print('not implemented yet unsorted data block reading') # to be implemented
        
        elif self['id'] in ('##SD', b'##SD'): # 
            if record.signalDataType==6:
                format='ISO8859'
            elif record.signalDataType==7:
                format='utf-8'
            elif record.signalDataType==8:
                format='<utf-16'
            elif record.signalDataType==9:
                format='>utf-16'
            pointer=0
            buf=[]
            nElement=0
            while pointer<self['length']-24:
                VLSDLen=unpack('I', fid.read(4))[0] # length of data
                buf.append(fid.read(VLSDLen).decode(format).rstrip('\x00'))
                pointer+=VLSDLen+4
                nElement+=1
            if nElement>1:
                buf=equalizeStringLength(buf)
                self['data']=fromrecords(buf, names=str(record.name))
            else: # only one data
                self['data']=buf
        
        elif self['id'] in ('##DZ', b'##DZ'): # zipped data block
            self['dz_org_block_type']=self.mdfblockreadCHAR(fid, 2)
            self['dz_zip_type']=self.mdfblockread(fid, UINT8, 1)
            if zip_type is not None: # HLBlock used
                self['dz_zip_type']=zip_type
            self['dz_reserved']=self.mdfblockreadBYTE(fid, 1)
            self['dz_zip_parameter']=self.mdfblockread(fid, UINT32, 1)
            self['dz_org_data_length']=self.mdfblockread(fid, UINT64, 1)
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
                M = self['dz_org_data_length']//N
                temp=array(memoryview(self['data'][:M*N]))
                tail=array(memoryview(self['data'][M*N:]))
                temp=temp.reshape(N, M).T.ravel()
                if len(tail)>0:
                    temp=append(temp, tail)
                self['data']=temp.tostring()
            if channelList is None and sortedFlag: # reads all blocks if sorted block and no channelList defined
                record.numberOfRecords = self['dz_org_data_length'] // record.recordLength
                self['data']=fromstring(self['data'] , dtype = record.numpyDataRecordFormat, shape = record.numberOfRecords , names=record.dataRecordName)
            elif channelList is not None and sortedFlag: # sorted data but channel list requested
                print('not implemented yet sorted compressed data block reading with channelList') # to be implemented
                
            else: # unsorted reading
                # reads only the channels using offset functions, channel by channel. Not yet ready
                buf={}
                position=0
                recordIdCFormat=record[list(record.keys())[0]]['record'].recordIDCFormat
                recordIDsize=record[list(record.keys())[0]]['record'].recordIDsize
                VLSDStruct=Struct('I')
                # initialise data structure
                for recordID in record:
                    for channelName in record[recordID]['record'].channelNames:
                        buf[channelName]=[]
                # read data
                while position<len(self['data']):
                    recordID=recordIdCFormat.unpack(self['data'][position:position+recordIDsize])[0]
                    if not record[recordID]['record'].Flags & 0b1: # not VLSD CG)
                        temp=record.readRecord(recordID, self['data'][position:position+record[recordID]['record'].recordLength+1], channelList)
                        position += record[recordID]['record'].recordLength+1
                        for channelName in temp:
                            buf[channelName].append(temp[channelName])
                    else: # VLSD CG
                        position+=recordIDsize
                        VLSDLen=VLSDStruct.unpack(self['data'][position:position+4])[0] # VLSD length
                        position+=4
                        temp=self['data'][position:position+VLSDLen-1]
                        if record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType==6:
                            temp=temp.decode('ISO8859')
                        elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType==7:
                            temp=temp.decode('utf-8')
                        elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType==8:
                            temp=temp.decode('<utf-16')
                        elif record[recordID]['record'].VLSD_CG[recordID]['channel'].signalDataType==9:
                            temp=temp.decode('>utf-16')
                        buf[record[recordID]['record'].VLSD_CG[recordID]['channelName']].append(temp)
                        position += VLSDLen
                # convert list to array
                for chan in buf:
                    buf[chan]=array(buf[chan])
                self['data']=buf

def equalizeStringLength(buf):
    maxlen=len(max(buf, key=len))
    for i in range(len(buf)): # resize string to same length, numpy constrain
        buf[i]+=' '*(maxlen-len(buf[i]))
    return buf

def append_field(rec, name, arr, numpy_dtype=None):
    arr = asarray(arr)
    if numpy_dtype is None:
        numpy_dtype = arr.dtype
    newdtype = dtype(rec.dtype.descr + [(name, numpy_dtype)])
    newrec = empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec
    
def change_field_name(arr, old_name, new_name):
    names=list(arr.dtype.names)
    for n in range(len(names)):
        if names[n]==old_name:
            break
    names[n]=new_name
    names=tuple(names)
    arr.dtype.names=names
    return arr

class DATA(dict):
    def __init__(self, fid, pointer):
        self.fid=fid
        self.pointerTodata=pointer
        self.type='sorted'
    def addRecord(self, record):
        self[record.recordID]={}
        self[record.recordID]['record']=record
        # detects VLSD CG
        for recordID in self[record.recordID]['record'].VLSD_CG:
            self[recordID]['record'].VLSD_CG=self[record.recordID]['record'].VLSD_CG
    def read(self, channelList, zip=None):
        if len(self)==1: #sorted dataGroup
            recordID=list(self.keys())[0]
            record=self[recordID]['record']
            self[recordID]['data']=self.loadSorted( record, zip=None, nameList=channelList)
            for cn in record.VLSD: # VLSD channels
                if channelList is None or record[cn].name in channelList:
                    temp=DATA(self.fid,  record[cn].data)
                    temp=temp.loadSorted(record[cn], zip=None, nameList=channelList)
                    rec=self[recordID]['data']['data']
                    rec=change_field_name(rec, convertName(record[cn].name), convertName(record[cn].name)+'_offset')
                    rec=append_field(rec, convertName(record[cn].name), temp['data'])
                    self[recordID]['data']['data']=rec.view(recarray)
        else: # unsorted DataGroup
            self.type='unsorted'
            self['data']=self.loadUnsorted( zip=None, nameList=channelList)
            
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
                    temps['data']=DATABlock(self.fid, temps['dl_data'][0], record, zip, channelList=nameList)
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
                    temps['data']['data']=hstack(temps['data']['data']) # concatenate data list
                    temps['data']['data']=temps['data']['data'].view(recarray) # vstack output ndarray instead of recarray
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
    
    def loadUnsorted(self, zip=None, nameList=None):
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
                    temps['data']=DATABlock(self.fid, temps['dl_data'][0], self, nameList, sortedFlag=False)
                elif temps['dl_count']==0:
                    raise('empty datalist')
                else:
                    # define size of each block to be read
                    temps['data']={}
                    temps['data']['data']=[]
                    temps['data']['data']=[DATABlock(self.fid, temps['dl_data'][dl], self, zip, nameList, sortedFlag=False)['data'] for dl in range(temps['dl_count'])]
                    if temps['dl_dl_next']:
                         self.pointerTodata=temps['dl_dl_next']
                         temps['data']['data'].append(self.loadUnsorted(record, zip, nameList)['data'])
                    #Concatenate all data
                    temps['data']['data']=hstack(temps['data']['data']) # concatenate data list
                    temps['data']['data']=temps['data']['data'].view(recarray) # vstack output ndarray instead of recarray
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
            temps['data']=self.loadUnsorted(record, zip=temps['hl_zip_type'], nameList=nameList)#, record, zip=temps['hl_zip_type'], nameList=nameList)
        else:
            temps['data']=DATABlock(self.fid, self.pointerTodata, self, zip_type=zip, channelList=nameList, sortedFlag=False)
        return temps['data']
    def readRecord(self, recordID, buf, channelList=None):
        return self[recordID]['record'].readRecordBuf(buf, channelList)
        
class recordChannel():
    def __init__(self, info, dataGroup, channelGroup, channelNumber, recordIDsize):
        self.name=info['CNBlock'][dataGroup][channelGroup][channelNumber]['name']
        self.channelNumber=channelNumber
        self.signalDataType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data_type']
        self.bitCount = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_count']
        self.channelType = info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_type']
        self.dataFormat=arrayformat4( self.signalDataType, self.bitCount )
        if self.channelType in (1, 4, 5): # if VSLD or sync or max length channel
            self.data=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_data']
            if self.channelType==1: # if VSLD
                self.dataFormat=arrayformat4(0, self.bitCount )
        self.RecordFormat=((self.name, convertName(self.name)),  self.dataFormat)
        if not self.signalDataType in (13, 14): # processed by several channels
            if not self.channelType==1: # if VSLD
                self.Format=datatypeformat4( self.signalDataType, self.bitCount )
            else:
                self.Format=datatypeformat4( 0, self.bitCount )
            self.CFormat=Struct(self.Format)
        if self.bitCount%8==0:
            self.nBytes=self.bitCount // 8
        else:
            self.nBytes=self.bitCount // 8 + 1
        self.bitOffset=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_bit_offset']
        self.byteOffset=info['CNBlock'][dataGroup][channelGroup][channelNumber]['cn_byte_offset']
        self.posBeg=recordIDsize+self.byteOffset
        self.posEnd=recordIDsize+self.byteOffset+self.nBytes
        self.VLSD_CG_Flag=False
        
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
        self.recordIDCFormat=''
        self.dataGroup=dataGroup
        self.channelGroup=channelGroup
        self.numpyDataRecordFormat=[]
        self.dataRecordName=[]
        self.master={}
        self.Flags=0
        self.VLSD_CG={}
        self.VLSD=[]
        self.MLSD={}
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
                self.recordIDCFormat=Struct('B')
            elif self.recordIDsize==2:
                self.numpyDataRecordFormat.append( ( format, 'uint16' ) )
                self.recordIDCFormat='H'
            elif self.recordIDsize==3:
                self.numpyDataRecordFormat.append( ( format, 'uint32' ) )
                self.recordIDCFormat='I'
            elif self.recordIDsize==4:
                self.numpyDataRecordFormat.append( ( format, 'uint64' ) )
                self.recordIDCFormat='L'
        self.recordID=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_record_id']
        self.recordLength=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_data_bytes']
        self.recordLength+= info['CGBlock'][self.dataGroup][self.channelGroup]['cg_invalid_bytes']
        self.numberOfRecords=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_cycle_count']
        self.Flags=info['CGBlock'][self.dataGroup][self.channelGroup]['cg_flags']
        if 'MLSD' in info:
            self.MLSD=info['MLSD']

        for channelNumber in list(info['CNBlock'][self.dataGroup][self.channelGroup].keys()):
            channel=recordChannel(info, self.dataGroup, self.channelGroup, channelNumber, self.recordIDsize)
            if channel.channelType in (2, 3): # master channel found
                self.master['name']=channel.name
                self.master['number']=channelNumber
            if channel.channelType in (0, 1, 2, 4, 5): # not virtual channel
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
                if 'VLSD_CG' in info: # is there VLSD CG
                    for recordID in info['VLSD_CG']:# look for VLSD CG Channel
                        if info['VLSD_CG'][recordID]['cg_cn']==(self.channelGroup, channelNumber):
                            self.VLSD_CG[recordID]=info['VLSD_CG'][recordID]
                            self.VLSD_CG[recordID]['channel']=channel
                            self.VLSD_CG[recordID]['channelName']=channel.name
                            self[-1].VLSD_CG_Flag=True
                            break
                if channel.channelType==1:
                    self.VLSD.append(channel.channelNumber)
            elif channel.channelType in (3, 6): # virtual channel
                pass # channel calculated based on record index later in conversion function

    def readSortedRecord(self, fid, pointer, channelList=None):
        # reads record, only one channel group per datagroup
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
            
    def readRecordBuf(self, buf, channelList=None):
        # read stream of record bytes
        temp={}
        if channelList is None:
            channelList = self.channelNames
        for channel in self: # list of recordChannels from channelList
            if channel.name in channelList and not channel.VLSD_CG_Flag:
                temp[channel.name] = channel.CFormat.unpack(buf[channel.posBeg:channel.posEnd])[0]
        return temp # returns dictionary of channel with its corresponding values
            
class mdf4(dict):
    """ mdf file class
    It imports mdf files version 4.1
    To use : yop= mdfreader.mdf('FileName.mf4')"""
    
    def __init__( self, fileName = None, info=None,multiProc = False,  channelList=None, convertAfterRead=True):
        self.masterChannelList = {}
        self.multiProc = False # flag to control multiprocessing, default deactivate, giving priority to mdfconverter
        self.author=''
        self.organisation=''
        self.project=''
        self.subject=''
        self.comment=''
        self.time=''
        self.date=''
        self.convertAfterRead=True
        # clears class from previous reading and avoid to mess up
        self.clear()
        if fileName == None and info!=None:
            self.fileName=info.fileName
        else:
            self.fileName=fileName
        self.read4(self.fileName, info, multiProc, channelList)

    ## reads mdf file
    def read4( self, fileName=None, info = None, multiProc = False, channelList=None, convertAfterRead=True):
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
                        if info['CNBlock'][dataGroup][channelGroup][channel]['cn_type'] in (2, 3):
                            masterDataGroup[dataGroup]=info['CNBlock'][dataGroup][channelGroup][channel]['name']

                buf.read(channelList)

                # Convert channels to physical values
                OkBuf=len(buf)>0 and 'data' in buf[list(buf.keys())[0]] and buf[list(buf.keys())[0]]['data'] is not None
                if self.multiProc and OkBuf: 
                    proc.append( Process( target = processDataBlocks4, args = ( Q, buf, info, dataGroup, channelList, self.multiProc) ) )
                    proc[-1].start()
                elif OkBuf: # for debugging purpose, can switch off multiprocessing
                    L.update(processDataBlocks4( None, buf, info, dataGroup, channelList, self.multiProc))
                elif buf.type=='unsorted':
                    L=buf['data']['data']
                else:
                    print('No data in dataGroup '+ str(dataGroup))

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
                        L.pop(channelName, None) # free memory
                        if info['CNBlock'][dataGroup][channelGroup][channel]['cn_type']==4: # sync channel
                            # attach stream to be synchronised
                            self[channelName]['attachment']=ATBlock(fid, info['CNBlock'][dataGroup][channelGroup][channel]['cn_data'])
                        convType = info['CCBlock'][dataGroup][channelGroup][channel]['cc_type']
                        if convType in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10): # needs conversion
                            self[channelName]['conversion'] = {}
                            self[channelName]['conversion']['parameters'] = {}
                            self[channelName]['conversion']['type'] = convType
                            if 'cc_val' in info['CCBlock'][dataGroup][channelGroup][channel]:
                                self[channelName]['conversion']['parameters']['cc_val'] = info['CCBlock'][dataGroup][channelGroup][channel]['cc_val']
                            if 'cc_ref' in info['CCBlock'][dataGroup][channelGroup][channel]:
                                self[channelName]['conversion']['parameters']['cc_ref'] = info['CCBlock'][dataGroup][channelGroup][channel]['cc_ref']
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
        
        fid.close() # close file
        
        if self.convertAfterRead: 
            self.convertAllChannel4()
        #print( 'Finished in ' + str( time.clock() - inttime ) )

    def getChannelData4(self, channelName):
        # returns data of channel
        if channelName in self:
            return self.convert4(channelName)
        else:
            raise('Channel not in dictionary')
    
    def convert4(self, channelName):
        # returns data converted if necessary
        if 'conversion' in self[channelName]: # there is conversion property
            if self[channelName]['conversion']['type'] == 1:
                return linearConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'])
            elif self[channelName]['conversion']['type'] == 2:
                return rationalConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'])
            elif self[channelName]['conversion']['type'] == 3:
                return formulaConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_ref']['Comment'])
            elif self[channelName]['conversion']['type'] == 4:
                return valueToValueTableWInterpConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'])
            elif self[channelName]['conversion']['type'] == 5:
                return valueToValueTableWOInterpConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'])
            elif self[channelName]['conversion']['type'] == 6:
                return valueRangeToValueTableConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'])
            elif self[channelName]['conversion']['type'] == 7:
                return valueToTextConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'], self[channelName]['conversion']['parameters']['cc_ref'])
            elif self[channelName]['conversion']['type'] == 8:
                return valueRangeToTextConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'], self[channelName]['conversion']['parameters']['cc_ref'])
            elif self[channelName]['conversion']['type'] == 9:
                return textToValueConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_val'], self[channelName]['conversion']['parameters']['cc_ref'])
            elif self[channelName]['conversion']['type'] == 10:
                return textToTextConv(self[channelName]['data'], self[channelName]['conversion']['parameters']['cc_ref'])
            else:
                return self[channelName]['data']
        else:
            return self[channelName]['data']
    
    def convertChannel4(self, channelName):
        self[channelName]['data'] = self.convert4(channelName)
        if 'conversion' in self[channelName]:
            self[channelName].pop('conversion')
    
    def convertAllChannel4(self):
        for channel in self:
            self.convertChannel4(channel)

def convertName(channelName):
    if PythonVersion<3:
        channelIdentifier=channelName.encode('utf-8')+'_title'
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
        dataType = 'S'+str(numberOfBits//8) # not directly processed
    elif signalDataType == 7: # UTF-8
        dataType = 'U' # not directly processed
    elif signalDataType in (8, 9): # UTF-16
        dataType = 'U' # not directly processed
    elif signalDataType == 10: # bytes array
        dataType = 'V'+str(numberOfBits//8) # not directly processed
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
        dataType = str(numberOfBits//8)+'s'
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
                    temp = buf[recordID]['data']['data'].__getattribute__( convertName(recordName) ) # extract channel vector
                else :# virtual channel 
                    temp=arange(buf[recordID]['record'].numberOfRecords)
                    
                if info['CNBlock'][dataGroup][channelGroup][channel]['cn_sync_type'] in (2, 3, 4): # master channel 
                    channelName = 'master' + str( dataGroup )

                # Process concatenated bits inside uint8
                if not chan.bitCount//8.0==chan.bitCount/8.0: # if channel data do not use complete bytes
                    mask = int(pow(2, chan.bitCount+1)-1) # masks isBitUnit8
                    if chan.signalDataType in (0,1, 2, 3): # integers
                        temp =  right_shift(temp,  chan.bitOffset)
                        temp =  bitwise_and(temp,  mask )
                        L[channelName] = temp
                    else: # should not happen
                        print('bit count and offset not applied to correct data type')
                        L[channelName] = temp
                else: #data using bull bytes
                    L[channelName] = temp
                
                # MLSD
                if chan.channelType==5:
                    pass #print('MLSD masking needed')
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
    try:
        from sympy import lambdify, symbols
    except:
        print('Please install sympy to convert channel ')
    X=symbols('X')
    expr=lambdify(X, formula, 'numpy')
    return expr(vect) 

def valueToValueTableWOInterpConv(vect, cc_val):
    val_count=2*int(len(cc_val) /2)
    intVal = [cc_val[i][0] for i in range(0, val_count, 2)]
    physVal = [cc_val[i][0] for i in range(1, val_count, 2)]
    if all( diff( intVal ) > 0 ):
        try:
            from scipy import interpolate
        except:
            raise('Please install scipy to convert channel')
        f = interpolate.interp1d( intVal, physVal , kind='nearest', bounds_error=False) # nearest
        return f(vect) # fill with Nan out of bounds while should be bounds
    else:
        print(( 'X values for interpolation of channel are not increasing' ))

def valueToValueTableWInterpConv(vect, cc_val):
    val_count=2*int(len(cc_val) /2)
    intVal = [cc_val[i][0] for i in range(0, val_count, 2)]
    physVal = [cc_val[i][0] for i in range(1, val_count, 2)]
    if all( diff( intVal ) > 0 ):
        return interp( vect, intVal, physVal ) # with interpolation
    else:
        print(( 'X values for interpolation of channel are not increasing' ))

def valueRangeToValueTableConv(vect, cc_val):
    val_count=int(len(cc_val) /3)
    key_min = [cc_val[i][0] for i in range(0, 3*val_count+1, 3)]
    key_max = [cc_val[i][0] for i in range(1, 3*val_count+1, 3)]
    value = [cc_val[i][0] for i in range(2, 3*val_count+1, 3)]
    # look up in range keys
    for Lindex in range(len(vect )):
        key_index=0 # default index if not found
        for i in range(val_count):
            if  key_min[i] < vect[Lindex] < key_max[i]:
                key_index=i
                break
        vect[Lindex]=value[key_index]
    return vect
def valueToTextConv(vect, cc_val, cc_ref):
    val_count=int(len(cc_val))
    temp=[]
    for Lindex in range(len(vect )):
        key_index=val_count # default index if not found
        for i in range(val_count):
            if  vect[Lindex] == cc_val[i]:
                key_index=i
                break
        temp.append(cc_ref[key_index])
    return temp
    
def valueRangeToTextConv(vect, cc_val, cc_ref):
    val_count=int(len(cc_val) /2)
    key_min = [cc_val[i][0] for i in range(0, 2*val_count, 2)]
    key_max = [cc_val[i][0] for i in range(1, 2*val_count, 2)]
    # look up in range keys
    temp=[]
    for Lindex in range(len(vect)):
        key_index=val_count # default index if not found
        for i in range(val_count):
            if  key_min[i] < vect[Lindex] < key_max[i]:
                key_index=i
                break
        temp.append(cc_ref[key_index])
    return temp
    
def textToValueConv(vect, cc_val, cc_ref):
    ref_count=len(cc_ref)
    temp=[]
    for Lindex in range(len(vect)):
        key_index=ref_count # default index if not found
        for i in range(ref_count):
            if  vect[Lindex] == cc_ref[i]:
                key_index=i
                break
        temp.append(cc_val[key_index])
    return temp

def textToTextConv(vect, cc_ref):
    ref_count=len(cc_ref)-2
    for Lindex in range(len(vect)):
        key_index=ref_count+1 # default index if not found
        for i in range(0, ref_count, 2):
            if  vect[Lindex] == cc_ref[i]:
                key_index=i
                break
        vect[Lindex]=cc_ref[key_index]
    return vect
