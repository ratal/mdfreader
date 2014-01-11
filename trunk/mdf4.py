# -*- coding: utf-8 -*-
""" Measured Data Format blocks paser for version 4.x
Created on Thu Dec 10 12:57:28 2014

:Author: `Aymeric Rateau <http://code.google.com/p/mdfreader/>`__

"""
import numpy
from math import log, exp
from sys import platform
from mdfinfo4 import info4, DATA
channelTime={'time','TIME'} # possible channel time writing, adjust if necessary

class mdf4(dict):
    """ mdf file class
    It imports mdf files version 4.1
    To use : yop= mdfreader.mdf('FileName.mf4')"""
    
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
            info = info4( self.fileName,  None )

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
        for dg in info['DGBlock'].keys():
            if not info['DGBlock']['dg_data']==0: # data exist
                #reads data block
                buf=DATA(fid, info['DGBlock']['dg_data'])
                
                numpyDataRecordFormat = []
                dataRecordName = []
                # defines record ID
                if info['DGBlock']['dg_data']==0: # no record ID
                    pass
                elif info['DGBlock']['dg_data']==1:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint8' ) )
                elif info['DGBlock']['dg_data']==2:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint16' ) )
                elif info['DGBlock']['dg_data']==3:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint32' ) )
                elif info['DGBlock']['dg_data']==4:
                    numpyDataRecordFormat.append( ( 'RecordID', 'uint64' ) )
                
                # defines data record for each channel group 
                if info['dg_rec_id_size']==0: # sorted data, no record ID, only one channel group
                    pass
                else: #unsorted data
                    for cg in info['CGBlock'][dg].keys():
                        recordId=info['CGBlock'][dg][cg]['cg_record_id']
                        for cn in info['CNBlock'][dg][cg][cn].keys():
            
                # reads all records according to their recordId
                recordSize=info['CGBlock'][dg][cg]['cg_dataBytes']
                
                # converts channel group records into channels
                buf = numpy.core.records.fromrecords( buf, dtype = numpyDataRecordFormat[recordId], shape = numberOfRecords[recordId] , names=dataRecordName[recordId])
