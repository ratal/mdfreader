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
        'master' : master channel of channel (time, crank angle, etc.)
        'description' : Description of channel

    Method plot('channelName') will plot the channel versus corresponding time
    Method resample(samplingTime) will resample all channels by interpolation and make only one master channel

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
from math import ceil
from mdf3reader import mdf3
from mdf4reader import mdf4
import numpy

from sys import version_info
from os.path import splitext
PythonVersion=version_info
PythonVersion=PythonVersion[0]

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
            nameList=channelNameList.listChannels3(self.fileName)
        else:
            from mdfinfo4 import info4
            channelNameList=info4()
            nameList=channelNameList.listChannels4(self.fileName)
        return nameList

class mdf( mdf3,  mdf4 ):
    """ mdf file class
    To use : yop= mdfreader.mdf('FileName.dat')
    Some additional useful methods:
    Resample : yop.resample(SamplingRate_in_secs)
    plot a specific channel : yop.plot('ChannelName')
    export to csv file : yop.exportCSV() , specific filename can be input
    export to netcdf : yop.exportNetCDF() """
    
    def __init__( self, fileName = None ):
        self.fileName = None
        self.VersionNumber=None
        self.masterChannelList = {}
        self.author=''
        self.organisation=''
        self.project=''
        self.subject=''
        self.comment=''
        self.time=''
        self.date=''
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
        
        print(self.fileName)
        
        # read file blocks
        info=mdfinfo(self.fileName)
        
        self.VersionNumber=info['IDBlock']['id_ver']
        if self.VersionNumber<400: # up to version 3.x not compatible with version 4.x
            self.read3(self.fileName, info, multiProc, channelList)
        else: #MDF version 4.x
            self.read4(self.fileName, info, multiProc, channelList)
    
    def write(self, fileName=None):
        #write mdf
        if self.fileName is not None:
            self.fileName = fileName
        else:
            print('Please provide filename')
        self.write3(fileName)
            
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
                    # plot using matplotlib the channel versus master channel
                    if len(list(self.masterChannelList.keys()))==1: # Resampled signals
                        masterName = self.masterChannelList.keys()
                        if masterName == '': # resampled channels, only one time channel most probably called 'master'
                            masterName ='master'
                        if masterName in list(self.keys()): # time channel properly defined
                            plt.plot( self[masterName]['data'], self[channelName]['data'] )
                            plt.xlabel( masterName + ' [' + self[masterName]['unit'] + ']' )
                        else: # no time channel found
                            plt.plot( self[channelName]['data'] )
                    else: # not resampled
                        if self[channelName]['master'] in list(self.keys()): # master channel is proper channel name
                            plt.plot( self[self[channelName]['master']]['data'], self[channelName]['data'] )
                            plt.xlabel( self[channelName]['master'] + ' [' + self[self[channelName]['master']]['unit'] + ']' )
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

    def resample( self, samplingTime = 0.1, masterChannelName=None ):
        """ Resamples mdf channels inside object for one single time"""
        # resample all channels to one sampling time vector
        if len(list(self.masterChannelList.keys()))>1: # Not yet resampled
            channelNames = list(self.keys())
            minTime = maxTime = []
            if masterChannelName is None:
                masterChannelName='master'
            self[masterChannelName] = {}
            unit = ''
            for master in list(self.masterChannelList.keys()):
                if master in self and len( self[master]['data'] ) > 5: # consider groups having minimum size 
                    minTime.append( self[master]['data'][0] )
                    maxTime.append( self[master]['data'][len( self[master]['data'] ) - 1] )
                    if self[master]['unit'] != '':
                        unit = self[master]['unit']
            self[masterChannelName]['data'] = numpy.arange( min( minTime ),max( maxTime ),samplingTime )
            self[masterChannelName]['unit'] = unit
            self[masterChannelName]['description'] = 'Unique master channel'

            # Interpolate channels
            timevect=[]
            for Name in channelNames:
                try:
                    if Name not in list(self.masterChannelList.keys()):
                        timevect = self[self[Name]['master']]['data']
                        if not self[Name]['data'].dtype.kind in ('S', 'U'): # if channel not array of string
                            self[Name]['data'] = numpy.interp( self[masterChannelName]['data'], timevect, self[Name]['data'] )
                            if masterChannelName in self[Name]:
                                del self[Name][masterChannelName]
                except:
                    if len( timevect ) != len( self[Name]['data'] ):
                        print(( Name + ' and time channel ' + self[Name][masterChannelName] + ' do not have same length' ))
                    elif not numpy.all( numpy.diff( timevect ) > 0 ):
                        print(( Name + ' has non regularly increasing time channel ' + self[Name][masterChannelName] ))
            # remove time channels in masterChannelList
            for ind in list(self.masterChannelList.keys()):
                del self[ind]
            self.masterChannelList = {} # empty dict
            self.masterChannelList[masterChannelName] = list(self.keys()) 
        else:
            pass

    def exportToCSV( self, filename = None, sampling = 0.1 ):
        # Export mdf file into csv
        # If no name defined, it will use original mdf name and path
        # By default, sampling is 0.1sec but can be changed
        import csv
        self.resample( sampling )
        if filename == None:
            filename = splitext(self.fileName)[0]
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
            filename = splitext(self.fileName)[0]
            filename = filename + '.nc'
        f = NetCDF.NetCDFFile( filename, 'w' )
        setattr( f, 'Author',(self.author))
        setattr( f, 'Date', (self.date))
        setattr( f, 'Time', (self.time))
        setattr( f, 'Organization', (self.organisation))
        setattr( f, 'ProjectName', (self.project))
        setattr( f, 'Subject', (self.subject))
        setattr( f, 'Comment', (self.comment))
        # Create dimensions having name of all time channels
        for time in list(self.masterChannelList.keys()):
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
            if len( list(self.masterChannelList.keys()) ) == 1: # mdf resampled
                var[name] = f.createVariable( CleanedName, type, ( list(self.masterChannelList.keys())[0], ) )
            else: # not resampled
                var[name] = f.createVariable( CleanedName, type, ( str(self[name]['master']), ) )
            # Create attributes
            setattr( var[name], 'title', CleanedName )
            setattr( var[name], 'units', self[name]['unit'])
            setattr( var[name], 'Description', self[name]['description'])
            if name in list(self.masterChannelList.keys()):
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
            filename = splitext(self.fileName)[0]
            filename = filename + '.hdf'
        f = h5py.File( filename, 'w' ) # create hdf5 file
        filegroup=f.create_group(os.path.basename(filename)) # create group in root associated to file
        attr = h5py.AttributeManager( filegroup ) # adds group/file metadata
        attr.create('Author',(self.author))
        attr.create('Date', (self.date))
        attr.create('Time', (self.time))
        attr.create('Organization', (self.organisation))
        attr.create('ProjectName', (self.project))
        attr.create('Subject', (self.subject))
        attr.create('Comment', (self.comment))
        if len( list(self.masterChannelList.keys()) ) > 1:
            # if several time groups of channels, not resampled
            groups = {}
            ngroups = 0
            grp = {}
            for channel in list(self.keys()):
                if self[channel]['master'] not in list(groups.keys()):
                    # create new time group
                    ngroups += 1
                    groups[self[channel]['master'] ] = ngroups
                    grp[ngroups] = filegroup.create_group( self[channel]['master'] )
                dset = grp[groups[self[channel]['master'] ]].create_dataset( channel, data = self[channel]['data'] )
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
            filename = splitext(self.fileName)[0]
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
            filename = filename = splitext(self.fileName)[0]
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
            filename = splitext(self.fileName)[0]
            filename = filename + '.xlsx'
        channels=list(self.keys())
        maxRows=max([len(self[channel]['data']) for channel in list(self.keys())]) # find max column length
        maxCols=len(list(self.keys())) # number of columns
        print('Creating Excel sheet')
        if len( list(self.masterChannelList.keys()) ) > 1: # not resampled data, can be long, writing cell by cell !
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
            if channel not in channelList and not 'master' in channel and channel not in list(self.masterChannelList.keys()) :
                # avoid to remove time channels otherwise problems with resample
                removeChannels.append(channel)
        if not len(removeChannels)==0:
            [self.masterChannelList[self[channel]['master']].remove(channel) for channel in removeChannels]
            [self.pop(channel) for channel in removeChannels]

    def copy(self):
        # copy a mdf class
        yop=mdf()
        yop.multiProc=self.multiProc
        yop.fileName=self.fileName
        yop.masterChannelList=self.masterChannelList
        for channel in list(self.keys()):
            yop[channel]=self[channel]
        return yop

    def mergeMdf(self, mdfClass):
        # merges data of 2 mdf classes
        # both must have been ressampled, otherwise, impossible to know time vector to match
        # create union of both lists
        unionedList=list(mdfClass.keys()) and list(self.keys())
        if not len(list(self.masterChannelList.keys()))==1:
            print( 'Data not resampled')
            raise
        initialTimeSize=len(self['master']['data'])
        for channel in unionedList:
            if channel in mdfClass and channel in self: # channel exists in both class
                if not channel=='master':
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
                refill=numpy.empty(len(mdfClass['master']['data']))
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
        datetimeInfo=numpy.datetime64(self.date.replace(':','-')+' '+self.time)
        originalKeys=self.keys()   
        # first find groups of channels sharing same time vector
        # initialise frame description
        #frame={}
        #for group in list(self.masterChannelList.keys()):
         #   frame[group]=[]
        # create list of channel for each time group
        #[frame[self[channel]['master']].append(channel) for channel in self.keys()]
        # create frames
        for group in list(self.masterChannelList.keys()):
            temp={}
            # convert time channel into timedelta
            time=datetimeInfo+numpy.array(self[group]['data']*10E6, dtype='timedelta64[us]')
            for channel in self.masterChannelList[group]:
                temp[channel]=pd.Series(self[channel]['data'], index=time)
            self[group+'_group']=pd.DataFrame(temp)
            self[group+'_group'].pop(group) # delete time channel, no need anymore
        # clean rest of self from data and time channel information
        [self[channel].pop('data') for channel in originalKeys]
        [self[channel].pop('master') for channel in originalKeys]
        self.masterGroups=[] # save time groups name in list
        [self.masterGroups.append(group+'_group') for group in list(self.masterChannelList.keys())]
        self.masterChannelList={}

if __name__ == "__main__":
    try:
        from multiprocessing import freeze_support
        freeze_support()
    except:
        None
    mdf()

