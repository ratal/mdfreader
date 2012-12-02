# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""
from PyQt4.QtGui import QMainWindow, QFileDialog, QAction
from PyQt4.QtCore import pyqtSignature, SIGNAL

from Ui_mdfreaderui import Ui_MainWindow
from io import open
from multiprocessing import Process,Pool,cpu_count
from mdfreader import mdfinfo,mdf

from sys import version_info
PythonVersion=version_info
PythonVersion=PythonVersion[0]
MultiProc=True # multiproc switch, for debug purpose

class MainWindow(QMainWindow, Ui_MainWindow, QFileDialog):
    """
    Class documentation goes here.
    """
    def __init__(self, parent = None):
        """
        Constructor
        """
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.fileNames=[] # files to convert
        self.mdfClass=mdf() # instance of mdf
        self.mdfinfoClass=mdfinfo() # instance of mdfinfo
        self.convertSelection='Matlab' # by default Matlab conversion is selected
        self.MergeFileBool=False # by default
        self.labFileName=[] # .lab file name
        self.actionPlotSelectedChannel = QAction("Plot", self.SelectedChannelList) # context menu to allow plot of channel
        self.SelectedChannelList.addAction(self.actionPlotSelectedChannel )
        self.connect(self.actionPlotSelectedChannel, SIGNAL("triggered()"), self.plotSelected)
        self.actionPlotChannel = QAction("Plot", self.channelList) # context menu to allow plot of channel
        self.channelList.addAction(self.actionPlotChannel)
        self.connect(self.actionPlotChannel, SIGNAL("triggered()"), self.plot)
        self.actionFileRemove= QAction("Delete", self.FileList) # context menu to remove selected file from list
        self.FileList.addAction(self.actionFileRemove)
        self.connect(self.actionFileRemove, SIGNAL("triggered()"), self.FileRemove)
    
    @pyqtSignature("")
    def on_browse_clicked(self):
        """
        Will open a dialog to browse for files
        """
        self.fileNames=QFileDialog.getOpenFileNames(self, "Select Measurement Files",filter=("MDF file (*.dat *.mdf)"))
        if not len(self.fileNames)==0:
            self.FileList.addItems(self.fileNames)
            self.mdfinfoClass.__init__()
            self.cleanChannelList()
            self.cleanSelectedChannelList()
            ChannelList=convertChannelList(self.mdfinfoClass.listChannels(str(self.fileNames[0])))
            self.SelectedChannelList.addItems(ChannelList)
    
    def cleanSelectedChannelList(self):
        # remove all items from list
        self.SelectedChannelList.clear()
        [self.SelectedChannelList.takeItem(0) for i in range(self.SelectedChannelList.count())]
        
    def cleanChannelList(self):
        # remove all items from list
        self.channelList.clear()
        [self.channelList.takeItem(0) for i in range(self.channelList.count())]

    @pyqtSignature("")
    def on_Convert_clicked(self):
        """
       Will convert mdf files into selected format
        """
        # create list of channels to be converted for all files
        channelList=[]
        [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
        channelList=list(set(channelList)) # remove duplicates
        # Process all mdf files recursively
        if self.FileList.count()>0: # not empty list
            ncpu=cpu_count() # to still have response from PC
            if ncpu<1:
                ncpu=1
            pool = Pool(processes=ncpu)
            if not self.MergeFileBool: # export all files separatly
                convertFlag=True
                convertSelection=self.convertSelection
                resampleValue=float(self.resampleValue.text())
                #resample if requested
                if self.resample.checkState():
                    if not len(self.resampleValue.text())==0:
                        resampleFlag=True
                    else:
                        print('Empty field for resampling')
                        raise 
                else:
                    resampleFlag=False
                args=[(str(self.FileList.takeItem(0).text()),channelList,resampleFlag,resampleValue,convertFlag,convertSelection) for i in range(self.FileList.count())]
                if MultiProc:
                    result=pool.map_async(processMDFstar,args)
                    result.get() # waits until finished
                else:
                    result=list(map(processMDFstar,args)) 
                self.cleanChannelList()
            elif self.FileList.count(): # Stack files data if min 2 files in list
                # import first file
                if len(self.resampleValue.text())==0:
                    print('Wrong value for resampling')
                    raise 
                convertFlag=False
                convertSelection=self.convertSelection
                resampleValue=float(self.resampleValue.text())
                resampleFlag=True # always resample when merging
                fileName=str(self.FileList.item(0).text()) # Uses first file name for the converted file
                # list filenames
                args=[(str(self.FileList.takeItem(0).text()),channelList,resampleFlag,resampleValue,convertFlag,convertSelection) for i in range(self.FileList.count())]
                if MultiProc:
                    res=pool.map_async(processMDFstar,args)
                    result=res.get()
                else:
                    result=map(processMDFstar,args) # no multiproc, for debug
                # Merge results
                self.mdfClass.__init__() # clear memory
                self.mdfClass.fileName=fileName
                self.mdfClass.multiProc=False
                buffer=self.mdfClass.copy()
                res=result.pop(0)
                self.mdfClass.update(res[0])
                self.mdfClass.timeChannelList=res[1]
                for res in result: # Merge
                    buffer.__init__()
                    buffer.update(res[0])
                    buffer.timeChannelList=res[1]
                    self.mdfClass.mergeMdf(buffer)
                # Export
                if self.convertSelection=='Matlab':
                    self.mdfClass.exportToMatlab()
                elif self.convertSelection=='csv':
                    self.mdfClass.exportToCSV()
                elif self.convertSelection=='netcdf':
                    self.mdfClass.exportToNetCDF()
                elif self.convertSelection=='hdf5':
                    self.mdfClass.exportToHDF5()
                elif self.convertSelection=='excel':
                    self.mdfClass.exportToExcel()
                elif self.convertSelection=='excel2010':
                    self.mdfClass.exportToXlsx()
                self.cleanChannelList()
                #self.cleanSelectedChannelList()
                self.mdfClass.__init__() # clear memory
    
    @pyqtSignature("QListWidgetItem*")
    def on_FileList_itemClicked(self, item):
        """
        If user click on file list
        """
        # Refresh list of channels from selected file
        self.mdfinfoClass.__init__()
        #self.mdfinfoClass.readinfo(item)
        self.cleanChannelList()
        ChannelList=convertChannelList(self.mdfinfoClass.listChannels(str(item.text())))
        self.channelList.addItems(ChannelList)
        self.mdfinfoClass.__init__() # clean object to free memory
    
    @pyqtSignature("bool")
    def on_matlab_clicked(self, checked):
        """
        Selects Matlab conversion
        """
        self.convertSelection='Matlab'
    
    @pyqtSignature("bool")
    def on_netcdf_clicked(self, checked):
        """
        Selects netcdf conversion.
        """
        self.convertSelection='netcdf'
    
    @pyqtSignature("bool")
    def on_hdf5_clicked(self, checked):
        """
        Selects hdf5 conversion.
        """
        self.convertSelection='hdf5'
    
    @pyqtSignature("bool")
    def on_csv_clicked(self, checked):
        """
        Selects csv conversion.
        """
        self.convertSelection='csv'
    
    @pyqtSignature("bool")
    def on_excel_clicked(self, checked):
        """
        Selects excel conversion.
        """
        self.convertSelection='excel'
 
    @pyqtSignature("bool")
    def on_excel2010_clicked(self, checked):
        """
        Selects excel conversion.
        """
        self.convertSelection='excel2010'

    @pyqtSignature("")
    def on_LabFileBrowse_clicked(self):
        """
        selects lab file from browser.
        """
        self.labFileName=QFileDialog.getOpenFileName(self, "Select Lab Files", filter=("Lab file (*.lab)"))
        if not len(self.labFileName)==0:
            self.LabFile.del_() # clear linedit
            self.LabFile.insert(str(self.labFileName)) # replace linedit field by browsed file name
            # read lab file
            labfile=open(str(self.labFileName), 'r')
            self.labChannelList=[]
            ine = labfile.readline() # read first line [lab]
            while 1:
                line = labfile.readline()
                if not line:
                    break
                self.labChannelList.append(line.replace('\n',''))
            self.cleanSelectedChannelList() # Clear Selected file list
            self.SelectedChannelList.addItems(self.labChannelList)
    
    def plot(self):
        #Finds selected file and read it
        selectedFile=self.FileList.selectedItems()
        self.mdfClass.__init__(str(selectedFile[0].text())) # read file
        # list items selected in listWidget
        Channels=self.channelList.selectedItems()
        selectedChannels=[]
        [selectedChannels.append(str(Channels[i].text())) for i in range(len(Channels))]
        # plot channels
        self.mdfClass.plot(selectedChannels)
        
    def plotSelected(self):
        # plots channels from selected list
        selectedFile=self.FileList.selectedItems()
        if not len(selectedFile)==0:
            self.mdfClass.__init__(str(selectedFile[0].text())) # read file
        else:
            self.mdfClass.__init__(str(self.FileList[0].text())) # read file
        # list items selected in listWidget
        Channels=self.SelectedChannelList.selectedItems()
        selectedChannels=[]
        [selectedChannels.append(str(Channels[i].text())) for i in range(len(Channels))]
        # plot channels
        self.mdfClass.plot(selectedChannels)
    def FileRemove(self):
        # removes selected file
        selectionList=self.FileList.selectedItems()
        [self.FileList.takeItem(self.FileList.row(selectionList[i])) for i in range(len(selectionList))]
        
    def on_SelectedChannelList_dropEvent(self):
        # avoids to have duplicates in list when channel is dropped
        channelList=[] 
        [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
        channelList=list(set( channelList)) #  removeDuplicates
        self.SelectedChannelList.clear()
        self.SelectedChannelList.addItems(channelList)
    
    @pyqtSignature("bool")
    def on_MergeFiles_toggled(self, checked):
        """
        Slot documentation goes here.
        """
        # toggle flag to merge files
        self.MergeFileBool= not self.MergeFileBool
        if self.MergeFileBool:
            self.resample.setCheckState(2)

def processMDF(fileName,channelist,resampleFlag,resampleValue,convertFlag,convertSelection):
    # Will process file according to defined options
    yop=mdf()
    yop.multiProc=False # already multiprocessed
    yop.read(fileName,channelList=channelist)
    if resampleFlag:
        yop.resample(resampleValue)
    if convertFlag:
        if convertSelection=='Matlab':
            yop.exportToMatlab()
        elif convertSelection=='csv':
            yop.exportToCSV()
        elif convertSelection=='netcdf':
            yop.exportToNetCDF()
        elif convertSelection=='hdf5':
            yop.exportToHDF5()
        elif convertSelection=='excel':
            yop.exportToExcel()
        elif convertSelection=='excel2010':
            yop.exportToXlsx()
    yopPicklable={} # picklable dict and not object
    for channel in list(yop.keys()):
        yopPicklable[channel]=yop[channel]
    return [yopPicklable,yop.timeChannelList]

def processMDFstar(args):
    try:
        return processMDF(*args)
    except Exception, e:
        print('Error, following file might be corrupted : '+args[0]) # Shows fileName and parameters to help finding corrupted files
        print('Please re-try by removing this file from the list and restart mdfconverter to kill processes and clean memory')
        raise(e) # produce error and stops all processes

def convertChannelList(channelList):
    if PythonVersion<3:
        return [str(name) for name in channelList]
    else:
        return [(name).decode('ascii', 'replace') for name in channelList]
