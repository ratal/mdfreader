# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""

from PyQt4.QtGui import QMainWindow, QFileDialog, QAction
from PyQt4.QtCore import pyqtSignature, SIGNAL, QStringList

from Ui_mdfreaderui import Ui_MainWindow
import io
#from mdfreader import mdf, mdfinfo

class MainWindow(QMainWindow, Ui_MainWindow, QFileDialog):
    """
    Class documentation goes here.
    """
    def __init__(self, parent = None, mdfClass=None, mdfinfoClass=None):
        """
        Constructor
        """
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.fileNames=[] # files to convert
        self.mdfClass=mdfClass # instance of mdf
        self.mdfinfoClass=mdfinfoClass # instance of mdfinfo
        self.convertSelection='Matlab' # by default Matlab conversion is selected
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
        self.fileNames=QFileDialog.getOpenFileNames(self, "Select Measurement Files", filter=("MDF file (*.dat)"))
        if not self.fileNames.isEmpty():
            self.FileList.addItems(self.fileNames)
            self.mdfinfoClass.__init__()
            self.mdfinfoClass.listChannels(str(self.fileNames[0]))
            self.cleanChannelList()
            self.cleanSelectedChannelList()
            self.SelectedChannelList.addItems(self.mdfinfoClass.channelNameList)
    
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
        # Process all mdf files recursively
        if self.FileList.count()>0: # not empty list
            # create list of channels to be converted for all files
            channelList=QStringList([]) # pass by QStringList to use removeDuplicates
            [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
            channelList.removeDuplicates()
            for i in range(self.FileList.count()):
                # read first file of list and removes it from list
                self.mdfClass.__init__()
                self.mdfClass.read(str(self.FileList.takeItem(0).text()), multiProc = True, channelList=channelList)
                self.show()
                #resample if requested
                if self.resample.checkState():
                    if not self.resampleValue.text().isEmpty():
                        self.mdfClass.resample(float(self.resampleValue.text()))
                        print('Resampled successfully'+str(self.fileNames[i]))
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
            self.cleanChannelList()
            #self.cleanSelectedChannelList()
    
    @pyqtSignature("QListWidgetItem*")
    def on_FileList_itemClicked(self, item):
        """
        If user click on file list
        """
        # Refresh list of channels from selected file
        self.mdfinfoClass.__init__()
        self.mdfinfoClass.listChannels(str(item.text())) # read file
        #self.mdfinfoClass.readinfo(item)
        self.cleanChannelList()
        self.channelList.addItems(self.mdfinfoClass.channelNameList)
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
        self.convertSlection='netcdf'
    
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
    
    @pyqtSignature("")
    def on_LabFileBrowse_clicked(self):
        """
        selects lab file from browser.
        """
        self.labFileName=QFileDialog.getOpenFileName(self, "Select Lab Files", filter=("Lab file (*.lab)"))
        if not self.labFileName.isEmpty():
            self.LabFile.del_() # clear linedit
            self.LabFile.insert(str(self.labFileName)) # replace linedit field by browsed file name
            # read lab file
            labfile=io.open(str(self.labFileName), 'r')
            self.labChannelList=[]
            ine = labfile.readline() # read first line [lab]
            while 1:
                line = labfile.readline()
                if not line:
                    break
                self.labChannelList.append(line)
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
        channelList=QStringList([]) # pass by QStringList to use removeDuplicates
        [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
        channelList.removeDuplicates()
        self.SelectedChannelList.clear()
        self.SelectedChannelList.addItems(channelList)
        
