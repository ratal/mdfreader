# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""

from PyQt4.QtGui import QMainWindow, QFileDialog
from PyQt4.QtCore import pyqtSignature

from Ui_mdfreaderui import Ui_MainWindow
import mdfreader

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
        self.fileNames=[]
        self.mdf=mdfreader.mdf()
        self.mdfinfo=mdfreader.mdfinfo()
        self.convertSelection='Matlab'
    
    @pyqtSignature("")
    def on_browse_clicked(self):
        """
        Will open a dialog to browse for files
        """
        self.fileNames=QFileDialog.getOpenFileNames(self, "Select Files", filter=("MDF file (*.dat)"))
        if not self.fileNames.isEmpty():
            self.FileList.clear()
            self.FileList.addItems(self.fileNames)
            self.mdfinfo.readinfo(str(self.fileNames[0]))
            self.channelList.clear()
            self.channelList.addItems(self.mdfinfo.channelNameList)
    
    @pyqtSignature("")
    def on_Convert_clicked(self):
        """
       Will convert mdf files into select format
        """
        # Process all mdf files recursively
        for i in range(len(self.fileNames)):
            self.mdf.fileName=str(self.fileNames[i])
            self.mdf.read(self.mdf.fileName)
            #resample if requested
            if self.resample.checkState():
                if not self.resampleValue.isEmpty():
                    self.mdf.resample(self.resampleValue)
            if self.convertSelection=='Matlab':
                self.mdf.exportToMatlab()
            elif self.convertSelection=='csv':
                self.mdf.exportToCSV()
            elif self.convertSelection=='netcdf':
                self.mdf.exportToNetCDF()
            elif self.convertSelection=='hdf5':
                self.mdf.exportToHDF5()
            elif self.convertSelection=='excel':
                self.mdf.exportToExcel()
    
    @pyqtSignature("QListWidgetItem*")
    def on_FileList_itemClicked(self, item):
        """
        If user click on file list
        """
        # Refresh list of channels from selected file
        self.mdfinfo.readinfo(item)
        self.channelList.clear()
        self.channelList.addItems(self.mdfinfo.channelNameList)
    
    @pyqtSignature("bool")
    def on_matlab_clicked(self, checked):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        self.convertSelection='Matlab'
    
    @pyqtSignature("bool")
    def on_netcdf_clicked(self, checked):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        self.convertSlection='netcdf'
    
    @pyqtSignature("bool")
    def on_hdf5_clicked(self, checked):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        self.convertSelection='hdf5'
    
    @pyqtSignature("bool")
    def on_csv_clicked(self, checked):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        self.convertSelection='csv'
    
    @pyqtSignature("bool")
    def on_excel_clicked(self, checked):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        self.convertSelection='excel'
    
    @pyqtSignature("bool")
    def on_resample_toggled(self, checked):
        """
        Slot documentation goes here.
        """
        # Activates or deactivates input for resampling time
        self.resampleValue.enabled()
