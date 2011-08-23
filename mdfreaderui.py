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
        self.mdfinfoClass=mdfreader.mdfinfo()
        self.convertSelection='Matlab'
        self.labFileName=[]
    
    @pyqtSignature("")
    def on_browse_clicked(self):
        """
        Will open a dialog to browse for files
        """
        self.fileNames=QFileDialog.getOpenFileNames(self, "Select Measurement Files", filter=("MDF file (*.dat)"))
        if not self.fileNames.isEmpty():
            self.FileList.clear()
            self.FileList.addItems(self.fileNames)
            self.mdfinfoClass.readinfo(str(self.fileNames[0]))
            self.channelList.clear()
            self.channelList.addItems(self.mdfinfoClass.channelNameList)
    
    @pyqtSignature("")
    def on_Convert_clicked(self):
        """
       Will convert mdf files into selected format
        """
        # Process all mdf files recursively
        mdf=mdfreader.mdf()
        for i in range(len(self.fileNames)):
            mdf.fileName=str(self.fileNames[i])
            mdf.read(mdf.fileName)
            #resample if requested
            if self.resample.checkState():
                if not self.resampleValue.text().isEmpty():
                    mdf.resample(float(self.resampleValue.text()))
            if self.convertSelection=='Matlab':
                mdf.exportToMatlab(mdf.fileName)
            elif self.convertSelection=='csv':
                mdf.exportToCSV(mdf.fileName)
            elif self.convertSelection=='netcdf':
                mdf.exportToNetCDF(mdf.fileName)
            elif self.convertSelection=='hdf5':
                mdf.exportToHDF5(mdf.fileName)
            elif self.convertSelection=='excel':
                mdf.exportToExcel(mdf.fileName)
    
    @pyqtSignature("QListWidgetItem*")
    def on_FileList_itemClicked(self, item):
        """
        If user click on file list
        """
        # Refresh list of channels from selected file
        self.mdfinfoClass.readinfo(item)
        self.channelList.clear()
        self.channelList.addItems(self.mdfinfoClass.channelNameList)
    
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
            self.LabFile.del_()
            self.LabFile.insert(str(self.labFileName))
