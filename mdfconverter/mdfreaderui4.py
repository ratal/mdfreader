# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""

from sys import version_info, path
from os.path import dirname, abspath
root = dirname(abspath(__file__))
path.append(root)

from io import open
from multiprocessing import Pool, cpu_count
from mdfreader import *

from PyQt4.QtGui import QMainWindow, QFileDialog, QAction
from PyQt4.QtCore import pyqtSignature, SIGNAL
from Ui_mdfreaderui4 import Ui_MainWindow

PythonVersion = version_info
PythonVersion = PythonVersion[0]
MultiProc = True  # multiprocess switch, for debug purpose put False


class MainWindow(QMainWindow, Ui_MainWindow, QFileDialog):

    """
    Class documentation goes here.
    """

    def __init__(self, parent=None):
        """
        Constructor
        """
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.fileNames = []  # files to convert
        self.mdfClass = mdf()  # instance of mdf
        self.mdfinfoClass = mdfinfo()  # instance of mdfinfo
        self.convertSelection = 'Matlab'  # by default Matlab conversion is selected
        self.MergeFiles = False  # by default
        self.labFileName = []  # .lab file name
        self.defaultPath = None  # default path to open for browsing files
        self.actionPlotSelectedChannel = QAction("Plot", self.SelectedChannelList)  # context menu to allow plot of channel
        self.SelectedChannelList.addAction(self.actionPlotSelectedChannel)
        self.connect(self.actionPlotSelectedChannel, SIGNAL("triggered()"), self.plotSelected)
        self.actionPlotChannel = QAction("Plot", self.channelList)  # context menu to allow plot of channel
        self.channelList.addAction(self.actionPlotChannel)
        self.connect(self.actionPlotChannel, SIGNAL("triggered()"), self.plot)
        self.actionFileRemove = QAction("Delete", self.FileList)  # context menu to remove selected file from list
        self.FileList.addAction(self.actionFileRemove)
        self.connect(self.actionFileRemove, SIGNAL("triggered()"), self.FileRemove)

    @pyqtSignature("")
    def on_browse_clicked(self):
        """
        Will open a dialog to browse for files
        """
        if self.defaultPath is None:
            self.fileNames = QFileDialog.getOpenFileNames(self, "Select Measurement Files", filter=("MDF file (*.dat *.mdf *.mf4 *.mfx *.mfxz)"))
            self.defaultPath = dirname(str(self.fileNames[0]))
        else:
            self.fileNames = QFileDialog.getOpenFileNames(self, "Select Measurement Files", self.defaultPath, filter=("MDF file (*.dat *.mdf *.mf4 *.mfx *.mfxz)"))
        if not len(self.fileNames) == 0:
            self.FileList.addItems(self.fileNames)
            self.mdfinfoClass.__init__()
            self.cleanChannelList()
            self.cleanSelectedChannelList()
            ChannelList = convertChannelList(self.mdfinfoClass.listChannels(str(self.fileNames[0])))
            self.SelectedChannelList.addItems(ChannelList)
            self.FileList.setItemSelected(self.FileList.item(0), True)

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
        channelList = []
        [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
        channelList = list(set(channelList))  # remove duplicates
        # Process all mdf files recursively
        if self.FileList.count() > 0:  # not empty list
            ncpu = cpu_count()  # to still have response from PC
            if ncpu < 1:
                ncpu = 1
            pool = Pool(processes=ncpu)
            if self.MergeFiles or self.FileList.count() < 2:  # export all files separately, inverted bool
                convertFlag = True
                convertSelection = self.convertSelection
                resampleValue = float(self.resampleValue.text())
                # re-sample if requested
                if self.resample.checkState():
                    if not len(self.resampleValue.text()) == 0:
                        resampleFlag = True
                    else:
                        raise Exception('Empty field for resampling')
                else:
                    resampleFlag = False
                args = [(str(self.FileList.takeItem(0).text()), channelList, resampleFlag, resampleValue, convertFlag, convertSelection) for i in range(self.FileList.count())]
                if MultiProc:
                    result = pool.map_async(processMDFstar, args)
                    result.get()  # waits until finished
                else:
                    result = list(map(processMDFstar, args))
                self.cleanChannelList()
            elif self.FileList.count() >= 2:  # Stack files data if min 2 files in list
                # import first file
                fileNameList = []
                if len(self.resampleValue.text()) == 0:
                    raise Exception('Wrong value for re-sampling')
                convertFlag = False
                convertSelection = self.convertSelection
                resampleValue = float(self.resampleValue.text())
                resampleFlag = True  # always re-sample when merging
                fileName = str(self.FileList.item(0).text())  # Uses first file name for the converted file
                # list filenames
                args = [(str(self.FileList.takeItem(0).text()), channelList, resampleFlag, resampleValue, convertFlag, convertSelection) for i in range(self.FileList.count())]
                if MultiProc:
                    res = pool.map_async(processMDFstar, args)
                    result = res.get()
                else:
                    result = map(processMDFstar, args)  # no multiprocess, for debug
                # Merge results
                self.mdfClass.__init__()  # clear memory
                self.mdfClass.fileName = fileName  # First filename will be used for exported file name
                self.mdfClass.multiProc = False  # do not use multiproc inside mdfreader while already using from mdfreaderui level
                buffer = self.mdfClass.copy()  # create/copy empty class in buffer
                res = result.pop(0)  # extract first file data from processed list
                self.mdfClass.update(res[0])  # initialize mdfclass wih first file data
                self.mdfClass.masterChannelList = res[1]  # initialize masterChannelList
                fileNameList.append(res[2])  # record merged file in list
                for res in result:  # Merge
                    buffer.__init__()  # clean buffer class
                    buffer.update(res[0])  # assigns next class to buffer
                    buffer.masterChannelList = res[1]
                    fileNameList.append(res[2])
                    self.mdfClass.mergeMdf(buffer)  # merge buffer to merged class mdfClass
                # Export
                if self.convertSelection == 'Matlab':
                    self.mdfClass.exportToMatlab()
                elif self.convertSelection == 'csv':
                    self.mdfClass.exportToCSV()
                elif self.convertSelection == 'netcdf':
                    self.mdfClass.exportToNetCDF()
                elif self.convertSelection == 'hdf5':
                    self.mdfClass.exportToHDF5()
                elif self.convertSelection == 'excel':
                    self.mdfClass.exportToExcel()
                elif self.convertSelection == 'excel2010':
                    self.mdfClass.exportToXlsx()
                elif self.convertSelection == 'mdf3':
                    self.mdfClass.write(fileName + '_new')
                self.cleanChannelList()
                print('File list merged :')
                for file in fileNameList:  # prints files merged for checking
                    print(file)
                self.mdfClass.__init__()  # clear memory

    @pyqtSignature("QListWidgetItem*")
    def on_FileList_itemClicked(self, item):
        """
        If user click on file list
        """
        # Refresh list of channels from selected file
        self.mdfinfoClass.__init__()
        # self.mdfinfoClass.readinfo(item)
        self.cleanChannelList()
        ChannelList = convertChannelList(self.mdfinfoClass.listChannels(str(item.text())))
        self.channelList.addItems(ChannelList)
        self.mdfinfoClass.__init__()  # clean object to free memory

    @pyqtSignature("bool")
    def on_matlab_clicked(self, checked):
        """
        Selects Matlab conversion
        """
        self.convertSelection = 'Matlab'

    @pyqtSignature("bool")
    def on_netcdf_clicked(self, checked):
        """
        Selects netcdf conversion.
        """
        self.convertSelection = 'netcdf'

    @pyqtSignature("bool")
    def on_hdf5_clicked(self, checked):
        """
        Selects hdf5 conversion.
        """
        self.convertSelection = 'hdf5'

    @pyqtSignature("bool")
    def on_csv_clicked(self, checked):
        """
        Selects csv conversion.
        """
        self.convertSelection = 'csv'

    @pyqtSignature("bool")
    def on_excel_clicked(self, checked):
        """
        Selects excel conversion.
        """
        self.convertSelection = 'excel'

    @pyqtSignature("bool")
    def on_excel2010_clicked(self, checked):
        """
        Selects excel conversion.
        """
        self.convertSelection = 'excel2010'

    @pyqtSignature("bool")
    def on_mdf3_clicked(self, checked):
        """
        Selects MDF3.3 conversion.
        """
        self.convertSelection = 'mdf3'

    @pyqtSignature("")
    def on_LabFileBrowse_clicked(self):
        """
        selects lab file from browser.
        """
        self.labFileName = QFileDialog.getOpenFileName(self, "Select Lab Files", filter=("Lab file (*.lab)"))
        if not len(self.labFileName) == 0:
            self.LabFile.del_()  # clear linedit
            self.LabFile.insert(str(self.labFileName))  # replace linedit field by browsed file name
            # read lab file
            labfile = open(str(self.labFileName), 'r')
            self.labChannelList = []
            line = labfile.readline()  # read first line [lab]
            while True:
                line = labfile.readline()
                if not line:
                    break
                self.labChannelList.append(line.replace('\n', ''))
            self.cleanSelectedChannelList()  # Clear Selected file list
            self.SelectedChannelList.addItems(self.labChannelList)

    def plot(self):
        # Finds selected file and read it
        selectedFile = self.FileList.selectedItems()
        self.mdfClass.__init__(str(selectedFile[0].text()))  # read file
        # list items selected in listWidget
        Channels = self.channelList.selectedItems()
        selectedChannels = []
        [selectedChannels.append(str(Channels[i].text())) for i in range(len(Channels))]
        # plot channels
        self.mdfClass.plot(selectedChannels)

    def plotSelected(self):
        # plots channels from selected list
        selectedFile = self.FileList.selectedItems()
        if not len(selectedFile) == 0:
            self.mdfClass.__init__(str(selectedFile[0].text()))  # read file
        else:
            self.mdfClass.__init__(str(self.FileList[0].text()))  # read file
        # list items selected in listWidget
        Channels = self.SelectedChannelList.selectedItems()
        selectedChannels = []
        [selectedChannels.append(str(Channels[i].text())) for i in range(len(Channels))]
        # plot channels
        self.mdfClass.plot(selectedChannels)

    def FileRemove(self):
        # removes selected file
        selectionList = self.FileList.selectedItems()
        [self.FileList.takeItem(self.FileList.row(selectionList[i])) for i in range(len(selectionList))]

    def on_SelectedChannelList_dropEvent(self):
        # avoids to have duplicates in list when channel is dropped
        channelList = []
        [channelList.append(str(self.SelectedChannelList.item(i).text())) for i in range(self.SelectedChannelList.count())]
        channelList = list(set(channelList))  # removeDuplicates
        self.SelectedChannelList.clear()
        self.SelectedChannelList.addItems(channelList)

    @pyqtSignature("bool")
    def on_MergeFiles_toggled(self, checked):
        """
        Slot documentation goes here.
        """
        # toggle flag to merge files
        self.MergeFiles = not self.MergeFiles
        if self.MergeFiles:
            self.resample.setCheckState(2)


def processMDF(fileName, channelist, resampleFlag, resampleValue, convertFlag, convertSelection):
    # Will process file according to defined options
    yop = mdf()
    yop.multiProc = False  # already multiprocessed
    yop.convertAfterRead = True
    yop.read(fileName)  # reads complete file
    yop.keepChannels(channelist)  # removes unnecessary channels
    if resampleFlag:
        yop.resample(resampleValue)
    if convertFlag:
        if convertSelection == 'Matlab':
            yop.exportToMatlab()
        elif convertSelection == 'csv':
            yop.exportToCSV()
        elif convertSelection == 'netcdf':
            yop.exportToNetCDF()
        elif convertSelection == 'hdf5':
            yop.exportToHDF5()
        elif convertSelection == 'excel':
            yop.exportToExcel()
        elif convertSelection == 'excel2010':
            yop.exportToXlsx()
        elif convertSelection == 'mdf3':
            yop.write(fileName + '_new')
    yopPicklable = {}  # picklable dict and not object
    for channel in list(yop.keys()):
        yopPicklable[channel] = yop[channel]
    return [yopPicklable, yop.masterChannelList, yop.fileName]


def processMDFstar(args):
    try:
        return processMDF(*args)
    except:
        print('Error, following file might be corrupted : ' + args[0])  # Shows fileName and parameters to help finding corrupted files
        raise Exception('Please re-try by removing this file from the list and restart mdfconverter to kill processes and clean memory')


def convertChannelList(channelList):
    if PythonVersion < 3:
        return [str(name) for name in channelList]
    else:
        return [(name) for name in channelList]
