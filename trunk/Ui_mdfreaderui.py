# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/ratal/workspace/mdfreader/mdfreaderui.ui'
#
# Created: Sat Apr 12 00:13:31 2014
#      by: PyQt4 UI code generator 4.9.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 653)
        self.TopLayout = QtGui.QWidget(MainWindow)
        self.TopLayout.setObjectName(_fromUtf8("TopLayout"))
        self.gridLayout = QtGui.QGridLayout(self.TopLayout)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_2 = QtGui.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label = QtGui.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 1, 1, 1)
        self.label_3 = QtGui.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 0, 2, 1, 1)
        self.Lists = QtGui.QSplitter(self.TopLayout)
        self.Lists.setOrientation(QtCore.Qt.Horizontal)
        self.Lists.setObjectName(_fromUtf8("Lists"))
        self.FileList = QtGui.QListWidget(self.Lists)
        self.FileList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.FileList.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.FileList.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.FileList.setDragEnabled(True)
        self.FileList.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.FileList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.FileList.setAlternatingRowColors(True)
        self.FileList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.FileList.setProperty("isWrapping", False)
        self.FileList.setResizeMode(QtGui.QListView.Adjust)
        self.FileList.setObjectName(_fromUtf8("FileList"))
        self.channelList = QtGui.QListWidget(self.Lists)
        self.channelList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.channelList.setAcceptDrops(True)
        self.channelList.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.channelList.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.channelList.setDragEnabled(True)
        self.channelList.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.channelList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.channelList.setAlternatingRowColors(True)
        self.channelList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.channelList.setObjectName(_fromUtf8("channelList"))
        self.SelectedChannelList = QtGui.QListWidget(self.Lists)
        self.SelectedChannelList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.SelectedChannelList.setAcceptDrops(True)
        self.SelectedChannelList.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.SelectedChannelList.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)
        self.SelectedChannelList.setDragEnabled(True)
        self.SelectedChannelList.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.SelectedChannelList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.SelectedChannelList.setAlternatingRowColors(True)
        self.SelectedChannelList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.SelectedChannelList.setSelectionBehavior(QtGui.QAbstractItemView.SelectItems)
        self.SelectedChannelList.setObjectName(_fromUtf8("SelectedChannelList"))
        self.gridLayout.addWidget(self.Lists, 1, 0, 4, 3)
        self.browse = QtGui.QPushButton(self.TopLayout)
        self.browse.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.browse.setAutoDefault(False)
        self.browse.setDefault(False)
        self.browse.setObjectName(_fromUtf8("browse"))
        self.gridLayout.addWidget(self.browse, 1, 3, 1, 1)
        self.Options = QtGui.QSplitter(self.TopLayout)
        self.Options.setOrientation(QtCore.Qt.Vertical)
        self.Options.setObjectName(_fromUtf8("Options"))
        self.splitter = QtGui.QSplitter(self.Options)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName(_fromUtf8("splitter"))
        self.verticalLayoutWidget = QtGui.QWidget(self.splitter)
        self.verticalLayoutWidget.setObjectName(_fromUtf8("verticalLayoutWidget"))
        self.ConvertSelect = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.ConvertSelect.setMargin(0)
        self.ConvertSelect.setObjectName(_fromUtf8("ConvertSelect"))
        self.matlab = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.matlab.setEnabled(True)
        self.matlab.setChecked(True)
        self.matlab.setObjectName(_fromUtf8("matlab"))
        self.ConvertSelect.addWidget(self.matlab)
        self.netcdf = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.netcdf.setObjectName(_fromUtf8("netcdf"))
        self.ConvertSelect.addWidget(self.netcdf)
        self.hdf5 = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.hdf5.setChecked(False)
        self.hdf5.setObjectName(_fromUtf8("hdf5"))
        self.ConvertSelect.addWidget(self.hdf5)
        self.csv = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.csv.setObjectName(_fromUtf8("csv"))
        self.ConvertSelect.addWidget(self.csv)
        self.excel = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.excel.setObjectName(_fromUtf8("excel"))
        self.ConvertSelect.addWidget(self.excel)
        self.excel2010 = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.excel2010.setObjectName(_fromUtf8("excel2010"))
        self.ConvertSelect.addWidget(self.excel2010)
        self.mdf3 = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.mdf3.setObjectName(_fromUtf8("mdf3"))
        self.ConvertSelect.addWidget(self.mdf3)
        self.horizontalLayoutWidget = QtGui.QWidget(self.splitter)
        self.horizontalLayoutWidget.setObjectName(_fromUtf8("horizontalLayoutWidget"))
        self.Resample = QtGui.QVBoxLayout(self.horizontalLayoutWidget)
        self.Resample.setMargin(0)
        self.Resample.setObjectName(_fromUtf8("Resample"))
        self.resample = QtGui.QCheckBox(self.horizontalLayoutWidget)
        self.resample.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.resample.setChecked(True)
        self.resample.setObjectName(_fromUtf8("resample"))
        self.Resample.addWidget(self.resample)
        self.resampleValue = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.resampleValue.setEnabled(True)
        self.resampleValue.setInputMethodHints(QtCore.Qt.ImhFormattedNumbersOnly|QtCore.Qt.ImhPreferNumbers)
        self.resampleValue.setObjectName(_fromUtf8("resampleValue"))
        self.Resample.addWidget(self.resampleValue)
        self.layoutWidget = QtGui.QWidget(self.splitter)
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.LabFile_2 = QtGui.QVBoxLayout(self.layoutWidget)
        self.LabFile_2.setMargin(0)
        self.LabFile_2.setObjectName(_fromUtf8("LabFile_2"))
        self.LabFileBrowse = QtGui.QPushButton(self.layoutWidget)
        self.LabFileBrowse.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.LabFileBrowse.setObjectName(_fromUtf8("LabFileBrowse"))
        self.LabFile_2.addWidget(self.LabFileBrowse)
        self.LabFile = QtGui.QLineEdit(self.layoutWidget)
        self.LabFile.setEnabled(True)
        self.LabFile.setObjectName(_fromUtf8("LabFile"))
        self.LabFile_2.addWidget(self.LabFile)
        self.MergeFiles = QtGui.QCheckBox(self.Options)
        self.MergeFiles.setEnabled(True)
        self.MergeFiles.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.MergeFiles.setChecked(True)
        self.MergeFiles.setObjectName(_fromUtf8("MergeFiles"))
        self.gridLayout.addWidget(self.Options, 2, 3, 1, 1)
        self.Convert = QtGui.QPushButton(self.TopLayout)
        self.Convert.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium))
        self.Convert.setObjectName(_fromUtf8("Convert"))
        self.gridLayout.addWidget(self.Convert, 3, 3, 1, 1)
        spacerItem = QtGui.QSpacerItem(113, 20, QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 4, 3, 1, 1)
        MainWindow.setCentralWidget(self.TopLayout)
        self.label_2.setBuddy(self.FileList)
        self.label.setBuddy(self.channelList)
        self.label_3.setBuddy(self.SelectedChannelList)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MDF Converter", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "File List", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "File Channel List", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Selected Channel List", None, QtGui.QApplication.UnicodeUTF8))
        self.FileList.setWhatsThis(QtGui.QApplication.translate("MainWindow", "File list to be converted", None, QtGui.QApplication.UnicodeUTF8))
        self.FileList.setSortingEnabled(False)
        self.channelList.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Channel list inside the selected file", None, QtGui.QApplication.UnicodeUTF8))
        self.channelList.setSortingEnabled(True)
        self.SelectedChannelList.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Selected channel list to be exported", None, QtGui.QApplication.UnicodeUTF8))
        self.SelectedChannelList.setSortingEnabled(True)
        self.browse.setToolTip(QtGui.QApplication.translate("MainWindow", "Click and select MDF file for conversion", None, QtGui.QApplication.UnicodeUTF8))
        self.browse.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Click to browse for MDF files to be converted", None, QtGui.QApplication.UnicodeUTF8))
        self.browse.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.matlab.setText(QtGui.QApplication.translate("MainWindow", "Matlab .mat", None, QtGui.QApplication.UnicodeUTF8))
        self.netcdf.setText(QtGui.QApplication.translate("MainWindow", "netCDF", None, QtGui.QApplication.UnicodeUTF8))
        self.hdf5.setText(QtGui.QApplication.translate("MainWindow", "HDF5", None, QtGui.QApplication.UnicodeUTF8))
        self.csv.setText(QtGui.QApplication.translate("MainWindow", "CSV", None, QtGui.QApplication.UnicodeUTF8))
        self.excel.setText(QtGui.QApplication.translate("MainWindow", "Excel 2003", None, QtGui.QApplication.UnicodeUTF8))
        self.excel2010.setText(QtGui.QApplication.translate("MainWindow", "Excel 2010", None, QtGui.QApplication.UnicodeUTF8))
        self.mdf3.setText(QtGui.QApplication.translate("MainWindow", "MDF3.3", None, QtGui.QApplication.UnicodeUTF8))
        self.resample.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Click to resample according to below sampling time in seconds", None, QtGui.QApplication.UnicodeUTF8))
        self.resample.setText(QtGui.QApplication.translate("MainWindow", "Resample", None, QtGui.QApplication.UnicodeUTF8))
        self.resampleValue.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Resampling time in seconds", None, QtGui.QApplication.UnicodeUTF8))
        self.resampleValue.setText(QtGui.QApplication.translate("MainWindow", "1.0", None, QtGui.QApplication.UnicodeUTF8))
        self.LabFileBrowse.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Click to selected file containing channel list", None, QtGui.QApplication.UnicodeUTF8))
        self.LabFileBrowse.setText(QtGui.QApplication.translate("MainWindow", "Lab file", None, QtGui.QApplication.UnicodeUTF8))
        self.LabFile.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Chosen lab file", None, QtGui.QApplication.UnicodeUTF8))
        self.MergeFiles.setText(QtGui.QApplication.translate("MainWindow", "MergeFiles", None, QtGui.QApplication.UnicodeUTF8))
        self.Convert.setWhatsThis(QtGui.QApplication.translate("MainWindow", "Click to convert all files according your selected options", None, QtGui.QApplication.UnicodeUTF8))
        self.Convert.setText(QtGui.QApplication.translate("MainWindow", "Convert", None, QtGui.QApplication.UnicodeUTF8))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

