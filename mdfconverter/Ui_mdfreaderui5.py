# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'mdfreaderui.ui'
#
# Created by: PyQt5 UI code generator 5.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 653)
        self.TopLayout = QtWidgets.QWidget(MainWindow)
        self.TopLayout.setObjectName("TopLayout")
        self.gridLayout = QtWidgets.QGridLayout(self.TopLayout)
        self.gridLayout.setObjectName("gridLayout")
        self.label_2 = QtWidgets.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.TopLayout)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 0, 2, 1, 1)
        self.Lists = QtWidgets.QSplitter(self.TopLayout)
        self.Lists.setOrientation(QtCore.Qt.Horizontal)
        self.Lists.setObjectName("Lists")
        self.FileList = QtWidgets.QListWidget(self.Lists)
        self.FileList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.FileList.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.FileList.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.FileList.setDragEnabled(True)
        self.FileList.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.FileList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.FileList.setAlternatingRowColors(True)
        self.FileList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.FileList.setProperty("isWrapping", False)
        self.FileList.setResizeMode(QtWidgets.QListView.Adjust)
        self.FileList.setObjectName("FileList")
        self.channelList = QtWidgets.QListWidget(self.Lists)
        self.channelList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.channelList.setAcceptDrops(True)
        self.channelList.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.channelList.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.channelList.setDragEnabled(True)
        self.channelList.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.channelList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.channelList.setAlternatingRowColors(True)
        self.channelList.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.channelList.setObjectName("channelList")
        self.SelectedChannelList = QtWidgets.QListWidget(self.Lists)
        self.SelectedChannelList.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.SelectedChannelList.setAcceptDrops(True)
        self.SelectedChannelList.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.SelectedChannelList.setEditTriggers(
            QtWidgets.QAbstractItemView.NoEditTriggers
        )
        self.SelectedChannelList.setDragEnabled(True)
        self.SelectedChannelList.setDragDropMode(QtWidgets.QAbstractItemView.DragDrop)
        self.SelectedChannelList.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.SelectedChannelList.setAlternatingRowColors(True)
        self.SelectedChannelList.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        self.SelectedChannelList.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectItems
        )
        self.SelectedChannelList.setObjectName("SelectedChannelList")
        self.gridLayout.addWidget(self.Lists, 1, 0, 4, 3)
        self.browse = QtWidgets.QPushButton(self.TopLayout)
        self.browse.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.browse.setAutoDefault(False)
        self.browse.setDefault(False)
        self.browse.setObjectName("browse")
        self.gridLayout.addWidget(self.browse, 1, 3, 1, 1)
        self.Options = QtWidgets.QSplitter(self.TopLayout)
        self.Options.setOrientation(QtCore.Qt.Vertical)
        self.Options.setObjectName("Options")
        self.splitter = QtWidgets.QSplitter(self.Options)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName("splitter")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.splitter)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.ConvertSelect = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.ConvertSelect.setContentsMargins(0, 0, 0, 0)
        self.ConvertSelect.setObjectName("ConvertSelect")
        self.matlab = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.matlab.setEnabled(True)
        self.matlab.setChecked(True)
        self.matlab.setObjectName("matlab")
        self.ConvertSelect.addWidget(self.matlab)
        self.netcdf = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.netcdf.setObjectName("netcdf")
        self.ConvertSelect.addWidget(self.netcdf)
        self.hdf5 = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.hdf5.setChecked(False)
        self.hdf5.setObjectName("hdf5")
        self.ConvertSelect.addWidget(self.hdf5)
        self.csv = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.csv.setObjectName("csv")
        self.ConvertSelect.addWidget(self.csv)
        self.excel = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.excel.setObjectName("excel")
        self.ConvertSelect.addWidget(self.excel)
        self.excel2010 = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.excel2010.setObjectName("excel2010")
        self.ConvertSelect.addWidget(self.excel2010)
        self.mdf3 = QtWidgets.QRadioButton(self.verticalLayoutWidget)
        self.mdf3.setObjectName("mdf3")
        self.ConvertSelect.addWidget(self.mdf3)
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.splitter)
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.Resample = QtWidgets.QVBoxLayout(self.horizontalLayoutWidget)
        self.Resample.setContentsMargins(0, 0, 0, 0)
        self.Resample.setObjectName("Resample")
        self.resample = QtWidgets.QCheckBox(self.horizontalLayoutWidget)
        self.resample.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.resample.setChecked(True)
        self.resample.setObjectName("resample")
        self.Resample.addWidget(self.resample)
        self.resampleValue = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.resampleValue.setEnabled(True)
        self.resampleValue.setInputMethodHints(
            QtCore.Qt.ImhFormattedNumbersOnly | QtCore.Qt.ImhPreferNumbers
        )
        self.resampleValue.setObjectName("resampleValue")
        self.Resample.addWidget(self.resampleValue)
        self.layoutWidget = QtWidgets.QWidget(self.splitter)
        self.layoutWidget.setObjectName("layoutWidget")
        self.LabFile_2 = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.LabFile_2.setContentsMargins(0, 0, 0, 0)
        self.LabFile_2.setObjectName("LabFile_2")
        self.LabFileBrowse = QtWidgets.QPushButton(self.layoutWidget)
        self.LabFileBrowse.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.LabFileBrowse.setObjectName("LabFileBrowse")
        self.LabFile_2.addWidget(self.LabFileBrowse)
        self.LabFile = QtWidgets.QLineEdit(self.layoutWidget)
        self.LabFile.setEnabled(True)
        self.LabFile.setObjectName("LabFile")
        self.LabFile_2.addWidget(self.LabFile)
        self.MergeFiles = QtWidgets.QCheckBox(self.Options)
        self.MergeFiles.setEnabled(True)
        self.MergeFiles.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.MergeFiles.setChecked(True)
        self.MergeFiles.setObjectName("MergeFiles")
        self.gridLayout.addWidget(self.Options, 2, 3, 1, 1)
        self.Convert = QtWidgets.QPushButton(self.TopLayout)
        self.Convert.setLocale(
            QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.Belgium)
        )
        self.Convert.setObjectName("Convert")
        self.gridLayout.addWidget(self.Convert, 3, 3, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(
            113, 20, QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum
        )
        self.gridLayout.addItem(spacerItem, 4, 3, 1, 1)
        MainWindow.setCentralWidget(self.TopLayout)
        self.label_2.setBuddy(self.FileList)
        self.label.setBuddy(self.channelList)
        self.label_3.setBuddy(self.SelectedChannelList)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MDF Converter"))
        self.label_2.setText(_translate("MainWindow", "File List"))
        self.label.setText(_translate("MainWindow", "File Channel List"))
        self.label_3.setText(_translate("MainWindow", "Selected Channel List"))
        self.FileList.setWhatsThis(
            _translate("MainWindow", "File list to be converted")
        )
        self.FileList.setSortingEnabled(False)
        self.channelList.setWhatsThis(
            _translate("MainWindow", "Channel list inside the selected file")
        )
        self.channelList.setSortingEnabled(True)
        self.SelectedChannelList.setWhatsThis(
            _translate("MainWindow", "Selected channel list to be exported")
        )
        self.SelectedChannelList.setSortingEnabled(True)
        self.browse.setToolTip(
            _translate("MainWindow", "Click and select MDF file for conversion")
        )
        self.browse.setWhatsThis(
            _translate("MainWindow", "Click to browse for MDF files to be converted")
        )
        self.browse.setText(_translate("MainWindow", "Browse"))
        self.matlab.setText(_translate("MainWindow", "Matlab .mat"))
        self.netcdf.setText(_translate("MainWindow", "netCDF"))
        self.hdf5.setText(_translate("MainWindow", "HDF5"))
        self.csv.setText(_translate("MainWindow", "CSV"))
        self.excel.setText(_translate("MainWindow", "Excel 2003"))
        self.excel2010.setText(_translate("MainWindow", "Excel 2010"))
        self.mdf3.setText(_translate("MainWindow", "MDF3.3"))
        self.resample.setWhatsThis(
            _translate(
                "MainWindow",
                "Click to resample according to below sampling time in seconds",
            )
        )
        self.resample.setText(_translate("MainWindow", "Resample"))
        self.resampleValue.setWhatsThis(
            _translate("MainWindow", "Resampling time in seconds")
        )
        self.resampleValue.setText(_translate("MainWindow", "1.0"))
        self.LabFileBrowse.setWhatsThis(
            _translate("MainWindow", "Click to selected file containing channel list")
        )
        self.LabFileBrowse.setText(_translate("MainWindow", "Lab file"))
        self.LabFile.setWhatsThis(_translate("MainWindow", "Chosen lab file"))
        self.MergeFiles.setText(_translate("MainWindow", "MergeFiles"))
        self.Convert.setWhatsThis(
            _translate(
                "MainWindow",
                "Click to convert all files according your selected options",
            )
        )
        self.Convert.setText(_translate("MainWindow", "Convert"))


if __name__ == "__main__":
    import sys

    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
