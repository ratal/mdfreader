# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '/home/ratal/workspace/mdfreader/mdfreaderui.ui'
#
# Created: Thu Aug 18 22:24:25 2011
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(800, 600)
        self.centralWidget = QtGui.QWidget(MainWindow)
        self.centralWidget.setObjectName("centralWidget")
        self.Convert = QtGui.QPushButton(self.centralWidget)
        self.Convert.setGeometry(QtCore.QRect(660, 530, 97, 27))
        self.Convert.setObjectName("Convert")
        self.FileList = QtGui.QListView(self.centralWidget)
        self.FileList.setGeometry(QtCore.QRect(20, 50, 291, 511))
        self.FileList.setObjectName("FileList")
        self.channelList = QtGui.QListView(self.centralWidget)
        self.channelList.setGeometry(QtCore.QRect(340, 50, 271, 511))
        self.channelList.setObjectName("channelList")
        self.browse = QtGui.QCommandLinkButton(self.centralWidget)
        self.browse.setGeometry(QtCore.QRect(30, 10, 101, 31))
        self.browse.setToolTip("")
        self.browse.setAutoDefault(False)
        self.browse.setDefault(False)
        self.browse.setObjectName("browse")
        self.splitter = QtGui.QSplitter(self.centralWidget)
        self.splitter.setGeometry(QtCore.QRect(630, 10, 151, 165))
        self.splitter.setOrientation(QtCore.Qt.Vertical)
        self.splitter.setObjectName("splitter")
        self.verticalLayoutWidget = QtGui.QWidget(self.splitter)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtGui.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.matlab = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.matlab.setEnabled(True)
        self.matlab.setObjectName("matlab")
        self.verticalLayout.addWidget(self.matlab)
        self.netcdf = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.netcdf.setObjectName("netcdf")
        self.verticalLayout.addWidget(self.netcdf)
        self.hdf5 = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.hdf5.setObjectName("hdf5")
        self.verticalLayout.addWidget(self.hdf5)
        self.csv = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.csv.setObjectName("csv")
        self.verticalLayout.addWidget(self.csv)
        self.excel = QtGui.QRadioButton(self.verticalLayoutWidget)
        self.excel.setObjectName("excel")
        self.verticalLayout.addWidget(self.excel)
        self.horizontalLayoutWidget = QtGui.QWidget(self.splitter)
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtGui.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.resample = QtGui.QCheckBox(self.horizontalLayoutWidget)
        self.resample.setObjectName("resample")
        self.horizontalLayout.addWidget(self.resample)
        self.lineEdit = QtGui.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit.setEnabled(False)
        self.lineEdit.setObjectName("lineEdit")
        self.horizontalLayout.addWidget(self.lineEdit)
        self.FileList_2 = QtGui.QPushButton(self.centralWidget)
        self.FileList_2.setGeometry(QtCore.QRect(630, 210, 105, 24))
        self.FileList_2.setObjectName("FileList_2")
        self.FileListEdit = QtGui.QLineEdit(self.centralWidget)
        self.FileListEdit.setEnabled(False)
        self.FileListEdit.setGeometry(QtCore.QRect(630, 260, 113, 25))
        self.FileListEdit.setObjectName("FileListEdit")
        MainWindow.setCentralWidget(self.centralWidget)
        self.actionBrowse = QtGui.QAction(MainWindow)
        self.actionBrowse.setObjectName("actionBrowse")

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.resample, QtCore.SIGNAL("clicked()"), self.lineEdit.show)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MDF Converter", None, QtGui.QApplication.UnicodeUTF8))
        self.Convert.setText(QtGui.QApplication.translate("MainWindow", "Convert", None, QtGui.QApplication.UnicodeUTF8))
        self.browse.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.matlab.setText(QtGui.QApplication.translate("MainWindow", "Matlab .mat", None, QtGui.QApplication.UnicodeUTF8))
        self.netcdf.setText(QtGui.QApplication.translate("MainWindow", "netCDF", None, QtGui.QApplication.UnicodeUTF8))
        self.hdf5.setText(QtGui.QApplication.translate("MainWindow", "HDF5", None, QtGui.QApplication.UnicodeUTF8))
        self.csv.setText(QtGui.QApplication.translate("MainWindow", "CSV", None, QtGui.QApplication.UnicodeUTF8))
        self.excel.setText(QtGui.QApplication.translate("MainWindow", "Excel", None, QtGui.QApplication.UnicodeUTF8))
        self.resample.setText(QtGui.QApplication.translate("MainWindow", "Resample", None, QtGui.QApplication.UnicodeUTF8))
        self.FileList_2.setText(QtGui.QApplication.translate("MainWindow", "File List", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBrowse.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBrowse.setToolTip(QtGui.QApplication.translate("MainWindow", "Open MDF file", None, QtGui.QApplication.UnicodeUTF8))
        self.actionBrowse.setShortcut(QtGui.QApplication.translate("MainWindow", "Ctrl+O", None, QtGui.QApplication.UnicodeUTF8))


if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

