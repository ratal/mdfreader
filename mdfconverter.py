from PyQt4 import QtCore, QtGui
from mdfreaderui import MainWindow
from sys import argv, exit, platform
from multiprocessing import freeze_support

if __name__ == "__main__":
    freeze_support()
    app = QtGui.QApplication( argv )
    myapp = MainWindow(  )
    myapp.show()
    exit( app.exec_() )