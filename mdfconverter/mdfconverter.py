from PyQt4 import QtGui
try:
    from .mdfreaderui import MainWindow
except:
    from mdfreaderui import MainWindow
from sys import argv, exit
from multiprocessing import freeze_support

def main():
    freeze_support()
    app = QtGui.QApplication(argv)
    myapp = MainWindow()
    myapp.show()
    exit(app.exec_())

if __name__ == "__main__":
    main()
