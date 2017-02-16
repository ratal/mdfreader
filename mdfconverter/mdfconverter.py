from PyQt4 import QtGui

from sys import argv, exit, path
from os.path import dirname, abspath
root = dirname(abspath(__file__))
path.append(root)
from mdfreaderui import MainWindow
from multiprocessing import freeze_support

def main():
    freeze_support()
    app = QtGui.QApplication(argv)
    myapp = MainWindow()
    myapp.show()
    exit(app.exec_())

if __name__ == "__main__":
    main()
