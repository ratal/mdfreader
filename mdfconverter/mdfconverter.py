from sys import argv, exit, path
from os.path import dirname, abspath
root = dirname(abspath(__file__))
path.append(root)

from multiprocessing import freeze_support

from PyQt5.QtWidgets import QApplication
from mdfreaderui5 import MainWindow


def main():
    freeze_support()
    app = QApplication(argv)
    myapp = MainWindow()
    myapp.show()
    exit(app.exec_())


if __name__ == "__main__":
    main()
