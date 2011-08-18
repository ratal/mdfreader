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
    
    @pyqtSignature("")
    def on_browse_clicked(self):
        """
        Slot documentation goes here.
        """
        self.file = str(QFileDialog.getOpenFileName(self, "Select File", ("MDF file (*.dat)")))
        self.mdf=mdfreader.mdf(self.file)
    
    @pyqtSignature("")
    def on_Convert_clicked(self):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        raise NotImplementedError
