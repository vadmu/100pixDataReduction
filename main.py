#!/usr/bin/python

import sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np

# comment
# comment2

class Plot2dWidget(pg.GraphicsLayoutWidget):
    def __init__(self, *args):
        pg.GraphicsLayoutWidget.__init__(self, *args)
        self.view = self.addViewBox()
        self.view.setAspectLocked(True)
        self.img = pg.ImageItem(border='w')
        self.view.addItem(self.img)
        
    def setImage(self, img):
        self.img.setImage(img)



class MainWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)

        self.resize(800, 600)
        self.setWindowTitle('mainwindow')
        self.plot = Plot2dWidget()
        self.openFileButton = QtGui.QPushButton("Open File")
        self.button = QtGui.QPushButton("Plot")
                
        
        self.layout = QtGui.QHBoxLayout()
        self.layout2 = QtGui.QVBoxLayout()
        self.layout2.addWidget(self.openFileButton)
        self.layout2.addWidget(self.button)
        self.layout.addLayout(self.layout2)
        self.layout.addWidget(self.plot)

        self.setLayout(self.layout)
        self.openFileButton.clicked.connect(self.openFile)
        self.button.clicked.connect(self.openFile2)
        
    def openFile(self):
        image = np.random.randint(1024, size=(1000, 500))
        self.plot.setImage(image)
        
   
    def openFile2(self):
        image = np.random.randint(1024, size=(500, 500))
        self.plot.setImage(image)    
        
#        textEdit = QtGui.QTextEdit()
#        self.setCentralWidget(textEdit)
#
#        exit = QtGui.QAction('Exit', self)
#        exit.setShortcut('Ctrl+Q')
#        exit.setStatusTip('Exit application')
#        self.connect(exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))
#
#        self.statusBar()
#
#        menubar = self.menuBar()
#        file = menubar.addMenu('&File')
#        file.addAction(exit)
#
#        toolbar = self.addToolBar('Exit')
#        toolbar.addAction(exit)

app = QtGui.QApplication(sys.argv)
main = MainWindow()
main.show()
sys.exit(app.exec_())

#self.toolbar = self.addToolBar('Exit')
#self.toolbar.addAction(self.exit)