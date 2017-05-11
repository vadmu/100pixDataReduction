#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import h5py

class coolWidget(QtGui.QGroupBox):
    setActiveSignal = QtCore.pyqtSignal()
    enableNextSignal = QtCore.pyqtSignal()
    def __init__(self):
        QtGui.QGroupBox.__init__(self)
        self.setMouseTracking(True)        
        self.setStyleSheet(
        """
        QGroupBox{border: 2px solid #777; border-radius: 10px;}
        *[active = 'true']{background:#afa;}
        *[hover = 'true']{border: 2px solid #f00;}   
        """)
        self.setMinimumWidth(50)
        self.setProperty("active", False)
        self.setProperty("hover", False)
        self.setEnabled(False)
        
    def enterEvent(self,event):
        if self.isEnabled():
            self.setProperty("hover", True)
            self.updateStyle()

    def leaveEvent(self,event):
        if self.isEnabled():
            self.setProperty("hover", False)
            self.updateStyle()
            
    def updateStyle(self):
        self.setStyle(self.style())        
                          

class LoadFileWidget(coolWidget):
    def __init__(self):
        coolWidget.__init__(self)
        self.openFileButton = QtGui.QPushButton("Open File") 
        self.openFileButton.clicked.connect(self.openFile)
        self.lineEditFile = QtGui.QLineEdit()
        self.lineScan = QtGui.QSpinBox()
        self.lineScan.setMaximum(99999)
        self.line100pix = QtGui.QLineEdit()
        self.lineEnergy = QtGui.QLineEdit()
        self.lineI0 = QtGui.QLineEdit()
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(self.openFileButton, 0, 0)
        self.layout.addWidget(QtGui.QLabel("Scan #"), 1, 0)
        self.layout.addWidget(QtGui.QLabel("100pix data"), 2, 0)
        self.layout.addWidget(QtGui.QLabel("Energy"), 3, 0)
        self.layout.addWidget(QtGui.QLabel("I0"), 4, 0)
        self.layout.addWidget(self.lineEditFile, 0, 1)
        self.layout.addWidget(self.lineScan, 1, 1)
        self.layout.addWidget(self.line100pix, 2, 1)
        self.layout.addWidget(self.lineEnergy, 3, 1)
        self.layout.addWidget(self.lineI0, 4, 1)
        self.setLayout(self.layout)  
        self.setEnabled(True)
        self.setProperty("active", True)
        self.lineScan.editingFinished.connect(self.changeScanNumber)
        self.lineEditFile.editingFinished.connect(self.searchDataPaths)
        self.scansList = []
        self.lastScanIndex = 0
        self.lineScan.setSpecialValueText(" ");
        self.line100pix.setText("data/sis3302")
        self.lineEnergy.setText("data/energy_all")
        self.lineI0.setText("data/i0")
    
    def mousePressEvent(self,event):
        if self.isEnabled():
            self.setProperty("active", True)            
            self.setActiveSignal.emit()
                 
    def openFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '/home')
        self.lineEditFile.setText(fname) 
        self.searchDataPaths()
                    
    def changeScanNumber(self):
        i = self.lineScan.value()
        if self.scansList:
            if not (i in self.scansList):
                if (i < self.scansList[self.lastScanIndex]):
                    self.lastScanIndex = (self.lastScanIndex - 1) % len(self.scansList)
                    self.lineScan.setValue(self.scansList[self.lastScanIndex])
                    self.searchDataPaths()
                elif (i >self.scansList[self.lastScanIndex]):
                    self.lastScanIndex = (self.lastScanIndex + 1) % len(self.scansList)
                    self.lineScan.setValue(self.scansList[self.lastScanIndex])
                    self.searchDataPaths()
            else:
                self.lastScanIndex = self.scansList.index(i)
                self.searchDataPaths()
        else:
            self.lineScan.setValue(0)
            
    def searchDataPaths(self):
        fname = self.lineEditFile.text()
        if os.path.exists(fname):
            if fname[-4:] == '.nxs':
                f = h5py.File(fname, "r")
                self.scansList = []
                for i in f.keys():
                    if i[:4] == 'scan':
                        self.scansList.append(int(i[4:]))
                        
                if self.scansList:
                    self.lineScan.setValue(self.scansList[self.lastScanIndex%(len(self.scansList))])
                    self.data100pix = self.getDataByPath(f[i], self.line100pix.text())
                    self.dataEnergy = self.getDataByPath(f[i], self.lineEnergy.text())
                    self.dataI0 = self.getDataByPath(f[i], self.lineI0.text())
                    self.enableNextSignal.emit()
                    print fname, "loaded"
                
    def getDataByPath(self, item, path):        
        for j in path.split("/"):
            item = item[j]
        return item
        
        
class CalibrationWidget(coolWidget):
    def __init__(self):
        coolWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(QtGui.QPushButton("Aling Peaks"), 0, 0)
        self.layout.addWidget(QtGui.QLabel("Noise Level"), 0, 1)
        self.setLayout(self.layout) 
        
    def mousePressEvent(self,event):
        if self.isEnabled():
            self.setProperty("active", True)
            self.setActiveSignal.emit()

        
class RejectionWidget(coolWidget):
    def __init__(self):
        coolWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(QtGui.QPushButton("Reject Outliers"), 0, 0)
        self.layout.addWidget(QtGui.QLabel("?"), 0, 1)
        self.setLayout(self.layout)  

        
class ROIWidget(coolWidget):
    def __init__(self):
        coolWidget.__init__(self)
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(QtGui.QPushButton("Apply"), 0, 0)
        self.layout.addWidget(QtGui.QLabel("ROI1"), 0, 1)
        self.setLayout(self.layout)          

        
class XASWidget(coolWidget):
    def __init__(self):
        coolWidget.__init__(self) 
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(QtGui.QPushButton("Save XAS"), 0, 0)
        self.layout.addWidget(QtGui.QLabel("?"), 0, 1)
        self.setLayout(self.layout)  
        
        
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

        self.setGeometry(200, 50, 1200, 800)
        self.setWindowTitle('Data Reduction 100pixHPGe')
        self.plot = Plot2dWidget() 
        self.loadFileWidget = LoadFileWidget()
        self.calibrationWidget = CalibrationWidget()
        self.rejectionWidget = RejectionWidget()
        self.roiWidget = ROIWidget()
        self.xasWidget = XASWidget()                    
        self.layoutV = QtGui.QVBoxLayout()        
        self.layoutV.addWidget(self.loadFileWidget)
        self.layoutV.addWidget(self.calibrationWidget)
        self.layoutV.addWidget(self.rejectionWidget)
        self.layoutV.addWidget(self.roiWidget)
        self.layoutV.addWidget(self.xasWidget)
        self.layoutH = QtGui.QHBoxLayout()
        self.layoutH.addLayout(self.layoutV)
        self.layoutH.addWidget(self.plot)     
        self.setLayout(self.layoutH)
        self.activeTabIndex = 0
        self.tabsList = [self.loadFileWidget, self.calibrationWidget, self.rejectionWidget, self.roiWidget, self.xasWidget]
        for tab in self.tabsList:
            tab.setActiveSignal.connect(self.setActiveTab)
            tab.enableNextSignal.connect(self.enableNext)

        
    def openFile(self):
        image = np.random.randint(1024, size=(1000, 500))
        self.plot.setImage(image)
        
    def setActiveTab(self):
        newTabIndex = self.tabsList.index(self.sender())
        if newTabIndex != self.activeTabIndex:
            self.tabsList[self.activeTabIndex].setProperty("active", False)
            self.tabsList[self.activeTabIndex].updateStyle()
            self.tabsList[newTabIndex].setProperty("active", True)
            self.tabsList[newTabIndex].updateStyle()
            self.activeTabIndex = newTabIndex
    
    def enableNext(self):
        newTabIndex = self.tabsList.index(self.sender()) + 1
        if newTabIndex < len(self.tabsList):
            self.tabsList[newTabIndex].setEnabled(True)            
                
            
app = QtGui.QApplication(sys.argv)
main = MainWindow()
main.show()
sys.exit(app.exec_())

