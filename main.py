#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import h5py
from functions import CalibratePeaks

class coolWidget(QtGui.QGroupBox):
    setActiveSignal = QtCore.pyqtSignal()
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
    fileLoaded = QtCore.pyqtSignal()
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
        fname = str(self.lineEditFile.text())
        if os.path.exists(fname):
            if fname[-4:] == '.nxs':
                f = h5py.File(fname, "r")
                self.scansList = []
                for i in f.keys():
                    if i[:4] == 'scan':
                        self.scansList.append(int(i[4:]))
                        
                if self.scansList:
                    scan = [i for i in f.keys() if i[:4] == 'scan'][self.lastScanIndex%(len(self.scansList))]
                    self.lineScan.setValue(self.scansList[self.lastScanIndex%(len(self.scansList))])
                    self.data100pix = self.getDataByPath(f[scan], str(self.line100pix.text()))
                    self.dataEnergy = self.getDataByPath(f[scan], str(self.lineEnergy.text()))
                    self.dataI0 = self.getDataByPath(f[scan], str(self.lineI0.text()))
                    self.fileLoaded.emit()
                    print fname, "- scan %d"%self.lineScan.value(), "loaded"
                
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
        
        
class Plot2dWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        pg.setConfigOption('background', 'w')
        self.view = pg.GraphicsView()
        self.vb = pg.ViewBox(border='w')
#        self.img = pg.ImageItem()
#        self.view.setAspectLocked(True)
        self.view.setCentralItem(self.vb)
        self.histWidget = pg.HistogramLUTWidget()
        self.histWidget.item.gradient.loadPreset('bipolar')
        self.energySpinBox = QtGui.QSpinBox()
        self.energySpinBox.setMinimum(1)
        self.energySpinBox.setMaximum(1)
        self.energySpinBox.valueChanged.connect(self.setImage)
        
#        self.view.addItem(self.img)
        self.layout = QtGui.QGridLayout() 
        self.layout.addWidget(self.view,0, 0)
        self.layout.addWidget(self.histWidget,0 ,1)
        self.layout.addWidget(self.energySpinBox,1 ,1)
        self.layout.setSpacing(0)
        self.setLayout(self.layout)
    
    def setData(self, data):
        self.data = data
        self.img = pg.ImageItem(self.data[0].T, border='w')
        self.vb.addItem(self.img)
#        self.vb.autoRange()
        self.energySpinBox.setMaximum(len(self.data))
        self.histWidget.setImageItem(self.img)
        self.histWidget.item.fillHistogram(False)
        self.minLevel = 0
        self.maxLevel = max([d.max() for d in data])
        self.histWidget.item.setHistogramRange(0, self.maxLevel)
        self.setImage(1)        
        self.histWidget.item.sigLevelChangeFinished.connect(self.changeLevels)
    
    def setImage(self, index):
        self.img.setImage(self.data[index - 1].T, autoLevels=False, levels = [self.minLevel, self.maxLevel])
#        self.histWidget.setImageItem(self.img)
#        self.histWidget.item.imageChanged()

    def changeLevels(self):
        self.minLevel, self.maxLevel = self.histWidget.item.getLevels()
        self.img.setLevels([self.minLevel, self.maxLevel])
        

        
class MainWindow(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.setGeometry(200, 50, 1200, 800)
        self.setWindowTitle('Data Reduction 100pixHPGe')
        self.plotWidget = Plot2dWidget() 
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
        self.layoutH.addWidget(self.plotWidget)     
        self.setLayout(self.layoutH)
        self.activeTabIndex = 0
        self.tabsList = [self.loadFileWidget, self.calibrationWidget, self.rejectionWidget, self.roiWidget, self.xasWidget]
        
        self.loadFileWidget.fileLoaded.connect(self.onFileLoad)
        
        self.plotWidget.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        for tab in self.tabsList:
            tab.setMinimumWidth(500)
            tab.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum))
            tab.setActiveSignal.connect(self.setActiveTab)
        
    def setActiveTab(self):
        newTabIndex = self.tabsList.index(self.sender())
        if newTabIndex != self.activeTabIndex:
            self.tabsList[self.activeTabIndex].setProperty("active", False)
            self.tabsList[self.activeTabIndex].updateStyle()
            self.tabsList[newTabIndex].setProperty("active", True)
            self.tabsList[newTabIndex].updateStyle()
            self.activeTabIndex = newTabIndex
            
    def onFileLoad(self):
        self.calibrationWidget.setEnabled(True)
        self.rejectionWidget.setEnabled(True)
        self.roiWidget.setEnabled(True)
        self.plotWidget.setData(self.loadFileWidget.data100pix)       
                  
            
app = QtGui.QApplication(sys.argv)
main = MainWindow()
main.show()
sys.exit(app.exec_())

