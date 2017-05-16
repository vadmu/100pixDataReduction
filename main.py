#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import h5py
from functions import peaks_detection

class coolWidget(QtGui.QGroupBox):
    sigSetActive = QtCore.pyqtSignal()
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
    sigFileLoaded = QtCore.pyqtSignal()
    def __init__(self):
        coolWidget.__init__(self)
        self.openFileButton = QtGui.QPushButton("Open File")         
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

        self.scansList = []
        self.lastScanIndex = 0
        self.lineScan.setSpecialValueText(" ");
        self.line100pix.setText("data/sis3302")
        self.lineEnergy.setText("data/energy_all")
        self.lineI0.setText("data/i0")
        
        self.openFileButton.clicked.connect(self.openFile)
        self.lineScan.setKeyboardTracking(False)
        self.lineScan.valueChanged.connect(self.changeScanNumber)
#        self.lineEditFile.returnPressed.connect(self.searchDataPaths)
        self.lineEditFile.editingFinished.connect(self.searchDataPaths)
    
    def mousePressEvent(self,event):
        if self.isEnabled():
            self.setProperty("active", True)            
            self.sigSetActive.emit()
                 
    def openFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open File', '/home')
        self.lineEditFile.setText(fname) 
        self.searchDataPaths()
                    
    def changeScanNumber(self):
        i = self.lineScan.value()
        if self.scansList:
            if not (i in self.scansList):
                if (i < self.scansList[self.lastScanIndex]):
                    self.lastScanIndex = (self.lastScanIndex - 1) % len(self.scansList)
                    self.lineScan.blockSignals(True) # otherwise file loads twice...
                    self.lineScan.setValue(self.scansList[self.lastScanIndex])
                    self.lineScan.blockSignals(False)
                    self.searchDataPaths()
                elif (i >self.scansList[self.lastScanIndex]):
                    self.lastScanIndex = (self.lastScanIndex + 1) % len(self.scansList)
                    self.lineScan.blockSignals(True) # otherwise file loads twice...
                    self.lineScan.setValue(self.scansList[self.lastScanIndex])
                    self.lineScan.blockSignals(False)
                    self.searchDataPaths()
            else:
                self.lastScanIndex = self.scansList.index(i)
                self.searchDataPaths()
        else:
            self.lineScan.blockSignals(True) # otherwise file loads twice...
            self.lineScan.setValue(0)
            self.lineScan.blockSignals(False)
            
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
                    self.lineScan.blockSignals(True) # otherwise file loads twice...
                    self.lineScan.setValue(self.scansList[self.lastScanIndex%(len(self.scansList))])
                    self.lineScan.blockSignals(False)
                    self.data100pix = self.getDataByPath(f[scan], str(self.line100pix.text()))
                    self.dataEnergy = self.getDataByPath(f[scan], str(self.lineEnergy.text()))
                    self.dataI0 = self.getDataByPath(f[scan], str(self.lineI0.text()))
                    self.sigFileLoaded.emit()
                    print fname, "- scan %d"%self.lineScan.value(), "loaded"
                
    def getDataByPath(self, item, path):        
        for j in path.split("/"):
            item = item[j]
        return item
        
        
class CalibrationWidget(coolWidget):
    sigAlingPeaks = QtCore.pyqtSignal()
    sigAlingAll = QtCore.pyqtSignal()
    def __init__(self):
        coolWidget.__init__(self)
        self.minV = 0.01
        self.maxV = 50.
        self.step = 0.01
        self.value = 5.00
        
        self.openFileButton = QtGui.QPushButton("Open Calibration")         
        self.lineEditFile = QtGui.QLineEdit()
        self.slider = QtGui.QSlider()
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.slider.setMinimum(round(self.minV/self.step))
        self.slider.setMaximum(round(self.maxV/self.step))
        self.lineEdit = QtGui.QLineEdit(str(self.value))
        self.slider.setValue(self.value/self.step)
        self.lineEdit.setValidator(QtGui.QDoubleValidator())                
        self.peaksButton = QtGui.QPushButton("Aling Peaks")
        self.alignAllButton = QtGui.QPushButton("Aling All")
        self.peaksButton.clicked.connect(self.sigAlingPeaks)
        self.alignAllButton.clicked.connect(self.sigAlingAll)
        
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(self.openFileButton, 0, 0)
        self.layout.addWidget(self.lineEditFile, 0, 1)
        self.layout.addWidget(QtGui.QLabel("Noise Level"), 1, 0)
        self.layout.addWidget(self.lineEdit, 1, 1)
        self.layout.addWidget(self.slider, 2, 1)
        self.layout.addWidget(self.peaksButton, 2, 0)     
        self.layout.addWidget(self.alignAllButton, 3, 0)  
        self.setLayout(self.layout) 
        
        self.openFileButton.clicked.connect(self.openCalibrationFile)
        self.lineEdit.editingFinished.connect(self.lineEditValueChanged)
        self.slider.valueChanged.connect(self.sliderValueChanged)
        
    def mousePressEvent(self,event):
        if self.isEnabled():
            self.setProperty("active", True)
            self.sigSetActive.emit()
            
    def lineEditValueChanged(self):
        self.value = float(str(self.lineEdit.text()))
        self.slider.setValue(round(self.value/self.step))
        
    def sliderValueChanged(self, v):
        self.value = v*self.step
        self.lineEdit.setText(str(self.value))
        
    def openCalibrationFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open Calibration File', '/home')
        self.lineEditFile.setText(fname) 
        
        

        
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
        self.indexEnergyCurrent = 0
#        self.img = pg.ImageItem()
#        self.view.setAspectLocked(True)
        self.view.setCentralItem(self.vb)
        self.histWidget = pg.HistogramLUTWidget()
        self.histWidget.item.gradient.loadPreset('bipolar')
        self.energySpinBox = QtGui.QSpinBox()
        self.energySpinBox.setMinimum(0)
        self.energySpinBox.setMaximum(0)
        self.slider = QtGui.QSlider()
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)
        self.slider.setOrientation(QtCore.Qt.Horizontal)
        self.energyLabel = QtGui.QLabel("Energy:  %7.2f eV"%0)
        self.slider.valueChanged.connect(self.energySpinBox.setValue)
        self.energySpinBox.valueChanged.connect(self.slider.setValue)
        self.energySpinBox.valueChanged.connect(self.setEnergyLabelText)

        
#        self.view.addItem(self.img)
        self.layout = QtGui.QGridLayout() 
        self.layoutSelection = QtGui.QHBoxLayout() 
        self.layout.addWidget(self.view,0, 0)
        self.layout.addWidget(self.histWidget,0 ,1)
        self.layoutSelection.addWidget(self.energySpinBox)        
        self.layoutSelection.addWidget(self.energyLabel)
        self.layoutSelection.addWidget(self.slider)
        self.layout.addLayout(self.layoutSelection, 1 ,0, 1, 2)
        
        self.layout.setSpacing(0)
        self.setLayout(self.layout)
        self.showRaw = True
        self.peakind = []
    
    def setData(self, data, energy):
        self.data = np.array(data) #makes a copy
        self.dataAlinged = np.array(data)
        self.energy = energy
        self.img = pg.ImageItem(self.data[0].T, border='w')
        self.vb.addItem(self.img)
        self.vb.autoRange(padding = 0)
#        self.vb.autoRange()
        self.energySpinBox.setMinimum(1)
        self.energySpinBox.setMaximum(len(self.data))
        self.slider.setMinimum(1)
        self.slider.setMaximum(len(self.data))
        self.pixels = np.arange(len(self.data[0]))
#        self.energySpinBox.setValue(1)
        self.histWidget.setImageItem(self.img)
        self.histWidget.item.fillHistogram(False)
        self.minLevel = 0
        self.maxLevel = max([d.max() for d in data])
        self.histWidget.item.setHistogramRange(0, self.maxLevel)
        self.setImage(1)        
        self.histWidget.item.sigLevelChangeFinished.connect(self.changeLevels)
        self.energySpinBox.valueChanged.connect(self.setImage)

    
    def setImage(self, index):
        if self.showRaw:
            self.img.setImage(self.data[index - 1].T, autoLevels=False, levels = [self.minLevel, self.maxLevel])
            self.indexEnergyCurrent = index - 1
        else:
            self.img.setImage(self.dataAlinged[index - 1].T, autoLevels=False, levels = [self.minLevel, self.maxLevel])
            self.indexEnergyCurrent = index - 1
#        self.histWidget.setImageItem(self.img)
#        self.histWidget.item.imageChanged()

    def changeLevels(self):
        self.minLevel, self.maxLevel = self.histWidget.item.getLevels()
        self.img.setLevels([self.minLevel, self.maxLevel])
        
    def alingPeaks(self):
        snr = self.sender().value
        self.peakind = []
#        self.progress = QtGui.QProgressDialog("Alinging Peaks", "Cancel", 0, len(self.pixels), self)
#        self.progress.setWindowTitle(" ")
#        self.progress.setMinimumSize(400, 100)
#        self.progress.setWindowModality(QtCore.Qt.WindowModal)
        for pix in self.pixels:
            spectrum = self.data[self.indexEnergyCurrent][pix]
#            self.progress.setValue(pix)
#            if self.progress.wasCanceled(): break
            if sum(spectrum) > 0:
                ind, noise = peaks_detection(spectrum.astype(np.float), np.arange(1, 10), snr)
                if ind:
                    self.dataAlinged[self.indexEnergyCurrent][pix] = (np.roll(spectrum, 1500 - ind[-1]))
                    self.peakind.append(ind[-1])
                else:
                    self.peakind.append(0)
            else:
                self.peakind.append(0)
                
#        self.progress.setValue(len(self.pixels))
        self.showRaw = False
        self.setImage(self.indexEnergyCurrent + 1)

    def alingAll(self):
        if self.peakind:
#            self.progress = QtGui.QProgressDialog("Alinging All Channels", "Cancel", 0, len(self.energy), self)
#            self.progress.setWindowTitle(" ")
#            self.progress.setMinimumSize(400, 100)
#            self.progress.setWindowModality(QtCore.Qt.WindowModal)
            for ei in range(len(self.energy)):
#                self.progress.setValue(ei)
#                if self.progress.wasCanceled(): break
                for pix in self.pixels:
                    spectrum = self.data[ei][pix]    
                    self.dataAlinged[ei][pix] = (np.roll(spectrum, 1500 - self.peakind[pix]))
                    
#            self.progress.setValue(len(self.energy))
            self.showRaw = False
            self.setImage(self.indexEnergyCurrent + 1)    
        print "ok"
        
    def setEnergyLabelText(self, v):
        self.energyLabel.setText("Energy:  %7.2f eV"%self.energy[v-1])
        

        
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
        
        self.loadFileWidget.sigFileLoaded.connect(self.onFileLoad)
        self.calibrationWidget.sigAlingPeaks.connect(self.plotWidget.alingPeaks)
        self.calibrationWidget.sigAlingAll.connect(self.plotWidget.alingAll)
        
        self.plotWidget.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        for tab in self.tabsList:
            tab.setMinimumWidth(500)
            tab.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum))
            tab.sigSetActive.connect(self.setActiveTab)
        
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
        self.plotWidget.setData(self.loadFileWidget.data100pix, self.loadFileWidget.dataEnergy)       
                  
            
app = QtGui.QApplication(sys.argv)
main = MainWindow()
main.show()
sys.exit(app.exec_())

