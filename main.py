#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import h5py
from functions import findPeaks, relatePix, align


########################################
#  Logger
########################################
import logging

class QtHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
    def emit(self, record):
        record = self.format(record)
        if record: XStream.stdout().write('%s\n'%record)

class XStream(QtCore.QObject):
    _stdout = None
    _stderr = None
    messageWritten = QtCore.pyqtSignal(str)
    def flush( self ):
        pass
    def fileno( self ):
        return -1
    def write( self, msg ):
        if ( not self.signalsBlocked() ):
            self.messageWritten.emit(unicode(msg))
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
        return XStream._stderr        
        
logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(levelname)-8s [%(asctime)-15s] %(message)s')
formatter.datefmt = '%Y-%m-%d %H:%M:%S'
logging.basicConfig(format ='%(levelname)-8s [%(asctime)-15s] %(message)s', 
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level = logging.INFO, 
                    filename = 'info.log')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(formatter)
logger.addHandler(console)

handlerStatus = QtHandler()
handlerStatus.setLevel(logging.INFO)
handlerStatus.setFormatter(formatter)
logger.addHandler(handlerStatus)


########################################
#  Data operations 
########################################
class Data(QtCore.QObject):
    updateState = QtCore.pyqtSignal()
    fileLoaded = QtCore.pyqtSignal()
#    energyChanged = QtCore.pyqtSignal(int)
    def __init__(self):
        QtCore.QObject.__init__(self)
        self.arrayRaw = []
        self.arrayAlinged = []
        self.arrayCalibrated = []
        self.arrayEnergy = []
        self.arrayI0 = []  
        self.xas = []
        self.roi = [0, 1]
        self.scanList = []
        self.dimPixels = 0
        self.dimChannels = 0
        self.energyIndex = 0
        self.pixRej = [False]*100
        self.pixOff = [False]*100
        self.rej1 = 0
        self.rej2 = 0
        self.ae = 1.
        self.be = 0.
        self.ch2ea = np.zeros(100)
        self.ch2eb = np.zeros(100)
        self.a=np.zeros(100)
        self.b=np.zeros(100)
        self.specEnergy = []
        self.specSum = []
    
    def alingPeaks(self, ct=700, pn=2, md=50, th=0.3):
        #cut - cuting point in channels (offset in the spectra)
        #peak_nr - number of peaks
        #md - minimum distance between peaks
        #th - threshold of noise, 1:max, 0:min
        index=0        
        for pix in range(self.dimPixels):         
            suma=sum(self.arrayRaw[self.energyIndex][pix])
            index=index+1
            if suma>0:
                break
        # Pixel with non-zero counts as reference spectrum:
        spectrum_ref=self.arrayRaw[self.energyIndex][index-1]
        # finds peaks in the reference spectrum:
        IR=findPeaks(spectrum_ref, pn, md, th, ct)
# TODO: remove dependence on N = 100     	
        self.a=np.zeros(100)
        self.b=np.zeros(100)
        # Finds peaks in each spectrum and relates the peak positions to the reference spectrum:
        for pix in range(self.dimPixels):
            spectrum = self.arrayRaw[self.energyIndex][pix]
            suma=sum(spectrum)
            if suma>0:
                I=findPeaks(spectrum, pn, md, th, ct)
#                logger.info(str(pix) + " " + str(suma) + " " + str(I))
                if IR.shape==I.shape:
                    self.a[pix], self.b[pix]=relatePix(I,IR)
            else:
                self.pixOff[pix] = True
            # Alignes the spectra:
            self.arrayAlinged[self.energyIndex][pix]=align(pix, spectrum, self.dimChannels, self.a, self.b) 
        logger.info("Data at %.2f eV aligned with parameters: %d %d %d %.2f"%(self.arrayEnergy[self.energyIndex], ct, pn, md, th)) 
        self.updateState.emit()

                
    def alignAll(self):
        for ei in range(len(self.arrayEnergy)):
            for pix in range(self.dimPixels):
                spectrum = self.arrayRaw[ei][pix]
                self.arrayAlinged[ei][pix]=align(pix, spectrum, self.dimChannels, self.a, self.b)
        logger.info("Data aligned") 
        self.calculateSum()
        self.updateState.emit()        
    
    def getDataByPath(self, item, path):        
            for j in path.split("/"):
                item = item[j]
            return item    
            
    def loadDataFileNXS(self, fname, scanIndex = 0, path = ["data/sis3302", "data/energy_all", "data/i0"]):
        self.scanList = []
        if os.path.exists(fname):
            # supposed that *.nxs file has "scan0000" format in the root  
            try:
                if fname[-4:] == '.nxs':
                    f = h5py.File(fname, "r")
                    self.scanList = []
                    for i in f.keys():
                        if i[:4] == 'scan':
                            self.scanList.append(i)
                    if self.scanList:
                        scan = self.scanList[scanIndex]
                        self.arrayRaw = self.getDataByPath(f[scan], path[0])
                        self.arrayAlinged = np.array(self.arrayRaw) # a copy
                        self.arrayEnergy = self.getDataByPath(f[scan], path[1])
                        self.arrayI0 = self.getDataByPath(f[scan], path[2])
                        self.dimPixels = len(self.arrayRaw[0])
                        self.dimChannels = len(self.arrayRaw[0][0])
                        logger.info("Data file loaded: %s - %s"%(fname, scan))
                        self.xas = np.zeros_like(self.arrayEnergy)
            except:
                logger.error("Cannot load the file: %s - %s"%(fname, scan))
        else:
            logger.warning("No such file: %s"%fname) 
        self.fileLoaded.emit()
    
    def loadDataFileDAT(self, fname):
        pass    
    
    def loadCalibrationFile(self, fname):
        self.ch2ea, self.ch2eb = np.loadtxt(fname, unpack = True) 
        found = False
        for i in range(len(self.ch2ea)):         
            if self.ch2ea[i]!=0 and self.ch2eb[i]!=0:
                if not found:
                    self.ae = self.ch2ea[i]
                    self.be = self.ch2eb[i]
                    found = True
            else:
                self.pixOff[i] = True
                
        self.a = self.ch2ea/self.ae * [not i for i in self.pixOff]
        self.b = (self.ch2eb - self.be)/self.ae * [not i for i in self.pixOff]     
        self.alignAll()
        logger.info("Calibration file loaded: %s"%fname)
#        self.updateState.emit() 
    
    def makeCalibrationFile(self, fname):
        np.savetxt(fname, np.column_stack((self.ch2ea, self.ch2eb)), fmt = '%10.3f')
        logger.info("Calibration file saved: %s"%fname)
        self.updateState.emit() 

    def calibrate(self, ch1, e1, ch2, e2):
        try:
            self.ae = (e2-e1)/(ch2-ch1)
        except:
            self.ae = 1.            
        self.be = e1 - ch1*self.ae    
        self.ch2ea = self.a * self.ae * [not i for i in self.pixOff]
        self.ch2eb = (self.ae*self.b + self.be) * [not i for i in self.pixOff]
        logger.info("Data calibrated: a=%.3f  b=%.3f"%(self.ae , self.be))
        self.updateState.emit() 
            
    def reject(self, start, end):
#        r1=self.sender().r1 #initial point in the region
#        r2=self.sender().r2 #final point in the region
#        sigma=self.sender().sigma #sigma of shot noise
#        max_counts=self.sender().max_counts #maximum number of counts as a threshold of good pixel
#        ei=self.sender().ei #incident energy
#
#        # Non zero pixels:
#        suma=sum(self.dataAlinged[ei][:,], axis=1)
#        nonzeroDataAligned=self.dataAlinged[ei][suma>max_counts,]
#        nonzero_pix=nonzeroDataAligned.shape[0]
#
#        # Normalized spectra:
#        max_pix=self.dataAlinged[ei][:,].max(1) #maximum at eaxh pixel
#        data_norm=(self.dataAlinged[ei][:].T/max_pix).T # normalized spectra (n,m)
#        data_norm[np.isnan(data_norm)] = 0 #gets rid of Nan due to divisions by zero (pixels with zero counts)
#
#        # Reject outliers:
#        mean_spectrum=np.mean(data_norm, axis=0)
#        shot_noise=sqrt(mean_spectrum) 
#
##        goodpix=np.ones(100, dtype=np.int8) # initialize pixels as good pixels (ones)
#
#        for pix in self.pixels:
#            if sum(self.dataAlinged[ei][pix,])>max_counts:
#                diff=abs(data_norm[pix,]-mean_spectrum) #difference between spectrum and mean spectrum
##                f=plt.subplot(10,10,pix)
##                plt.plot(mean_spectrum[r1:r2], 'r')
##                plt.plot(data_norm[pix,r1:r2], 'b')
##                plt.plot(diff[r1:r2], 'g')
##                plt.plot(sigma*shot_noise[r1:r2], 'm')
##                plt.xticks([], [])
##                plt.yticks([], [])
##                plt.title("pix: "+str(pix))
#                ## If the difference is bigger than the shot noise pixel is turned off (zero)
#                if sum(diff[r1:r2]>sigma*shot_noise[r1:r2])>0: 
##                    f.patch.set_facecolor('gray')
##                    goodpix[pix]=0.
#            else:
#                diff=np.zeros(2048)
#                goodpix[pix]=0.
        
        
        self.updateState.emit() 
        
    def calculateSum(self):
        self.specSum = np.zeros((len(self.arrayEnergy), self.dimChannels))
        self.specEnergy = np.arange(self.dimChannels) * self.ae + self.be
        self.roi = [self.specEnergy[-1]*0.5, self.specEnergy[-1]*0.7]
        for ei in range(len(self.arrayEnergy)):
            for pix in range(self.dimPixels):
                self.specSum[ei] += self.arrayAlinged[ei][pix]
#        self.updateState.emit() 
    def setROI(self, roi1, roi2):
        self.roi = [roi1, roi2]
        self.updateState.emit() 

    def calculateXAS(self):
        i1 = int((self.roi[0] - self.be) /self.ae)
        i2 = int((self.roi[1] - self.be) /self.ae) + 1
        self.xas = np.zeros_like(self.arrayEnergy)
        for ei in range(len(self.arrayEnergy)):
            for pix in range(self.dimPixels):
                self.xas[ei] += float(sum(self.arrayAlinged[ei][pix][i1:i2])) / self.arrayI0[ei]
        self.updateState.emit() 
        logger.info("XAS calculated: %.2f "%(sum(self.xas)))
        
    def saveXASFile(self, fname):
        np.savetxt(fname, np.column_stack((self.arrayEnergy, self.xas)))
        logger.info("XAS file saved: %s"%fname)
        self.updateState.emit() 
         
         
# a global instance
data = Data()
# data.loadDataFileNXS("test_Y_long.nxs", 3)
# data.alingPeaks()

    
########################################
#  Settings Widgets
########################################                
        
class LoadFileWidget(QtGui.QTabWidget):
    def __init__(self):
        QtGui.QTabWidget.__init__(self)
        self.nxsWidget = QtGui.QWidget()
        self.datWidget = QtGui.QWidget()
        self.openFileButton1 = QtGui.QPushButton("Open File")     
        self.openFileButton2 = QtGui.QPushButton("Open File") 
        self.lineEditFile1 = QtGui.QLineEdit()
        self.lineEditFile2 = QtGui.QLineEdit()
        self.lineScan = QtGui.QComboBox()
        self.line100pix = QtGui.QLineEdit()
        self.lineEnergy = QtGui.QLineEdit()
        self.lineI0 = QtGui.QLineEdit()
        self.layout1 = QtGui.QGridLayout()
        self.layout1.addWidget(self.openFileButton1, 0, 0)
        self.layout1.addWidget(QtGui.QLabel("Scan #"), 1, 0)
        self.layout1.addWidget(QtGui.QLabel("100pix data"), 2, 0)
        self.layout1.addWidget(QtGui.QLabel("Energy"), 3, 0)
        self.layout1.addWidget(QtGui.QLabel("I0"), 4, 0)
        self.layout1.addWidget(self.lineEditFile1, 0, 1)
        self.layout1.addWidget(self.lineScan, 1, 1)
        self.layout1.addWidget(self.line100pix, 2, 1)
        self.layout1.addWidget(self.lineEnergy, 3, 1)
        self.layout1.addWidget(self.lineI0, 4, 1)
        self.layout2 = QtGui.QGridLayout()
        self.layout2.addWidget(self.openFileButton2, 0, 0)
        self.layout2.addWidget(self.lineEditFile2, 0, 1)
        
        self.nxsWidget.setLayout(self.layout1)  
        self.datWidget.setLayout(self.layout2)
        self.addTab(self.nxsWidget, "NXS File")
        self.addTab(self.datWidget, "ASCII File")
        
        self.lastScanIndex = 0
        self.line100pix.setText("data/sis3302")
        self.lineEnergy.setText("data/energy_all")
        self.lineI0.setText("data/i0")
        
        self.lineScan.currentIndexChanged.connect(self.changeScanNumber)
        self.openFileButton1.clicked.connect(self.openFileDialog1)
        self.openFileButton2.clicked.connect(self.openFileDialog2)
        self.lineEditFile1.editingFinished.connect(self.openFile1)
        self.lineEditFile2.editingFinished.connect(self.openFile2)
        # self.line100pix.editingFinished.connect(self.openFile)
        # self.lineEnergy.editingFinished.connect(self.openFile)
        # self.lineI0.editingFinished.connect(self.openFile)
        
        

    def openFileDialog1(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        self.lineEditFile1.setText(fname) 
        self.openFile1()
    
    def openFileDialog2(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open File')
        self.lineEditFile2.setText(fname) 
        self.openFile2()
        
    def openFile1(self, i = 0):
        self.lineScan.currentIndexChanged.disconnect(self.changeScanNumber)
        self.lineScan.clear() 
        fname = str(self.lineEditFile1.text())
        path0 = str(self.line100pix.text())
        path1 = str(self.lineEnergy.text())
        path2 = str(self.lineI0.text())
        data.loadDataFileNXS(str(fname), i, [path0, path1, path2])           
        for scan in data.scanList:
            self.lineScan.addItem(scan)
        self.lineScan.currentIndexChanged.connect(self.changeScanNumber)
    
    def openFile2(self):
        pass
        
    def changeScanNumber(self, i):
        fname = str(self.lineEditFile1.text())
        path0 = str(self.line100pix.text())
        path1 = str(self.lineEnergy.text())
        path2 = str(self.lineI0.text())
        data.loadDataFileNXS(str(fname), i, [path0, path1, path2])    
        

class CalibrationWidget(QtGui.QTabWidget):
    def __init__(self):
        QtGui.QTabWidget.__init__(self)
        self.alWidget = QtGui.QWidget()
        self.calWidget = QtGui.QWidget()
        self.step = 0.01
        self.lineEditCT = QtGui.QLineEdit("700")        
        self.lineEditMD = QtGui.QLineEdit("50")
        self.lineEditPN = QtGui.QLineEdit("2")
        self.lineEditPN.setDisabled(True)
        self.lineEditTH = QtGui.QLineEdit("0.3")               
        self.alignCurButton = QtGui.QPushButton("Aling Current")
        self.alignAllButton = QtGui.QPushButton("Aling All")
        self.alignCurButton.clicked.connect(self.alingCur)
        self.alignAllButton.clicked.connect(self.alingAll)
        
        self.lineEditPos1 = QtGui.QLineEdit("500")
        self.lineEditPos2 = QtGui.QLineEdit("1500")
        self.lineEditEnergy1 = QtGui.QLineEdit("5000")
        self.lineEditEnergy2 = QtGui.QLineEdit("15000")
        self.openFileButton = QtGui.QPushButton("Open Calibration")         
        self.lineEditFile = QtGui.QLineEdit()
        self.calibrateButton = QtGui.QPushButton("Calibrate") 
        self.saveCalibrationButton = QtGui.QPushButton("Save Calibration") 
        
        self.layout1 = QtGui.QGridLayout()
        self.layout1.addWidget(QtGui.QLabel("Offset:"), 0, 0)
        self.layout1.addWidget(self.lineEditCT, 0, 1)
        self.layout1.addWidget(QtGui.QLabel("Noise Threshold:"), 1, 0)	
        self.layout1.addWidget(self.lineEditTH, 1, 1)
        self.layout1.addWidget(QtGui.QLabel("Minimum Distance:"), 2, 0)
        self.layout1.addWidget(self.lineEditMD, 2, 1)
        self.layout1.addWidget(QtGui.QLabel("# of Peaks:"), 3, 0)
        self.layout1.addWidget(self.lineEditPN, 3, 1)
        self.layout1.addWidget(self.alignCurButton, 4, 1)     
        self.layout1.addWidget(self.alignAllButton, 4, 0) 	
        self.layout2 = QtGui.QGridLayout()
        self.layout2.addWidget(self.openFileButton, 0, 0)
        self.layout2.addWidget(self.lineEditFile, 0, 1, 1, 3)
        self.layout2.addWidget(QtGui.QLabel("Line 1 pos:"), 1, 0)
        self.layout2.addWidget(QtGui.QLabel("Line 2 pos:"), 2, 0)
        self.layout2.addWidget(self.lineEditPos1, 1, 1)
        self.layout2.addWidget(self.lineEditPos2, 2, 1)
        self.layout2.addWidget(QtGui.QLabel("e1:"), 1, 2)
        self.layout2.addWidget(QtGui.QLabel("e2:"), 2, 2)
        self.layout2.addWidget(self.lineEditEnergy1, 1, 3)
        self.layout2.addWidget(self.lineEditEnergy2, 2, 3)
        self.layout2.addWidget(self.calibrateButton, 3, 3)
        self.layout2.addWidget(self.saveCalibrationButton, 3, 0)
        
        self.alWidget.setLayout(self.layout1) 
        self.calWidget.setLayout(self.layout2)                 
        self.addTab(self.alWidget, "Alignment")
        self.addTab(self.calWidget, "Calibration")
        
        self.openFileButton.clicked.connect(self.openCalibrationFile)
        self.saveCalibrationButton.clicked.connect(self.saveCalibrationFile)
        self.calibrateButton.clicked.connect(self.calibrate)

        
    def openCalibrationFile(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open Calibration File')        
        data.loadCalibrationFile(str(fname))
        self.lineEditFile.setText(fname) 
        
    def saveCalibrationFile(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Calibration File', "calibration.dat")
        data.makeCalibrationFile(str(fname))
        self.lineEditFile.setText(fname)
        
    def alingCur(self):        
        data.alingPeaks(ct=int(self.lineEditCT.text()),                      
                        pn=int(self.lineEditPN.text()), 
                        md=int(self.lineEditMD.text()), 
                        th=float(self.lineEditTH.text()))
                        
    def alingAll(self): 
        data.alignAll()
        
    def calibrate(self):
        data.calibrate(int(self.lineEditPos1.text()),
                       float(self.lineEditEnergy1.text()),
                       int(self.lineEditPos2.text()),
                       float(self.lineEditEnergy2.text()))
        
    def setValuePos1(self, i):
        self.setCurrentIndex(1)
        self.lineEditPos1.setText(str(i))
        self.lineEditEnergy1.setText(str(int(i*data.ae + data.be)))
        
    def setValuePos2(self, i):
        self.setCurrentIndex(1)
        self.lineEditPos2.setText(str(i))
        self.lineEditEnergy2.setText(str(int(i*data.ae + data.be)))
        
        
class RejectionWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self) 
        self.rejectButton = QtGui.QPushButton("Reject Outliers")  
        self.rejectButton.clicked.connect(self.reject)
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(self.rejectButton, 0, 1) 
        self.setLayout(self.layout) 
    
    def reject(self):
        data.reject(500 , 1500)
    
                     
class ROIWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.roiB = QtGui.QLineEdit("400") 
        self.roiE = QtGui.QLineEdit("600")        
        self.btnApply = QtGui.QPushButton("Apply")           
        self.layout = QtGui.QGridLayout()
        self.layout.addWidget(self.btnApply, 0, 3)
        self.layout.addWidget(QtGui.QLabel("ROI"), 0, 0)
        self.layout.addWidget(self.roiB, 0, 1)
        self.layout.addWidget(self.roiE, 0, 2)
        self.setLayout(self.layout) 
        self.roiB.editingFinished.connect(self.setROI)
        self.roiE.editingFinished.connect(self.setROI)
        self.btnApply.clicked.connect(data.calculateXAS)
    
    def updateROItext(self):
        self.roiB.setText("%7.2f"%data.roi[0])
        self.roiE.setText("%7.2f"%data.roi[1])
        
    def setROI(self):
        roi1 = float(self.roiB.text())
        roi2 = float(self.roiE.text())
        data.setROI(roi1, roi2)

        
class XASWidget(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self) 
        self.layout = QtGui.QGridLayout()
        self.btnSave = QtGui.QPushButton("Save XAS")
        self.layout.addWidget(self.btnSave, 0, 0)
        self.setLayout(self.layout)  
        self.btnSave.clicked.connect(self.saveXASFile)
        
    def saveXASFile(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Calibration File', "xas.dat")
        data.saveXASFile(str(fname))
#        self.lineEditFile.setText(fname)

        
########################################
#  Plot Widgets
########################################          

class PlotMap(QtGui.QWidget):
    pos1changed = QtCore.pyqtSignal(int)
    pos2changed = QtCore.pyqtSignal(int)
    def __init__(self):
        QtGui.QWidget.__init__(self)
        pg.setConfigOption('background', 'w')
        self.view = pg.GraphicsView()
        self.vb = pg.ViewBox(border='w')
        self.vb.setLimits(minXRange = 1, maxXRange = 2048, xMin = 0, xMax = 2048,
                          minYRange = 1, maxYRange = 100, yMin = 0, yMax = 100)
        self.img = None
#        self.indexEnergyCurrent = 0
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
#        self.energySpinBox.valueChanged.connect(data.setEnergyIndex)

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
        
        data.fileLoaded.connect(self.onFileLoad)
        data.updateState.connect(self.updateImage)
#        data.energyChanged.connect(self.energySpinBox.setValue)
        
        self.isoLine1 = pg.InfiniteLine(angle=90, movable=True, bounds = [1, 2046], pen=pg.mkPen(color=(150, 100, 150), width=2))
        self.isoLine2 = pg.InfiniteLine(angle=90, movable=True, bounds = [1, 2046], pen=pg.mkPen(color=(150, 200, 0), width=2))
        self.vb.addItem(self.isoLine1)
        self.vb.addItem(self.isoLine2)
        self.isoLine1.setValue(500)
        self.isoLine2.setValue(1500)
        self.isoLine1.setZValue(1001) # bring iso line above
        self.isoLine2.setZValue(1000) # bring iso line above
        self.isoLine1.sigPositionChanged.connect(self.onLine1PosChanged)
        self.isoLine2.sigPositionChanged.connect(self.onLine2PosChanged)
   
    def onLine1PosChanged(self):
#        logger.warning(self.isoLine1.pos())
        self.pos1changed.emit(self.isoLine1.pos().x())
        
    def onLine2PosChanged(self):
        self.pos2changed.emit(self.isoLine2.pos().x())
        
    def setEnergyLabelText(self, v):
        self.energyLabel.setText("Energy:  %7.2f eV"%data.arrayEnergy[v-1])
        
    def onFileLoad(self):
        if self.img: 
            self.vb.removeItem(self.img)
        self.img = pg.ImageItem(data.arrayAlinged[0].T, border='w')
        self.minLevel = 0
        self.maxLevel = max([d.max() for d in data.arrayAlinged])
        self.vb.addItem(self.img)
        self.vb.autoRange(padding = 0)
        self.energySpinBox.setMinimum(1)
        self.energySpinBox.setMaximum(len(data.arrayAlinged))
        self.energySpinBox.setValue(data.energyIndex +1)
        self.slider.setMinimum(1)
        self.slider.setMaximum(len(data.arrayAlinged))
        self.pixels = np.arange(data.dimPixels)
        self.histWidget.setImageItem(self.img)
        self.histWidget.item.fillHistogram(False)
        self.histWidget.item.setHistogramRange(0, self.maxLevel)
        self.histWidget.item.sigLevelChangeFinished.connect(self.changeLevels)
        self.energySpinBox.valueChanged.connect(self.updateImage)
       
        
    def updateImage(self):
        data.energyIndex = self.energySpinBox.value() - 1
        self.img.setImage(data.arrayAlinged[data.energyIndex].T, autoLevels=False, levels = [self.minLevel, self.maxLevel])

    def changeLevels(self):
        self.minLevel, self.maxLevel = self.histWidget.item.getLevels()
        self.img.setLevels([self.minLevel, self.maxLevel])    
        
class PlotCalibration(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.win = pg.TableWidget(editable = True, sortable=False)
        self.layout = QtGui.QVBoxLayout() 
        self.layout.addWidget(self.win)
        self.setLayout(self.layout)
        self.win.setData([])
        data.updateState.connect(self.replot)
        self.win.setFormat("%.3f")
        
    def replot(self):  
        self.win.setData(np.vstack((data.ch2ea, data.ch2eb)).T)

        
class PlotRejection(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self) 
        self.win = pg.GraphicsLayoutWidget()
        self.plots = []
# TODO: remove dependence on N = 100        
        for i in range(5):
            for j in range(10):
                pl = self.win.addPlot(enableMenu = False)
                pl.showAxis("left", show=False)
                pl.showAxis("bottom", show=False)
                pl.setMouseEnabled(x=False, y=False)
                pl.vb.setBackgroundColor((0, 0, 0, 255))
                x = np.cos(np.linspace(0, 2*np.pi, 1000))
                y = np.sin(np.linspace(0, 4*np.pi, 1000))
                pl.plot(x, y)
                self.plots.append(pl)                
            self.win.nextRow() 
        
        self.layout = QtGui.QVBoxLayout() 
        self.layout.addWidget(self.win)
        self.setLayout(self.layout)
        data.updateState.connect(self.replot)
        
    def mousePressEvent(self, event):
#        logger.warning("clicked %d %d"%(event.globalX(), event.globalY()))
#        logger.warning(self.plots[0].viewGeometry())
        for i, pl in enumerate(self.plots):
            rect = pl.viewGeometry()
            if (rect.x() < event.globalX() < rect.x() + rect.width() ) and \
                (rect.y() < event.globalY() < rect.y() + rect.height()):
                data.pixRej[i] =  not data.pixRej[i]  
                self.replot()
                
    def replot(self):   
        for i, pl in enumerate(self.plots):         
            if data.pixOff[i]:
                pl.vb.setBackgroundColor((0, 0, 0, 255))
            elif data.pixRej[i]:
                pl.vb.setBackgroundColor((120, 0, 0, 255))
            else:
                pl.vb.setBackgroundColor((0, 120, 0, 255))
                        
class PlotROI(QtGui.QWidget):
    roiChanged = QtCore.pyqtSignal()
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.pl = pg.PlotWidget()
        self.roi = pg.LinearRegionItem([-0.5, 0.5], bounds=[-1, 1], movable=True)
        self.pl.getViewBox().setBackgroundColor((0, 0, 0, 255))
        x = np.cos(np.linspace(0, 2*np.pi, 1000))
        y = np.sin(np.linspace(0, 4*np.pi, 1000))
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

               
        self.layout = QtGui.QVBoxLayout() 
        self.layoutSelection = QtGui.QHBoxLayout() 
        self.layout.addWidget(self.pl)
        self.layoutSelection.addWidget(self.energySpinBox)        
        self.layoutSelection.addWidget(self.energyLabel)
        self.layoutSelection.addWidget(self.slider)
        self.layout.addLayout(self.layoutSelection)        
        self.layout.setSpacing(0)
        self.setLayout(self.layout)
        self.pl.plot(x, y)
        self.pl.addItem(self.roi )
        
        data.fileLoaded.connect(self.onFileLoad)
        data.updateState.connect(self.updateImage)
        
    def setEnergyLabelText(self, v):
        self.energyLabel.setText("Energy:  %7.2f eV"%data.arrayEnergy[v-1])
        
    def onFileLoad(self):
#        if self.img: 
#            self.vb.removeItem(self.img)
#        self.img = pg.ImageItem(data.arrayAlinged[0].T, border='w')
#        self.vb.addItem(self.img)
#        self.vb.autoRange(padding = 0)
#        self.minLevel = 0
#        self.maxLevel = max([d.max() for d in data.arrayAlinged])
        self.energySpinBox.setMinimum(1)
        self.energySpinBox.setMaximum(len(data.arrayAlinged))
        self.slider.setMinimum(1)
        self.slider.setMaximum(len(data.arrayAlinged))
#        self.pixels = np.arange(data.dimPixels)
#        self.histWidget.setImageItem(self.img)
#        self.histWidget.item.fillHistogram(False)

#        self.histWidget.item.setHistogramRange(0, self.maxLevel)
#        self.histWidget.item.sigLevelChangeFinished.connect(self.changeLevels)
        self.energySpinBox.valueChanged.connect(self.updateImage)
        self.roi.sigRegionChanged.connect(self.onChangeROI)
        
    def onChangeROI(self):
        data.roi = self.roi.getRegion()
        self.roiChanged.emit()
        
        
    def updateImage(self):

        data.energyIndex = self.energySpinBox.value() - 1
        if len(data.specEnergy):
            self.pl.clear()
            self.pl.plot(data.specEnergy, data.specSum[data.energyIndex])
            self.roi.setBounds([data.specEnergy[0], data.specEnergy[-1]])
            self.roi.setRegion(data.roi)
            self.pl.addItem(self.roi )
                
class PlotXAS(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.pl = pg.PlotWidget()
        self.pl.getViewBox().setBackgroundColor((50, 0, 50, 255))
        x = np.cos(np.linspace(0, 2*np.pi, 1000))
        y = np.sin(np.linspace(0, 4*np.pi, 1000))
        self.pl.plot(x, y)
        self.layout = QtGui.QVBoxLayout() 
        self.layout.addWidget(self.pl)
        self.setLayout(self.layout)
        data.updateState.connect(self.replot)
        
    def replot(self):        
        self.pl.clear()
        self.pl.plot(data.arrayEnergy, data.xas)
   
            
class PlotTabWidget(QtGui.QTabWidget):
    def __init__(self):
        QtGui.QTabWidget.__init__(self)
        self.rawWidget = PlotMap()
        self.calWidget = PlotCalibration()
        self.rejWidget = PlotRejection()
        self.roiWidget = PlotROI()
        self.xasWidget = PlotXAS()
        
        self.addTab(self.rawWidget, "Data Raw/Aligned")
        self.addTab(self.calWidget, "Calibration")
        self.addTab(self.rejWidget, "Rejection")
        self.addTab(self.roiWidget, "ROI")
        self.addTab(self.xasWidget, "XAS Spectrum")
        

        

        
########################################
#  Main Window
########################################         
class MainWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setGeometry(50, 50, 1000, 800)
        self.setWindowTitle('Data Reduction 100pixHPGe')
        self.loadFileWidget = LoadFileWidget()
        self.calibrationWidget = CalibrationWidget()
        self.rejectionWidget = RejectionWidget()
        self.roiWidget = ROIWidget()
        self.xasWidget = XASWidget()          
        self.plotWidget = PlotTabWidget()   

        self.plotWidget.rawWidget.pos1changed.connect(self.calibrationWidget.setValuePos1)
        self.plotWidget.rawWidget.pos2changed.connect(self.calibrationWidget.setValuePos2)
        self.plotWidget.roiWidget.roiChanged.connect(self.roiWidget.updateROItext)
        
        self.layoutV = QtGui.QVBoxLayout()  
        self.layoutV.addWidget(self.loadFileWidget)
        self.layoutV.addWidget(self.calibrationWidget)
        self.layoutV.addWidget(self.rejectionWidget)
        self.layoutV.addWidget(self.roiWidget)
        self.layoutV.addWidget(self.xasWidget)
        self.layoutH = QtGui.QHBoxLayout()
        self.layoutH.addLayout(self.layoutV)
        self.layoutH.addWidget(self.plotWidget)   
        self.cw = QtGui.QWidget()        
        self.cw.setLayout(self.layoutH)
        self.setCentralWidget(self.cw)
        self.activeTabIndex = 0
        self.tabsList = [self.loadFileWidget, self.calibrationWidget, self.rejectionWidget, self.roiWidget, self.xasWidget]
        self.plotWidget.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Minimum))
        for tab in self.tabsList:
            tab.setMinimumWidth(400)
            tab.setSizePolicy(QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Minimum))
        self.statusBar().showMessage("Ready")
        XStream.stdout().messageWritten.connect(self.statusBar().showMessage)
                  
            
app = QtGui.QApplication(sys.argv)
main = MainWindow()  
main.show()
sys.exit(app.exec_())

