import os
import FileIO.loadData
import CustomWidgets.plot2
import Analysis.kSpace
import Analysis.thetaSpace
from functools import partial

import sys
from time import time
import numpy as np
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from PyQt5.QtWidgets import QMainWindow, qApp, QApplication, QFileDialog, QInputDialog
from PyQt5.QtGui import QTransform
from pathlib import Path


import arpes_gui



def showEnergyFunc(self,ui):
    ui.bindingEnergyDisplay.setText('Energy='+"%.4f" % self.axis2Position[ui.topGraphWidget.currentIndex]+'eV')

def loadFileDiag(self,ui):
    home_dir = str(Path.home())
    if self.lastTimeDir!='':
        home_dir=self.lastTimeDir
    fname = QFileDialog.getOpenFileName(self,'Open file', home_dir)
    self.lastTimeDir=os.path.dirname(fname[0])
    self.rawData=FileIO.loadData.load(self,ui,fname[0])
    #set energy axis
    self.fermi_pos=Analysis.thetaSpace.findMaxSlope(np.sum(self.rawData[len(self.rawData)//2],0))
    print(self.axis2Position[self.fermi_pos])
    self.axis2Position=self.axis2Position-self.axis2Position[self.fermi_pos]
    dataShape=np.shape(self.rawData)
    print(dataShape)
    trans=QTransform()
    trans.translate(self.axis1Position[0],self.axis2Position[0])
    trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
    #ui=self.findChildren(QWidget,name='centralwidget')[0]
    #self.centralWidget().topGraphWidget.setImage(self.rawData,transform=trans)
    ui.topGraphWidget.setImage(self.rawData,transform=trans)
    ui.topGraphWidget.view.invertY(False)
    ui.topGraphWidget.view.getViewBox().setAspectLocked(lock=False)
    #self.topGraphWidget.view.getViewBox().scaleBy(x=1.,y=0.1)
    ui.topGraphWidget.view.setLabel('bottom','θ','°')
    ui.topGraphWidget.view.setLabel('left','E_b','eV')
    self.displayData=self.rawData
    
def convertFileDiag(self,ui):
    home_dir = str(Path.home())
    fname = QFileDialog.getOpenFileName(self, 'Open file', home_dir)
    FileIO.loadData.TXT2NPY(self,fname[0])
    
def normalizeFunc(self,ui):
    self.rawData=Analysis.thetaSpace.normalizeData(self.rawData)[0]
    trans=QTransform()
    trans.translate(self.axis1Position[0],self.axis2Position[0])
    self.trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
    ui.topGraphWidget.setImage(self.rawData,transform=trans)

def maskFunc(self,ui):
    text, ok = QInputDialog.getText(self, 'Mask Range?','How many lines do you want to remove from the left and right part of the data')
    if ok:
        self.rawData=Analysis.thetaSpace.maskData(self.rawData,int(text))
        trans=QTransform()
        trans.translate(self.axis1Position[0],self.axis2Position[0])
        trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
        ui.topGraphWidget.setImage(self.rawData,transform=trans)
    else:
        return None

def fixFermiSurfFunc(self,ui):
    for i in range(len(self.axis1Position)):
        if self.axis1Position[i]*self.axis1Position[i+1] <0.:
            break
    
    Analysis.thetaSpace.quickCalFL(self.rawData, self.fermi_pos,i)
    trans=QTransform()
    trans.translate(self.axis1Position[0],self.axis2Position[0])
    trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
    ui.topGraphWidget.setImage(self.rawData,transform=trans)

def findBandPosFunc(self,ui):
    xyRange=self.topGraphWidget.view.getViewBox().viewRange()
    for i in range(len(self.kxPosition)):
        if self.kxPosition[i]<xyRange[0][0]:
            xMin=i
        if self.kxPosition[i]<xyRange[0][1]:
            xMax=i

    for i in range(len(self.kzPosition)):
        if self.kzPosition[i]<xyRange[1][0]:
            yMin=i
        if self.kzPosition[i]<xyRange[1][1]:
            yMax=i

    for i in range(len(self.axis2Position)):
        if self.axis2Position[i]<-0.5:
            tMin=i
        if self.axis2Position[i]<0.015:
            tMax=i  

    self.bandPos=np.vstack((self.bandPos,Analysis.kSpace.findBandPos(self.displayData,self.axis2Position,self.kxPosition,self.kzPosition,[tMin,tMax,xMin,xMax,yMin,yMax])))
    
    #self.bandPos2=Analysis.kSpace.findBandPos(self.displayData,self.axis2Position,self.kxPosition,self.kzPosition,[tMin,tMax,len(self.kxPosition)-xMax,len(self.kxPosition)-xMin,yMin,yMax])
    
    self.bandDisplay = gl.GLViewWidget()
    self.plotUI = CustomWidgets.plot2.Ui_Form()
    self.plotUI.setupUi(self.bandDisplay)
    self.bandDisplay.show()

    g = gl.GLGridItem()
    self.bandDisplay.addItem(g)
    
        
    sp1 = gl.GLScatterPlotItem(pos=self.bandPos, size=0.025, color=[1.0,0.,0.,0.75], pxMode=False)
    #sp2 = gl.GLScatterPlotItem(pos=self.bandPos2, size=0.025, color=[1.0,0.,0.,0.75], pxMode=False)
    sp1.translate(0.,-self.kzPosition[int((yMax+yMin)/2)],0.)
    #sp2.translate(0.,-self.kzPosition[int((yMax+yMin)/2)],0.)
    self.bandDisplay.addItem(sp1)
    #self.bandDisplay.addItem(sp2)

def kxkyConvertFunc(self,ui):
    timeStart=time()
    self.polarPosition=np.zeros(np.size(self.rawData,0))
    deltaPolar=(float(ui.polarEndBox.text())-float(ui.polarStartBox.text()))/(np.size(self.rawData,0)-1)
    for i in range(np.size(self.rawData,0)):
        self.polarPosition[i]=float(ui.polarStartBox.text())+i*deltaPolar
    
    self.kSpaceData=Analysis.kSpace.ToKyKxSpace(self.rawData, self.axis1Position, self.axis2Position, self.polarPosition, int(ui.energyBox.text()), int(ui.v0Box.text()))
    
    if ui.radioBtn21.isChecked():
        self.displayData=Analysis.kSpace.dataNoneInterpolation(self.kSpaceData[0])
    elif ui.radioBtn22.isChecked():
        self.displayData=Analysis.kSpace.dataLinearInterpolation(self.kSpaceData[0])
    else:
        self.displayData=self.kSpaceData[0]
    
    self.kxPosition=self.kSpaceData[1]
    self.kyPosition=self.kSpaceData[2]
    
    ui.topGraphWidget.view.setLabel('bottom','kx','1/Å')
    ui.topGraphWidget.view.setLabel('left','ky','1/Å')
    trans=QTransform()
    trans.translate(self.kxPosition[0],self.kyPosition[0])
    trans.scale((self.kxPosition[len(self.kxPosition)-1]-self.kxPosition[0])/len(self.kxPosition),(self.kyPosition[len(self.kyPosition)-1]-self.kyPosition[0])/len(self.kyPosition))
    ui.topGraphWidget.setImage(self.displayData,transform=trans)
    print(time()-timeStart)

    self.bandPos=np.zeros((0,3))

def kxkzConvertFunc(self,ui):
    self.photonEnergyPosition=np.zeros(np.size(self.rawData,0))
    deltaE=(int(ui.energyStartBox.text())-int(ui.energyEndBox.text()))/(np.size(self.rawData,0)-1)
    for i in range(np.size(self.rawData,0)):
        self.photonEnergyPosition[i]=int(ui.energyStartBox.text())-i*deltaE
    
    self.kSpaceData=Analysis.kSpace.ToKzKxSpace(self.rawData, self.axis1Position, self.axis2Position, self.photonEnergyPosition, int(ui.polarBox.text()), int(ui.v0Box.text()))
    
    if ui.radioBtn21.isChecked():
        self.displayData=Analysis.kSpace.dataNoneInterpolation(self.kSpaceData[0])
    elif ui.radioBtn22.isChecked():
        self.displayData=Analysis.kSpace.dataLinearInterpolation(self.kSpaceData[0])
    else:
        self.displayData=self.kSpaceData[0]
    
    self.kxPosition=self.kSpaceData[1]
    self.kzPosition=self.kSpaceData[2]
    
    ui.topGraphWidget.view.setLabel('bottom','kx','1/Å')
    ui.topGraphWidget.view.setLabel('left','kz','1/Å')
    trans=QTransform()
    trans.translate(self.kxPosition[0],self.kzPosition[0])
    trans.scale((self.kxPosition[len(self.kxPosition)-1]-self.kxPosition[0])/len(self.kxPosition),(self.kzPosition[len(self.kzPosition)-1]-self.kzPosition[0])/len(self.kzPosition))
    self.topGraphWidget.setImage(self.displayData,transform=trans)
    

    self.bandPos=np.zeros((0,3))

def minGradFunc(self,ui):
    self.displayData=Analysis.kSpace.minimumGradient(self.displayData)
    #self.trans.reset()
    # self.trans.translate(self.kxPosition[0],self.kyPosition[0])
    # self.trans.scale((self.kxPosition[len(self.kxPosition)-1]-self.kxPosition[0])/len(self.kxPosition),(self.kyPosition[len(self.kyPosition)-1]-self.kyPosition[0])/len(self.kyPosition))
    ui.topGraphWidget.setImage(self.displayData)



def gaussianBlurFunc(self,ui):
    self.displayData=Analysis.kSpace.gaussian_blur(self.displayData) 
    ui.topGraphWidget.setImage(self.displayData)

        
        
def main():
    app = QApplication(sys.argv)
    mainWindow = QMainWindow()
    ui=arpes_gui.Ui_MainWindow()
    ui.setupUi(mainWindow)
    
    mainWindow.show()

    mainWindow.lastTimeDir=''
    
    ui.topGraphWidget.view=pg.PlotItem()
    ui.topGraphWidget.view.setAspectLocked(True)
    ui.topGraphWidget.ui.graphicsView.setCentralItem(ui.topGraphWidget.view)
    ui.topGraphWidget.view.addItem(ui.topGraphWidget.roi)
    ui.topGraphWidget.view.addItem(ui.topGraphWidget.normRoi)
    ui.topGraphWidget.view.addItem(ui.topGraphWidget.imageItem)
    ui.topGraphWidget.view.register(ui.topGraphWidget.name)

    ui.topGraphWidget.sigTimeChanged.connect(partial(showEnergyFunc,mainWindow,ui))
    ui.exitAct.setStatusTip('Exit application')
    ui.exitAct.triggered.connect(qApp.quit)
    ui.loadFileAct.triggered.connect(partial(loadFileDiag,mainWindow,ui))
    ui.txt2npyAct.triggered.connect(partial(convertFileDiag,mainWindow,ui))
    ui.normalAct.triggered.connect(partial(normalizeFunc,mainWindow,ui))
    ui.maskAct.triggered.connect(partial(maskFunc,mainWindow,ui))
    ui.fixFermiAct.triggered.connect(partial(fixFermiSurfFunc,mainWindow,ui))
    ui.findBandPosAct.triggered.connect(partial(findBandPosFunc,mainWindow,ui))
    ui.minGradAct.triggered.connect(partial(minGradFunc,mainWindow,ui))
    ui.kxkyBtn.clicked.connect(partial(kxkyConvertFunc,mainWindow,ui))
    ui.kxkzBtn.clicked.connect(partial(kxkzConvertFunc,mainWindow,ui))
    ui.gaussianBlurAct.triggered.connect(partial(gaussianBlurFunc,mainWindow,ui))

    #default select tab 1.
    ui.tabWidget.setCurrentIndex(0)

    #set main window in the center of screen
    screen = app.primaryScreen()
    size = screen.size()
    w = mainWindow.geometry().width()
    h = mainWindow.geometry().height()
    x = (size.width() - w) // 2
    y = (size.height() - h) // 2
    mainWindow.move(x, y)

    

    

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

