from hashlib import new
from locale import currency
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
from PyQt5.QtWidgets import QMainWindow, qApp, QApplication, QFileDialog, QInputDialog ,QDialog
from PyQt5.QtGui import QTransform
from pathlib import Path


import GUI.arpes_main
import GUI.image_translation

class MainWindow(QMainWindow, GUI.arpes_main.Ui_MainWindow):
    def showEnergyFunc(self):
        self.bindingEnergyDisplay.setText('Energy='+"%.4f" % self.axis2Position[self.topGraphWidget.currentIndex]+'eV')

    #generate QTransform for the image based on currentStatus
    def generateTransform(self):
        trans=QTransform()
        if self.currentStatus=='Ready':
            trans.translate(self.axis1Position[0],self.axis2Position[0])
            trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
            return trans
        elif self.currentStatus=='kxky':
            trans.translate(self.kxPosition[0],self.kyPosition[0])
            trans.scale((self.kxPosition[len(self.kxPosition)-1]-self.kxPosition[0])/len(self.kxPosition),(self.kyPosition[len(self.kyPosition)-1]-self.kyPosition[0])/len(self.kyPosition))
            return trans
        elif self.currentStatus=='kxkz':
            trans.translate(self.kxPosition[0],self.kzPosition[0])
            trans.scale((self.kxPosition[len(self.kxPosition)-1]-self.kxPosition[0])/len(self.kxPosition),(self.kzPosition[len(self.kzPosition)-1]-self.kzPosition[0])/len(self.kzPosition))
            return trans

    def loadFileDiag(self):
        home_dir = str(Path.home())
        if self.lastTimeDir!='':
            home_dir=self.lastTimeDir
        fname = QFileDialog.getOpenFileName(self,'Open file', home_dir)
        self.lastTimeDir=os.path.dirname(fname[0])
        self.rawData=FileIO.loadData.load(self,fname[0])
        #set energy axis
        self.fermi_pos=Analysis.thetaSpace.findMaxSlope(np.sum(self.rawData[len(self.rawData)//2],0))
        print(self.axis2Position[self.fermi_pos])
        self.axis2Position=self.axis2Position-self.axis2Position[self.fermi_pos]
        dataShape=np.shape(self.rawData)
        print(dataShape)
        self.currentStatus='Ready'
        trans=QTransform()
        trans.translate(self.axis1Position[0],self.axis2Position[0])
        trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
        self.topGraphWidget.setImage(self.rawData,transform=trans)
        #self.topGraphWidget.view.invertY(False)
        self.topGraphWidget.view.getViewBox().setAspectLocked(lock=False)
        self.topGraphWidget.view.setLabel('bottom','θ','°')
        self.topGraphWidget.view.setLabel('left','E_b','eV')
        self.displayData=self.rawData
        
    def convertFileDiag(self):
        home_dir = str(Path.home())
        fname = QFileDialog.getOpenFileName(self, 'Open file', home_dir)
        FileIO.loadData.TXT2NPY(self,fname[0])
        
    def normalizeFunc(self):
        self.rawData=Analysis.thetaSpace.normalizeData(self.rawData)[0]
        self.topGraphWidget.setImage(self.rawData,transform=self.generateTransform())

    def maskFunc(self):
        text, ok = QInputDialog.getText(self, 'Mask Range?','How many lines do you want to remove from the left and right part of the data')
        if ok:
            self.rawData=Analysis.thetaSpace.maskData(self.rawData,int(text))
            self.topGraphWidget.setImage(self.rawData,transform=self.generateTransform())
        else:
            return None

    def fixFermiSurfFunc(self):
        for i in range(len(self.axis1Position)):
            if self.axis1Position[i]*self.axis1Position[i+1] <0.:
                break
        #Ask for the temperature and use the temperature to calculate the fermi surface
        text, ok = QInputDialog.getText(self, 'Fix Fermi Surface?','What is the temperature?')
        if ok:
            Analysis.thetaSpace.slowCalFL(self.rawData, self.fermi_pos,i,self.axis2Position,float(text))
            self.topGraphWidget.setImage(self.rawData,transform=self.generateTransform())

    def findBandPosFunc(self):
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

    def kxkyConvertFunc(self):
        timeStart=time()
        self.polarPosition=np.zeros(np.size(self.rawData,0))
        deltaPolar=(float(self.polarEndBox.text())-float(self.polarStartBox.text()))/(np.size(self.rawData,0)-1)
        for i in range(np.size(self.rawData,0)):
            self.polarPosition[i]=float(self.polarStartBox.text())+i*deltaPolar
        
        self.kSpaceData=Analysis.kSpace.ToKyKxSpace(self.rawData, self.axis1Position, self.axis2Position, self.polarPosition, int(self.energyBox.text()), int(self.v0Box.text()))
        
        if self.radioBtn21.isChecked():
            self.displayData=Analysis.kSpace.dataNoneInterpolation(self.kSpaceData[0])
        elif self.radioBtn22.isChecked():
            self.displayData=Analysis.kSpace.dataLinearInterpolation(self.kSpaceData[0])
        else:
            self.displayData=self.kSpaceData[0]
        
        self.kxPosition=self.kSpaceData[1]
        self.kyPosition=self.kSpaceData[2]
        self.topGraphWidget.view.getViewBox().setAspectLocked(lock=True)
        self.topGraphWidget.view.setLabel('bottom','kx','1/Å')
        self.topGraphWidget.view.setLabel('left','ky','1/Å')
        self.currentStatus='kxky'
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.bandPos=np.zeros((0,3))

    def kxkzConvertFunc(self):
        self.photonEnergyPosition=np.zeros(np.size(self.rawData,0))
        deltaE=(int(self.energyStartBox.text())-int(self.energyEndBox.text()))/(np.size(self.rawData,0)-1)
        for i in range(np.size(self.rawData,0)):
            self.photonEnergyPosition[i]=int(self.energyStartBox.text())-i*deltaE
        
        self.kSpaceData=Analysis.kSpace.ToKzKxSpace(self.rawData, self.axis1Position, self.axis2Position, self.photonEnergyPosition, int(self.polarBox.text()), int(self.v0Box.text()))
        
        if self.radioBtn21.isChecked():
            self.displayData=Analysis.kSpace.dataNoneInterpolation(self.kSpaceData[0])
        elif self.radioBtn22.isChecked():
            self.displayData=Analysis.kSpace.dataLinearInterpolation(self.kSpaceData[0])
        else:
            self.displayData=self.kSpaceData[0]
        
        self.kxPosition=self.kSpaceData[1]
        self.kzPosition=self.kSpaceData[2]
        self.topGraphWidget.view.getViewBox().setAspectLocked(lock=True)
        self.topGraphWidget.view.setLabel('bottom','kx','1/Å')
        self.topGraphWidget.view.setLabel('left','kz','1/Å')
        self.currentStatus='kxkz'
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())


    def minGradFunc(self):
        self.displayData=Analysis.kSpace.minimumGradient(self.displayData)
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())

    def gaussianBlurFunc(self):
        self.displayData=Analysis.kSpace.gaussian_blur(self.displayData) 
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())


    def integrateAlongEaxisFunc(self):
        #pop a window to get the integration range
        text,ok=QInputDialog.getText(self,'Integration Range','Enter the range of integration along the energy axis (eV) (e.g. 0.01)')
        if ok:
            fRange=float(text)
            iRange=int(fRange/abs(self.axis2Position[1]-self.axis2Position[0]))
            print(iRange)
            self.displayData=Analysis.kSpace.integrateAlongEaxis(self.displayData,iRange)
            #reset scale and display new data
            self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())
            

    def openAlignFermiSurfFunc(self):
        self.transformList=[]
        ch=Dialog_translate(self,self.displayData,self.kxPosition,self.kyPosition)
        ch.open()

    def applyAlignFermiSurface(self):
        axis1Zero=np.where(self.kxPosition>0)[0][0]
        axis2Zero=np.where(self.kyPosition>0)[0][0]
        print(self.transformList)
        for i in range(len(self.transformList)):
            if self.transformList[i][0]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][0],axis=1)
            elif self.transformList[i][1]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][1],axis=2)
            elif self.transformList[i][2]!=0:
                self.displayData=Analysis.thetaSpace.rotate3D(self.displayData,self.transformList[i][2],axis1Zero,axis2Zero)

        #display the new data
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())

    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        self.lastTimeDir=''
        self.topGraphWidget.view=pg.PlotItem()
        self.topGraphWidget.view.setAspectLocked(True)
        self.topGraphWidget.ui.graphicsView.setCentralItem(self.topGraphWidget.view)
        self.topGraphWidget.view.addItem(self.topGraphWidget.roi)
        self.topGraphWidget.view.addItem(self.topGraphWidget.normRoi)
        self.topGraphWidget.view.addItem(self.topGraphWidget.imageItem)
        self.topGraphWidget.view.register(self.topGraphWidget.name)

        self.topGraphWidget.sigTimeChanged.connect(self.showEnergyFunc)
        self.exitAct.setStatusTip('Exit application')
        self.exitAct.triggered.connect(qApp.quit)
        self.loadFileAct.triggered.connect(self.loadFileDiag)
        self.txt2npyAct.triggered.connect(self.convertFileDiag)
        self.normalAct.triggered.connect(self.normalizeFunc)
        self.maskAct.triggered.connect(self.maskFunc)
        self.fixFermiAct.triggered.connect(self.fixFermiSurfFunc)
        self.findBandPosAct.triggered.connect(self.findBandPosFunc)
        self.minGradAct.triggered.connect(self.minGradFunc)
        self.kxkyBtn.clicked.connect(self.kxkyConvertFunc)
        self.kxkzBtn.clicked.connect(self.kxkzConvertFunc)
        self.gaussianBlurAct.triggered.connect(self.gaussianBlurFunc)
        self.integrateAlongEaxisAct.triggered.connect(self.integrateAlongEaxisFunc)
        self.alignFermiSufAct.triggered.connect(self.openAlignFermiSurfFunc)

        #default select tab 1.
        self.tabWidget.setCurrentIndex(0)
        self.v0Box.setText('13')

        #set main window in the center of screen
        screen = qApp.primaryScreen()
        size = screen.size()
        w = self.geometry().width()
        h = self.geometry().height()
        x = (size.width() - w) // 2
        y = (size.height() - h) // 2
        self.move(x, y)

class Dialog_translate(QDialog, GUI.image_translation.Ui_Dialog):
    def generateTransform(self):
        trans=QTransform()
        trans.translate(self.axis1Position[0],self.axis2Position[0])
        trans.scale((self.axis1Position[len(self.axis1Position)-1]-self.axis1Position[0])/len(self.axis1Position),(self.axis2Position[len(self.axis2Position)-1]-self.axis2Position[0])/len(self.axis2Position))
        return trans

    def loadImageViewStatus(self):
        result=np.zeros(3)
        result[0]=self.leftGraphWidget.currentIndex
        minMaxLevels=self.leftGraphWidget.getHistogramWidget().getLevels()
        result[1]=minMaxLevels[0]
        result[2]=minMaxLevels[1]
        return result

    def restoreImageViewStatus(self,status):
        self.leftGraphWidget.setCurrentIndex(int(status[0]))
        self.leftGraphWidget.getHistogramWidget().setLevels(status[1],status[2])

    def moveRightFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaX=int(self.deltaXBox.text())
        #roll all elements in 2D matrix self.displayData up by deltaX
        if np.sum(self.transformList[-1]*np.array([deltaX,0,0]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([deltaX,0,0])
        else:
            self.transformList.append(np.array([deltaX,0,0]))
        self.displayData=np.roll(self.displayData,deltaX,axis=1)
        #set new display data
        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)
        

    def moveLeftFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaX=int(self.deltaXBox.text())
        if np.sum(self.transformList[-1]*np.array([-deltaX,0,0]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([-deltaX,0,0])
        else:
            self.transformList.append(np.array([-deltaX,0,0]))
        self.displayData=np.roll(self.displayData,-deltaX,axis=1)
        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)

    def moveUpFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaY=int(self.deltaXBox.text())
        if np.sum(self.transformList[-1]*np.array([0,deltaY,0]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([0,deltaY,0])
        else:
            self.transformList.append(np.array([0,deltaY,0]))
        self.displayData=np.roll(self.displayData,deltaY,axis=2)
        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)
    
    def moveDownFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaY=int(self.deltaXBox.text())
        if np.sum(self.transformList[-1]*np.array([0,-deltaY,0]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([0,-deltaY,0])
        else:
            self.transformList.append(np.array([0,-deltaY,0]))
        self.displayData=np.roll(self.displayData,-deltaY,axis=2)
        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)

    def clockwiseFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaDeg=float(self.deltaDegBox.text())*np.pi/180.
        if np.sum(self.transformList[-1]*np.array([0,0,deltaDeg]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([0,0,deltaDeg])
        else:
            self.transformList.append(np.array([0,0,deltaDeg]))

        self.displayData=self.rawData
        for i in range(len(self.transformList)):
            if self.transformList[i][0]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][0],axis=1)
            elif self.transformList[i][1]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][1],axis=2)
            elif self.transformList[i][2]!=0:
                self.displayData=Analysis.thetaSpace.rotate3D(self.displayData,self.transformList[i][2],self.axis1Zero,self.axis2Zero)

        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)

    def counterClockwiseFunc(self):
        currentStatus=self.loadImageViewStatus()
        deltaDeg=float(self.deltaDegBox.text())*np.pi/180.
        if np.sum(self.transformList[-1]*np.array([0,0,-deltaDeg]))!=0:
            self.transformList[-1]=self.transformList[-1]+np.array([0,0,-deltaDeg])
        else:
            self.transformList.append(np.array([0,0,-deltaDeg]))

        self.displayData=self.rawData
        for i in range(len(self.transformList)):
            if self.transformList[i][0]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][0],axis=1)
            elif self.transformList[i][1]!=0:
                self.displayData=np.roll(self.displayData,self.transformList[i][1],axis=2)
            elif self.transformList[i][2]!=0:
                self.displayData=Analysis.thetaSpace.rotate3D(self.displayData,self.transformList[i][2],self.axis1Zero,self.axis2Zero)

        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        self.restoreImageViewStatus(currentStatus)

    def addRefLineTool(self):
        #open a QInputDialog to get the reference line parameters
        text,ok=QInputDialog.getText(self,'Reference Line','Enter the reference line parameters in the form of "x1,y1,x2,y2,..."')
        if ok:
            #store in a array
            refLine=np.array(text.split(','),dtype=float)
            print(refLine)
            #check if the array is of the correct size
            if len(refLine)%2==0:
                imv_v=self.leftGraphWidget.getView()
                refPlot=pg.PlotDataItem(x=refLine[0::2],y=refLine[1::2],pen=pg.mkPen(color='r',width=1))
                imv_v.addItem(refPlot)

  


    def acceptFunc(self):
        print('Accept change')
        self.parent().transformList=self.transformList
        self.parent().applyAlignFermiSurface()
        self.close()
    
    def __init__(self, parent=None, displayData=None, axis1Position=None, axis2Position=None):
        super(Dialog_translate, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle('Align Image')

        self.axis1Position=axis1Position
        self.axis2Position=axis2Position
        self.displayData=displayData
        self.rawData=displayData
        self.transformList=list()
        self.transformList.append(np.array([0,0,0]))
        #find zero elements in axis1 and axis2
        self.axis1Zero=np.where(self.axis1Position>0)[0][0]
        self.axis2Zero=np.where(self.axis2Position>0)[0][0]


        #initialize pg.ImageView with view=pg.PlotItem()
        self.leftGraphWidget.view=pg.PlotItem()
        self.leftGraphWidget.view.setAspectLocked(True)
        self.leftGraphWidget.ui.graphicsView.setCentralItem(self.leftGraphWidget.view)
        self.leftGraphWidget.view.addItem(self.leftGraphWidget.roi)
        self.leftGraphWidget.view.addItem(self.leftGraphWidget.normRoi)
        self.leftGraphWidget.view.addItem(self.leftGraphWidget.imageItem)
        self.leftGraphWidget.view.register(self.leftGraphWidget.name)

        self.leftGraphWidget.setImage(self.displayData,transform=self.generateTransform())
        imv_v=self.leftGraphWidget.getView()
        ref_x=pg.PlotCurveItem(x=[self.axis1Position[0],self.axis1Position[len(self.axis1Position)-1]],y=[0,0],pen=pg.mkPen('r',width=1))
        ref_y=pg.PlotCurveItem(x=[0,0],y=[self.axis2Position[0],self.axis2Position[len(self.axis2Position)-1]],pen=pg.mkPen('r',width=1))
        imv_v.addItem(ref_x)
        imv_v.addItem(ref_y)


        
        self.deltaXBox.setText('5')
        self.deltaDegBox.setText('0.5')
        self.upButton.clicked.connect(self.moveUpFunc)
        self.downButton.clicked.connect(self.moveDownFunc)
        self.leftButton.clicked.connect(self.moveLeftFunc)
        self.rightButton.clicked.connect(self.moveRightFunc)
        self.clockwiseButton.clicked.connect(self.clockwiseFunc)
        self.counterClockwiseButton.clicked.connect(self.counterClockwiseFunc)
        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.acceptFunc)
        self.toolButton.clicked.connect(self.addRefLineTool)
        
    def open(self):
        self.show()


        
def main():
    app = QApplication(sys.argv)
    mainWindow = MainWindow()
    
    mainWindow.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

