from ast import Num
from hashlib import new
from locale import currency
import os

from matplotlib.pyplot import connect, sci
import scipy.optimize
import FileIO
import CustomWidgets.plot2
import Analysis
from functools import partial

import sys
from time import time
import numpy as np
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from PyQt5.QtWidgets import QMainWindow, qApp, QApplication, QFileDialog, QInputDialog ,QDialog
from PyQt5.QtGui import QTransform
from pathlib import Path


import GUI

class MainWindow(QMainWindow, GUI.arpes_main.Ui_MainWindow):
    def showEnergyFunc(self):
        if self.currentStatus=='kxky':
            self.bindingEnergyDisplay.setText('Energy='+"%.4f" % self.axis2Position[self.topGraphWidget.currentIndex]+'eV     kz='+"%.4f" % self.kSpaceData[4]+'Å')
            try:
                self.topGraphWidget.getHistogramWidget().setLevels(self.blackWhiteLeveslRatio[0]*np.mean(self.displayData[self.topGraphWidget.currentIndex]),self.blackWhiteLeveslRatio[1]*np.mean(self.displayData[self.topGraphWidget.currentIndex]))
            except AttributeError:
                pass
    def getBlackWhiteLevelsFunc(self):
        try:
            if np.mean(self.displayData[self.topGraphWidget.currentIndex]!=0):
                self.blackWhiteLeveslRatio=self.topGraphWidget.getHistogramWidget().getLevels()/np.mean(self.displayData[self.topGraphWidget.currentIndex])
        except AttributeError:
            pass
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
        myFileDialog=QFileDialog()
        myFileDialog.setFileMode(QFileDialog.ExistingFiles)
        fname = myFileDialog.getOpenFileNames(self,'Open file', home_dir)
        if len(fname[0])==1:
            self.lastTimeDir=os.path.dirname(fname[0][0])
            self.rawData=FileIO.loadData.load(self,fname[0][0])
            #set energy axis
            self.fermi_pos=Analysis.thetaSpace.findMaxSlope(np.sum(self.rawData[len(self.rawData)//2],0))
            print(self.axis2Position[self.fermi_pos])
            self.axis2Position=self.axis2Position-self.axis2Position[self.fermi_pos]
        elif len(fname[0])>1:
            self.lastTimeDir=os.path.dirname(fname[0][0])
            for i in range(len(fname[0])):
                if i==0:
                    rawData1=FileIO.loadData.load(self,fname[0][i])
                    tempAxis1=self.axis1Position.copy()
                    tempAxis2=self.axis2Position.copy()
                    self.rawData=np.zeros((len(fname[0]),len(rawData1),len(rawData1[0])))
                    self.rawData[i]=rawData1.copy()
                    result,err=scipy.optimize.curve_fit(partial(Analysis.thetaSpace.fermiDiracDis,20,0.01),self.axis2Position,np.sum(self.rawData[i],0),p0=[0,np.mean(self.axis2Position),0.5*np.mean(np.sum(self.rawData[i],0))],maxfev=10000)
                    self.fermi_pos=np.where(self.axis2Position>result[1])[0][0]
                else:
                    rawData1=FileIO.loadData.load(self,fname[0][i])
                    self.rawData[i]=rawData1.copy()
                    result,err=scipy.optimize.curve_fit(partial(Analysis.thetaSpace.fermiDiracDis,20,0.01),self.axis2Position,np.sum(self.rawData[i],0),p0=[0,np.mean(self.axis2Position),0.5*np.mean(np.sum(self.rawData[i],0))],maxfev=10000)
                    c_fermi_pos=np.where(self.axis2Position>result[1])[0][0]
                    self.rawData[i]=np.roll(self.rawData[i],-self.fermi_pos+c_fermi_pos,axis=1)
                
                self.axis1Position=tempAxis1
                self.axis2Position=tempAxis2
                #set energy axis
                self.axis2Position=self.axis2Position-self.axis2Position[self.fermi_pos]
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
        
    def saveScanGeometryFunc(self):
        if self.transformList!=None:
            fname = QFileDialog.getSaveFileName(self, 'Save file', self.lastTimeDir)
            np.save(fname[0],self.transformList, allow_pickle=True)
    
    def loadScanGeometryFunc(self):
        home_dir = str(Path.home())
        if self.lastTimeDir!='':
            home_dir=self.lastTimeDir
        fname = QFileDialog.getOpenFileName(self, 'Open file', home_dir)
        self.transformList=np.load(fname[0], allow_pickle=True)
        self.applyAlignFermiSurface()
    

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
            self.rawData=Analysis.thetaSpace.slowCalFL(self.rawData, self.fermi_pos,i,self.axis2Position,float(text))
            newAxis2Length=len(self.rawData[0][0])
            deltaE=self.axis2Position[1]-self.axis2Position[0]
            newAxis2=np.linspace(self.axis2Position[0]-deltaE*(newAxis2Length-len(self.axis2Position)),self.axis2Position[len(self.axis2Position)-1],newAxis2Length)
            self.axis2Position=newAxis2
            self.topGraphWidget.setImage(self.rawData,transform=self.generateTransform())

    def findBandPosFunc(self):
        #broken, do not use
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
        ch=Dialog_translate(self,self.displayData,self.kxPosition,self.kyPosition,self.axis2Position)
        ch.open()

    def applyAlignFermiSurface(self):
        axis1Zero=np.where(self.kxPosition>0)[0][0]
        axis2Zero=np.where(self.kyPosition>0)[0][0]
        print(self.transformList)
        for i in range(len(self.transformList)):
            if self.transformList[i][0]!=0:
                self.displayData=np.roll(self.displayData,int(self.transformList[i][0]),axis=1)
            elif self.transformList[i][1]!=0:
                self.displayData=np.roll(self.displayData,int(self.transformList[i][1]),axis=2)
            elif self.transformList[i][2]!=0:
                self.displayData=Analysis.thetaSpace.rotate3D(self.displayData,self.transformList[i][2],axis1Zero,axis2Zero)

        #display the new data
        self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())

    def openGetCutFunc(self):
        self.cutInfo=list()
        ch=Dialog_getCut(self)
        ch.open()

    def generateCut(self):
        print(self.cutInfo)
        cutData=np.zeros((len(self.kxPosition),len(self.axis2Position)))
        kLength=np.sqrt(np.sum(np.square(self.cutInfo[0]-self.cutInfo[1])))
        newKaxis=np.linspace(0,kLength,len(cutData))
        pointSeries=np.zeros((len(cutData),2))
        pointSeries[:,0]=np.linspace(self.cutInfo[0][0],self.cutInfo[1][0],len(cutData))
        pointSeries[:,1]=np.linspace(self.cutInfo[0][1],self.cutInfo[1][1],len(cutData))
        deltaN=int(self.cutInfo[2]/(self.kxPosition[1]-self.kxPosition[0]))
        theta=np.arctan((self.cutInfo[1]-self.cutInfo[0])[1]/(self.cutInfo[1]-self.cutInfo[0])[0])
        distanceArray=np.sum(np.square(pointSeries),axis=1)
        minIndex=np.argmin(distanceArray)
        newKaxis=newKaxis-newKaxis[minIndex]
        for i in range(len(cutData)):
            x0=np.where(self.kxPosition>pointSeries[i,0])[0][0]
            y0=np.where(self.kyPosition>pointSeries[i,1])[0][0]

            for h in range(deltaN):
                deltaX=round(-h*np.sin(theta))
                deltaY=round(h*np.cos(theta))
                if x0+deltaX>=0 and x0+deltaX<len(self.kxPosition):
                    if y0+deltaY>=0 and y0+deltaY<len(self.kyPosition):
                        cutData[i,:]=self.displayData[:,x0+deltaX,y0+deltaY]

        ch=Dialog_translate(self,cutData,newKaxis,self.axis2Position)
        ch.leftGraphWidget.view.getViewBox().setAspectLocked(lock=False)
        ch.open()

        
        

    def symmetrizeFunc(self):
        #pop a window to get the N-fold symmetry
        text,ok=QInputDialog.getText(self,'N-fold Symmetry','Enter the N-fold symmetry (e.g. 2)') 
        if ok:
            if text =='mirror':
                x0=np.where(self.kxPosition>0)[0][0]
                newDisplayData=self.displayData.copy()
                newDisplayData=np.roll(np.flip(newDisplayData,axis=1),2*x0-len(newDisplayData[0]),axis=1)
                newDisplayData+=self.displayData
                self.displayData=newDisplayData
                self.topGraphWidget.setImage(self.displayData,transform=self.generateTransform())

            else:
                n=int(text)
                newDisplayData=self.displayData.copy()
                x0=np.where(self.kxPosition>0)[0][0]
                y0=np.where(self.kyPosition>0)[0][0]
                for i in range(n):
                    newDisplayData+=Analysis.thetaSpace.rotate3D(self.displayData,float(i)*2.*np.pi/float(n),x0,y0)
                self.displayData=newDisplayData
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
        self.topGraphWidget.getHistogramWidget().sigLevelsChanged.connect(self.getBlackWhiteLevelsFunc)
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
        self.alignFermiSufAct_2.triggered.connect(self.openAlignFermiSurfFunc)
        self.symmetrizeAct.triggered.connect(self.symmetrizeFunc)
        self.getCutAct.triggered.connect(self.openGetCutFunc)
        self.actionLoad_Existed_Result.triggered.connect(self.loadScanGeometryFunc)
        self.actionSave_Current_Result.triggered.connect(self.saveScanGeometryFunc)

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
    def indexChangedFunc(self):
        try:
            self.leftGraphWidget.getHistogramWidget().setLevels(self.blackWhiteLeveslRatio[0]*np.mean(self.displayData[self.leftGraphWidget.currentIndex]),self.blackWhiteLeveslRatio[1]*np.mean(self.displayData[self.leftGraphWidget.currentIndex]))
        except AttributeError:
            pass
        self.setWindowTitle('Energy='+"%.4f" % self.energyPosition[self.leftGraphWidget.currentIndex]+'eV')
        

        imv_v=self.leftGraphWidget.getView()
        #if existed self.refPlot, then remove it.
        if hasattr(self,'bands_energy'):
            if hasattr(self,'refBandsPlot'):
                imv_v.removeItem(self.refBandsPlot)
                del self.refBandsPlot
            contourPoints=Analysis.kSpace.getBandsContourList(self.bands_energy,self.bands_kx,self.bands_ky,self.energyPosition[self.leftGraphWidget.currentIndex]+float(self.bandsRefEnergy.text()))
            if len(contourPoints) >0:
                self.refBandsPlot=pg.PlotDataItem(x=contourPoints[:,0],y=contourPoints[:,1],connect=np.zeros(len(contourPoints),dtype=np.int32),symbol='o',symbolPen=pg.mkPen(color='r'),symbolSize=1,pxMode=True)
                imv_v.addItem(self.refBandsPlot)



    def getBlackWhiteLevelsFunc(self):
        if np.mean(self.displayData[self.leftGraphWidget.currentIndex]!=0):
            self.blackWhiteLeveslRatio=self.leftGraphWidget.getHistogramWidget().getLevels()/np.mean(self.displayData[self.leftGraphWidget.currentIndex])


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

    def addBandsTool(self):
        #open a FileDialog to get the band file
        bandFileName=QFileDialog.getOpenFileName(self,'Open Band File','.','*.gnu')
        if bandFileName[0]!='':
            bandsData=np.loadtxt(bandFileName[0])
            N=60
            #interpolate the data to N points
            self.bands_energy,self.bands_kx,self.bands_ky=Analysis.kSpace.upscaleKgrid(bandsData,N)

    def acceptFunc(self):
        print('Accept change')
        self.parent().transformList=self.transformList
        self.parent().applyAlignFermiSurface()
        self.close()
    
    def __init__(self, parent=None, displayData=None, axis1Position=None, axis2Position=None, energyPosition=None):
        super(Dialog_translate, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle('Align Image')

        self.axis1Position=axis1Position
        self.axis2Position=axis2Position
        self.energyPosition=energyPosition
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
        #ref_x=pg.PlotCurveItem(x=[self.axis1Position[0],self.axis1Position[len(self.axis1Position)-1]],y=[0,0],pen=pg.mkPen('r',width=1))
        #ref_y=pg.PlotCurveItem(x=[0,0],y=[self.axis2Position[0],self.axis2Position[len(self.axis2Position)-1]],pen=pg.mkPen('r',width=1))
        #imv_v.addItem(ref_x)
        #imv_v.addItem(ref_y)


        
        self.deltaXBox.setText('5')
        self.deltaDegBox.setText('0.5')
        self.leftGraphWidget.sigTimeChanged.connect(self.indexChangedFunc)
        self.leftGraphWidget.getHistogramWidget().sigLevelsChanged.connect(self.getBlackWhiteLevelsFunc)
        self.upButton.clicked.connect(self.moveUpFunc)
        self.downButton.clicked.connect(self.moveDownFunc)
        self.leftButton.clicked.connect(self.moveLeftFunc)
        self.rightButton.clicked.connect(self.moveRightFunc)
        self.clockwiseButton.clicked.connect(self.clockwiseFunc)
        self.counterClockwiseButton.clicked.connect(self.counterClockwiseFunc)
        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.acceptFunc)
        self.toolButton.clicked.connect(self.addRefLineTool)
        self.addBandsButton.clicked.connect(self.addBandsTool)
        
    def open(self):
        self.show()


class Dialog_getCut(QDialog, GUI.getCut_Dialog.Ui_Dialog):
    def __init__(self, parent=None):
        super(Dialog_getCut, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle('Get Cut Info')

        self.deltaKBox.setText('0.1')

        self.buttonBox.rejected.connect(self.close)
        self.buttonBox.accepted.connect(self.acceptFunc)
    
    def acceptFunc(self):
        print('Accept change')
        self.parent().cutInfo=list()
        self.parent().cutInfo.append(np.fromstring(self.startPointBox.text(),dtype=float,sep=','))
        self.parent().cutInfo.append(np.fromstring(self.endPointBox.text(),dtype=float,sep=','))
        self.parent().cutInfo.append(float(self.deltaKBox.text()))
        self.parent().generateCut()
        self.close()

    def open(self):
        return self.show()
        
def main():
    app = QApplication(sys.argv)
    mainWindow = MainWindow()
    
    mainWindow.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

