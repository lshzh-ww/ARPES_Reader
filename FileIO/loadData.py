import scipy.io
import numpy
from PyQt5.QtWidgets import QInputDialog,QFileDialog
from pathlib import Path
import zipfile
import configparser
import io

def load(self,filename):
    fileType=filename[len(filename)-3:]
    if fileType=='npy':
        return npyFile(self,filename)
    elif fileType=='mat':
        return matFile(self,filename)
    elif fileType=='txt':
        return txtFile(self,filename)
    elif fileType=='zip':
        return zipFile(self,filename)   


def matFile(self,filename):
    data=scipy.io.loadmat(filename)
    allKeys=data.keys()
    text, ok = QInputDialog.getText(self, 'Input Dialog',str(allKeys))
    if ok:
        M = data.get(str(text))
    else:
        return None
    M=numpy.swapaxes(M,0,2)
    return M


#For DA30 mode at BNL
def zipFile(self,filename):
    da30File=zipfile.ZipFile(filename)
    infoWrapper=io.TextIOWrapper(da30File.open('Spectrum_XPS.ini'),encoding='utf-8')
    config = configparser.ConfigParser()
    config.read_file(infoWrapper)
    self.axis1Position=numpy.linspace(float(config['spectrum']['heightoffset']),float(config['spectrum']['heightoffset'])+(float(config['spectrum']['height'])-1)*float(config['spectrum']['heightdelta']),int(config['spectrum']['height']))
    self.axis2Position=numpy.linspace(float(config['spectrum']['widthoffset']),float(config['spectrum']['widthoffset'])+(float(config['spectrum']['width'])-1)*float(config['spectrum']['widthdelta']),int(config['spectrum']['width']))
    self.polarStartBox.setText(config['spectrum']['depthoffset'])
    endPolar=float(config['spectrum']['depthoffset'])+(float(config['spectrum']['depth'])-1)*float(config['spectrum']['depthdelta'])
    self.polarEndBox.setText(str(endPolar))
    M=numpy.frombuffer(da30File.read('Spectrum_XPS.bin'),dtype=numpy.float32)
    M=numpy.reshape(M,(int(config['spectrum']['depth']),int(config['spectrum']['height']),int(config['spectrum']['width'])))

    
    infoWrapper2=io.TextIOWrapper(da30File.open('XPS.ini'),encoding='utf-8')
    config2 = configparser.ConfigParser(strict=False)
    config2.read_file(infoWrapper2)
    self.energyBox.setText(config2['SES']['Excitation Energy'])
    M=M.copy()
    da30File.close()
    print(M.shape)
    #M=numpy.swapaxes(M,0,2)
    return M

def txtFile(self,filename):
    data=numpy.loadtxt(filename)
    text, ok = QInputDialog.getText(self, 'Input Dialog','data nums, data length and data width?')
    dataShape=numpy.fromstring(text,dtype=numpy.int,sep=' ')
    data=data.reshape(dataShape[0],dataShape[1],dataShape[2])
    data=numpy.swapaxes(data,1,2)
    self.axis1Position=data[0,1:,0]
    self.axis2Position=data[0,0,1:]
    M=data[:,1:,1:]
    return M
    
def npyFile(self,filename):
    data=numpy.load(filename)
    print(data.shape)
    #self.axis1Position=data[0,1:,0]
    #self.axis2Position=data[0,0,1:]
    self.axis1Position=numpy.linspace(-16.853389+60*0.048221,-16.853389+738*0.048221,679)
    self.axis2Position=numpy.linspace(54.426453,54.426453+609*0.001635,610)
    #return data[:,1:,1:]
    return data[:,:,:]


def TXT2NPY(self,filename):
    data=numpy.loadtxt(filename)
    text, ok = QInputDialog.getText(self, 'Input Dialog','data nums, data length and data width?')
    dataShape=numpy.fromstring(text,dtype=numpy.int,sep=' ')
    data=data.reshape(dataShape[0],dataShape[1],dataShape[2])
    data=numpy.swapaxes(data,1,2)
    home_dir = str(Path.home())
    saveFilename = QFileDialog.getSaveFileName(self, 'Save File', home_dir)
    print(saveFilename[0])
    numpy.save(saveFilename[0],data)
    return data


