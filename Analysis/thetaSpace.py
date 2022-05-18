from multiprocessing.dummy import active_children
from matplotlib.pyplot import sci
import numpy
from numba import jit
import scipy.interpolate
import scipy.optimize
import scipy.ndimage
from functools import partial

def findMaxSlope(intensity):
    deltaH=4
    slopeArray=numpy.roll(intensity,deltaH)-numpy.roll(intensity,-deltaH)
    slopeArray[0:deltaH]=0.
    slopeArray[len(slopeArray)-deltaH:]=0.
    peak_list=list()
    for i in range(1,len(slopeArray)-2):
        if slopeArray[i]>slopeArray[i-1] and slopeArray[i]>slopeArray[i+1]:
            peak_list.append(i)
    for i in range(len(peak_list)-1):
        last_one=peak_list.pop()
        if slopeArray[peak_list[-1]] > 60.*abs(last_one):
            break
    return peak_list[-1]

def fermiDiracDis(Temp,sigma_wd,x,A,u,b):
    #8.61733e-5 is kb in eV/K
    try:
        return scipy.ndimage.gaussian_filter((A*x+b)/(numpy.exp((x-u)/(Temp*8.61733e-5))+1.),sigma_wd)
    except RuntimeWarning:
        return 0.

#fitting data by the fermi dirac distribution to get the zero point
def fitFermiDirac(energyAxis,intensityArray,Temp):
    # take convolution width as 10 meV
    convWidth=0.01
    sigma_wd=convWidth/abs(energyAxis[1]-energyAxis[0])

    # initialize parameters
    a_0=0.
    u_0=0.
    b_0=0.5*numpy.mean(intensityArray)
    # fitting
    try:
        result,err=scipy.optimize.curve_fit(partial(fermiDiracDis,Temp,sigma_wd),energyAxis,intensityArray,[a_0,u_0,b_0],maxfev=10000)
    except RuntimeError:
        return numpy.NaN
    return result[1]

@jit(nopython=True)
def maskData(data,N):
    indexs=len(data)
    rows=len(data[0])
    lines=len(data[0,0])
    
    for i in range(indexs):
        for j in range(rows):
            for k in range(lines):
                if j<N:
                    data[i,j,k]=0
                elif j>(rows-N):
                    data[i,j,k]=0
    
    return data

def normalizeData(data):
    rawIntensity=numpy.zeros(numpy.size(data,0))
    for i in range(numpy.size(data,0)):
        rawIntensity[i]=numpy.mean(data[i])

    for i in range(numpy.size(data,0)):
        data[i]=data[i]/rawIntensity[i]

    return data,rawIntensity

def simpleParabolaFunc(x,a,c):
    return a*x*x+c

def simpleLinearFunc(x,a,b):
    return a*x+b

def quickCalFL(data, fLPos, zeroThetaX):
    N=9
    fN=float(N)
    xAxisLen=float(len(data[0]))
    bins_array=numpy.zeros((N,len(data[0][0])))
    meta_para_fit=numpy.zeros((2,len(data)))
    for index in range(len(data)):
        for i in range(N):
            bins_array[i]=numpy.sum(data[index,int((float(i)-1.)*xAxisLen/fN):int(float(i)*xAxisLen/fN),:],0)
        
        xPos_list=list()
        yPos_list=list()
        for i in range(N):
            if numpy.sum(bins_array[i]) > 0.75* numpy.sum(bins_array[N//2]) and i!=0 and i!=N-1:
                yPos_list.append(findMaxSlope(bins_array[i]))
                xPos_list.append(int((float(i)-0.5)*xAxisLen/fN)-zeroThetaX)
        para=scipy.optimize.curve_fit(simpleParabolaFunc,xPos_list,yPos_list)
        #print(para[0])
        meta_para_fit[0,index]=para[0][0]
        meta_para_fit[1,index]=para[0][1]

    para=scipy.optimize.curve_fit(simpleLinearFunc,range(len(data)-10),meta_para_fit[0,:-10])
    meta_para_fit[0]=simpleLinearFunc(range(len(data)),para[0][0],para[0][1])
    para=scipy.optimize.curve_fit(simpleParabolaFunc,range(len(data)-10),meta_para_fit[1,:-10])
    meta_para_fit[1]=simpleParabolaFunc(range(len(data)),para[0][0],para[0][1])

    for index in range(len(data)):    
        for i in range(len(data[0])):
            step=int(simpleParabolaFunc(i-zeroThetaX,meta_para_fit[0,index],meta_para_fit[1,index]))
            data[index,i,:]=numpy.roll(data[index,i],-step+fLPos)
    
    return data

def slowCalFL(data,fLPos,zeroThetaX,energyAxis,Temp):
    N=15
    fN=float(N)
    #find location of the first value larger than -1 in energyAxis and the first value larger than 1 in energyAxis
    fittingRegion=numpy.zeros(2,dtype=int)
    fittingRegion[0]=numpy.where(energyAxis>-0.3)[0][0]
    #fittingRegion[1]=numpy.where(energyAxis>1.)[0][0]
    fittingRegion[1]=len(energyAxis)

    xAxisLen=float(len(data[0]))
    bins_array=numpy.zeros((N,len(data[0][0])))
    meta_para_fit=numpy.zeros((2,len(data)))
    for index in range(len(data)):
        for i in range(N):
            bins_array[i]=numpy.sum(data[index,int((float(i)-1.)*xAxisLen/fN):int(float(i)*xAxisLen/fN),:],0)
        
        xPos_list=list()
        yPos_list=list()
        for i in range(N):
            if numpy.sum(bins_array[i]) > 0.5* numpy.sum(bins_array[N//2]) and i!=0 and i!=N-1:
                #fit the data within fitting region
                result=fitFermiDirac(energyAxis[fittingRegion[0]:fittingRegion[1]],bins_array[i,fittingRegion[0]:fittingRegion[1]],Temp)
                if result!=numpy.NaN:
                    yPos_list.append(numpy.where(energyAxis>result)[0][0])
                    xPos_list.append(int((float(i)-0.5)*xAxisLen/fN)-zeroThetaX)
                else:
                    print("Filled at index:",index,'with value:',i)

        #paint fit result on the data by set the corresponding value to max value
        #max_brightness=numpy.max(data[index])
        #for i in range(len(xPos_list)):
        #    data[index,int(xPos_list[i]+zeroThetaX),yPos_list[i]]=max_brightness
            


        para=scipy.optimize.curve_fit(simpleParabolaFunc,xPos_list,yPos_list)
        meta_para_fit[0,index]=para[0][0]
        meta_para_fit[1,index]=para[0][1]

    para=scipy.optimize.curve_fit(simpleLinearFunc,range(len(data)-10),meta_para_fit[0,:-10])
    meta_para_fit[0]=simpleLinearFunc(range(len(data)),para[0][0],para[0][1])
    para=scipy.optimize.curve_fit(simpleParabolaFunc,range(len(data)-10),meta_para_fit[1,:-10])
    meta_para_fit[1]=simpleParabolaFunc(range(len(data)),para[0][0],para[0][1])

    for index in range(len(data)):    
        for i in range(len(data[0])):
            step=int(simpleParabolaFunc(i-zeroThetaX,meta_para_fit[0,index],meta_para_fit[1,index]))
            data[index,i,:]=numpy.roll(data[index,i],-step+fLPos)
    
    return meta_para_fit
    


# rotate a 2D array respect to a given position
def rotate2D(data,angle,x0,y0):
    xAxisLen=len(data)
    yAxisLen=len(data[0])
    newData=numpy.zeros((xAxisLen,yAxisLen))
    for i in range(xAxisLen):
        for j in range(yAxisLen):
            x_prime=x0+round(-(j-y0)*numpy.sin(angle)+(i-x0)*numpy.cos(angle))
            y_prime=y0+round((j-y0)*numpy.cos(angle)+(i-x0)*numpy.sin(angle))
            if x_prime>=0 and x_prime<xAxisLen and y_prime>=0 and y_prime<yAxisLen:
                newData[i,j]=data[x_prime,y_prime]
    
    return newData
    
@jit(nopython=True)
def rotate3D(data,angle,x0,y0):
    indexLen=len(data)
    xAxisLen=len(data[0])
    yAxisLen=len(data[0][0])
    newData=numpy.zeros((indexLen,xAxisLen,yAxisLen))
    for i in range(xAxisLen):
        for j in range(yAxisLen):
            x_prime=x0+round(-(j-y0)*numpy.sin(angle)+(i-x0)*numpy.cos(angle))
            y_prime=y0+round((j-y0)*numpy.cos(angle)+(i-x0)*numpy.sin(angle))
            if x_prime>=0 and x_prime<xAxisLen and y_prime>=0 and y_prime<yAxisLen:
                newData[:,i,j]=data[:,x_prime,y_prime]
    
    return newData