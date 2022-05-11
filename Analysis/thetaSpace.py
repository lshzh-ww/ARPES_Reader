import numpy
from numba import jit
import scipy.interpolate
import scipy.optimize

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

def fermiDiracDis(x,A,u):
    Temp=10.
    return A/(numpy.exp((x-u)/(Temp*8.61733e-5))+1.)

#@jit
def fixFermiZeroPoint(data,energyAxis,fitXStart,fitXEnd,fitYStart,fitYEnd,Temp):
    #points for get fermi arc
    N=11
    energyArray=energyAxis[fitYStart:fitYEnd]
    splineData=numpy.zeros((2,N),dtype=numpy.float64)
    xIndex=numpy.arange(fitXStart,fitXEnd+1,1)
    for i in range(len(data)):
        for j in range(N):
            xPos=int(j*(fitXEnd-fitXStart)/(N-1)+fitXStart)
            intensityArray=data[i,xPos,fitYStart:fitYEnd]
            try:
                result=scipy.optimize.curve_fit(fermiDiracDis,energyArray,intensityArray,[1,0])
            except RuntimeError:
                splineData[0,j]=xPos
            else:   
                splineData[0,j]=xPos
                delta=0
                for k in range(len(energyAxis)):
                    if energyAxis[k]>result[0][1] and energyAxis[k]<0.:
                        delta=delta+1
                    elif energyAxis[k]<result[0][1] and energyAxis[k]>0.:
                        delta=delta-1
                splineData[1,j]=delta
        print(splineData)
        tck=scipy.interpolate.splrep(splineData[0],splineData[1],s=0)
        diffArray=scipy.interpolate.splev(xIndex,tck)
        for j in range(len(xIndex)):
            data[i,xIndex[j]]=numpy.roll(data[i,xIndex[j]],int(diffArray[j]))

    return data


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


# rotate a 2D array respect to a given position
def rotate2D(data,angle,x0,y0):
    xAxisLen=len(data[0])
    yAxisLen=len(data)
    x0=int(x0)
    y0=int(y0)
    angle=numpy.radians(angle)
    newData=numpy.zeros((yAxisLen,xAxisLen))
    for i in range(xAxisLen):
        for j in range(yAxisLen):
            newData[j,i]=data[y0-yAxisLen//2+int(round(j*numpy.sin(angle)+i*numpy.cos(angle))),x0-xAxisLen//2+int(round(j*numpy.cos(angle)-i*numpy.sin(angle)))]
    return newData