import numpy
import math
import scipy.interpolate
import scipy.optimize
import scipy.ndimage
import scipy.spatial.transform.rotation

from numba import jit

#constant is sqrt(2*m*e)/(h bar) *10e-10
kConstant=0.5123167
pi=3.1415926

@jit(nopython=True)
def theta2kMap(thax,azi0,tilt0,pol0,bindingE,hv,innerPot):
    pi=3.1415926
    ec = 1.6e-19
    me=9.109e-31
    hbar=1.0546e-34
    eul=numpy.array([azi0,tilt0,pol0])
    eul*=pi/180
    electronE=bindingE+hv-5.
    rotm = numpy.array([[numpy.cos(eul[0])*numpy.cos(eul[1]),       numpy.sin(eul[1])*numpy.sin(eul[2])*numpy.cos(eul[0]) - numpy.sin(eul[0])*numpy.cos(eul[2]),      numpy.sin(eul[1])*numpy.cos(eul[0])*numpy.cos(eul[2]) + numpy.sin(eul[0])*numpy.sin(eul[2])],
                        [numpy.sin(eul[0])*numpy.cos(eul[1]),       numpy.sin(eul[0])*numpy.sin(eul[1])*numpy.sin(eul[2]) + numpy.cos(eul[0])*numpy.cos(eul[2]),      numpy.sin(eul[0])*numpy.sin(eul[1])*numpy.cos(eul[2]) - numpy.sin(eul[2])*numpy.cos(eul[0])],
                        [-numpy.sin(eul[1]),                        numpy.sin(eul[2])*numpy.cos(eul[1]),                                                              numpy.cos(eul[1])*numpy.cos(eul[2])]])
    kMap = numpy.zeros((len(thax),3))
    kMap[:,0]=numpy.sin(thax*pi/180.)
    kMap[:,2]=numpy.cos(thax*pi/180.)
    for i in range(len(thax)):
        kMap[i]=numpy.dot(rotm,kMap[i])
        kMap[i,0:2]=kMap[i,0:2]*numpy.sqrt(electronE*ec*2*me)/hbar*1e-10
        kMap[i,2]=(numpy.sqrt(((electronE*(kMap[i,2]**2)+innerPot)*ec)*2*me)/hbar)*1e-10
    return kMap


@jit(nopython=True)
def ToKzKxSpace(data, axis1Position, axis2Position, PhotonEnergyPosition, polar,v0):
    xRange=len(axis1Position)//2
    yRange=len(axis1Position)//2
    
    displayData=numpy.zeros((len(axis2Position),xRange,yRange),dtype=numpy.float32)
    kxMax=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+PhotonEnergyPosition[0])*math.sin(axis1Position[len(axis1Position)-1]*pi/180.0)*1.1
    kxMin=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+PhotonEnergyPosition[0])*math.sin(axis1Position[0]*pi/180.0)*1.1
    kzMax=kConstant*math.sqrt((axis2Position[len(axis2Position)-1]+PhotonEnergyPosition[0])*math.pow(math.cos(polar*pi/180.0),2)+v0)*1.1
    kzMin=kConstant*math.sqrt((axis2Position[0]+PhotonEnergyPosition[len(PhotonEnergyPosition)-1])*math.pow(math.cos(polar*pi/180.0)*math.cos(axis1Position[len(axis1Position)-1]*pi/180.0),2)+v0)*0.9


    for k in range(len(axis2Position)):
        for i in range(len(PhotonEnergyPosition)):
            for j in range(len(axis1Position)):
                kx=kConstant*math.sqrt(axis2Position[k]+PhotonEnergyPosition[i])*math.sin(axis1Position[j]*pi/180.0)
                kz=kConstant*math.sqrt((axis2Position[k]+PhotonEnergyPosition[i])*math.pow(math.cos(polar*pi/180.0)*math.cos(axis1Position[j]*pi/180.0),2)+v0)
                indexX=int(xRange*(kx-kxMin)/(kxMax-kxMin))
                indexY=int(yRange*(kz-kzMin)/(kzMax-kzMin))
                displayData[k,indexX,indexY]=data[i,j,k]


    print('done')
    print(kzMax,kzMin,kxMax,kxMin)
    kxPosition=numpy.arange(kxMin,kxMax,(kxMax-kxMin)/xRange)
    kzPosition=numpy.arange(kzMin,kzMax,(kzMax-kzMin)/yRange)
    energyPosition=axis2Position

    return displayData,kxPosition,kzPosition,energyPosition

@jit(nopython=True)
def ToKyKxSpace(data, axis1Position, axis2Position, polarPosition, Energy,v0):
    kxMax=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+Energy)*math.sin(axis1Position[len(axis1Position)-1]*pi/180.0)*1.1
    kxMin=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+Energy)*math.sin(axis1Position[0]*pi/180.0)*1.1
    kyMax=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+Energy)*math.sin(polarPosition[len(polarPosition)-1]*pi/180.0)*1.1
    kyMin=kConstant*math.sqrt(axis2Position[len(axis2Position)-1]+Energy)*math.sin(polarPosition[0]*pi/180.0)*1.1
    print(kxMax,kxMin,kyMax,kyMin)
    kzAvg=0.
    xRange=len(axis1Position)//2
    #yRange=len(axis1Position)//2
    yRange=int(xRange*(kyMax-kyMin)/(kxMax-kxMin))
    displayData=numpy.zeros((len(axis2Position),xRange,yRange),dtype=numpy.float32)

    for k in range(len(axis2Position)):
        for i in range(len(polarPosition)):
            kMap=theta2kMap(axis1Position,0,0,polarPosition[i],axis2Position[k],Energy,v0)
            for j in range(len(axis1Position)):
                kzAvg+=kMap[j,2]
                kx=kMap[j,0]
                ky=kMap[j,1]
                indexX=int(xRange*(kx-kxMin)/(kxMax-kxMin))
                indexY=int(yRange*(ky-kyMin)/(kyMax-kyMin))
                displayData[k,indexX,indexY]=data[i,j,k]

    kzAvg=kzAvg/len(axis1Position)/len(axis2Position)/len(polarPosition)
    print('done')
    print(kyMax,kyMin,kxMax,kxMin)
    kxPosition=numpy.arange(kxMin,kxMax,(kxMax-kxMin)/xRange)
    kzPosition=numpy.arange(kyMin,kyMax,(kyMax-kyMin)/yRange)
    energyPosition=axis2Position

    return displayData,kxPosition,kzPosition,energyPosition,kzAvg

@jit(nopython=True)
def dataNoneInterpolation(displayData):
    indexs=len(displayData)
    rows=len(displayData[0])
    lines=len(displayData[0,0])
    
    y1=0
    y2=0
    for i in range(indexs):
        for j in range(rows):
            for k in range(lines):
                if displayData[i,j,k]!=0:
                    if y1!=0:
                        y2=k
                        for l in range(y1,y2):
                            displayData[i,j,l]=displayData[i,j,y1]
                        y1=y2
                    else:
                        y1=k
            y1=0
            y2=0
    return displayData

@jit(nopython=True)
def dataLinearInterpolation(displayData):
    indexs=len(displayData)
    rows=len(displayData[0])
    lines=len(displayData[0,0])
    
    y1=0
    y2=0
    for i in range(indexs):
        for j in range(rows):
            for k in range(lines):
                if displayData[i,j,k]!=0:
                    if y1!=0:
                        y2=k
                        for l in range(y1,y2):
                            displayData[i,j,l]=displayData[i,j,y1]+(l-y1)*(displayData[i,j,y2]-displayData[i,j,y1])/(y2-y1)
                        y1=y2
                    else:
                        y1=k
            y1=0
            y2=0
    return displayData

@jit
def dataSplineInterpolation(displayData):
    indexs=len(displayData)
    rows=len(displayData[0])
    lines=len(displayData[0,0])
    
    splineData=numpy.zeros((2,lines),dtype=numpy.float32)
    interIndex=numpy.arange(0,lines,1)
    count=0
    for i in range(indexs):
        for j in range(rows):
            for k in range(lines):
                if displayData[i,j,k]!=0:
                    splineData[0,count]=k
                    splineData[1,count]=displayData[i,j,k]
                    count=count+1
            if count>3:
                tck=scipy.interpolate.splrep(splineData[0,0:count],splineData[1,0:count],s=0)
                displayData[i,j,int(splineData[0,0]):int(splineData[0,count-1])]=scipy.interpolate.splev(interIndex[int(splineData[0,0]):int(splineData[0,count-1])],tck)
            count=0
        print(i,'/',indexs)
    return displayData
    
def lorentzFunc(x,x0,I,gamma):
    return I/((x-x0)*(x-x0)+gamma*gamma)

def findBandPos(displayData,tPosition,xPosition,yPosition,range3D):
    tMin=range3D[0]
    tMax=range3D[1]
    xMin=range3D[2]
    xMax=range3D[3]
    yMin=range3D[4]
    yMax=range3D[5]
    xArray=xPosition[xMin:xMax]
    bandPoints=[]
    for i in numpy.arange(tMin,tMax,10):
        for k in numpy.arange(yMin,yMax,5):
            intensityArray=displayData[i,xMin:xMax,k]
            try:
                result=scipy.optimize.curve_fit(lorentzFunc,xArray,intensityArray,[xPosition[int((xMin+xMax)/2)],1,1])
            except RuntimeError:
                print(i,k)
            else:
                if result[0][0]<xPosition[xMax] and result[0][0]>xPosition[xMin]:
                    bandPoints.append([result[0][0],yPosition[k],tPosition[i]])

    return numpy.array(bandPoints)

@jit(nopython=True)
def minimumGradient(displayData):
    N=2
    indexs=len(displayData)
    rows=len(displayData[0])
    lines=len(displayData[0,0])

    gradientMap=numpy.zeros((indexs,rows,lines),dtype=numpy.float32)

    for i in range(N,indexs-N):
        for j in range(N,rows-N):
            for k in range(N,lines-N):
                temp=displayData[i,j,k]
                cache=0.

                cache+=(temp-displayData[i,j+N,k+N])**2
                cache+=(temp-displayData[i,j-N,k-N])**2
                cache+=(temp-displayData[i,j-N,k+N])**2
                cache+=(temp-displayData[i,j+N,k-N])**2
                
                cache/=2.

                cache+=(temp-displayData[i,j+N,k])**2
                cache+=(temp-displayData[i,j-N,k])**2
                cache+=(temp-displayData[i,j,k+N])**2
                cache+=(temp-displayData[i,j,k-N])**2

                if cache!=0:
                    gradientMap[i,j,k]=1./math.sqrt(cache)
            
    return gradientMap

def gaussian_blur(displayData):
    for index in range(len(displayData)):
        displayData[index]=scipy.ndimage.gaussian_filter(displayData[index],sigma=3)
    return displayData

@jit(nopython=True)
def integrateAlongEaxis(displayData,intRange):
    newDisplayData=numpy.zeros((len(displayData),len(displayData[0]),len(displayData[0,0])),dtype=numpy.float32)
    for i in range(intRange,len(displayData)-intRange):
        for j in range(len(displayData[0])):
            for k in range(len(displayData[0,0])):
                newDisplayData[i,j,k]=numpy.sum(displayData[i-intRange:i+intRange,j,k])

    return newDisplayData


def upscaleKgrid(qe_bands_contour,N):
    """Input data shuuld be qe bands calculation with kgrid type tpiba_c, read by np.loadtxt(). Grid should be a square. Will interpolate to NxN grid"""

    #find the first repeated point in qe_bands_contour[:,0]
    for i in range(len(qe_bands_contour[:,0])):
        if qe_bands_contour[i,0]==qe_bands_contour[i+1,0]:
            data_n=i
            break

    #create original kgrid position list
    original_points=numpy.zeros(((data_n+1)**2,2),dtype=numpy.float32)
    for i in range((data_n+1)):
        for j in range(data_n+1):
            original_points[i*(data_n+1)+j,0]=qe_bands_contour[i,0]
            original_points[i*(data_n+1)+j,1]=qe_bands_contour[j,0]        
    
    dense_y=numpy.linspace(qe_bands_contour[0,0],qe_bands_contour[data_n,0],N)
    dense_x=numpy.linspace(qe_bands_contour[0,0],qe_bands_contour[data_n,0],N)
    dense_x,dense_y=numpy.meshgrid(dense_x,dense_y)

    bands_num=int(len(qe_bands_contour[:,0])/((data_n+1)**2))
    bands_interpolation=numpy.zeros((bands_num,N,N),dtype=numpy.float32)

    for i in range(bands_num):
        bands_interpolation[i]=scipy.interpolate.griddata(original_points,qe_bands_contour[i*(data_n+1)**2:(i+1)*(data_n+1)**2,1],(dense_x,dense_y),method='cubic')

    return bands_interpolation,dense_x,dense_y

@jit(nopython=True)
def getBandsContourList(bands_interpolation,dense_x,dense_y,ref_E):
    list_kxky=list()
    delta_kx=abs(dense_x[0,1]-dense_x[0,0])
    delta_ky=abs(dense_y[1,0]-dense_y[0,0])
    for bands_i in range(len(bands_interpolation)):
        grid_z1=bands_interpolation[bands_i]
        if numpy.max(grid_z1)<ref_E or numpy.min(grid_z1)>ref_E:
            continue

        for i in range(1,len(dense_x)-1):
            for j in range(1,len(dense_x[0])-1):
                kx_0=dense_x[i,j]
                ky_0=dense_y[i,j]
                if grid_z1[i,j]>ref_E:
                    if grid_z1[i-1,j]<ref_E:
                        list_kxky.append((kx_0-delta_kx*(grid_z1[i,j]-ref_E)/(grid_z1[i,j]-grid_z1[i-1,j]),ky_0))
                    if grid_z1[i+1,j]<ref_E:
                        list_kxky.append((kx_0+delta_kx*(grid_z1[i,j]-ref_E)/(grid_z1[i,j]-grid_z1[i+1,j]),ky_0))
                    if grid_z1[i,j-1]<ref_E:
                        list_kxky.append((kx_0,ky_0-delta_ky*(grid_z1[i,j]-ref_E)/(grid_z1[i,j]-grid_z1[i,j-1])))
                    if grid_z1[i,j+1]<ref_E:
                        list_kxky.append((kx_0,ky_0+delta_ky*(grid_z1[i,j]-ref_E)/(grid_z1[i,j]-grid_z1[i,j+1])))

    contourArray=numpy.array(list_kxky)
    contourArray=numpy.append(contourArray,contourArray*numpy.array([1,-1]),axis=0)
    contourArray=numpy.append(contourArray,contourArray*numpy.array([-1,1]),axis=0)
    return contourArray