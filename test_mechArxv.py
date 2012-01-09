import os
import glob
from mesh import *
from meshUtils import *
import numpy as np

meshDirPath = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/'
arxvInfo = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes.info'
arxvData = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes.db'

ver = np.zeros([3], dtype = np.int32)
arxvHdrDtype = np.dtype( [ ('info' , np.int64  , 5),
                           ('parms', np.float64, 5)
                         ] # 64 bytes
                       )

ver[0] = 0
ver[1] = 0
ver[2] = 0  
#--------------------------------------------------------------------
dbInfoFObj = file(arxvInfo, 'wb')
dbDataFObj = file(arxvData, 'wb')

files = []
for infile in glob.glob( os.path.join(meshDirPath, '*id-0000*') ):
    files.append(infile)

nMeshes = np.zeros( 1, dtype = np.int32 )
nMeshes[0] = len(files)

print 'found %s' % (nMeshes)

# writing the version number to the arxive header file
ver.tofile( dbInfoFObj)
    
hdrArr = np.ndarray( nMeshes, dtype = arxvHdrDtype )

# constructing the database info header file and 
# dumping the meshes data into big files
for i in np.arange(nMeshes):
    
    fName = files[i]
    print fName
    
    m = mesh(fName)
    mData = m.getData()

    hdrArr[i]['info'][0] = i # mesh number 
    hdrArr[i]['info'][1] = 0 # data file index
    hdrArr[i]['info'][2] = mData['hdr']['nSteps']  
    hdrArr[i]['info'][3] = mData['hdr']['nSpecs']
    hdrArr[i]['info'][4] = dbDataFObj.tell()  # offset from the start of file
     

    hdrArr[i]['parms'][0] = mData['hdr']['G0']
    hdrArr[i]['parms'][1] = mData['hdr']['nGas']
    hdrArr[i]['parms'][2] = mData['hdr']['gammaMech']
    hdrArr[i]['parms'][3] = np.float64(0.0)
    hdrArr[i]['parms'][4] = np.float64(0.0)

    mData.tofile( dbDataFObj )

dbDataFObj.close()

nMeshes.tofile( dbInfoFObj )
hdrArr.tofile( dbInfoFObj )
dbInfoFObj.close()

#-----------------------------------------------------------------------------


dbInfoFObj = file(arxvInfo, 'rb')
dbDataFObj = file(arxvData, 'rb')

# re-reading the header
verR = np.fromfile( dbInfoFObj, dtype = (np.int32, 3), count = 1)
verR = verR[0]

nMeshesR = np.fromfile( dbInfoFObj, dtype = np.int32, count = 1)
nMeshesR = nMeshesR[0]

hdrArrR = np.fromfile( dbInfoFObj, dtype = arxvHdrDtype, count = nMeshesR)

print verR
print nMeshesR
print hdrArrR['info'] - hdrArr['info']
dbInfoFObj.close()

mDummy = mesh()

# re-reading the data into a list
dbDataFObj = file(arxvData, 'rb')
for i in np.arange(nMeshesR):
    print i
    
    nSteps = hdrArrR[i]['info'][2]
    nSpecs = hdrArrR[i]['info'][3]
    
    mDtype = mDummy.constructMeshDtype(nSpecs, nSteps)
    m = np.fromfile( dbDataFObj, dtype = mDtype, count = 1)

    print m['hdr']['G0']        - hdrArrR[i]['parms'][0]
    print m['hdr']['nGas']      - hdrArrR[i]['parms'][1]
    print m['hdr']['gammaMech'] - hdrArrR[i]['parms'][2]
    print '---------------------------------'


