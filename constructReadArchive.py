from mesh import *
from meshUtils import *
import numpy as np
from time import *

home = '/home/mher'

#runDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-0.5/'
#runDirPath =  home + '/ism/runs/oneSided/uniformSweep2-foo/'
#runDirPath =  home + '/ism/runs/oneSided/testOneSidedPDRGrid4/'
#runDirPath =  home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/'
#runDirPath =  home + '/ism/runs/oneSided/dynamicMeshTest1/'
#runDirPath =  home + '/ism/runs/oneSided/foo/'
runDirPath =  home + '/ism/runs/oneSided/dynamicMeshTest1-subset/'

# constructing the archive
t0 = time()
arxvW = meshArxv( dirPath = runDirPath )
arxvW.construct( writeDb = True )
print 'time constructing %f' % (time() - t0)


# reading the archive 
t0 = time()
arxvR = meshArxv( dirPath = runDirPath )
arxvR.readDb( )
print 'time reading %f' % (time() - t0)

# checking the accuracy
diff = 0.0
for i in np.arange(arxvR.nMeshes):
    mr = arxvR.meshes[i]
    diff += (mr['hdr']['G0']        - arxvR.infoAll[i]['parms'][0])   
    diff += (mr['hdr']['nGas']      - arxvR.infoAll[i]['parms'][1])   
    diff += (mr['hdr']['gammaMech'] - arxvR.infoAll[i]['parms'][2])
    
if diff == 0.0:
    print 'database reading check passed'
else:      
    print 'database might be corrupt'
