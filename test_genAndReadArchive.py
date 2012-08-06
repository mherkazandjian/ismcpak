from mesh import *
from meshUtils import *
import numpy as np
from time import *

#runDirPath =  '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/'
runDirPath = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/'

# constructing the archive
t0 = time()
arxvW = meshArxv(  )
arxvW.construct( runDirPath )
print 'time constructing %f' % (time() - t0)

# reading the archive 
t0 = time()
arxvR = meshArxv(  )
arxvR.readDb( runDirPath )
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