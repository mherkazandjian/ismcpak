from meshUtils import *
from time import *

home = '/home/mher'

#runDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-0.5/'
#runDirPath =  home + '/ism/runs/oneSided/uniformSweep2-foo/'
#runDirPath =  home + '/ism/runs/oneSided/testOneSidedPDRGrid4/'
#runDirPath =  home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/'
#runDirPath =  home + '/ism/runs/oneSided/dynamicMeshTest1/'
#runDirPath =  home + '/ism/runs/oneSided/foo/'
#runDirPath =  home + '/ism/runs/oneSided/dynamicMeshTest1-subset/'
#runDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-1.0/'
#runDirPath =  home + '/ism/runs/oneSided/surfaceGrid-z-0.1/'
#runDirPath =  home + '/ism/runs/oneSided/singleModels-z-2.0/' 
#runDirPath =  home + '/ism/runs/oneSided/sph-db-z-0.2/' 
#runDirPath =  home + '/ism/runs/oneSided/dynamicMesh-z-2.0/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-sweep/'
runDirPath =  home + '/ism/runs/oneSided/dynamicMesh-z-1.0-750-mw-CR/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1_solar/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10_solar/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-100_solar/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1000_solar/'
#runDirPath =  home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10000_solar/'

#runDirPath =  home + '/ism/runs/oneSided/surfaceGrid-z-2.0-high-res-no-gmech/'
#runDirPath =  home + '/ism/runs/oneSided/singleModels-z-2.0/' 

# constructing the archive
t0 = time()
arxvW = meshArxv(dirPath = runDirPath)
arxvW.construct( meshNamePrefix = 'mesh', writeDb = True )
print 'time constructing %f' % (time() - t0)


# reading the archive 
t0 = time()
arxvR = meshArxv( dirPath = runDirPath )
arxvR.readDb( )
arxvR.checkIntegrity()
print 'time reading %f' % (time() - t0)