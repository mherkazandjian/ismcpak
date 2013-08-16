from mesh import *
from meshUtils import *
import numpy as np
from time import *

def constructAndReadArxiv( runDirPath ):

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


constructAndReadArxiv( '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1_solar/' )
constructAndReadArxiv( '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10_solar/' )
constructAndReadArxiv( '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-100_solar/' )
constructAndReadArxiv( '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-1000_solar/' )
constructAndReadArxiv( '/home/mher/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-10000_solar/' )


