from meshUtils import *
from time import *
import os

home = os.environ['HOME']

runDirPath =  home + '/ism/runs/tests/dynamicGrid/'

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
