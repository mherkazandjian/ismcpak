"""
<keywords>
example, pdr, amuse, read, write, construct, database, meshes, archive
</keywords>
<description>
Example script that reads meshes writen by e.g run_dynamicGmechGrid.py
and converts the individual meshes into a single database. For testing
and demonstration purposes the written database is read back.

run using the command:

    $ cd $ISMCPAK/test
    $ python constructReadArchive.py

</description>
"""
from meshUtils import *
from time import *
import os

runDirPath =  '../../data/oneSided/surfaceGrid-z-1.0-no-gmech/'

# constructing the archive
t0 = time()
arxvW = meshArxv(dirPath=runDirPath)
arxvW.construct(meshNamePrefix='mesh', writeDb=True)
print 'time constructing %f' % (time() - t0)


# reading the archive
t0 = time()
arxvR = meshArxv( dirPath = runDirPath )
arxvR.readDb( )
arxvR.checkIntegrity()
print 'number of meshes in the database = %d' % arxvR.nMeshes
print 'time reading %f' % (time() - t0)
