from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *


db1Path = '/home/mher/ism/runs/oneSided/uniformSweepNew-1/'
db2Path = '/home/mher/ism/runs/oneSided/uniformSweepNew-2/'
db3Path = '/home/mher/ism/runs/oneSided/uniformSweepNew-1and2/'

# making the first database
arxv1 = meshArxv()
arxv1.construct( db1Path, writeDb = True )

# making the second database
arxv2 = meshArxv()
arxv2.construct( db2Path, writeDb = True )

# merging both
arxv1.mergeDbs(db2Path, db3Path)

# re-reading the merged database
t0 = time()
arxv3 = meshArxv(  )
arxv3.readDb( db3Path )
print 'time reading %f' % (time() - t0)
arxv3.checkIntegrity()


print 'done'
