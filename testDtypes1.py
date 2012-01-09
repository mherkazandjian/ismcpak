import os
import glob
from mesh import *
from meshUtils import *
import pickle
from numpy import *
import numpy as np
from chemicalNetwork import *
from enumSpecies import *
from mesh import *
from string import *
import os

meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/mesh.dat-id-002920'
abunSpecFname = 'data/species.inp'

specFileUsed = enumSpecies(abunSpecFname, offsetFromZero=1 )
eSpcs = specFileUsed.genPythonClass(offset=0)

m = mesh(meshFname);

#print m
#m.plot(eSpcs)
#pyl.show()

fName = '/home/mher/ism/foo.npy'

data = m.getData()
dt = data.dtype
print data['hdr']

data.tofile(fName)
fName = '/home/mher/ism/foo.npy'
np.fromfile(fName, dtype=dt)
print data['hdr'] # not sure if it is loading right!!!

# loading is not working
data = m.getData()
dt = data.dtype
print data['hdr']

np.save(fName, data)
np.load(fName) 