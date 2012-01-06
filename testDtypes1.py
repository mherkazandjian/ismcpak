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

meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/mesh.dat-id-002920'
abunSpecFname = 'data/species.inp'

specFileUsed = enumSpecies(abunSpecFname, offsetFromZero=1 )
eSpcs = specFileUsed.genPythonClass(offset=0)

m = mesh(meshFname);

#print m
m.plot(eSpcs)
pyl.show()

data = m.getData()
m.saveData('/home/mher/ism/foo.sav')

