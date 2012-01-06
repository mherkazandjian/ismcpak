from numpy import *

from time import *

from chemicalNetwork import *
from enumSpecies import *
from mesh import *

meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/mesh.dat-id-002920'
abunSpecFname = 'data/species.inp'

specFileUsed = enumSpecies(abunSpecFname, offsetFromZero=1 )
eSpcs = specFileUsed.genPythonClass(offset=0)

m = mesh(meshFname);

m.plot(eSpcs)
pyl.show()