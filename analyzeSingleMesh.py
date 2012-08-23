from numpy import *

from time import *
import pylab as pyl
from chemicalNetwork import *
from enumSpecies import *
from mesh import *
from ismUtils import *
from utils import fetchNestedDtypeValue

#------------------mesh data file-------------------------------------
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000-0.1.7'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat'
#meshFname      = '/data1/mher/ism/runs/oneSided/testOneSidedPDRGrid/meshes/mesh.dat-id-000000'
meshFname      = '/data1/mher/ism/runs/oneSided/runSingleMesh/meshes/mesh.dat-10000percent'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat-adaptive-0.1.10'
#meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/meshes/mesh.dat-id-000119'
metallicity   = 1.0
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

# elements and basic species from which all the other species are made
import baseSpecies
baseSpec = baseSpecies.baseSpecies()

#------------------------------------------------------------------
net = chemicalNetwork(rxnFile, baseSpec, UMISTVER = 'umist99')
# reading the species to be removed from a file
net.removeSpecies( underAbunFile = underAbunFile )
net.removeSpecies( species = removeManual )
# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = specNumFile)

# setting up the PDR slab
m = mesh(meshFname, net, metallicity)
TMean1, nColls, NCOLVG1, = m.getRadexParameters('H2', 'CO', 2*0.01)
print TMean1, nColls, NCOLVG1

m = mesh( meshFname, net, metallicity )
m.plot()

pyl.show()
print 'done'

# 378.37764026 1182.21211321 8.05680609553e+18