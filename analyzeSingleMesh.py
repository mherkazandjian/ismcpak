from numpy import *

from time import *
import pylab as pyl
from chemicalNetwork import *
from enumSpecies import *
from mesh import *
from ismUtils import *

#------------------mesh data file-------------------------------------
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000-0.1.7'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat'
meshFname      = '/data1/mher/ism/runs/oneSided/testOneSidedPDRGrid/meshes/mesh.dat-id-000000'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat-adaptive-0.1.10'
#meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/meshes/mesh.dat-id-000119'
metallicity   = 1.0
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

# elements and basic species from which all the other species are made
baseSpec = [  specie('CRPHOT', specType = -1, charge=0 , init=1),
              specie('PHOTON', specType = -1, charge=0 , init=1),
              specie('CRP'   , specType = -1, charge=0 , init=1),
              specie('PAH'   , specType = 0 , charge=0 , init=1),
              specie('H2V'   , specType = 0 , charge=0 , init=1, comp = [ ['H',2] ]),
              specie('13C'   , specType = 0 , charge=0 , init=1),
              specie('Na'    , specType = 0 , charge=0 , init=1),
              specie('Mg'    , specType = 0 , charge=0 , init=1),
              specie('Si'    , specType = 0 , charge=0 , init=1),
              specie('Cl'    , specType = 0 , charge=0 , init=1),
              specie('Fe'    , specType = 0 , charge=0 , init=1),
              specie('He'    , specType = 0 , charge=0 , init=1),
              specie('H'     , specType = 0 , charge=0 , init=1),
              specie('M'     , specType = -1, charge=0 , init=1),
              specie('C'     , specType = 0 , charge=0 , init=1),
              specie('N'     , specType = 0 , charge=0 , init=1),
              specie('O'     , specType = 0 , charge=0 , init=1),
              specie('P'     , specType = 0 , charge=0 , init=1),
              specie('S'     , specType = 0 , charge=0 , init=1),
              specie('e-'    , specType = 0 , charge=-1, init=1) ]
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

m.plot()
fig2 = pyl.figure()
"""

pyl.plot(m.data['therm']['heating'],'r')
pyl.hold(True)
pyl.plot(m.data['therm']['cooling'],'ro')
"""

pyl.semilogy(m.data['state']['Av'], m.data['state']['gasT'],'r')

pyl.show()
print 'done'

# 378.37764026 1182.21211321 8.05680609553e+18