from numpy import *

from time import *
import pylab as pyl
from chemicalNetwork import *
from enumSpecies import *
from mesh import *

HOME = '/home/mher'
#------------------mesh data file-------------------------------------
metallicity   = 1.0
meshFname     = HOME + '/ism/runs/oneSided/single_mesh/meshes/mesh.dat-id-000000'
#-----------------chemical network parameters------------------------
rxnFile       = HOME + '/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = HOME + '/ism/code/ismcpak/data/species.inp'
underAbunFile = HOME + '/ism/code/ismcpak/data/underabundant.inp'
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
net.remove_species( underAbunFile = underAbunFile )
net.remove_species( species = removeManual )
# reading the species number and their corresponding indies and abundances from ascii files
net.assign_numbers_to_species(fileName = specNumFile)

m = mesh(meshFname, chemNet=net, metallicity=metallicity)
m.plot(plot_Av_range=[0,30])

pyl.show()
print 'done'

