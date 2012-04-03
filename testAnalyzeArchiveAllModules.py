from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *

# runs
#   uniformSweep2-z-2/
#   uniformSweep2-highRes-z-1.0
#   uniformSweep2-z-0.5
#   
#   uniformSweep2-z-1.0-no-mech
#   uniformSweep2-z-2-no-mech
#   uniformSweep2-z-0.5-no-mech
# 
#---------------------------Archive parameters-----------------------
#runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'
runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'


abunSpecFname = '/home/mher/ism/code/ismcpak/data/species.inp'
lgammaMechSec = -20.0
metallicity   = 2.0
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
specAbunFile  = '/home/mher/ism/code/ismcpak/data/abun.out'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

gridsRes = 20

T        = 800.0
Av       = 20.0
nDens    = 10**5.5
G0       = 10**5
zeta     = 5e-17
albedo   = 0.6

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
# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv(  )
arxv.readDb( runDirPath )
arxv.checkIntegrity()
print 'time reading %f' % (time() - t0)
#------------------------------------------------------------------
# read and setting up the chemical network used in the 
t0 = time()
# settin up the orignial netowrk
net = chemicalNetwork(rxnFile, baseSpec, UMISTVER = 'umist99')
# reading the species to be removed from a file
net.removeSpecies( underAbunFile = underAbunFile )
net.removeSpecies( species = removeManual )
# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = specNumFile)
arxv.setChemicalNetwork(net) # assiginig the chemical network to the archive
#-------------------------------------------------------------------
# plotting stuff

arxv.plotGrid(gridsRes, lgammaMechSec )

print 'done'
