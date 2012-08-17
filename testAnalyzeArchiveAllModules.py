from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *

#---------------------------Archive parameters-----------------------
#runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/'
runDirPath    = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/'
#runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'
gridsRes      = 10
lgammaMechSec = -30.0
metallicity   = 1.0
plotRangenG0  = [[0,6],[0,6]]

radexParms    = { 'radexPath'         : '/home/mher/ism/code/radex/Radex/bin/radex',  
                  'molDataDirPath'    : '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                  'specStr'           : 'CO',
#                  'xH2_Min'          : 2*0.0000000001
                  'xH2_Min'           : -1.0,
                  'collisionPartners' : ['H2','H+','H'],
#                  'collisionPartners' : ['H2','H','H+','e-']
#                  'collisionPartners' : ['H2']
#                  'collisionPartners' : ['H2','H+','e-','H']
                  'plotTransitionInGrid' : 0
                }
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']


# elements and basic species from which all the other species are made
import baseSpecies
baseSpec = baseSpecies.baseSpecies()

#------------------------------------------------------------------
# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv( metallicity = metallicity )
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

pyl.ioff()
arxv.plotGrid(gridsRes, lgammaMechSec, radexParms, ranges = plotRangenG0)

pyl.show()
"""
#arxv.saveGridsToFiles(gridsRes, lgammaMechSec, radexParms)
pyl.show()
"""
print 'done'
