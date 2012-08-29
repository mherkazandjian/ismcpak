from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *

#---------------------------Archive parameters-----------------------
runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweepNew-1and2/'
metallicity   = 1.0

#quantity       = ['state', 'gasT']
quantity       = ['therm', 'heating']
#quantity       = ['fineStructureCoolingComponents', 'C+', 'rate', '1-0']
plotRange_nG0  = [[0,6],[0,6]]
slabIdx        = 0
res            = [4, 4]
adaptive       = False   # fixed section in mechanical heating
lgammaMechSec  = -30.0
#adaptive       = True    # adaptive grid in percent of surface heating
#lgammaMechSec  = 0.001

log10            = True

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

arxv.showGrid(quantity = quantity,
              slabIdx  = slabIdx,
              ranges   = plotRange_nG0,
              res      = res,
              zSec     = lgammaMechSec,
              adaptive = adaptive,
              log10    = log10)

"""
surfHeatinFac = 1.0
# defining the points in the uniform 2D grid
ranges = plotRange_nG0                                                                                                                                                                                                 
grid_x, grid_y = np.mgrid[ranges[0][0]:ranges[0][1]:complex(0,res[0]),
                          ranges[1][0]:ranges[1][1]:complex(0,res[1])]
nPts = np.product(grid_x.shape)
xNew = grid_x.reshape(nPts)
yNew = grid_y.reshape(nPts)
zNew = xNew.copy()

zNew[:]  = -30.0
dataNew  = np.array( [xNew, yNew, zNew] ).T
f2       = arxv.construct3DInterpolationFunction(quantity = ['therm', 'heating'], slabIdx  = 0, log10 = True)
zNew[:]  = np.log10(surfHeatinFac * (10.0**(f2(dataNew))) )
print zNew
"""

print 'done'