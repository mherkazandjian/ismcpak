from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#=====================================================================================
#=========================================PARAMETERS==================================
#=====================================================================================

home = '/home/mher'
#---------------------------Archive parameters-----------------------
# database to analyze
#runDirPath    = home + '/ism/runs/oneSided/dynamicMeshTest1/'
runDirPath    = home + '/ism/runs/oneSided/uniformSweepNew-1and2/'
#runDirPath    = home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/'

# reference database
runDirPath2  = home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/' 

#quantity       = ['state', 'gasT']
#quantity       = ['therm', 'heating']
quantity       = ['fineStructureCoolingComponents', 'C+', 'rate', '1-0']

qx             = ['hdr', 'nGas']
qy             = ['hdr', 'G0']

plotRanges  = [[0,6],[0,6],[-16,-30]]
slabIdx        = 0
metallicity    = 1.0
res            = [100, 100] # resolution of the grid
relativeGmech  = False  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                       # False => 3rd dim is gMech 
zSec           = -24  # section in the 3D dimension to be used for generating 
                        # grids. usuall this is the log10 of the mechanical heating
                        # it can be used as the ratio of mechanical heating to the
                        # surface heating(gMech = 0)
log10          = True # use the log of the quantity in interpolating and plotting

radexParms    = { 'radexPath'         : home + '/ism/code/radex/Radex/bin/radex',  
                  'molDataDirPath'    : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
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
rxnFile       = home + '/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = home + '/ism/code/ismcpak/data/species.inp'
underAbunFile = home + '/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']
#=====================================================================================
#=====================================================================================
#=====================================================================================

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


if relativeGmech:
    # interpolation function as a function of n,G0,gMech/gSurface(gmech=0)
    f = arxv.getInterpFunctionGmechToSurfaceHeating(quantity, referenceDatabasePath = runDirPath2, slabIdx = slabIdx, log10 = log10)
else:
    f = None


arxv.showGrid(quantity = quantity,
              slabIdx  = slabIdx,
              ranges   = plotRanges,
              res      = res,
              zSec     = zSec,
              log10    = log10,
              fInterp  = f)

arxv.plot_3D_grid_point(ranges = plotRanges)

pyl.show()
print 'done'