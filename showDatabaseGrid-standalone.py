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

#---------------------------Archive parameters-----------------------
# database to analyze
runDirPath    = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/'
#runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweepNew-1and2/'
#runDirPath    = '/home/mher/ism/runs/oneSided/surfaceGridHighRes-z-1.0/'

# reference database
runDirPath2  = '/home/mher/ism/runs/oneSided/uniformSweepNew-1and2/' 
minGmech     = -30.0
metallicity  = 1.0

#quantity       = ['state', 'gasT']
quantity       = ['therm', 'heating']
#quantity       = ['fineStructureCoolingComponents', 'C+', 'rate', '1-0']
plotRange_nG0  = [[0,6],[0,6]]
slabIdx        = -1
res            = [100, 100]
lgammaMechSec  = -30.0     ###;;; check with minGmech

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

# reading the reference archive
print 'setting up the reference archive'
t0 = time()
arxvRef = meshArxv( metallicity = metallicity )
arxvRef.readDb( runDirPath2 )
arxvRef.checkIntegrity()
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

"""
arxv.showGrid(quantity = quantity,
              slabIdx  = slabIdx,
              ranges   = plotRange_nG0,
              res      = res,
              zSec     = lgammaMechSec,
              log10    = log10)

sys.exit()
"""

x = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'nGas']) )
y = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'G0']) )
z = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'gammaMech']) )
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x,y,z, 'o')
ax.set_xlim(-1, 7)
ax.set_ylim(-1, 7)
ax.set_zlim(-50, -10)
ax.set_xlabel('nGas')
ax.set_ylabel('G0')
ax.set_zlabel('gMech')
pyl.show()
print 'done'
sys.exit()

# getting the surface heating of the meshes in the database from the reference database
x = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'nGas']) )
y = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'G0']) )
gMechZero = x.copy()
gMechZero[:] = minGmech # lowest mechanical energy used (in log)
f         = arxvRef.construct3DInterpolationFunction(quantity = ['therm', 'heating'], slabIdx  = 0, log10 = True)
dataNew   = np.array( [y, x, gMechZero] ).T #### swap the x with y
gammaSurf = f(dataNew)

z = np.log10(arxv.getQuantityFromAllMeshes( ['hdr', 'gammaMech']) )

print 'x      = ', x[0:20]
print 'y      = ', y[0:20]
print 'Gmech0 = ', gMechZero[0:20]
print 'GSurf  = ', gammaSurf[0:20]
print 'z      = ', z[0:5]

r = 10.0**z / 10.0**gammaSurf
print 'ratio  = ', r

values = arxv.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx)

if log10 != None and log10 == True:
    values[:] = np.log10(values[:])

data = np.array([x, y, r]).T  #3D

ti = time()
f = interpolate.LinearNDInterpolator(data, values) # getting the interpolation function     
tf = time()
print 'constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti)

arxv.showGrid(quantity = quantity,
              slabIdx  = slabIdx,
              ranges   = plotRange_nG0,
              res      = res,
              zSec     = 0.01,
              log10    = log10,
              fInterp  = f)
print 'done'