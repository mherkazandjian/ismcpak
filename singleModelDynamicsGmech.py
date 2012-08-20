from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *
#---------------------------------------------------------------------------

runDirPath    = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/'

gridsRes      = 10
lgammaMechSec = -30.0
metallicity   = 1.0
plotRangenG0  = [[0,6],[0,6]]
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

# getting the basic species defined in baseSpecies.py
import baseSpecies
baseSpecs = baseSpecies.baseSpecies()

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
net = chemicalNetwork(rxnFile, baseSpecs, UMISTVER = 'umist99')
# reading the species to be removed from a file
net.removeSpecies( underAbunFile = underAbunFile )
net.removeSpecies( species = removeManual )
# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = specNumFile)
arxv.setChemicalNetwork(net) # assiginig the chemical network to the archive
#-------------------------------------------------------------------
# plotting stuff

#pyl.ioff()
#arxv.plotGrid(gridsRes, lgammaMechSec, radexParms)

#tGas = arxv.computeSurfaceTemperatureGrid(res = 20, ranges = plotRanges)
#tGas.imshow()

fig   = pyl.figure(0, (6,6))
imAxs = pyl.axes([0.2, 0.2, 0.5, 0.5])
proc  = arxv.computeHeatingCoolingGrid(slabIndex = 0, meshInds = None, whichThermal = 'heating', whichProcess = 'photo', res = 10, ranges = plotRangenG0)
im    = proc.imshow()
tickV =  [-26, -24, -22, -20, -18, -16, -14]
cntr  = proc.plotContour( levels = tickV, colors = 'black' )
pyl.clabel(cntr, inline=1, fontsize=10)
pyl.title('heating|cooling process')
cbarAxs = pyl.axes([0.1, 0.8, 0.8, 0.1 ])
pyl.title('log10 of process')
cbar11 = pyl.colorbar(im, cax=cbarAxs, ax=pyl.gca(), orientation = 'horizontal')
pyl.show()

##########################################################################
##########################################################################
################TESTIN INTERPOLATION# TESTIN INTERPOLATION################
################TESTIN INTERPOLATION# TESTIN INTERPOLATION################
################TESTIN INTERPOLATION# TESTIN INTERPOLATION################
##########################################################################
##########################################################################
from scipy import interpolate
import numpy as np
import time

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

nPts = np.prod(proc.shape)
x = proc.cntrd[0].reshape(nPts)
y = proc.cntrd[1].reshape(nPts)
data = np.array([x,y]).T # the shape should be (nPts, nDim), thats why i take the transpose                                                                                                                                                
z = proc.reshape(nPts)

ti = time.time()
f = interpolate.LinearNDInterpolator(data, z) # getting the interpolation function                                                                                                                                                           
tf = time.time()
print 'constructed he interpolation function from %d points in %f seconds' % (nPts, tf-ti)

ti = time.time()
fNear = interpolate.NearestNDInterpolator(data, z) # getting the interpolation function
tf = time.time()
print 'constructed he interpolation (nearest neighbour) function from %d points in %f seconds' % (nPts, tf-ti)

nPts    = 100000 # number of points from the interpolation function will be constructed                                                                                                                                                 
nInterp = 100000   # number of points to interpolate over                                                                                                                                                                                      

# generating new points where interpolation will be done for new values                                                                                                                                                                      
# (new reigon is half the size inside the old reigon to avoid nan's for this simple example)                                                                                                                                                 
xNew = (np.random.rand(nInterp)*(6.0 - 0.0) + 0.0)
yNew = (np.random.rand(nInterp)*(6.0 - 0.0) + 0.0)
dataNew = np.array([xNew,yNew]).T # the shape should be (nInterp, nDim), thats why i take the transpose                                                                                                                                    

ti = time.time()
zNew = f(dataNew)
tf = time.time()
print 'interpolated %d points in %f seconds at a rate of %e pts/sec' % (nInterp, tf-ti, nInterp / (tf-ti))

ti = time.time()
zNew1 = fNear(dataNew)
tf = time.time()
print 'interpolated (nearest neighbour) %d points in %f seconds at a rate of %e pts/sec' % (nInterp, tf-ti, nInterp / (tf-ti))


# defining the points in the uniform 2D grid                                                                                                                                                                                                 
grid_x, grid_y = np.mgrid[0.0:6.0:100j, 0.0:6.0:100j]

# interpolating with different methods                                                                                                                                                                                                       
from scipy.interpolate import griddata
grid_z0 = griddata(dataNew, zNew, (grid_x, grid_y), method='linear')

import matplotlib.pyplot as plt
plt.subplot(111)
plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
plt.title('Linear')
plt.show()

print 'done'
