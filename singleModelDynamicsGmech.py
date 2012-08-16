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


"""
#arxv.saveGridsToFiles(gridsRes, lgammaMechSec, radexParms)
pyl.show()
"""
print 'done'
