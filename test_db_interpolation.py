#test interpolation of data of meshes from the DB within the DB bounds and outside  
#-----------------------------------------------------------------------------------------
import time
import sys
import os

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate, spatial
import numpy
import pylab

import meshUtils
from mylib.utils.misc  import xselect
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

home = '/home/mher'
pdrDatabaseDirPath =  home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/'
#-----------------------------------------------------------------------------
#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = pdrDatabaseDirPath, readDb = True, set_defaults = True)

fLinear = arxvPDR.construct3DInterpolationFunction(quantity = ['state', 'gasT'], slabIdx  = 0, log10 = True, interpolator='linear')

gMechNew  = numpy.linspace(-50.0, -10.0, 100.0)
nNew  = numpy.ones(gMechNew.shape)*(2.3)
G0New = numpy.ones(gMechNew.shape)*(1.8)

dataNew = numpy.array([nNew, G0New, gMechNew]).T
values = arxvPDR.getQuantityFromAllMeshes(['state', 'gasT'], slabIdx = 0)
values = numpy.log10(values)

ti = time.time()
tLinear = fLinear(dataNew)
tf = time.time()

print tLinear

data = numpy.array([arxvPDR.grid_x, arxvPDR.grid_y, arxvPDR.grid_z]).T  #3D coordinates
tree=spatial.cKDTree(data, leafsize=1)

indsNan = numpy.isnan(tLinear).nonzero()[0]
ptsNan = dataNew[indsNan,:]
r = tree.query(ptsNan)

pylab.figure()
pylab.plot(gMechNew, tLinear)
pylab.plot(gMechNew, tLinear,'x')
tLinear[indsNan] = values[r[1]]
pylab.plot(gMechNew[indsNan], tLinear[indsNan],'ko')
pylab.ylim([-1, 5.0])

for ind in r[1]:
    print data[ind,:], values[ind]
    
pylab.show()
#tNew = np.reshape(tNew, grid_x.shape)
#
#return tNew, grid_x, grid_y

"""
#getting an interpolated quantity from the PDR meshes for the states of the SPH gas
#particles.
dataInterp = arxvPDR.grid_interpolator(gas.rho, gas.fuvheat, gas.dethdt,
                                       quantity = ['state', 'gasT'], slabIdx = 0,)
"""

#plotting the gmech vs n_gas of the sph particles
"""
pylab.loglog(n_gas_cgs[0::100], g_mech[0::100],'.')
pylab.xlim([1, 1e6])
pylab.ylim([1e-26, 1e-16])
pylab.ylabel(r'$\Gamma$ [erg cm$^{-3}$ s$^{-1}]')
pylab.xlabel(r'n [cm$^{-3}$]')
"""
#info = xselect(x = numpy.log10(n_gas_cgs), y = numpy.log10(gasT) )
"""
n_gas_pdr = 10.0**arxvPDR.grid_x
gasT_pdr = arxvPDR.getQuantityFromAllMeshes(['state', 'gasT'], slabIdx = 0)

pylab.loglog(n_gas_cgs, gasT, 'r.')
pylab.loglog(n_gas_pdr, gasT_pdr, 'b.')
pylab.show()
"""

print 'done'


