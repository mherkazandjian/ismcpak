# in this test file, we load a database of PDR meshes and extract information
# at conditions determined by the SPH simulations
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate
import numpy
import pylab

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants
from mylib.utils.misc  import xselect

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

home = '/home/mher'
#path to fi snapshot
filename = home + "/ism/runs/galaxies/test1/test/test.000000"
#PDR database file
pdrDatabaseDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-1.0/'
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/gMech_sph.eps'
#-----------------------------------------------------------------------------

#loading the sph simulation data
gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

#plotting the temperature vs the gas density
meanmwt=1.3|units.amu
n_gas_cgs = (gas.rho/meanmwt).value_in( units.cm **-3 )
gasT = gas.temperat.value_in(units.K)
G0 = 6.54*gas.fuvheat.value_in( units.none )
g_mech = gas.dethdt.as_quantity_in( units.erg / (units.g * units.s ) )
# converting the mechanical heating rate from per unit mass (of H gas)     
# to per unit volume (see notesISM.odt)
g_mech = g_mech.value_in( g_mech.unit ) # gMech in erg / (g s )
g_mech = g_mech * 1.6474e-24 * n_gas_cgs # gMech in erg / (cm^3 s) 

#-------------------------------------------------
width  = 3.0    #figure width (non normalized) 
height = 3.4  #figure heigh (non normalized)
as_rat = width/height #aspect ratio of the figure

ax_xs  = 0.25 #axses x start (normalized)
ax_ys  = 0.18 #axses y start (normalized)
ax_sz  = 0.73 #axses size (normalized)

fig    = pylab.figure(0, figsize = (width, height) )
fig.set_facecolor('white')
ax1    = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])

#plotting the gmech vs n_gas of the sph particles
x = numpy.log10(n_gas_cgs[0::100])
y = numpy.log10(g_mech[0::100])

ax1.plot(x, y,'.', markersize=1)

#plotting the boundary curve
indsBoundaryCurve = [19735, 19662,  3511,  4001, 11167, 13883,   157, 10948,  6004, 8977, 19920, 
                     66, 13913,  6037,  8996,  6019, 2,    31, 3896,   109,  8995, 19935,  6223,  
                     8987,    33, 19759,  4101, 34,  5991,  1408,  3757,  6339, 10691, 10645,] 
#ax1.plot(x[indsBoundaryCurve], y[indsBoundaryCurve], 'ko')
                     
                               
ax1.set_xlim([0, 6])
ax1.set_ylim([-26, -21])
ax1.set_ylabel(r'$\Gamma$ [erg cm$^{-3}$ s$^{-1}$]', size = 'large')
ax1.set_xlabel(r'n [cm$^{-3}$]', size = 'large')

fig.savefig(imageSavePath)

pylab.show()
