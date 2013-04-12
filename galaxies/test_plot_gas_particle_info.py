from matplotlib import pylab as pyl
import numpy as np

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants

import meshUtils
 
home = '/home/mher'

fnameSave = home + '/ism/runs/galaxies/test1/test/xyz.npz'


# reading the full snapshot
if True:
    filename = home + "/ism/runs/galaxies/test1/test/test.000000"
    gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)
 
    #print gas     
    gMech   =  gas.dethdt
    fuvheat =  gas.fuvheat
    nGas    =  gas.rho
    
    print gMech[0], fuvheat[0], nGas[0]

arxv = meshUtils.meshArxv()    
x,y,z = arxv.grid_interpolator(gas.rho, gas.fuvheat, gas.dethdt)

pyl.plot(x, z,'.')
pyl.xlabel('$log_{10} n_{gas} ( cm^{-3} )$')
pyl.ylabel('$log_{10} \Gamma_{mech} ( erg.cm^{-3}.s^{-1) }$')
pyl.xlim( 0, 4)
pyl.ylim( -30, -21)
pyl.show()

#pyl.plot(numpy.log10(x), numpy.log10(y),'.')
#pyl.show()



#extracting the temperatures
tKin = gas[:].temperat.value_in(units.K)

#extracting the masses of the gas particles
m = gas[:].mass.value_in(units.MSun)

#extracting the positions in pc
px, py, pz = (gas[:].x.value_in(units.parsec), gas[:].y.value_in(units.parsec), gas[:].z.value_in(units.parsec))

#extracting the h2 fraction
h2frac = gas[:].h2frac.value_in(units.none)
"""

#saving some of the snapshot information
"""
numpy.savez(fnameSave, lx = lx, ly = ly, lz = lz,
                       px = px, py = py, pz = pz,
                       m = m,
                       tKin = tKin,
                       h2frac = h2frac
            )
"""
print 'done'