import time, os

import matplotlib
matplotlib.use('Qt4Agg')
from pylab import *

from numpy import repeat, nan, arange, where, hstack, ones

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units

from galaxies import fi_utils

fig, axs = subplots(2, 3, figsize=(12,8))

#===========================================================================================================
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params =\
 {
   'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',    # the path of the dir containing the simulation
   'snap_index' : 4,
   'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
   'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
               'sph':{
                      'min_log_n_use'  : -3.0,      
                      'min_log_G0_use' : -3.0,
                      'min_log_gm_use' : -50.0,
                     },
               #the size of the box to be displayed (particles outside the range are discarded)
               'box_size' : [-20, 20] | units.kpc,
               },
 }
#===========================================================================================================
 

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

# setting the filename
suffix = '%06d' % params['snap_index']
snapName = 'fiout.%s' % suffix 
filename = params['rundir'] + '/firun/' + snapName 
    
#loading the sph simulation data
print 'loading snapshot %s : ' % filename
gas_fi, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

print 'gas fraction = ', gas_fi.mass.sum()/(gas_fi.mass.sum() + dark.mass.sum())

# getting the gas particles in cgs units
gas = fi_utils.convert_units_to_pdr_units(gas_fi, metallicity)

print 'done reading fi snapshot : %s' % filename
print 'number of sph particles in snapshot = %d' %  len(gas)

gas.plot_PDF_CDF(fit_func_rng=[1e-2, 1e4], in_ax=axs, color='r')

#===========================
# saving the attributes of the gas particle to an ascii file
import numpy

data = numpy.array(
                   [
                    gas.x,
                    gas.y,
                    gas.z,
                    gas.radius,
                    gas.vdisp,
                    gas.mass,
                    gas.n,
                    gas.G0,
                    gas.gmech,
                    gas.id,
                    gas.n,
#                    gas.Pe,
                    gas.T,
                    gas.vx,
                    gas.vy,
                    gas.vz,
                   ]
                  )

numpy.savetxt('/home/mher/snapshot.out', data)

#===========================
show()
print 'done'