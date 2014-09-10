import time

import matplotlib
matplotlib.use('Qt4Agg')

from numpy import *

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.datamodel import Particles
from amuse.units import units

from galaxies import fi_utils
#===========================================================================================================
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',    # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext', # the path of the dir containing the simulation
          'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                            },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-10, 10] | units.kpc,
                      },
          }

npp           = 5             # number of particles to be sampled from each SPH particle
n_min_sample  = 1e2            # particles with densities greater than this are sampled
fit_func_rng  = [1e-2, 1e+3]   # the range of densities used in constructing the function used for the sampling
save_sampled  = False  
snap_index    = 4
#===========================================================================================================
 
#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

# setting the filename
suffix = '%06d' % snap_index
snapName = 'fiout.%s' % suffix 
filename = params['rundir'] + '/firun/' + snapName 
    
#loading the sph simulation data
print 'loading snapshot %s : ' % filename
gas_fi, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)


# getting the gas particles in cgs units
gas = fi_utils.convert_units_to_pdr_units(gas_fi, metallicity)

print 'done reading fi snapshot : %s' % filename
print 'number of sph particles in snapshot = %d' %  len(gas)

###########################################################################################################

R = sqrt(gas.x*gas.x + gas.y*gas.y)

# selecting the original gas particles
gas_selection_amuse = gas[where( (R > 1.0)*(R < 2.0))]
gas_selection = fi_utils.gas_set(len(gas_selection_amuse)) 
gas_selection.copy_attr_from_amuse_set(gas_selection_amuse, ['n'])


## extending the densities to include higher ones
n_s, w_s, gas_gt, w_gt, gas_lt = gas_selection.sample_higher_densities(npp = npp, 
                                                                       n_min_sample = n_min_sample,
                                                                       fit_func_rng=fit_func_rng)



print 'done'