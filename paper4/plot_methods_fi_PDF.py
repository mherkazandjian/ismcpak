import matplotlib
matplotlib.use('Qt4Agg')

import numpy

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units

from galaxies import fi_utils
from paper4 import paper_plots
#===========================================================================================================
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',        # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100', # the path of the dir containing the simulation
          'use_sampled_set' : True,
          'snap_index'  : 4,
          'imres'       : 100,                                                  # resolution of the maps to be produced imres x imres
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

save_figs     = True
fig_paths     = {
                 'fig1' : '/home/mher/ism/docs/paper04/src/figs/methods/PDF_T_n.eps',
                 'fig2' : '/home/mher/ism/docs/paper04/src/figs/methods/PDF_n_fit.eps',
                 }
#===========================================================================================================
 
## extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

## setting the filename
suffix = '%06d' % params['snap_index']
snapName = 'fiout.%s' % suffix 
filename = params['rundir'] + '/firun/' + snapName 

## appending the extention if the parameters instruct using a sampled snapshot
if 'use_sampled_set' in params and params['use_sampled_set'] == True:
    filename += '.ext.npz'

## loading the sph simulation data
print 'loading snapshot %s : ' % filename
if 'use_sampled_set' in params and params['use_sampled_set'] == True:
    
    gas, attr_read = fi_utils.load_gas_particle_info(filename)
    
else:
    gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

    ## getting a new particle set for the gas particles with the attributes we will be using later
    ## and converting the units to the same as those of the PDR models    
    gas = fi_utils.convert_units_to_pdr_units(gas, metallicity)


print 'done reading fi snapshot : %s' % filename
print 'number of sph particles in snapshot = %d' %  len(gas)

####################################################################################### 
##########################making the plots############################################# 
####################################################################################### 

## making the plot for the original distribution
paper_plots.plot_methods_fig(gas, save_figs, fig_paths, params) 

print 'done'