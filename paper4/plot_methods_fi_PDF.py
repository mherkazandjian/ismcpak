import time

import matplotlib
matplotlib.use('Qt4Agg')

import numpy

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.datamodel import Particles
from amuse.units import units

from galaxies import fi_utils
from paper4 import paper_plots
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

npp           = 100             # number of particles to be sampled from each SPH particle
n_min_sample  = 1e2            # particles with densities greater than this are sampled
fit_func_rng  = [1e-2, 1e+3]   # the range of densities used in constructing the function used for the sampling
snap_index    = 4
save_figs     = True
fig_paths     = {
                 'fig1' : '/home/mher/ism/docs/paper04/src/figs/methods/PDF_T_n.eps',
                 'fig2' : '/home/mher/ism/docs/paper04/src/figs/methods/PDF_n_fit.eps',
                 }
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

## extending the densities to include higher ones
n_s, w_s, gas_gt, w_gt, gas_lt = gas.sample_higher_densities(npp = npp, 
                                                             n_min_sample = n_min_sample,
                                                             fit_func_rng=fit_func_rng)

## making the new particle set
gas_s = Particles(len(n_s))

## attributes to be assigined to the sampled particles
attr_list = ('Av', 'G0', 'Pe', 'T', 'gmech', 'mass', 'n', 'vdisp', 'vx', 'vy', 'vz', 'x', 'y', 'z',)

## creating a new particle set by repeating the attributes of the parent
## particles npp times (the new densities are assigned later)
for attr in attr_list:
    
    ## getting the data of the attribute of the gas particles which were sampled
    attr_gas_gt = getattr(gas_gt, attr)
    
    ## repeating the attributes npp times
    attr_gas_s = numpy.repeat(attr_gas_gt, npp)
    
    ## setting the repeated attributes to the sampled gas particle set
    setattr(gas_s, attr, attr_gas_s)

## setting the weight attribute to gas_s, gas_gt, gas_lt
gas_lt.weights = 1.0
gas_gt.weights = w_gt
gas_s.weights = w_s   

## assigning the sampled densities to the sampeld gas particle set
gas_s.n = n_s

## assigning the keys to the sampled particles to their parent's keys
gas_lt.children = numpy.ones(len(gas_lt))*numpy.nan
gas_lt.parent   = numpy.ones(len(gas_lt))*numpy.nan

gas_gt.children = numpy.arange(len(gas_gt))
gas_gt.parent   = numpy.ones(len(gas_gt))*numpy.nan

gas_s.children = numpy.ones(len(gas_s))*numpy.nan
gas_s.parent   = numpy.repeat(gas_gt.children, npp) 

## aggregating all the particles sets into one set and writing it to the disk
gas_all = fi_utils.gas_set(len(gas_lt) + len(gas_gt) + len(gas_s))

# adding the new attributes related to the sampled particles 
new_attr_list = attr_list + ('weights', 'children', 'parent')

## setting the new attributes as zero arrays to the new particle set (the new densities are assigned later)
for attr in new_attr_list:
    
    ## getting the data of the attribute of the gas particles which were sampled
    attr_data = numpy.hstack((getattr(gas_lt, attr),    
                              getattr(gas_gt, attr), 
                              getattr(gas_s, attr)
                              )
                            )
    
    ## setting the repeated attributes to the sampled gas particle set
    setattr(gas_all, attr, attr_data)


####################################################################################### 
##########################making the plots############################################# 
####################################################################################### 

## making the plot for the original distribution
#paper_plots.plot_methods_fig(gas_all, save_figs, fig_paths) 

print 'done'