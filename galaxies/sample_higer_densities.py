import time, os

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from numpy import repeat, nan, arange, where, hstack, ones

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units

from galaxies import fi_utils
#===========================================================================================================
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',    # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100', # the path of the dir containing the simulation
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

npp                     = 20              # number of particles to be sampled from each SPH particle
n_min_sample            = 100.0            # particles with densities greater than this are sampled
n_max_sample            = 1e6
fit_func_rng            = [1e-2, 1e+4]   # the range of densities used in constructing the function used for the sampling
save_sampled            = False

nbins_in_match_interval = 25
matching_tol            = 0.001 
nPasses                 = 1            
save_weights_func       = False
snap_index              = 4
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
                                                             n_max_sample = n_max_sample,
                                                             fit_func_rng=fit_func_rng,
                                                             plot=True)

## making the new particle set
gas_s = fi_utils.gas_set(len(n_s))
## attributes to be assigined to the sampled particles
attr_list = ('Av', 'G0', 'Pe', 'T', 'gmech', 'mass', 'n', 'vdisp', 
             'vx', 'vy', 'vz', 'x', 'y', 'z', 'id', 'radius')

## creating a new particle set by repeating the attributes of the parent
## particles npp times (the new densities are assigned later)
for attr in attr_list:
    
    ## getting the data of the attribute of the gas particles which were sampled
    attr_gas_gt = getattr(gas_gt, attr)
    
    ## repeating the attributes npp times
    attr_gas_s = repeat(attr_gas_gt, npp)
    
    ## setting the repeated attributes to the sampled gas particle set
    setattr(gas_s, attr, attr_gas_s)

## setting the weight attribute to gas_s, gas_gt, gas_lt
gas_lt.weights = 1.0
gas_gt.weights = w_gt
gas_s.weights = w_s   

## assigning the sampled densities to the sampeld gas particle set
gas_s.n = n_s

## assigning the keys to the sampled particles to their parent's keys
gas_lt.children = ones(len(gas_lt))*nan
gas_lt.parent   = ones(len(gas_lt))*nan

gas_gt.children = arange(len(gas_gt))
gas_gt.parent   = ones(len(gas_gt))*nan

gas_s.children = ones(len(gas_s))*nan
gas_s.parent   = repeat(gas_gt.children, npp) 

## setting the radii and the masses of the sampled particles (children and parents)
gas_gt.mass = gas_gt.mass * gas_gt.weights
gas_s.mass = gas_s.mass * gas_s.weights

gas_gt.radius = gas_gt.radius * (gas_gt.weights)**(1.0/3.0) 
gas_s.radius = gas_s.radius * (gas_s.weights * ((repeat(gas_gt.n, npp))/gas_s.n) )**(1.0/3.0) 

## checking if the parent/child assignment and weighting is done correctly
## might take lots of time
if False:
    t0 = time.time()
    for i, child in enumerate(gas_gt.children):
    
        inds = where(gas_s.parent == child)
        
        sum_weights = (gas_gt.weights[i] + gas_s.weights[inds].sum())
        if sum_weights != 1.0:
            print 'i = %d  child = %d, 1.0 - sum_weights = %e' % (i, child, 1.0 - sum_weights)
            raise ValueError('check order of the sampled particles')
    
        print i
    print 'time checking (seconds) = ', time.time() - t0


## aggregating all the particles sets into one set and writing it to the disk
gas_all = fi_utils.gas_set(len(gas_lt) + len(gas_gt) + len(gas_s))

# adding the new attributes related to the sampled particles 
new_attr_list = attr_list + ('weights', 'children', 'parent')

## setting the new attributes as zero arrays to the new particle set (the new densities are assigned later)
for attr in new_attr_list:
    
    ## getting the data of the attribute of the gas particles which were sampled
    attr_data = hstack((getattr(gas_lt, attr),
                        getattr(gas_gt, attr),  
                        getattr(gas_s, attr)
                        )
                       ) 
    
    ## setting the repeated attributes to the sampled gas particle set
    setattr(gas_all, attr, attr_data)


if save_sampled == True:
    
    print 'number of particles in each set:'
    print 'gas_lt : ', len(gas_lt)
    print 'gas_gt : ', len(gas_gt)
    print 'gas_s  : ', len(gas_s)
    print '------------------------------'
    print 'total : ',  len(gas_lt) + len(gas_gt) + len(gas_s)
      
    fname_ext = filename + '.ext'
    fi_utils.save_gas_particle_info(fname_ext, gas_all, new_attr_list)
    print 'saved snapshot with sampled data to \n\t %s' % fname_ext 

if save_weights_func == True:

    save_weights_func_path =  params['rundir'] + '/firun/' + 'weights_func.%06d' % snap_index
    
    w, inters, w_func = gas_all.match_weights(npp, fit_func_rng, n_min_sample, 
                                              [n_min_sample, gas_all.n.max()], 
                                              nbins_in_match_interval, matching_tol, 
                                              nPasses=nPasses,
                                              save_weights_info=save_weights_func_path)


pylab.show()
print 'done'