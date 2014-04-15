#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
from pylab import *

from amuse.units import units, nbody_system
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
from galaxies import fi_utils
import ismUtils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100',  # the path of the dir containing the simulation
          'imres' : 200,                                                 # resolution of the maps to be produced imres x imres
          'species' : [],  # ['CO']
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snap'    : 4,
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [0.01, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-10.0, 10.0] | units.kpc, #kpc
                      },

          'all_maps' : {  
                        'map_CO1-0': {
                                     'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                                     'v_rng'   : [-10.0, 4.0],
                                     'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                                     'as_log10': True,
                                     'func'    : numpy.average,
                                     'weights' : 'em_fluxKkms_CO1-0',
                                     },
                        'map_NH2'  : {
                                     'attr'    : 'pdr_NH2',
                                     'v_rng'   : [10.0, 30.0],
                                     'title'   : r'$N(H2)$', 
                                     'as_log10': True,
                                     'func'    : numpy.average,
                                     'weights' : 'em_fluxKkms_CO1-0',
                                     },
                        },
        'save_maps' : False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

snap_index = params['snap']
    
#path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snap_index + '.states.npz'  

#loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'], load_pdr=params['pdr_sph'])
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

#keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

#making the 2D histogram
print 'getting the spatial distrubutions'
hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
hist.info()
print 'done getting the spatial distributuions'

## keeping the gas particles within the region of interest
gas = gas[hist.inds_in]

### optional ###
## keeping only the original particles and discarding the sampled ones
#gas = gas[gas.get_inds_original_set()]
gas.set_radii(weighting='original-only', rundir=params['rundir'], snap_index=params['snap'])

n, vdisp, Pe, Av = gas.n, gas.vdisp, gas.Pe, gas.Av
r_original = gas.radius.copy()

## plotting stuff
subplot(221)
loglog(n[::100], vdisp[::100], '.')
xlabel('n')
ylabel('vdisp')

subplot(222)

loglog(n[::100], r_original[::100]*1000.0, 'b.', label='r_smooth', zorder=2)

#gas.set_radii_from_Av_and_using_original_weighting()
gas.set_radii_from_Av_by_conserving_total_area(rundir=params['rundir'], snap_index=params['snap'])
r_from_Av = gas.radius*1000.0
loglog(n[::100], r_from_Av[::100], 'r.', label='r_Av', zorder=1)

xlabel('n')
ylabel('length')
legend()


show()
