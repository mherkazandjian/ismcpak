#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units, nbody_system
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 25,                                                 # resolution of the maps to be produced imres x imres
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',      # the path to the dir containing the PDR database
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2/',          # the path to the dir containing the PDR database          
          'species' : ['CO'],
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(20, 20 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [3.0, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-2.0, 2.0] | units.kpc, #kpc
                      },

          'all_maps' : {  
                        'map_CO1-0': {
                                     'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                                     'v_rng'   : [-10.0, 4.0],
                                     'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                                     'as_log10': True,
                                     'func'    : numpy.sum,
                                     },
                        'map_NH2'  : {
                                     'attr'    : 'pdr_NH2',
                                     'v_rng'   : [10.0, 30.0],
                                     'title'   : r'$N(H2)$', 
                                     'as_log10': True,
                                     'func'    : numpy.sum,
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

def generate_maps(snap_index, params):
    
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

    #getting the maps for CO luminosity and H2 column density
    map_CO_info = params['all_maps']['map_CO1-0']
    map_CO_data = fi_utils.make_map(gas, hist, **map_CO_info)

    map_NH2_info = params['all_maps']['map_NH2']
    map_NH2_data = fi_utils.make_map(gas, hist, **map_NH2_info)
    
    map_X_factor = numpy.log10((10.0**map_CO_data) / 10.0**map_NH2_data)

    map_X_factor_plot = {'v_rng'   : [-30.0, -10.0],
                         'title'   : 'X_factor',   
                        } 
    fi_utils.plot_map(map_X_factor, params, map_X_factor_plot,
                      fi_utils.get_snapshot_time(snap_index, params), 
                      params,
                      snap_filename
                     )

    pylab.figure()
    pylab.plot(map_CO_data.flatten(), map_NH2_data.flatten(), '.')
    #
#

for snap in params['snaps']:    
    generate_maps(snap, params)

pylab.show()