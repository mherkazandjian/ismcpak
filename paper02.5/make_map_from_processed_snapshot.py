#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
from galaxies import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO', '13CO'],
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(4, 4 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 20000000.0],
                             'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                      },
          'all_maps' : {
                  'map1'   : {
                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                              'as_log10': True,
                              'func'    : numpy.mean,
                             },
                  'map2'   : {
                              'attr'    : 'em_fluxKkms_13CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'$f(L_{13CO(1-0} K.km.s-1))$', 
                              'as_log10': True,
                              'func'    : numpy.mean,
                             },
                        },
        'save_maps' : True,
        'save_image': False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

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

    #getting all the map
    for i, this_map in enumerate(params['all_maps']):
        
        this_map_info = params['all_maps'][this_map]
        
        map_data = fi_utils.make_map(gas, hist, **this_map_info)
        
        #saving the produced map file into the analysis dir
        filename = os.path.join(params['rundir'],'analysis', 'fiout.%06d.%s.npz' % (snap_index, this_map_info['attr']))
        
        if params['save_maps']:
            numpy.savez_compressed(filename, params=params, map_data=map_data)
            print 'saved map to file:\n\t\t\t%s' % filename
        
        fi_utils.plot_map(map_data, params, this_map_info, 
                          fi_utils.get_snapshot_time(snap_index, params), 
                          params,
                          snap_filename
                          )
    #
#
###############################################################################################################
#################################DWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARF#######################
#################################DWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARFDWARF#######################
###############################################################################################################
for snap in params['snaps']:    
    generate_maps(snap, params)

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO', '13CO'],
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(20, 20 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 20000000.0],
                             'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-2.5, 2.5] | units.kpc, #kpc
                      },
          'all_maps' : {
                  'map1'   : {
                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                              'as_log10': True,
                              'func'    : numpy.mean,
                             },
                  'map2'   : {
                              'attr'    : 'em_fluxKkms_13CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'$f(L_{13CO(1-0} K.km.s-1))$', 
                              'as_log10': True,
                              'func'    : numpy.mean,
                             },
                        },
        'save_maps' : True,
        'save_image': False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

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

    #getting all the map
    for i, this_map in enumerate(params['all_maps']):
        
        this_map_info = params['all_maps'][this_map]
        
        map_data = fi_utils.make_map(gas, hist, **this_map_info)
        
        #saving the produced map file into the analysis dir
        filename = os.path.join(params['rundir'],'analysis', 'fiout.%06d.%s.npz' % (snap_index, this_map_info['attr']))
        
        if params['save_maps']:
            numpy.savez_compressed(filename, params=params, map_data=map_data)
            print 'saved map to file:\n\t\t\t%s' % filename
        
        fi_utils.plot_map(map_data, params, this_map_info, 
                          fi_utils.get_snapshot_time(snap_index, params), 
                          params,
                          snap_filename
                          )
    #
#

for snap in params['snaps']:    
    generate_maps(snap, params)

pylab.show()

pylab.show()
