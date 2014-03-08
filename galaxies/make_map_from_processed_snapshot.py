#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
from numpy import arange
import pylab

from amuse.units import units
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100',  # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-test',  # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO'],#'CO'], #'13CO', 'HCN', 'HNC', 'HCO+'],
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
#          'weights' : 'original-only', #'by-number', #'matched',  #'original-only' ,#None ,#by-number          
          'weights' : 'by-number', #'matched',  #'original-only' ,#None ,#by-number          
#          'weights' : #'matched',  #'original-only' ,#None ,#by-number          
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
          'check'   : 'default',
#          'maps'   : {CO
#                      'attr' : 'n', #'mass', 'G0', 'gmech', 'Av'
#                      'v_rng': [-3.0, 4.0],
#                      'title': r'$f(\bar{n})$', 
#                      'log10': True
#                     },
#          'maps'   : {
#                      'attr' : 'em_fluxcgs_CO1-0', #'mass', 'G0', 'gmech', 'Av'
#                      'v_rng': [-10.0, -2.0],
#                      'title': r'$f(L_{CO(1-0})$', 
#                      'log10': True                      
#                     }, 
          'all_maps' : {


                        
#                  'map1'   : {
#                              'attr'    : 'em_fluxKkms_CO2-1', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
##                              'title'   : r'asadasdasd', #$f(L_{CO(10-9} K.km.s-1))$', 
#                              'as_log10': True,
##                              'func'    : numpy.mean,
#                              'func'    : numpy.average,
#                              'weights' : 'weights',
#                             },

                  'map1'   : {
                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'', 
                              'as_log10': True,
#                              'func'    : fi_utils.luminosity,
                              'func'    : fi_utils.mean_flux,
                             },


#                  'map2'   : {
#                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'weights',
#                             },
#                  'map1'   : {
#                              'attr'    : 'em_fluxKkms_HCN1-0', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
#                              'as_log10': True,
#                              #'func'    : numpy.mean,
#                              'func'    : numpy.average,
#                              'weights' : 'weights',
#                             },
#                  'map2'   : {
#                              'attr'    : 'em_fluxKkms_13CO1-0', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                              'title'   : r'$f(L_{13CO(1-0} K.km.s-1))$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'em_fluxKkms_CO1-0',
#                             },
#                  'map3'   : {
#                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'em_fluxKkms_CO1-0',
#                             },
#                  'map4'   : {
#                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
#                              'v_rng'   : [-10.0, 4.0],
#                             'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'em_fluxKkms_CO1-0',
#                             },
#                  'map5'   : {
#                              'attr'    : 'pdr_NCO',
#                              'v_rng'   : [10.0, 30.0],
#                              'title'   : r'$N(CO)$', 
#                              'as_log10': True,
#                              'func'    : numpy.sum,
#                             },
#                  'map6'   : {
#                              'attr'    : 'pdr_NH2',
#                              'v_rng'   : [10.0, 25.0],
#                              'title'   : r'$N(H2)$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'em_fluxKkms_CO1-0',                              
#                             },
#                  'map7'   : {
#                              'attr'    : 'pdr_N13CO',
#                              'v_rng'   : [10.0, 30.0],
#                              'title'   : r'$N(13CO)$', 
#                              'as_log10': True,
#                              'func'    : numpy.sum,
#                             },
#                  'map8'   : {
#                              'attr'    : 'G0',
#                              'v_rng'   : [-2.0, 4.0],
#                              'title'   : r'$G_0$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                             },
#                  'map9'   : {
#                              'attr'    : 'pdr_NH',
#                              'v_rng'   : [10.0, 25.0],
#                              'title'   : r'$N(H)$', 
#                              'as_log10': True,
#                              'func'    : numpy.mean,
#                              #'func'    : numpy.average,
#                              #'weights' : 'em_fluxKkms_CO1-0',                              
#                             },
          

                        },
        'save_maps' : False,
        'save_image': False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

'''
infos = [
         ['CO'  , arange(0,10,3)], 
         ['13CO', arange(0,10,3)], 
         ['HCN' , arange(0,6,2)], 
         ['HNC' , arange(0,6,2)], 
         ['HCO+', arange(0,6,2)], 
        ]  

from lineDict import lines

for i, info in enumerate(infos):
    print info
    for j, line in enumerate(info[1]):
        print line

        line_code = info[0] + '%d-%d' % (line+1, line)
        

        map_info = {
                    'attr'    : 'em_fluxKkms_' + line_code,
                    'v_rng'   : [-10.0, 4.0],
                    'title'   : r'Flux[%s] K km s$^{-1}$' % lines[line_code]['latex'],   
                    'as_log10': True,
                    'func'    : numpy.mean,
                    #'func'    : numpy.average,
                    #'weights' : 'weights',
                    }

        params['all_maps']['map%d%d' % (i,j)] = map_info
'''

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

def generate_maps(snap_index, params):
    
    ## path to processed fi snapshot  
    snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snap_index + '.states.npz'
    
    ## loading the processed sph simulation data with the emissions 
    logger.debug('loading proccessed snapshot %s : ' % snap_filename)
    gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'], 
                                                  load_pdr=params['pdr_sph'],
                                                  )
    
    logger.debug('done reading fi snapshot : %s' % snap_filename)
    logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))
    
    ## checking for weird particles and taking care of them
    gas.check_particles(params['check'], logger)

    ## setting the weights
    weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % snap_index
    gas.use_weights(weighting=params['weights'], weights_filename = weights_filename)
    
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
    return gas
#

for snap in params['snaps']:
    gas = generate_maps(snap, params)

pylab.show()
