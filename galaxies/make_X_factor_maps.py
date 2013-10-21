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
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
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

    #plotting the map of the X factor 
    map_NH2_info = params['all_maps']['map_NH2']
    map_NH2_data = fi_utils.make_map(gas, hist, **map_NH2_info)
    
    map_X_factor = numpy.log10((10.0**map_NH2_data/10.0**map_CO_data))

    map_X_factor_plot = {'v_rng'   : [18.0, 22.0],
                         'title'   : 'X_factor',   
                        } 
    
    fi_utils.plot_map(map_X_factor, params, map_X_factor_plot,
                      fi_utils.get_snapshot_time(snap_index, params), 
                      params,
                      snap_filename
                     )

    #plottng the X factor relationship    
    pylab.figure()

    LCO = map_CO_data.flatten()
    NH2 = map_NH2_data.flatten()
    
    pylab.plot(LCO[::10], NH2[::10], '.')
    x = numpy.linspace(-5.0, 5.0, 100)    
    pylab.plot(x, numpy.log10(2.0e20*(1.0 - 0.3)) + x,'--r')
    pylab.plot(x, numpy.log10(2.0e20) + x,'r')
    pylab.plot(x, numpy.log10(2.0e20*(1.0 + 0.3)) + x,'--r')
    pylab.xlim(-5,5)

    ################################################################################
    #checking how much of the emission is resulting from differnet section in L(CO) 
    ################################################################################
    LCO_sections = [ [1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, numpy.Inf]]
    
    LCO10_total = numpy.sum(10.0**map_CO_data)
    
    for LCO_section in LCO_sections: 
        
        print LCO_section
        
        inds_ij = numpy.where( (map_CO_data >= LCO_section[0])*(map_CO_data < LCO_section[1]) )
        
        if inds_ij[0].size == 0:
            continue
        
        inds_i, inds_j = inds_ij 
            
        print 'percent emission of selected particles', numpy.sum(10.0**map_CO_data[inds_ij])/LCO10_total
        
        gas_den_mean = 0.0
        n_part = 0L
        
        for n, i in enumerate(inds_i):
            
            j = inds_j[n]
            
            inds_gas_in_bin = hist.get_indicies([i, j])
        
            gas_den_mean += numpy.sum(gas[inds_gas_in_bin].n)
            n_part += inds_gas_in_bin.size
        
        gas_den_mean /= n_part
        
        #inds_gas = hist.get_indicies(i)
        print 'number of particles in this section of luminosity = %e, mean density = %e' % (n_part, gas_den_mean)
    #####################################
    
    return map_CO_data, map_NH2_data, hist, gas
    #
#

for snap in params['snaps']:    
    map_CO_data, map_NH2_data, hist, gas = generate_maps(snap, params)
    
pylab.show()