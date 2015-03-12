#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units, nbody_system
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
from galaxies import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

fig_save_path = '/home/mher/ism/docs/paper02.5/src/figs/x_factor.eps'
#fig_save_path = None

##########################################DISK GALAXY##################################################
params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 400,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO'],
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(4, 4 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [0.01, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                      },

          'all_maps' : {
                        'map_CO1-0': {
                                     'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                                     'v_rng'   : [-10.0, 4.0],
                                     'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                                     'as_log10': True,
                                     'func'    : numpy.mean,
                                     #'func'    : numpy.average,
                                     #'weights' : 'em_fluxKkms_CO1-0',                                     
                                     },
                        'map_NH2'  : {
                                     'attr'    : 'pdr_NH2',
                                     'v_rng'   : [10.0, 30.0],
                                     'title'   : r'$N(H2)$', 
                                     'as_log10': True,
                                     'func'    : numpy.mean,
                                     #'func'    : numpy.average,
                                     #'weights' : 'em_fluxKkms_CO1-0',                                     
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
    fig = pylab.figure(figsize=(4,4))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

    LCO = map_CO_data.flatten()
    NH2 = map_NH2_data.flatten()

    #plotting the line corresponding to the canonical X factor
    x = numpy.linspace(-5.0, 5.0, 100)    
    ax.plot(x, numpy.log10(2.0e20*(1.0 - 0.3)) + x,'--k', alpha=0.3)
    ax.plot(x, numpy.log10(2.0e20) + x,'k', alpha=0.3)
    ax.plot(x, numpy.log10(2.0e20*(1.0 + 0.3)) + x,'--k', alpha=0.3)

    #ax.plot([1.0, 1.0], [16, 25], '--g', alpha=0.3)
    
    #plotting NH2 vs CO luminosity
    ax.plot(LCO[::10], NH2[::10], '.', markersize=2)
    ax.set_xlabel(r'$\log_{10}$ [<W$_{{\rm CO}(1-0)}$> / K km s$^{-1}$]', size=10)
    ax.set_ylabel(r'$\log_{10}$ [<N(${\rm H}_2)$> / cm$^{-2}$]'   , size=10)
    ax.text(-2.5, 20,r'N(${\rm H}_2)$ = $10^{20}$ W$({{\rm CO}(1-0)})$', size=10, rotation=42.0)

    ##
    print 'getting the distribution as a function of L(CO)'
    hist_LCO = hist_nd(LCO.clip(-5.0,LCO.max()).reshape((1,LCO.size)), nbins=50.0, reverse_indicies=True, loc=True)
    hist_LCO.info()

    mean_NH2_fact = numpy.zeros(hist_LCO.f.shape, 'f8')
    for i in numpy.arange(hist_LCO.nBins[0]):
        
        inds_in_bin = hist_LCO.get_indicies([i])

        if inds_in_bin.size > 0:
            
            mean_NH2_fact[i] = numpy.log10(numpy.mean(10.0**NH2[inds_in_bin])) 
            
    print 'done getting the distributuion as a function of W(CO)'

    ax.plot(hist_LCO.f.cntrd.flatten(), mean_NH2_fact, 'g--', linewidth=2.5)
    #ax.plot([1, 1], [16, 25], 'k-.')
    ax.plot([0, 0], [16, 25], 'k-.')

    ax.set_ylim(16, 25)    
    ax.set_xlim(-8,3)

    ################################################################################
    #checking how much of the emission is resulting from differnet section in L(CO) 
    ################################################################################
    LCO_sections = [ [-1.0, 0.0], [0.0, 1.0], [2.0, 3.0], [3.0, 4.0], [4.0, numpy.Inf]]
    
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
    
    return map_CO_data, map_NH2_data, hist, gas, fig, ax
    #
#

for snap in params['snaps']:    
    map_CO_data_disk, map_NH2_data_disk, hist_disk, gas_disk, fig, ax = generate_maps(snap, params)
##########################################END PLOTTING STUFF FROM THE DISK GALAXY##################################################




##########################################DWARF GALAXY##################################################
params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 200,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO'],
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(20, 20 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [0.01, 29.9],  #sph particles with Av higher than this are clipped to this value                             
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
                                     #'func'    : numpy.mean,
                                     'func'    : numpy.average,
                                     'weights' : 'em_fluxKkms_CO1-0',                                     
                                     },
                        'map_NH2'  : {
                                     'attr'    : 'pdr_NH2',
                                     'v_rng'   : [10.0, 30.0],
                                     'title'   : r'$N(H2)$', 
                                     'as_log10': True,
                                     #'func'    : numpy.mean,
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
        
    return map_CO_data, map_NH2_data, hist, gas
    #
#

for snap in params['snaps']:    
    map_CO_data_dwarf, map_NH2_data_dwarf, hist_dwarf, gas_dwarf = generate_maps(snap, params)
##########################################END PLOTTING STUFF FROM THE DISK GALAXY##################################################

LCO = map_CO_data_dwarf.flatten()
NH2 = map_NH2_data_dwarf.flatten()
ax.plot(LCO[::2], NH2[::2], 'r.', markersize=4)

if fig_save_path != None:
    fig.savefig(fig_save_path)
print 'saved image file to :\n\t\t\t %s' % fig_save_path

pylab.show()