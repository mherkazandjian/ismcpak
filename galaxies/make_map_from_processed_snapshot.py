'''
 in this test file, we load a database of PDR meshes and extract information
 at conditions determined by the SPH simulations

  - exclude sph particles which have densities below which the LVG model codes can handel 
    (for example, one with n_gas < 1e-3))
  - exlcude particle with alpha > 1 (unless for very high densities (get the densitiy of 
    the specie whose  emission is to be computed and see if its abundance is too low, if 
    if is too low, just set the emission of the sph particle corresponding to that specie
    to zero.
  - abundances of CO for PDR modesls with n < 0.1 cm^-3 are ~ 1e-10 (which is too low compared to the
    high density modesl

'''
#########################################################################################################
import time
import sys
import os
import subprocess
import multiprocessing
from IPython.parallel import Client

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate
import numpy
from numpy import log10
import pylab
import logging

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants, nbody_system
import meshUtils
from mylib.utils.misc  import xselect, default_logger
from mylib.utils.histogram import hist_nd
from mylib.utils.interpolation import sectioned_4D_interpolator 
from fi_utils import parse_old_runinfo_file
import fi_utils
import lineDict
from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

snaps = numpy.arange(4, 4 + 1, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',      # the path to the dir containing the PDR database
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2/',          # the path to the dir containing the PDR database          
          'use_em' : True, 
          'species' : ['CO', '13CO'],
          
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [3.0, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                      },
#          'maps'   : {
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
          'maps'   : {
                      'attr' : 'em_fluxKkms_13CO1-0', #'mass', 'G0', 'gmech', 'Av'
                      'v_rng': [-10.0, 4.0],
                      'title': r'$f(L_{13CO(1-0} K.km.s-1)$', 
                      'log10': True                      
                     }, 
 
                      #'f_mean_n' : {'pos': [0,2], 'title': r'$f(\bar{n})$'  , 'v_rng': [-3.0, 4.0], 'log10': True},
                      #'f_mean_g0': {'pos': [0,3], 'title': r'$f(\bar{g_0})$', 'v_rng': [-3.0, 3.0], 'log10': True},
                      #--------
                      #'f_mean_gm'              : {'pos': [1,0], 'title': r'$f(\bar{\Gamma_m})$'                , 'v_rng': [-35.0, -22.0], 'log10': True},
                      #'f_mean_Av'              : {'pos': [1,1], 'title': r'$f(\bar{Av})$'                      , 'v_rng': [0.0  , 2.0] ,  'log10': True},
                      #'T_mean'                 : {'pos': [1,2], 'title': r'$f(\bar{T})$'                        , 'v_rng': [0.0  , 5.0] ,  'log10': True},
                      #'f_mean_em_C+-1-0'       : {'pos': [1,1], 'title': r'$f(L_{C^+ 158 \mu m})$'               , 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_no_gm_C+-1-0' : {'pos': [1,2], 'title': r'$f(L_{C^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_CO-7-6'       : {'pos': [1,3], 'title': r'$f(L_{O^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #--------
                      #'f_mean_em_CO1-0'       : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'                  , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_CO2-1'       : {'pos': [2,1], 'title': r'$f(L_{CO(1-0})$'                  , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_no_gm_CO-1-0' : {'pos': [2,1], 'title': r'$f(L_{CO(1-0)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #'f_mean_em_CO-3-2'       : {'pos': [2,2], 'title': r'$f(L_{CO(3-2})$'                  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #'f_mean_em_no_gm_CO-3-2' : {'pos': [2,3], 'title': r'$f(L_{CO(3-2)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #--------
                      #'f_mean_em_HCN-1-0'       : {'pos': [3,0], 'title': r'$f(L_{HCN(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
                      
        #'save_maps' : False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

def generate_maps(snapIndex, params):
    
    #path to processed fi snapshot  
    snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snapIndex + '.states.npz'  
    
    #loading the sph simulation data 
    logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
    gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'])
    
    logger.debug('done reading fi snapshot : %s' % snap_filename)
    logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))
    
    path = params['rundir'] + '/firun/runinfo'
    runinfo = parse_old_runinfo_file(path)

    snap_time = (float(runinfo['dtime'])*timeUnit.number) * (float(runinfo['noutbod']) * snap)

    #keeping gas particles within the specified ranges
    gas = fi_utils.select_particles(gas, params['ranges'])
    logger.debug('got the sph particles in the required ranges')
    logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

    #making the 2D histogram
    print 'getting the spatial distrubutions'
    
    hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
    hist.info()
    print 'done getting the spatial distributuions'

    #getting the map
    map_data = fi_utils.make_map(params['maps']['attr'], gas, hist, as_log10=params['maps']['log10'])
        
    #displaying all the maps in a single plot
    fig = pylab.figure(0, figsize=(8,8))
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    pylab.setp(ax, 
               'xlabel', 'x(kpc)', 'ylabel', 'y(kpc)', 'xlim', [bs_min, bs_max], 'ylim', [bs_min, bs_max],
               'aspect', 'equal', 'adjustable', 'datalim'
               )

    im = ax.imshow(map_data, 
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=params['maps']['v_rng'][0],  
                   vmax=params['maps']['v_rng'][1], 
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')
    ax.set_title(params['maps']['title'], size='large')
        
    pylab.colorbar(im, ax=ax, orientation='vertical')
        
    pylab.show()
#

for snap in snaps:    
    generate_maps(snap, params)

pylab.show()