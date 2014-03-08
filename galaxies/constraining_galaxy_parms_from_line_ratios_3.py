'''
In this routine, we compute the Xi2 statistic for line ratios by using available emission data from teh PDR
grids. Here the interpolated emissions are used...

THIS ROUTINE IS NOT COMPLETE, BUT THE BASICS ARE THERE, INTERPOALTING...ETC.. NO XI2 COMPUTATIONS YET

  - load pdr database
  - load computed line ratios as a function of beam size
  - load the emission for the specified line ratios to be included in the Xi2 statistic
  - compute the Xi2
  
'''
#########################################################################################################
import matplotlib
matplotlib.use('Qt4Agg')

import scipy
import numpy
import pylab
import time

from amuse.units import units

from mylib.utils.misc  import default_logger

import fi_utils
import meshUtils

#======================================================parameters=================================================

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',          # the path to the dir containing the PDR database
                    
          'use_em'  : True, 
          'use_pdr' : False,
          
          'snaps'  : numpy.arange(0, 20 + 1, 1),           
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200000.0],
                             'Av_clip'        :  [3.0, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-20.0, 20.0] | units.kpc, #kpc

                       #modes in the PDR arxv which are within those ranges will be used in constructing interpolation functions 
                      'interp'   : {'log_n'    : [-10.0,  10.0],
                                    'log_G0'   : [-2.0 ,  3.0 ],
                                    'log_gmech': [-50.0, -20.0],
                                    'Av'       : [ 3.0 ,  28.0],
                                    'Av_res'   : 1.0,   #resolution of constructing the interp funcs from the PDR data  
                                    },   #Av_res is used only for quantities interpolated upon from PDR meshes 
                      },
          
          'emission': {},
          
          'maps'   : {
                      'pdr_em_<q>fluxKkms</q>_CO1-0'  : {'log10': True},
                      'pdr_em_<q>fluxKkms</q>_CO2-1'  : {'pos': [2,1], 'title': r'$f(L_{CO(1-0}[K.km.s-1)$' , 'v_rng': [-10.0, 4.0] , 'log10': True},
                      'pdr_em_<q>fluxKkms</q>_CO6-4'  : {'pos': [2,1], 'title': r'$f(L_{CO(1-0}[K.km.s-1)$' , 'v_rng': [-10.0, 4.0] , 'log10': True},
                      'pdr_em_<q>fluxcgs</q>_13CO1-0' : {'pos': [2,2], 'title': r'$f(L_{13CO(1-0})$'        , 'v_rng': [-10.0, 4.0], 'log10': True},
                      },

          'interpolator' : scipy.interpolate.NearestNDInterpolator, 
          #'interpolator' : scipy.interpolate.LinearNDInterpolator, 
          'save_info'   : False,
          'save_species' : ['CO', '13CO'],
          }
#############################################################################################################

## setting up the logger object
logger = default_logger()

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)


## loading the molecular specie radex database
#arxvPDR.readDbsRadex(allDbs=True)
arxvPDR.readDbsRadex(species=['CO'], in_Av_rng=params['ranges']['interp']['Av'])

## the keys of the emission maps only
em_keys = fi_utils.get_emission_only_keys(params['maps'])

## making all the interpolation functions
em_interp_funcs, pdr_interp_funcs = fi_utils.make_interp_funcs_from_arxv(arxvPDR, em_keys, {}, params, logger)

## getting the emission of all the lines to be used in the minimization over a 4D grid



logn  = numpy.array([0.0])
logG0 = numpy.array([0.0])
logGM = numpy.array([-50.0])
AV    = numpy.array([10.0])
data_use = numpy.array([logn, logG0, logGM, AV]).T
em_a = em_interp_funcs['pdr_em_<q>fluxKkms</q>_CO2-1'](data_use)
em_b = em_interp_funcs['pdr_em_<q>fluxKkms</q>_CO1-0'](data_use)

print em_a/em_b 






