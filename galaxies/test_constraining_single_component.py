'''
    - compute the 2D histogram
    - from the histogram of the projected particles on the sky 
    - from the histogram compute the mean intensity per bin (pixel) for each line ( K.km.s^-1)
    - from this compute the luminosity within each pixel
    - 
    -
    -
    -
    -
    -
    -
'''
#########################################################################################################
import os, time, sys, pickle

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import scipy
from scipy.stats import chisqprob

from amuse.units import units

from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd
 
from galaxies import fi_utils
import meshUtils

import line_ratio_utils
from line_ratio_utils import ratio_sets as r
from line_ratio_utils import *
import constraining

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params = {
          ##################### parameters for making the mock maps #########################
          
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100',  # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'matched', #'by-number', #'by-number', #'matched',  #'original-only' ,#None ,#by-number          
           
          'snap_index': 4,
          'ranges'    : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                         'sph':{
                               'min_log_n_use'  : -3.0,
                               'min_log_G0_use' : -3.0,
                               'min_log_gm_use' : -50.0,
                               'Av_use'         :  [0.0, 20000000.0],
                               'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                               },
                      
                        #the size of the box to be displayed (particles outside the range are discarded)
                        'box_size' : [-8.0, 8.0] | units.kpc, 
                        },
          'check'   : 'default',          
          'em_unit'   : 'em_fluxKkms', #'em_fluxKkms', 'em_fluxcgs'
          'save_maps' : False,

          'error_bars'  : 0.2, 
          'obs_res'      : 21, #9
                    
          ##################### parameters for fitting for the maps #########################
          
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0/',   # the path to the dir containing the PDR database
          #'pdrDb' :  home + '/ism/runs/oneSided/dynamicMesh-z-1.0/',
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',  # the path to the dir containing the PDR database

          ### all the ratios up to J = 2-1
#          'line_ratios' : [
#                           'CO2-1/CO1-0',  
#
#                           '13CO2-1/13CO1-0', 
#
#                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', 
#                        ],

          ### all the ratios up to J = 4-3
          'line_ratios' : [
                           'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO5-4/CO1-0',  'CO6-5/CO1-0', 'CO7-6/CO1-0', 'CO8-7/CO1-0', 'CO9-8/CO1-0', 'CO10-9/CO1-0', 'CO11-10/CO1-0', 'CO12-11/CO1-0', 'CO13-12/CO1-0', 

                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',

                           '13CO1-0/CO1-0', '13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO5-4/CO1-0', '13CO6-5/CO1-0',

                           'HCO+1-0/CO1-0', 'HCO+4-3/CO1-0', 'HCO+7-6/CO1-0',

                           'HCN1-0/CO1-0', 'HCN3-2/CO1-0', 'HCN4-3/CO1-0',

                           'HNC1-0/CO1-0', 'HNC3-2/CO1-0',
                        ],

          ### all the ratios up to J = 6-5
#          'line_ratios' : [
#                           'CO2-1/CO1-0'    , 'CO3-2/CO1-0'    , 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 
#
#                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',
#
#                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', '13CO3-2/CO3-2', '13CO4-3/CO4-3', '13CO5-4/CO5-4', '13CO6-5/CO6-5',
#                       ],


          
          #'interpolator' : scipy.interpolate.NearestNDInterpolator, 
          'interpolator' : scipy.interpolate.LinearNDInterpolator,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

em_unit = params['em_unit'].replace('em_flux','')

## setting up the logger object
logger = default_logger()

## setting up the line ratio strings 
if 'line_ratios' in params:
    line_ratios = params['line_ratios']
elif 'lines' in params:
    line_ratios = line_ratio_utils.line_ratio_combinations(params['lines']['include'], params['lines']['combinations'])
 
## getting the lines involved in the ratios
lines = line_ratio_utils.all_lines_involved(line_ratios)

## setting up the line attribute array whose emission will be retried from the snapshot
attrs = {}
for i, line in enumerate(lines):
    
        attr = params['em_unit'] + '_' + line
        
        attrs[attr] = True
#

## the observed data
obs = line_ratio_utils.observations()

obs['CO1-0'] = {'fluxcgs':0.3*1e-16*1e3}
obs['CO2-1'] = {'fluxcgs':2.4*1e-16*1e3}
obs['CO3-2'] = {'fluxcgs':7.3*1e-16*1e3}
obs['CO4-3'] = {'fluxcgs':12.8*1e-16*1e3}
obs['CO5-4'] = {'fluxcgs':17.1*1e-16*1e3}
obs['CO6-5'] = {'fluxcgs':17.2*1e-16*1e3}
obs['CO7-6'] = {'fluxcgs':18.2*1e-16*1e3}
obs['CO8-7'] = {'fluxcgs':17.9*1e-16*1e3}
obs['CO9-8'] = {'fluxcgs':12.2*1e-16*1e3}
obs['CO10-9'] = {'fluxcgs':9.3*1e-16*1e3}
obs['CO11-10'] = {'fluxcgs':7.7*1e-16*1e3}
obs['CO12-11'] = {'fluxcgs':5.5*1e-16*1e3}
obs['CO13-12'] = {'fluxcgs':3.9*1e-16*1e3}

obs['13CO1-0'] = {'fluxcgs':0.02*1e-16*1e3}
obs['13CO2-1'] = {'fluxcgs':0.2*1e-16*1e3}
obs['13CO3-2'] = {'fluxcgs':0.8*1e-16*1e3}
obs['13CO5-4'] = {'fluxcgs':0.9*1e-16*1e3}
obs['13CO6-5'] = {'fluxcgs':0.7*1e-16*1e3}

obs['HCO+1-0'] = {'fluxcgs':0.006*1e-16*1e3}
obs['HCO+4-3'] = {'fluxcgs':0.2*1e-16*1e3}
obs['HCO+7-6'] = {'fluxcgs':0.5*1e-16*1e3}

obs['HCN1-0'] = {'fluxcgs':0.009*1e-16*1e3}
obs['HCN3-2'] = {'fluxcgs':0.2*1e-16*1e3}
obs['HCN4-3'] = {'fluxcgs':0.1*1e-16*1e3}

obs['HNC1-0'] = {'fluxcgs':0.008*1e-16*1e3}
obs['HNC3-2'] = {'fluxcgs':0.09*1e-16*1e3}

## setting the errors of the fluxes
obs.set_flux_errors(em_unit='fluxcgs', err_bar_percent=0.1)
 
## make line ratios from the mock luminosities
obs_ratios = line_ratio_utils.ratios()

## computing the line ratios that will be used in the minimization
obs_ratios.make_ratios(obs, params['line_ratios'], em_unit='fluxcgs')

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

model_em, grid_coords = arxvPDR.get_model_emission_from_involved_line_ratios(params['line_ratios'], 
                                                                             em_unit='Kkms')


## fitting
f = constraining.Xi2_line_ratios_single_component(obs_data = obs_ratios, 
                                                  model_data = model_em, 
                                                  model_parms = grid_coords,
                                                  line_ratios = params['line_ratios'],
                                                  )
f.compute_model_line_ratios()
f.compute_Xi2()
f.print_minima()

fig = pylab.figure(figsize=(16,10))

# plotting the line ratios in these axes
ax1 = fig.add_axes([0.05, 0.55, 0.25, 0.4])
ax2 = fig.add_axes([0.35, 0.55, 0.25, 0.4])

# plotting the line ratios
f.plot_results(fig = fig, ax = ax1)
#print '-------------------------'
f.plot_results(fig = fig, ax = ax2, no_gmech=True)

pylab.show()
