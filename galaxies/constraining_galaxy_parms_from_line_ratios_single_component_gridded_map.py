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
          'weights' : 'conserve-area',  #'conserve-area', 'matched',  #'original-only' ,#None ,#by-number          
           
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
#          'line_ratios' : [
#                           'CO2-1/CO1-0'    , 'CO3-2/CO1-0'    , 'CO4-3/CO1-0', 
#
#                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0',
#
#                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', '13CO3-2/CO3-2', '13CO4-3/CO4-3', 
#                        ],
          
          ### all the ratios up to J = 6-5
#          'line_ratios' : [
#                           'CO2-1/CO1-0'    , 'CO3-2/CO1-0'    , 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 
#
#                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',
#
#                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', '13CO3-2/CO3-2', '13CO4-3/CO4-3', '13CO5-4/CO5-4', '13CO6-5/CO6-5',
#                       ],

          ### other line ratios 1
          'line_ratios' : [
                           'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 'CO7-6/CO1-0', 
                           'CO8-7/CO1-0', 'CO9-8/CO1-0', 'CO10-9/CO1-0', 'CO11-10/CO1-0', 'CO12-11/CO1-0', 'CO13-12/CO1-0', 
                           'CO14-13/CO1-0', 'CO15-14/CO1-0', 

                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0', '13CO7-6/13CO1-0', 
                           '13CO8-7/13CO1-0', '13CO9-8/13CO1-0', '13CO10-9/13CO1-0', '13CO11-10/13CO1-0', '13CO12-11/13CO1-0', '13CO13-12/13CO1-0', 
                           '13CO14-13/13CO1-0', '13CO15-14/13CO1-0', 

                           '13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO4-3/CO1-0', '13CO5-4/CO1-0', '13CO6-5/CO1-0', '13CO7-6/CO1-0', 
                           '13CO8-7/CO1-0', '13CO9-8/CO1-0', '13CO10-9/CO1-0', '13CO11-10/CO1-0', '13CO12-11/CO1-0', '13CO13-12/CO1-0', 
                           '13CO14-13/CO1-0', '13CO15-14/CO1-0', 

                           'HCN1-0/CO1-0', 'HCN2-1/CO1-0', 'HCN3-2/CO1-0', 'HCN4-3/CO1-0', 'HCN5-4/CO1-0',  'HCN6-5/CO1-0', 'HCN7-6/CO1-0',

                           'HNC1-0/CO1-0', 'HNC2-1/CO1-0', 'HNC3-2/CO1-0', 'HNC4-3/CO1-0', 'HNC5-4/CO1-0',  'HNC6-5/CO1-0', 'HNC7-6/CO1-0',

                           'HCO+1-0/CO1-0', 'HCO+2-1/CO1-0', 'HCO+3-2/CO1-0', 'HCO+4-3/CO1-0', 'HCO+5-4/CO1-0',  'HCO+6-5/CO1-0', 'HCO+7-6/CO1-0',

#                           'SiO1-0/CO1-0', 'SiO2-1/CO1-0', 'SiO3-2/CO1-0', 'SiO4-3/CO1-0', 'SiO5-4/CO1-0',  'SiO6-5/CO1-0', 'SiO7-6/CO1-0',

#                           'CS1-0/CO1-0', 'CS2-1/CO1-0', 'CS3-2/CO1-0', 'CS4-3/CO1-0', 'CS5-4/CO1-0',  'CS6-5/CO1-0', 'CS7-6/CO1-0',

#                           'CN1-0/CO1-0', 'CN2-1/CO1-0', 'CN3-2/CO1-0', 'CN4-3/CO1-0', 'CN5-4/CO1-0',  'CN6-5/CO1-0', 'CN7-6/CO1-0',

                        ],

          
#          'lines'       : {
#                           'include'     : [
#                                            'CO1-0', 'CO2-1', 'CO3-2', 'CO4-3', 'CO5-4', 'CO6-5',                           
#                                            '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', '13CO5-4', '13CO6-5',                           
#                                           ],
#                           'combinations' :[
#                                            'CO/CO', '13CO/13CO', '13CO/CO'
#                                           ],
#                         },


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

## getting the unique transitions whose maps will be constructed
curve_attrs = attrs.keys()
print 'emissions to be extracted from the processed snapshopt'
for attr in curve_attrs: print '\t%s' % attr

## getting the species involved in those emissions
species = line_ratio_utils.species_involved(lines)
print 'Species invloved = ', species

## path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

## loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, species, load_pdr=params['pdr_sph'], load_only_em=curve_attrs)    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

## setting the radii based on the suggested weighting
gas.set_radii(weighting=params['weights'], rundir=params['rundir'], snap_index=params['snap_index'])

## checking for weird particles and taking care of them
gas.check_particles(params['check'], logger)

## keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

## makingt the 2d histogram of the gas particles
bs_min, bs_max = params['ranges']['box_size'].number

print 'getting the spatial distrubutions....'
hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx = bs_max, 
               nbins = params['imres'], reverse_indicies = True, loc = True)
hist.info()
print '\t\tdone getting the spatial distributuions'

## keeping the gas particles which are within the ranges of the histogram
gas = gas[hist.inds_in]

## getting the luminosity maps for each line
luminosity = {
              'lines' : numpy.array(lines, 'S'),
              'maps'  : {},
              'units' : {
                        'luminosity' : 'K.km.s^-1.kpc^2',
                        },
             }

print 'computing the luminosity from all the pixles in for each line map...'

print 'making a sample fine map (just for plotting purposes)...'
this_attr = attrs.keys()[0]
    
## the intensity map (intensity weight averaged intensity map)
#this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #func=numpy.average, weights=this_attr)
#this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.average, weights='weights')
this_map_luminosity = fi_utils.make_map(gas, hist, attr=this_attr, func=fi_utils.total_luminosity)
        
## getting the line code from the attribute name
line = this_attr.replace(params['em_unit']+'_','')
    
luminosity['maps'][line] = this_map_luminosity
    
print '\t\tfinished making the luminosity maps'

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

## setting up the object which grids the particles and fits
estimator = fi_utils.galaxy_gas_mass_estimator(luminosity_maps = luminosity, 
                                               gas = gas, 
                                               params = params, 
                                               hist = hist,
                                               line_ratios=line_ratios,
                                               arxvPDR=arxvPDR,
                                               em_unit=em_unit)

estimator.get_model_emission_from_pdr_arxv_involved_line_ratios(em_unit=em_unit) 

estimator.setup_observed_grid()

#estimator.estimate_mass_in_pixel(x_ind=4, y_ind=3)

#estimator.estimate_mass_in_all_pixels()

#getting the fits for all the pixels within a radia range
#pixels_info = estimator.estimate_mass_in_pixels_at_radius(r = 2)

#print 'total H2 mass          :', estimator.H2_mass_mesh.sum()
#print 'total H2 mass no gmech :', estimator.H2_mass_no_gmech_mesh.sum()