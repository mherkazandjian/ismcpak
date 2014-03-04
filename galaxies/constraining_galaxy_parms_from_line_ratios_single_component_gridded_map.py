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
import pylab

from amuse.units import units

from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd
 
from galaxies import fi_utils
import lineDict
import meshUtils

import line_ratio_utils
import mylib.units
from mylib.constants import M_SUN_SI
from mylib.utils.ndmesh import ndmesh
from mylib.utils.misc import scale

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {
          ##################### parameters for making the mock maps #########################
          
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext',  # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'by-number', #'by-number', #'by-number', #'matched',  #'original-only' ,#None ,#by-number          
           
          'snap_index': numpy.arange(4, 4 + 1, 1),
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
          'em_unit'   : 'em_fluxKkms',
          'lines'     : [
                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3', #'CO5-4'  , 'CO6-5'  ,
                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', #'13CO5-4', '13CO6-5',
                        ],
          'save_maps' : False,
                    
          ##################### parameters for fitting for the maps #########################
          
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',  # the path to the dir containing the PDR database
          #'pdrDb' :  home + '/ism/runs/oneSided/dynamicMesh-z-1.0/',

          ### all the ratios up to J = 2-1
#          'line_ratios' : [
#                           'CO2-1/CO1-0', #'CO3-2/CO1-0', 'CO4-3/CO1-0', #'CO5-4/CO1-0', 'CO6-5/CO1-0', 
#
#                           '13CO2-1/13CO1-0', #'13CO3-2/13CO1-0', '13CO4-3/13CO1-0', #'13CO5-4/13CO1-0', '13CO6-5/13CO1-0',
#
#                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', #'13CO3-2/CO3-2', '13CO4-3/CO4-3', #'13CO5-4/CO5-4', '13CO6-5/CO6-5',
#                        ],

          ### all the ratios up to J = 4-3
          'line_ratios' : [
                           'CO2-1/CO1-0'    , 'CO3-2/CO1-0'    , 'CO4-3/CO1-0', #'CO5-4/CO1-0', 'CO6-5/CO1-0', 

                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', #'13CO5-4/13CO1-0', '13CO6-5/13CO1-0',

                          '13CO1-0/CO1-0', '13CO2-1/CO2-1', '13CO3-2/CO3-2', '13CO4-3/CO4-3', #'13CO5-4/CO5-4', '13CO6-5/CO6-5',
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

          'error_bars'  : 0.2, 

          #'interpolator' : scipy.interpolate.NearestNDInterpolator, 
          'interpolator' : scipy.interpolate.LinearNDInterpolator, 
          'obs_res'      : 21,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

## setting up the logger object
logger = default_logger()

## setting up the line attribute array whose emission will be retried from the snapshot
attrs = {}
for i, line in enumerate(params['lines']):
    
        attr = params['em_unit'] + '_' + line
        
        attrs[attr] = True
#

## getting the unique transitions whose maps will be constructed
curve_attrs = attrs.keys()
print 'emissions to be extracted from the processed snapshopt'
for attr in curve_attrs: print '\t%s' % attr

## getting the species involved in those emissions
species = {}
for key in curve_attrs:
    species[lineDict.lines[key.replace(params['em_unit'] + '_', '')]['specStr']] = True
print 'Species invloved = ', species.keys()

## path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

## loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, species, load_pdr=params['pdr_sph'])    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

## setting the weights
weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % params['snap_index']
gas.use_weights(weighting=params['weights'], weights_filename = weights_filename)

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
              'lines' : numpy.array(params['lines'], 'S'),
              'maps'  : {},
              'units' : {
                        'luminosity' : 'K.km.s^-1.kpc^2',
                        },
             }

print 'computing the luminosity from all the pixles in for each line map...'

print 'making the maps of all the lines...'
for i, this_attr in enumerate(attrs):
    
    ## the intensity map (intensity weight averaged intensity map)
    #this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #func=numpy.average, weights=this_attr)
    this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.average, weights='weights')
    
    ## computing the luminsoty by mutiplying by the area of each pixel
    this_map_luminosity = this_map_intensity * hist.f.dl.prod()
    
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

## setting up the line ratio strings 
if 'line_ratios' in params:
    line_ratios = params['line_ratios']
elif 'lines' in params:
    line_ratios = line_ratio_utils.line_ratio_combinations(params['lines']['include'], params['lines']['combinations'])

## loading the emission info from all the models for all Avs (also check for the consistenscy of the 
## number of models...i.e same number of models for all the lines)

## make line ratios from the mock luminosities 
obs_mock_ratios_template = line_ratio_utils.ratios()

for line_ratio in line_ratios:

    line1, line2 = line_ratio_utils.lines_involved(line_ratio)
    
    v1, v2 = 1, 1
    
    obs_mock_ratios_template.make_ratios(
                                         {
                                           line1:{'fluxKkms': v1, 'err': params['error_bars']*v1}, 
                                           line2:{'fluxKkms': v2, 'err': params['error_bars']*v2}
                                         },
                                         ratios = [line_ratio],
                                         em_unit = 'fluxKkms'
                                        )

    print line_ratio, obs_mock_ratios_template[line_ratio]
    
obs_mock_ratios_template.species_and_codes()

                    
model_em = {}
for i, line in enumerate(obs_mock_ratios_template.codes):
    v, grid_coords = arxvPDR.get_emission_from_all_radex_dbs_for_Av_range(
                                                                          line = line, 
                                                                          Avs = 'all', 
                                                                          quantity = 'fluxKkms',
                                                                          keep_nans = True,
                                                                         )
    model_em[line] = 10.0**v
    print line, v.size, grid_coords.shape
    print '----------------------'

    ## some checks of the sizes
    if v.size != grid_coords.shape[0]:
        raise ValueError('number of elements in the emission values is different from the number of modesl.')
    
    if i == 0:
        nModels = v.size
    else:
        if nModels != v.size:
            raise ValueError('the number of elements for this line differes at least from that of one of the other lines')
#########

estimator = fi_utils.galaxy_gas_mass_estimator(luminosity_maps = luminosity, 
                                               gas = gas, 
                                               params = params, 
                                               hist = hist,
                                               line_ratios=line_ratios,
                                               arxvPDR=arxvPDR,)

estimator.get_model_emission_from_involved_line_ratios()

estimator.setup_observed_grid()

#estimator.estimate_mass_in_pixel(x_ind=4, y_ind=3)

#estimator.estimate_mass_in_all_pixels()

#getting the fits for all the pixels within a radia range
#pixels_info = estimator.estimate_mass_in_pixels_at_radius(r = 2)

#print 'total H2 mass          :', estimator.H2_mass_mesh.sum()
#print 'total H2 mass no gmech :', estimator.H2_mass_no_gmech_mesh.sum()








