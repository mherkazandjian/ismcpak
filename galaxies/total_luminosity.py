'''
plots the total luminosity of all the SPH particles.
'''
#########################################################################################################
import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units

from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd
 
from galaxies import fi_utils
import lineDict

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
          'weights' : 'original-only', #'by-number', #'by-number', #'by-number', #'matched',  #'original-only' ,#None ,#by-number          
          'obs_res'      : 21,
           
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
          'check'   : 'default',          
          'em_unit'   : 'em_fluxKkms',
          'lines'     : [
                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3', #'CO5-4'  , 'CO6-5'  ,
                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', #'13CO5-4', '13CO6-5',
                         'HCN1-0',
                         'HNC1-0',
                         'CS1-0',
                         'SiO1-0',
                         'HCO+1-0'
                        ],
          'save_maps' : False,
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

## checking for weird particles and taking care of them
gas.check_particles(params['check'], logger)

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

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

pylab.figure()

specsStrs =  ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS', 'SiO']
colors    =  ['k' , 'r'   , 'g'  , 'b'  , 'c'   , 'y' , 'm']

## plotting the total luminosity weighted by number of sampled points
weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % params['snap_index']
gas.use_weights(weighting='by-number', weights_filename = weights_filename)

sym       = '-'

gas.use_weights
for i, specStr in enumerate(specsStrs):
    
    print specStr
    
    x, y = gas.get_total_luminosity_ladder(specStr)
    
    pylab.semilogy(x, y*1e6, colors[i] + sym, label=specStr)

if True:
    ## plotting the total luminosity of the original points
    gas.use_weights(weighting='original-only')
        
    sym       = '--'
    
    
    for i, specStr in enumerate(specsStrs):
        
        print specStr
        
        x, y = gas.get_total_luminosity_ladder(specStr)
        
        pylab.semilogy(x, y*1e6, colors[i] + sym)
    

pylab.show()









