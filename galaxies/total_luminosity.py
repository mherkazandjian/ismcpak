'''
plots the total luminosity of all the SPH particles.
'''
#########################################################################################################
import os
import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units

from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd
import mylib.units 
from galaxies import fi_utils
import lineDict

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
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'original-only', #'by-number', #'matched',  #'original-only' ,#None ,#by-number          
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
          'em_unit'   : 'em_fluxKkms', # 'em_fluxcgs', 'em_fluxKkms'
#          'lines'     : ['CO1-0', '13CO1-0', 'HCN1-0', 'HNC1-0', 'CS1-0', 'SiO1-0', 'HCO+1-0'],
          'lines'     : ['CO1-0'],
          'save_data' : False,
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

## setting the radii and weights based on the suggested weighting
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

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

pylab.figure()

specsStrs =  ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS', 'SiO']
colors    =  ['k' , 'r'   , 'g'  , 'b'  , 'c'   , 'y' , 'm']

#specsStrs =  ['CO']
#colors    =  ['k']

#specsStrs =  ['HCN', 'HNC']
#colors    =  ['g'  , 'b'  ]

#specsStrs =  ['HCN', 'HNC', 'HCO+', 'CS', 'SiO']
#colors    =  ['g'  , 'b'  , 'c'   , 'y' , 'm']

print 'getting the luminsoty of the sampled matched set'

## plotting the total luminosity weighted by number of sampled points
gas.set_radii(weighting='matched', rundir=params['rundir'], snap_index=params['snap_index'])

sym       = '-'

em_unit = params['em_unit'].replace('em_flux', '')

lum_ladder = {} 

for i, specStr in enumerate(specsStrs):
    
    print specStr
    
    x, y = gas.get_total_luminosity_ladder(specStr, em_unit=em_unit)
    
    if em_unit == 'Kkms':
        y *= 1e6  # form K km /s kpc2 -> K km /s pc2
    elif em_unit == 'cgs':
            y = y * ((mylib.units.KPC2CM)**2) / mylib.constants.Lsun_erg_s # from ergs/cm2/s kpc2 -> Lsun
    else:
        raise ValueError('unknown unit %s' % em_unit)
    
    pylab.semilogy(x+1, y, colors[i] + sym, label=specStr)

    lum_ladder[specStr] = {'Jup': x+1, 'L': y}
    
print 'getting the luminsoty of the original set'
if True:
    ## plotting the total luminosity of the original points
    gas.set_radii(weighting='original-only', rundir=params['rundir'], snap_index=params['snap_index'])
    
    sym       = '--'
    
    lum_ladder_original = {}
    
    for i, specStr in enumerate(specsStrs):
        
        print specStr
        
        x, y = gas.get_total_luminosity_ladder(specStr, em_unit=em_unit)

        if em_unit == 'Kkms':
            y *= 1e6  # form K km /s kpc2 -> K km /s pc2
        elif em_unit == 'cgs':
            y = y * ((mylib.units.KPC2CM)**2) / mylib.constants.Lsun_erg_s # from ergs/cm2/s kpc2 -> Lsun
        else:
            raise ValueError('unknown unit %s' % em_unit)
        
        pylab.semilogy(x+1, y, colors[i] + sym)

        lum_ladder_original[specStr] = {'Jup': x+1, 'L': y}
    
pylab.legend(loc=0)

if em_unit == 'Kkms':
    ylabel = 'K km/s pc^2'
else:
    ylabel = ' L_sun'
    
pylab.ylabel(ylabel)
pylab.show()
