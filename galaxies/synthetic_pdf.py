#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab
import scipy

from numpy import *
from pylab import *

from amuse.units import units
from mylib.utils.misc  import default_logger 
from galaxies import fi_utils, coset9_sol_info
import meshUtils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params = {
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100',  # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : [],#'CO'], #'13CO', 'HCN', 'HNC', 'HCO+'],
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'original-only', #by-number', #'matched',  #'original-only' ,#None ,#by-number          
          'snap_index' : numpy.arange(4, 4 + 1, 1),
          'ranges'      : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
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
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

## path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'
    
## loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename)
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'])
 
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

## setting the radii and weights based on the suggested weighting
gas.set_radii(weighting=params['weights'], rundir=params['rundir'], snap_index=params['snap_index'])

## checking for weird particles and taking care of them
gas.check_particles(params['check'], logger)

#keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

gaso = gas[gas.get_inds_original_set()]

################################################################################################
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

F = arxvPDR.load_interp_func(info={'source':'radex'}, line='HCN1-0', quantity='fluxKkms')

########################

# radius density scaling function, r(kpc) = R(n[cm-3]) (kpc) 
r_func = coset9_sol_info.r_sph_kpc

## the probability density at that density (not in log scale
nPDF = coset9_sol_info.nPDF_sph.pdf
#nPDF = coset9_sol_info.nPDF_sph_test.pdf

nmin, nmax = -3, 6.0

fi_utils.luminosity_from_pdf(nPDF, r_func, gaso, gas, arxvPDR, F, params,  nmin, nmax)

print 'done'