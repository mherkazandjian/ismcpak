import time, os

import matplotlib
matplotlib.use('Qt4Agg')
from pylab import *

from numpy import repeat, nan, arange, where, hstack, ones
import numpy
import pylab

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units

from galaxies import fi_utils

#===========================================================================================================
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params =\
 {
   'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',    # the path of the dir containing the simulation
   'snap_index' : 4,
   'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
   'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
               'sph':{
                      'min_log_n_use'  : -3.0,      
                      'min_log_G0_use' : -3.0,
                      'min_log_gm_use' : -50.0,
                     },
               #the size of the box to be displayed (particles outside the range are discarded)
               'box_size' : [-20, 20] | units.kpc,
               },
 }
#===========================================================================================================

# load the full snapshot
#if True:
if False:
    #extracting/guessing the metallicity from the name of the directory of the run
    metallicity = fi_utils.guess_metallicity(params['rundir'])

    # setting the filename
    suffix = '%06d' % params['snap_index']
    snapName = 'fiout.%s' % suffix 
    filename = params['rundir'] + '/firun/' + snapName 

    #loading the sph simulation data
    print 'loading snapshot %s : ' % filename
    gas_fi, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

    print 'gas fraction = ', gas_fi.mass.sum()/(gas_fi.mass.sum() + dark.mass.sum())

    # getting the gas particles in cgs units
    gas = fi_utils.convert_units_to_pdr_units(gas_fi, metallicity)

    print 'done reading fi snapshot : %s' % filename
    print 'number of sph particles in snapshot = %d' %  len(gas)

    # selecting particles with density greater than a specified value
    inds = numpy.where( gas.n > 10 )
    gas_in = gas[inds]
    
    # saving the subset particles
    attr_list = [
                 'x', 'y', 'z', 'radius', 'vdisp', 'mass', 'n', 'G0',
                 'Av', 'gmech', 'id', 'Pe', 'T', 'vx', 'vy', 'vz',
                ]
    fi_utils.save_gas_particle_info('gas_subset.npz', gas_in, attr_list)
    

# load a subset of saved the particles read above
#if False:
if True:
    
    gas, attrs_loaded = fi_utils.load_gas_particle_info('gas_subset.npz', 
                                                           load_pdr=False)

    # selecting a set of particles within a spatial box
    
    inds = numpy.where( (gas.x > -660.0/1e3)*(gas.x < -560.0/1e3)*\
                        (gas.y > -400.0/1e3)*(gas.y < -350.0/1e3) )
    
    g = gas[inds]
    
    x,y,z,r = g.x*1e3, g.y*1e3, g.z*1e3, g.radius  # r is in pc
    G0, gmech, sigma = g.G0, g.gmech, g.vdisp 
    
    #pylab.plot(x,y,'.')
    #pylab.show()
    
    for p in g:
        
        # finding the relative distance of all the particles
        # relative to this particle
        d_rel = numpy.sqrt(
                           (p.x - g.x)**2 +
                           (p.y - g.y)**2 +
                           (p.z - g.z)**2
                          )
        
        # checking how many particles are within the smoothing 
        # length of the sph particle
        inds = numpy.where( d_rel < p.radius/1e3 )

        # getting the particles withing the smoothing length
        g_within_smooth = g[inds]
        
        #print p.radius, len(inds[0])
        
        #pylab.plot( p.x, p.y, 'mo')
        #pylab.plot( g.x, g.y, 'r.')
        #pylab.plot( g.x[inds], g.y[inds], 'b.')
        #pylab.show()
        
        
        print 'n    ', g_within_smooth.n.min(), g_within_smooth.n.max() 
        print 'Av   ', g_within_smooth.Av.mean(), g_within_smooth.Av.std()
        print 'G0   ', g_within_smooth.G0.mean(), g_within_smooth.G0.std()
        print 'gmech', g_within_smooth.gmech.mean(), g_within_smooth.gmech.std() 
        print 'vdisp', g_within_smooth.vdisp.mean(), g_within_smooth.vdisp.std() 
        
        print '-----------'
        
print 'done'
