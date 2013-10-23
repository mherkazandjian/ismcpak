import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
from numpy import exp, pi, sqrt

import scipy.ndimage

import pylab

import pyfits

from amuse.units import units
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
import fi_utils

import time

home = '/home/mher'

params = {
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'snap_index' : 4,

          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 20000000.0],
                             'Av_clip'        :  [3.0, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },

                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-10.0, 10.0] | units.kpc, #kpc
                     },
          
          #species data to be loaded
          'species': ['CO'],
          'pdr_sph': None,
        
          #cube info          
          'cube'   : {
                      'attr'     : 'em_fluxKkms_CO1-0'   , #the emission line to be made into a cube
                      'func'     : numpy.sum             , #function to be used to compute the total emission in the cell
                      'xy_rng'   : [-8.0, 8.0, -8.0, 8.0], #spatial bounds of the projected image (in kpc)
                      'im_res'   : 200,                    #spatial resolution of the cube in each dimension over the domain
                      'spec_rng' : [-100.0, 100.0],        #the range in the velocities in km/s
                      'spec_res' : 200,                   #number of velocity bins of the spectra to be constructed for each pixel
                      
                      'plot'     : {
                                    'as_log10'   : True,         #plot the map values in log10 scale
                                    'rng'        : [-10.0, 6.0], #range of value of the map
                                    'title'      : r'$f(L_{CO(1-0} K.km.s-1))$',
                                    #####info about the v channel maps
                                    'n_vsec_plt' : 20, #number of velocity channels maps
                                    'n_per_row'  : 5,  #number of maps per row 
                                   }, 
                      'l_width_mc': 1.0,   # Micro - Turbulance line width in km/s
                      'save_cube' : False, 
                      'save_fits' : True, 
                      },
         }
#############################################################################################################################

#setting up the logger object
logger = default_logger()

#path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

#getting the time of the snapshot
snap_time = fi_utils.get_snapshot_time(params['snap_index'], params)

#loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'], load_pdr=params['pdr_sph'])    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

#plotting the particles in the processed snapshot
fig_part = pylab.figure(figsize=(8,8))
ax = fig_part.add_axes([0.2, 0.2, 0.6, 0.6])
pylab.plot(gas.x[::100], gas.y[::100], '.', markersize=1)

#keeping gas particles within the specified ranges
gas_in_rng = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

#plotting the particles in the ranges specified by params['ranges']
pylab.plot(gas_in_rng.x[::100], gas_in_rng.y[::100], 'r+', markersize=1)

######################################### CONSTRUCTING EMISSION MAP ##########################################

#selecting the particles in the field of view of the cube
inds_in_map_spatial_ranges = numpy.where( 
                                         (gas_in_rng.x > params['cube']['xy_rng'][0])*(gas_in_rng.x < params['cube']['xy_rng'][1])*
                                         (gas_in_rng.y > params['cube']['xy_rng'][2])*(gas_in_rng.y < params['cube']['xy_rng'][3])
                                        )

gas_in_cube = gas_in_rng[inds_in_map_spatial_ranges]

#plotting the particles in field of view of the cube
pylab.plot(gas_in_cube.x[::100], gas_in_cube.y[::100], 'g+', markersize=1)


#binning the particles
hist = hist_nd(
               numpy.vstack((gas_in_cube.x, gas_in_cube.y)), 
               mn = [params['cube']['xy_rng'][0], params['cube']['xy_rng'][2]],   
               mx = [params['cube']['xy_rng'][1], params['cube']['xy_rng'][3]],
               nbins = params['cube']['im_res'], 
               reverse_indicies=True, 
               loc=True
              )

if len(gas_in_cube) != hist.data.shape[1]:
    raise ValueError('some gas particles have been excluded during binning')
 
pylab.draw()
pylab.show()

#declaring the array holding the emission
cube_map = numpy.zeros((params['cube']['im_res'],  params['cube']['im_res']), 'f8')

gas_em = getattr(gas_in_cube, params['cube']['attr'])

#looping over the bins of the 2D histogram of the x,y coordinates and computing 
#total emission in each pixel
for i in numpy.arange(params['cube']['im_res']):
        
    for j in numpy.arange(params['cube']['im_res']):
            
        inds_in_bin = hist.get_indicies([i,j])
            
        if inds_in_bin.size > 0:
            
            cube_map[i,j] = params['cube']['func'](gas_em[inds_in_bin])
        #
    #
#

fig_map = pylab.figure(figsize=(8,8))
ax_map = fig_map.add_axes([0.15, 0.085, 0.75, 0.75])

if params['cube']['plot']['as_log10'] == True:
    
    cube_map = numpy.log10(cube_map)
    
im = ax_map.imshow(cube_map.T, 
                   extent = params['cube']['xy_rng'],
                   vmin   = params['cube']['plot']['rng'][0], 
                   vmax   = params['cube']['plot']['rng'][1],
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

cbar_ax_map = fig_map.add_axes([0.2, 0.97, 0.6, 0.02]) 
cbar_ax_map.tick_params(axis='both', which='major', labelsize=30)

ax_map.set_title(params['cube']['plot']['title'], size=30, verticalalignment='bottom')
ax_map.tick_params(axis='both', which='major', labelsize=20)
        
pylab.colorbar(im, ax=ax_map, cax=cbar_ax_map, orientation='horizontal', 
                   ticks=numpy.linspace(params['cube']['plot']['rng'][0], params['cube']['plot']['rng'][1], 5))
    
pylab.figtext(0.01, 0.87, '%.2f' % snap_time + 'Gyr', 
              color='black', size='xx-large', weight='bold')

#plotting the pixel boundaries
ax_map.plot(hist.f.spos[0], hist.f.spos[1], 'w+', markersize=100)

ax_map.set_xlabel('x(kpc)', size='large')
ax_map.set_ylabel('y(kpc)', size='large')

######################################### CONSTRUCTING THE DATA CUBE ##########################################
# this can be easily parallelized by broadcasting the x,y,vlos, em, hist??, 
###############################################################################################################
data_cube = numpy.zeros((params['cube']['im_res'],  
                         params['cube']['im_res'], 
                         params['cube']['spec_res']), 'f8')

#getting the attributes of the gas which will be used to construct the data cube
x, y, vlos, em = gas_in_cube.x, gas_in_cube.y, -gas_in_cube.vz, getattr(gas_in_cube, params['cube']['attr'])

#parameters of a typical spectrum for the cube
v_min, v_max = params['cube']['spec_rng']
v_res = params['cube']['spec_res']

#the line width due to micro turbulance
v_width = params['cube']['l_width_mc']

#the velocity bins (.. todo:: xxx shift this to the centroid of the velcoty bins) 
v = numpy.linspace(v_min, v_max, v_res)

for i in numpy.arange(params['cube']['im_res']):
    
    print 'processing image row  i = %i' % i

    t0 = time.time()
    
    for j in numpy.arange(params['cube']['im_res']):
            
        inds_in_pixel = hist.get_indicies([i,j])
            
        if inds_in_bin.size > 0:
            
            n = inds_in_pixel.size

            #the spectrum of all the particles in the pixel 
            spect_pixel = numpy.zeros(v.size, 'f8')

            #constructing the spectrum for a pixel
            for l, ind in enumerate(inds_in_pixel):
                  
                em_p = em[ind]
                v_los_p = vlos[ind]
                
                #computing the shifted gaussian
                spect_this = exp( - (v - v_los_p)**2 / (2.0*v_width**2) )
                #scaling the gaussian 
                spect_this *= (em_p / (v_width * sqrt(2.0 * pi) ))
                
                spect_pixel += spect_this
            #
            dt = time.time() - t0
            
            #copying the spectrum of this pixel to the data cube
            data_cube[i, j, :] = spect_pixel[:]
        #
    #
    print 'constructed spectruma for this row of pixels in %.2e seconds (%.2f %% complete)' % (time.time() - t0, 100.0*numpy.double(i)/params['cube']['im_res'])

#
if params['cube']['save_cube'] == True:
    numpy.savez_compressed(os.path.join(params['rundir'], 'analysis', 'fiout.%06d.cube' % params['snap_index']) , 
                           [data_cube, params], 
                           names=['data_cube','params'])
if params['cube']['save_fits'] == True:
    cube_fits = data_cube.swapaxes(2,0)
    pyfits.writeto(os.path.join(params['rundir'], 'analysis', 'fiout.%06d.cube' % params['snap_index']) + '.fits',
                   cube_fits, 
                   clobber=True) 

######################################### Plotting the data cube ##########################################
# plotting the data cube in sections of velocity 
###########################################################################################################
fi_utils.plot_cube_sections(data_cube, params)

def plot_pixel_spectrum(i, j, data_cube, hist, params, zoom=None):
    '''
    '''
    
    if zoom == None:
        zoom = 1.0
        
    #----------------------------------
    #plotting the spectrum for a trial pixel of the map
    inds_in_pixel = hist.get_indicies([i, j])
    n = inds_in_pixel.size
    print 'number of sph particles in pixel = %d' % n  
    #ax_map.plot(x[inds_in_pixel], y[inds_in_pixel], 'y.', markersize=1, color='k')
    
    #constructing the spectrum corresponding to that pixel
    v_min, v_max = params['cube']['spec_rng']
    v_res = params['cube']['spec_res']
    
    #the velocity bins 
    v = numpy.linspace(v_min, v_max, v_res)
    
    #the spectrum of all the particles in the pixel 
    spect_pixel = data_cube[i,j,:]
    
    fig_spec = pylab.figure()
    ax_spec = fig_spec.add_axes([0.2, 0.2, 0.6, 0.6])
        
    #plotting the spectrum
    ax_spec.plot(v, spect_pixel, 'r')
    
    #rebinnign the spectrum to a lower resolution
    v_new           = scipy.ndimage.zoom(v          , zoom)
    spect_pixel_new = scipy.ndimage.zoom(spect_pixel, zoom)
    pylab.step(v_new, spect_pixel_new, 'b')
    pylab.plot(params['cube']['spec_rng'], [0.0,0.0], 'b--')
    
    pylab.xlabel(r'$v (km s^{-1})$')
    pylab.ylabel(r'T$_{\rm antenna}$')
        
    pylab.draw()
    pylab.show()    
#

    
pylab.draw()
pylab.show()


print 'done'