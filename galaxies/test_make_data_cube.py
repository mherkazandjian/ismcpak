import numpy
from numpy import exp, pi, sqrt

from amuse.datamodel import Particles
from mylib.utils.histogram import hist_nd
import pylab
import scipy.ndimage
import time

home = '/home/mher'

params = {
          #name of the file containing the post processed data
          #'filename' : '/home/mher/ism/runs/galaxies/coset2run4/coset-2-std/analysis/data.npz',
          'filename' : '/home/mher/ism/runs/galaxies/coset2run4/coset-9-sol/analysis/data.npz',
          
          #resolution of the map
          'imres' : 25,   # resolution of the maps to be produced imres x imres
          
          #spectrum info
          'spec_res' : 1000,             #number of velocity bins of the spectra to be constructed for each pixel
          'spec_rng' : [-300.0, 300.0],    #the range in the velocities
          
          #the size of the box to be displayed (particles outside the range are discarded)
          'box_size' : [-2.0, 2.0], #kpc

          #map value range
          'em_CO-1-0'  : {'title': r'$f(L_{CO(1-0})$', 'v_rng': [-10.0, -4.0], 'log10': True},
        
         }

#loading the data
data = numpy.load(params['filename'])

#assigning the contents to a particle set
gas = Particles(data['gas_x'].size)
gas.x = data['gas_x']
gas.y = data['gas_y']
gas.z = data['gas_z']
gas.vx = data['gas_vx']
gas.vy = data['gas_vy']
gas.vz = data['gas_vz']
gas.mass = data['gas_mass']
gas.em_CO_1_0 = data['em_CO_1_0']

#declaring some variables for convinience
bs_min, bs_max = params['box_size']
res = params['imres']

#making the map
hist = hist_nd(
               numpy.vstack((gas.x, gas.y)), 
               mn = bs_min,  
               mx = bs_max,
               nbins = res, 
               reverse_indicies=True, 
               loc=True
              )

if hist.data.shape[1] != len(gas):
    print 'warning: not all the gas particles are within the map boundaries'

#declaring the array holding the emissions
em_CO_1_0_map = numpy.zeros((res,res), 'f8')

#looping over the bins of the 2D histogram of the x,y coordinates and computing the averages of the maps
for i in numpy.arange(res):
        
    for j in numpy.arange(res):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            if inds_in_bin.size > 0:
                
                em_CO_1_0_map[i,j] = numpy.mean(gas.em_CO_1_0[inds_in_bin])

fig = pylab.figure(figsize=(8,8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

log_em_CO_1_0_map = numpy.log10(em_CO_1_0_map)
im = ax.imshow(log_em_CO_1_0_map.T, 
              extent=[bs_min, bs_max, bs_min, bs_max],
              vmin=params['em_CO-1-0']['v_rng'][0], 
              vmax=params['em_CO-1-0']['v_rng'][1], 
              interpolation='bessel', #intepolation used for imshow
              origin='lower')

#plotting the pixel boundaries
ax.plot(hist.f.spos[0], hist.f.spos[1], 'r+', markersize=100)

#plotting the spectrum for a trial pixel of the map
i, j = 12, 9 
inds_in_pixel = hist.get_indicies([i, j])
n = inds_in_pixel.size
print 'number of sph particles in pixel = %d' % n  
ax.plot(gas.x[inds_in_pixel], gas.y[inds_in_pixel], 'y.', markersize=1, color='k')


#constructing the spectrum corresponding to that pixel
v_min, v_max = params['spec_rng']
v_res = params['spec_res']

#the velocity bins 
v = numpy.linspace(v_min, v_max, v_res)

#the spectrum of all the particles in the pixel 
spect_pixel = numpy.zeros(v.size, 'f8')

v_width = 1.0

fig = pylab.figure()
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

#constructing the spectrum for a pixel
t0 = time.time()
for i, ind in enumerate(inds_in_pixel):  
    em_CO_1_0_this = gas.em_CO_1_0[ind]
    v_los_this = gas.vz[ind]
    print i, em_CO_1_0_this, v_los_this
    
    spect_this = exp( - (v - v_los_this)**2 / (2.0*v_width**2) )
    spect_this *= 1.0 / (v_width * sqrt(2.0 * pi) )
    spect_this *= em_CO_1_0_this
    
    spect_pixel += spect_this
    #pylab.plot(v, spect_this)
dt = time.time() - t0
print 'constructed spectrum in %.2e seconds' % (dt)
#plotting the spectrum
pylab.plot(v, spect_pixel, 'r')

#rebinnign the spectrum to a lower resolution
v_new           = scipy.ndimage.zoom(v          , 0.1)
spect_pixel_new = scipy.ndimage.zoom(spect_pixel, 0.1)
pylab.step(v_new, spect_pixel_new, 'b')
pylab.plot([-300.0, 300.0], [0.0,0.0], 'b')

pylab.xlabel(r'$v (km s^{-1})$')
pylab.ylabel(r'CO(1-0) spectrum')

pylab.show()
