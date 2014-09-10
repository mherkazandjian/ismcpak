import pyfits
import numpy
import pylab
from numpy import exp, zeros, sqrt, linspace, arange, meshgrid

fits_filename = '/home/mher/tmp/test_fits.fits'
nx, ny, nz    = 50, 50, 100  #data cube size (spatial_x_res, spatial_y_res, v_res)
x_rng         = [-2.0, 2.0]
y_rng         = [-3.0, 3.0] 
v_rng         = [-50.0, 50.0]
v_center      = 10.0 
v_width       = 10.0
s_width       = 0.3  #source width
###################################################################
cube = zeros((nx, ny, nz), 'f8')

#generating the signal (this will be assigned to all the pixel)
v = linspace(v_rng[0], v_rng[1], nz)
f = (1.0/(sqrt(2.0)*v_width))*exp(-(((v - v_center)**2)/v_width**2))

#plotting the spectrum (just to check if it is ok)
fig = pylab.figure()
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
ax.plot(v,f)

x = linspace(x_rng[0], x_rng[1], nx) 
y = linspace(y_rng[0], y_rng[1], ny)

x_grd, y_grd = meshgrid(x,y)
image = exp(-(x_grd**2 + y_grd**2)/2.0) 

for i in arange(nx):
    for j in arange(ny):        
        cube[i,j,:] = image[i,j]*f[:]
    #
#

cube = cube.swapaxes(2,0)
pyfits.writeto(fits_filename, cube, clobber=True) 

print 'done'