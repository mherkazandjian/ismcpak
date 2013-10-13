# plot the log of ngas, G0, and gmech of the sph particles in a 
# color coded grid
#==============================================================

import numpy as np
import pylab as pyl
import matplotlib.cm as cm

home = '/home/mher'

fnameLoad = home + '/ism/runs/galaxies/test1/test/xyz.npz'
levels = np.arange(-30, -20, 1)
ranges = ( 0, 6, 0, 6 )

f = np.load(fnameLoad)
lx, ly, lz = f['lx'], f['ly'], f['lz']
px, py, pz = f['px'], f['py'], f['pz']
m, tKin, h2frac = f['m'], f['tKin'], f['h2frac']

inds = np.where( (lx > 0.0)*(ly > 0.0) )
lx, ly, lz = lx[inds], ly[inds], lz[inds]

print 'n sph particles with this range = ', len(inds[0]) 

# defining the points in the uniform 2D grid
grid_x, grid_y = np.mgrid[ranges[0]:ranges[1]:100j, ranges[2]:ranges[3]:100j]

points = np.array([lx,ly]).T
values = lz

# interpolating with different methods
from scipy.interpolate import griddata
grid_z0 = griddata(points, values, (grid_x, grid_y), method='linear')

fig = pyl.figure(0, (6,6))
pyl.subplot(111)

pyl.plot(lx[::500], ly[::500], 'ko', alpha = 0.1, markersize = 5)
pyl.axis(ranges)
pyl.xlabel('$log_{10} n$')
pyl.ylabel('$log_{10} G_0$')

im  = pyl.imshow(grid_z0.T, interpolation = 'bilinear', origin = 'lower',
                cmap = cm.jet, extent = ranges )
pyl.xlabel('$log_{10} n$')
pyl.ylabel('$log_{10} G_0$')


#CS = pyl.contour(grid_z0, levels, origin = 'lower', extent = ranges, colors = 'black')
#pyl.clabel(CS)

# make a colorbar for the contour lines
CB = pyl.colorbar(im, shrink=0.8, extend='both')
CB.set_ticks( levels )
pyl.show()