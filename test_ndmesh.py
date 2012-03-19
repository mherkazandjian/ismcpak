import numpy as np
import pylab as pyl
from ndmesh import *
from time import *
import utils

# testing two 2 meshes
#--------------------------------------
nMin =  -2.0; nMax = 6.0;  nn  = 10 # density range # x axis (0)
GMin =   0.0; GMax = 10.0; nG0 = 5  # G0 range      # y axis (1)
#---------------------------------------

dn = np.float64((nMax - nMin)/nn)
dG = np.float64((GMax - GMin)/nG0)

x = ndmesh( (nn, nG0), dtype=np.float64 )
x.fill(0)

ranges = [ [nMin, nMax], [GMin, GMax] ]
x.setup( ranges )

fig, axs = pyl.subplots(1, 2, figsize = (12,6) )
ax1 = axs[0]; ax1_n = 121
ax2 = axs[1]; ax2_n = 122

pyl.subplot(ax1_n)
#fig = pyl.figure( figsize = (12,12) )
#ax = pyl.gca()
#fig.add_axes(ax)

pyl.xlim( xmin = nMin, xmax = nMax)
pyl.ylim( ymin = GMin, ymax = GMax)

#pyl.plot([1,2,3,34,5])

x.plotGrid(fig)
x.plotCntrd(fig)

#pyl.hold(True)


pyl.subplot(ax2_n)

xy = x.getCntrd()
xc = xy[0,:]; yc = xy[1,:]
x[:] = np.sin(xc + yc)

xpt = [0.0, 0.0, 0.5, 0.99]
ypt = [0.0, 0.5, 0.5, 0.5]
vpt = [2.1, 3.3, 4.6, 5.9]
xi =  utils.scale(xpt, 0, nn, 0, 1, integer=True)
yi =  utils.scale(ypt, 0, nG0, 0, 1, integer=True)
x[xi, yi] = vpt

ax2  = pyl.imshow( x.transpose(), cmap = pyl.cm.jet, origin = 'bottom', extent=[nMin, nMax, GMin, GMax] )

pyl.hold(True)
x.plotGrid(fig)
x.plotCntrd(fig)

pyl.contour(x, extent=[nMin, nMax, GMin, GMax])

pyl.show()
