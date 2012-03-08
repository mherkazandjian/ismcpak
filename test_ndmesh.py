import numpy as np
import pylab as pyl
from ndmesh import *
from time import *

# testing two 2 meshes
#--------------------------------------
nMin =  0.0; nMax = 6.0; dn = 0.25 # density range
GMin = -0.5; GMax = 6.0; dG = 0.25 #G0 range

gridCntrds = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]  # cell centroid coordinates
grid_spos  = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]  # cell starting coordinates
grid_epos  = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]  # cell ending coordinates

# centroids
gridCntrds[0,:] = gridCntrds[0,:] + dn/2.0
gridCntrds[1,:] = gridCntrds[1,:] + dG/2.0 

# start pos
grid_spos[0,:] = gridCntrds[0,:] - dn/2.0
grid_spos[1,:] = gridCntrds[1,:] - dG/2.0 

# end pos
grid_epos[0,:] = gridCntrds[0,:] + dn/2.0
grid_epos[1,:] = gridCntrds[1,:] + dG/2.0 

""" 
print  gridCntrds
print '----------------'
print  grid_spos
print '----------------'
print  grid_epos
print '----------------'
print gridCntrds[0,:,0] # row 1 centroids (fixed G0 varying n )
print gridCntrds[0,0,:] # col 1 centroids (fixed n  varying G0)
print '----------------'
print gridCntrds[0,:,1] # row 2 centroids (fixed G0 varying n )
print gridCntrds[0,2,:] # col  centroids (fixed n  varying G0)
"""

xc = grid_spos[0,:,:]
yc = grid_spos[1,:,:]
xe = grid_epos[0,:,:]
ye = grid_epos[1,:,:]

from matplotlib.lines import Line2D                        
from matplotlib.artist import setp
 
fig = pyl.figure()
ax = pyl.gca()
pyl.xlim( xmin = -2, xmax = 12)
pyl.ylim( ymin = -2, ymax = 12)
fig.add_axes(ax)


shape = np.array(gridCntrds.shape)
nDim   = shape[1]
nCells = np.prod( shape[1::] )
print "number of cells  = %d" % nCells

gc = np.resize(gridCntrds, (2, nCells))
gs = np.resize(grid_spos , (2, nCells))
ge = np.resize(grid_epos , (2, nCells))

t0 = time()
allLines = ()
for i in np.arange( nCells ):
    poss = gs[:,i] # [xs, ys, zs...]
    pose = ge[:,i] # [xe, ye, ze...]
    
    # extracting the starting and ending coords as individual variables for
    # clarity
    xs = poss[0]; xe = pose[0]
    ys = poss[1]; ye = pose[1]
    # defining the path for each reactangle
    xPathRect = [ xs, xe, xe, xs, xs ]
    yPathRect = [ ys, ys, ye, ye, ys ]
    # plotting the rectangle
    lines = pyl.plot( xPathRect , yPathRect )
    allLines = allLines + (lines,)
    
setp(allLines, linewidth=0.5, color='r')

print 'time rendering          : %f s' % (time() - t0)

pyl.show()
