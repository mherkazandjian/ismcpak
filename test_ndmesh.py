import numpy as np
import pylab as pyl
from ndmesh import *

# density range
nMin = 0.0
nMax = 6.0
dn   = 2.0

#G0 range
GMin = -0.5
GMax =  6.0
dG   =  2.0

gridCntrds = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]
grid_spos  = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]
grid_epos  = np.mgrid[ nMin:nMax:dn, GMin:GMax:dG ]

# centroids
gridCntrds[0,:] = gridCntrds[0,:] + dn/2.0
gridCntrds[1,:] = gridCntrds[1,:] + dG/2.0 

# start pos
grid_spos[0,:] = gridCntrds[0,:] - dn/2.0
grid_spos[1,:] = gridCntrds[1,:] - dG/2.0 

# end pos
grid_epos[0,:] = gridCntrds[0,:] + dn/2.0
grid_epos[1,:] = gridCntrds[1,:] + dG/2.0 
 
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
print gridCntrds[0,2,:] # col 3 centroids (fixed n  varying G0)