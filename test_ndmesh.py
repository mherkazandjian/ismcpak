import numpy as np
import pylab as pyl
from ndmesh import *

dims = (3, 5)
ndm = ndmesh( dims, dtype = np.int32)
ndm.set_ranges( [-1,1,5,10] )
#print ndm
#print ndm.getRanges()
ndm.misc()
ndm.centroids()
#print x


print ndm.dl
print ndm.ldim
