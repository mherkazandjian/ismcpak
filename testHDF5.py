import h5py
from numpy import *

# using "core" as the driver load everything into memeory
# it is written to disk optionally

f = h5py.File('/home/mher/foo.hdf5', 'w')

n = 10

x = zeros( (1,n), dtype=float64)
y = zeros( (1,n), dtype=float64)


dset = f.create_dataset('x', data=x)
dset = f.create_dataset('y', data=y)

f.close()
