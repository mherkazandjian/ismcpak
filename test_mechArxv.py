import os
import glob
from mesh import *
from meshUtils import *
import pickle

path = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/'

files = []
for infile in glob.glob( os.path.join(path, '*id-0000*') ):
    files.append(infile)
    
arxv = []
print 'found %s' % (len(files))
for file in files:
    print file
    m = mesh(file)
    arxv.append(m)

ma = meshArxv();
ma.setArxv(arxv)

mmm = mesh('/home/mher/mesh.dat-id-000444')

f = open('/home/mher/ism/foo.sav', 'wb')
pickle.dump(mmm , f)
f.close() 

"""
pkl_file = open('data.pkl', 'rb')

data1 = pickle.load(pkl_file)
pprint.pprint(data1)
"""

