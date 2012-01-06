import os
import glob
from mesh import *
from meshUtils import *
import pickle
from numpy import *
import numpy as np
from chemicalNetwork import *
from enumSpecies import *
from mesh import *

meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/meshes/mesh.dat-id-002920'
abunSpecFname = 'data/species.inp'

specFileUsed = enumSpecies(abunSpecFname, offsetFromZero=1 )
eSpcs = specFileUsed.genPythonClass(offset=0)

m = mesh(meshFname);

#m.plot(eSpcs)
#pyl.show()


hdrFormat =  [ 
                  ('version'  , np.int32  , 1),
                  ('G0'       , np.float64, 1),
                  ('nGas'     , np.float64, 1),
                  ('gammaMech', np.float64, 1),
                  ('nSteps'   , np.int32  , 1),
                  ('nSpecs'   , np.int32  , 1),
             ]
hdrDtype  = np.dtype( hdrFormat  )

#s = np.empty( 1, dtype = hdrDtype )
#print s


nSteps = 63
nSpecs = 309
dataFormat = [
                  ('gasT'  , np.float64, (nSteps, 1) ),
                  ('dustT' , np.float64, (nSteps, 1) ),
                  ('Av'    , np.float64, (nSteps, 1) ),
                  ('abun'  , np.float64, (nSpecs, nSteps) ),
             ]
dataDtype  = np.dtype( dataFormat  )

fooDtype = np.dtype( [ ('hdr', hdrDtype, 1), ('data', hdrDtype, 1) ]  )



s = np.empty( 1, dtype = fooDtype )
s['hdr']['version'] = 122
print s
print s['hdr']['version']

#desc = hdrDtype.descr
#print desc[3][1].__class__

desc = dataDtype.descr
print desc[3]






"""
meshFmt = [ 
            ('hdr' , np.dtype(hdrFormat) , 1),
            ('data', np.dtype(dataFormat), 1),
          ]

s = np.empty( [10,10], dtype=dtHdr  )
d = np.empty( [10,10], dtype=dtData )

meshDtype = np.dtype( meshFmt )

s = np.empty( 1, dtype = meshDtype )

print s.hdr
"""

"""

print len(s[0,0])
x = s[0,0]
x['version'] = 1
x['G0']      = 312

print x['version']
print x['G0']

hdrFmt = m.headerFormat()
print hdrFmt[0]

print d.dtype
#print s[0]'name']='asdad'
"""

"""
x = np.float64( (1,100) )
print s[0,0].fields

f = open('/home/mher/ism/foo.sav', 'wb')
pickle.dump(s , f)
f.close() 

a = float64(0)
print a.nbytes
"""
