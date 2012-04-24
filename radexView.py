"""
  plotting radex output for multiple LVG models 
"""

import subprocess
from time import *
import pylab as pyl
import numpy as np

from radex import *

Jall = np.arange(20) + 1

# path of the radex excutable
radexPath = '/home/mher/ism/code/radex/Radex/bin/radex'  
# parameters that will be passed to radex
inFile = { 'molData'                : 'co.dat'                              ,
           'outPath'                : 'foo'                              ,
           'freqRange'              : [0, 500000]                              ,
           'tKin'                   : None                           ,
           'collisionPartners'      : ['H2']                                  ,
           'nDensCollisionPartners' : [100.0]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e18                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 1                                       }

"""
inFile['molnDens'] = 1e23
parms = ( ('tKin', 10), 
          ('tKin', 20),
          ('tKin', 30),
          ('tKin', 40),
          ('tKin', 50),
          ('tKin', 60),
          ('tKin', 70),
          ('tKin', 80),
          ('tKin', 90), )
"""


inFile['tKin'] = 1000.0
parms = ( ('molnDens', 1e16), 
          ('molnDens', 1e17),
          ('molnDens', 1e18),
          ('molnDens', 1e19),
          ('molnDens', 1e20),
          ('molnDens', 1e21),
          ('molnDens', 1e22),
          ('molnDens', 1e23))


# creating the radex process instance
radexObj = radex(radexPath)
radexObj.setInFile( inFile )
fig, axs = radexObj.setupPlot( len(parms) )

t0 = time()

for i, p in enumerate( parms ):
    
    # re-setting the status
    radexObj.setDefaultStatus()
    
    print i, p[0], p[1]
    radexObj.setInFileParm(p[0], p[1]) # setting put the parameters
    radexObj.run(  checkInput = True ) # running

    if radexObj.getStatus() &  radexObj.FLAGS['SUCCESS']:    

        #print radexObj.getRawOutput()
        print 'number of iterations = %d' % radexObj.getNIter() 
        radexObj.plotModelInFigureColumn(Jall, axs[:,i], '%.1f' % p[1],   )

    else:
        print 'warnings'
        print '--------'
        print radexObj.warnings

    t1 = time()
    print 'time running radex = %f ' % (t1 - t0)    


radexObj.setLabels()
pyl.show()
