from radex import *
import subprocess
from time import *

# path of the radex excutable
radexPath = '/home/mher/ism/code/radex/Radex/bin/radex'  
# parameters that will be passed to radex
inFile = { 'molData'                : 'co.dat'                              ,
           'outPath'                : 'foo'                              ,
           'freqRange'              : [0, 50000]                              ,
           'tKin'                   : 20.0                                    ,
           'collisionPartners'      : ['H2']                                  ,
           'nDensCollisionPartners' : [1e2]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e24                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 1                                       }

# creating the radex process instance
radexObj = radex(radexPath)

t0 = time()
# setting put the parameters, running and parsing the output
radexObj.setInFile( inFile )
radexObj.run()
#print radexObj.getRawOutput()           
radexObj.parseOutput()
t1 = time()

print 'time running radex = %f ' % (t1 - t0)

print 'warnings'
print '--------'
print radexObj.warnings

print 'number of iterations = %d' % radexObj.getNIter() 
hcop10 = radexObj.getTransition(1)  # getting the info of the transiotion from 1->0
#print hcop10

# printing all transitions and fluxes
print 'header'
print '------'
print radexObj.outputHdr

for transition in radexObj.transitions:
    print transition['upper'], transition['lower'], transition['fluxcgs'] 
