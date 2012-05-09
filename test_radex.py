"""
driver script to run radex through python and read back the output
A single LVG run
"""

from radex import *
from time import *


# path of the radex excutable
radexPath      = '/home/mher/ism/code/radex/Radex/bin/radex'  
molDataDirPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'

# parameters that will be passed to radex single partner
"""
inFile = { 'specStr'                : 'CN'                              ,
           'freqRange'              : [0, 50000]                              ,
           'tKin'                   : 10.0                                    ,
           'collisionPartners'      : ['H2']                              ,
           'nDensCollisionPartners' : [1e3]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e14                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 1                                       }
"""

# parameters that will be passed to radex multiple partners
inFile = { 'specStr'                : 'CO'                                ,
           'freqRange'              : [0, 50000]                              ,
           'tKin'                   : 300.0                                    ,
           'collisionPartners'      : ['H2','H' , 'e-', 'H+', 'He']            ,
           'nDensCollisionPartners' : [1e1 , 1e2, 1e2, 1e1 , 1e3 ]            ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e18                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 1                                       }

# creating the radex process instance
radexObj = radex(radexPath, molDataDirPath)

t0 = time()
# setting put the parameters, running and parsing the output
radexObj.setInFile( inFile )

radexObj.run( checkInput = True, verbose = True)

if radexObj.getStatus() &  radexObj.FLAGS['SUCCESS'] :
    
        radexObj.parseOutput()
        t1 = time()

        print 'time running radex = %f ' % (t1 - t0)
    
        print 'number of iterations = %d' % radexObj.getNIter() 
        hcop10 = radexObj.getTransition(1)  # getting the info of the transiotion from 1->0
        #print hcop10
        
        # printing all transitions and fluxes
        print 'header'
        print '------'
        print radexObj.outputHdr
        
        for transition in radexObj.transitions:
            print transition['upper'], transition['lower'], transition['fluxcgs']
               
        radexObj.plotModel()
else:
    if radexObj.getStatus() &  radexObj.FLAGS['ITERWARN']:
        print 'did not converge'
        
        print 'warnings'
        print '--------'
        print radexObj.warnings
