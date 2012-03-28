from radex import *
import subprocess

radexPath = '/home/mher/ism/code/radex/Radex/bin/radex'
inFile = { 'molData'                : 'hco+.dat'                              ,
           'outPath'                : 'hco+.out'                              ,
           'freqRange'              : [0, 50000]                              ,
           'tKin'                   : 20.0                                    ,
           'collisionPartners'      : ['H2']                                  ,
           'nDensCollisionPartners' : [1e3]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e13                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 1                                       }

radexObj = radex(radexPath)
radexObj.setInFile( inFile )
#print radexObj.getInFile()

# writing the input file

radexObj.run()
radexObj.parseOutput()

print radexObj.rawOutput
print radexObj.outputHdr
for transition in radexObj.transitions:
    print transition['upper'], transition['lower'], transition['fluxcgs'] 
print 