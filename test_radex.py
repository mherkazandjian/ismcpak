from radex import *
import subprocess

radexPath = '/home/mher/ism/code/radex/Radex/bin/'
inFile = { 'molData'                : 'hco+.dat'                              ,
           'outPath'                : 'hco+.out'                              ,
           'freqRange'              : [0, 50000]                              ,
           'tKin'                   : 20.0                                    ,
           'collisionPartners'      : ['H2']                                  ,
           'nDensCollisionPartners' : [1e4]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e13                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 0                                       }

radexObj = radex(radexPath)
radexObj.setInFile( inFile )
print radexObj.getInFile()

# writing the input file
fInput = open( 'dummyInFile.inp', 'w')
fInput.write( radexObj.genInputFileContentAsStr() )
fInput.close()

command  = radexPath + 'radex' +  ' <  ' + 'dummyInFile.inp ' 
command += ' > /dev/null ' 
print command
os.system( command )

data = ''
outFile = open( 'hco+.out', 'r')
i = 0
for line in outFile:
#    print line
    if i >= 11 :
        data += line
    i += 1

outFile.close()

print data