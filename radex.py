import numpy as np
import pylab as pyl
import os 

class radex( ):
    """
    this class interfaces between radex. i.e run radex through python. 
    """

    def __init__(self, execPath):
        
        self.execPath  = execPath # absolute path to the radex executable
        self.inFile    = None     # dictionary holding all the input parameters
        self.nCollPart = None     # number of collision partners
         
    def setInFileParm(self, parm, value):
        self.inFile[parm] = value
    def getInFileParm(self, parm):
        return self.inFile[parm]

    def setInFile(self, inFile):
        self.inFile = inFile
    def getInFile(self):
        return self.inFile
    
    def genInputFileContentAsStr(self):
        
        self.nCollPart = len(self.inFile['collisionPartners'])
        
        strng = ''
        
        strng  += '%s\n' % self.inFile['molData']
        strng  += '%s\n' % self.inFile['outPath']
        strng  += '%d %d\n' % (self.inFile['freqRange'][0], self.inFile['freqRange'][1])
        strng  += '%f\n' % (self.inFile['tKin'])
        strng  += '%d\n' % self.nCollPart
        
        for i in np.arange( self.nCollPart ):
            strng += '%s\n' % self.inFile['collisionPartners'][i]
            strng += '%e\n' % self.inFile['nDensCollisionPartners'][i]

        strng += '%f\n' % self.inFile['tBack']
        strng += '%e\n' % self.inFile['molnDens']
        strng += '%f\n' % self.inFile['lineWidth']
        strng += '%d\n' % self.inFile['runAnother']
        
        return strng

"""    
inFile = { 'molData'                : '/data1/mher/ism/code/radex/Radex/data/',
           'outPath'                : None                                    ,
           'freqRange'              : [0, 0]                                  ,
           'tKin'                   : 20.0                                    ,
           'collisionPartners'      : ['H2']                                  ,
           'nDensCollisionPartners' : [1e4]                                   ,
           'tBack'                  : 2.73                                    ,
           'molnDens'               : 1e13                                    ,
           'lineWidth'              : 1.0                                     ,
           'runAnother'             : 0                                       }
"""