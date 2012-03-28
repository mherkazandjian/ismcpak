import numpy as np
import pylab as pyl
import subprocess
from string import *
import re

class radex( ):
    """
    this class interfaces between radex. i.e run radex through python.
    
    methods :
        __init__
        setInFileParm(parm, value)
        
        setInFile(inFile)
        getInFileParm(parm)
        getInFile()
        genInputFileContentAsStr()
        getOutpur()
        run()
        parseOutput()  
    """

    def __init__(self, execPath):
        
        self.execPath    = execPath # absolute path to the radex executable
        self.inFile      = None     # dictionary holding all the input parameters
        self.nCollPart   = None     # number of collision partners
        self.proccess    = None     # the subprocess.Popen object
        self.rawOutput   = None     # the output of the run
        self.lineInfo    = None     # line information from the output
        self.outputHdr   = None     # the header of the output, should be 
                                    # consistent with inFile
        self.transitions = None
        
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
        strng += '%d'   % self.inFile['runAnother']
        
        return strng

    # run radex with the set input 
    def run(self):
        radexInput = self.genInputFileContentAsStr()
        #print '###################'
        #print radexInput
        #print '###################'
        self.proccess = subprocess.Popen(self.execPath           , 
                                         stdin=subprocess.PIPE   ,  
                                         stdout=subprocess.PIPE  ,  
                                         stderr=subprocess.PIPE  )
        radexOutput = self.proccess.communicate(input=radexInput)[0]
        #print self.proccess.poll()
        #print '-----------------'
        #print radexOutput
        #print '-----------------'
        self.rawOutput = radexOutput 
        
    # run radex with the set input 
    def getOutput(self):
        return self.output

    def parseOutput(self):
        
        header = {}
        
        output = self.rawOutput
        lines  =  output.splitlines()
        nLines = len(lines)
        
        # getting the number of iterations
        lineSplt = lines[4].split()
        if lineSplt[0] != 'Warning:':
            header['nIter'] = np.int32(lineSplt[2])
            header['converged'] = 1
        else:
            header['nIter'] = np.int32(lineSplt[6])
            header['converged'] = 0

        lineNum = 5
        # getting the rest of the header
        lineSplt = lines[lineNum].split(':')
        header['version'] = lineSplt[1]
        lineNum += 1
        lineSplt = lines[lineNum].split(':')
        header['geometry'] = lineSplt[1]
        lineNum += 1
        lineSplt = lines[lineNum].split(':')
        header['molDataFile'] = lineSplt[1]
        lineNum += 1         
        lineSplt = lines[lineNum].split(':')
        header['tKin'] = np.float64(lineSplt[1])
        
        collPartnerDens = np.zeros(self.nCollPart) 
        for i in np.arange(self.nCollPart ):
            lineNum += 1
            lineSplt = lines[lineNum].split(':')
            collPartnerDens[i] = np.float64(lineSplt[1])
        header['collPartDens'] = collPartnerDens 

        lineNum += 1         
        lineSplt = lines[lineNum].split(':')
        header['tBack'] = np.float64(lineSplt[1])

        lineNum += 1         
        lineSplt = lines[lineNum].split(':')
        header['molnDens'] = np.float64(lineSplt[1])

        lineNum += 1         
        lineSplt = lines[lineNum].split(':')
        header['lineWidth'] = np.float64(lineSplt[1])

        self.outputHdr = header 
        
        transitions = []
        #--------------------------------------------------------------------------
        # parses a line containing the info of a transition line and returns a dict
        # of the into
        def parseLineData( line ):
            info = {}
            
            lineSplt = line.split()
            
            info['upper'   ] = np.int32(lineSplt[0]) 
            info['lower'   ] = np.int32(lineSplt[1]) 
            info['E_up'    ] = np.float64(lineSplt[2]) 
            info['Tex'     ] = np.float64(lineSplt[5]) 
            info['tau'     ] = np.float64(lineSplt[6]) 
            info['T_R'     ] = np.float64(lineSplt[7]) 
            info['pop_up'  ] = np.float64(lineSplt[8]) 
            info['pop_down'] = np.float64(lineSplt[9]) 
            info['fluxKkms'] = np.float64(lineSplt[10]) 
            info['fluxcgs' ] = np.float64(lineSplt[11]) 
            
            return info
        #--------------------------------------------------------------------------
            
        lineNum += 4        
        for line in lines[lineNum:nLines-1]:
            transition = parseLineData(line)
            transitions.append(transition)
            
        self.transitions = transitions
