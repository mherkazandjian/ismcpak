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
        getRawOutput()
        run()
        parseOutput()  
        tansitionInfo = getTransition( upper )
        warnings = getWarnings()
        nIter = getNIter()
        fig, axs = setupPlot(nx) # sets up the object to plot the radex output
    """

    def __init__(self, execPath):
        
        self.execPath    = execPath # absolute path to the radex executable
        self.inFile      = None     # dictionary holding all the input parameters
        self.nCollPart   = None     # number of collision partners
        self.proccess    = None     # the subprocess.Popen object
        self.rawOutput   = None     # the output of the run
        self.lineInfo    = None     # line information from the output
        self.warnings    = None     # the warnings dumped by radex
        self.nIter       = None     # number of iterations 
        self.outputHdr   = None     # the header of the output, should be 
                                    # consistent with inFile
        self.transitions = None

        #plotting attributes
        self.fig = None
        self.axs = None
        
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
    def getRawOutput(self):
        return self.rawOutput

    def parseOutput(self):
        
        output = self.rawOutput
        lines  =  output.splitlines()
        nLines = len(lines)
        
        # looking for the part of the output after the warning
        i = 4 # 5th line (skipping into)
        while True:
            if lines[i][0] == '*' and lines[i+1][0] == '*':
                break
            i += 1
        
        self.warnings = lines[4:i]  # extracting the warning
        #print self.warnings
        
        # collecting the header
        lineNumHdrStart = i
        while True:
            if lines[i][0] == '*' and lines[i+1][0] != '*':
                break
            i += 1
        lineNumHdrEnd = i+1
        
        self.outputHdr = lines[lineNumHdrStart:lineNumHdrEnd]
        #print  self.outputHdr

        lineNum = lineNumHdrEnd
        lineSplt = lines[lineNum].split()
        self.nIter = np.int32(lineSplt[3])          
        
        transitions = []
        #--------------------------------------------------------------------------
        # parses a line containing the info of a transition line and returns a dict
        # of the into
        def parseLineData( line ):
            info = {}
            
            lineSplt = line.split()
            #print lineSplt
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

        lineNum += 3        
        for line in lines[lineNum:nLines-1]:
            transition = parseLineData(line)
            transitions.append(transition)
            
        self.transitions = transitions

    def getTransition(self, upper):
        
        for transition in self.transitions:
            if transition['upper'] == upper :
                return transition
            
        #if upper value is not found, return Null
        return None
    
    def getWarnings(self):
        return self.warnings
    def getNIter(self):
        return self.nIter
    
    # sets up the object to plot the radex output ( nx rows and four coluns)
    def setupPlot(self, nx):
        # axs[0,0] is the one in the top left corner
        self.fig, self.axs = pyl.subplots(4, nx, sharex = False, sharey = False, figsize=(8,8) )

        pyl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9,
                            wspace=0.0, hspace=0.0)
        
        return (self.fig, self.axs)
    
    # plots the output of a certain model in a certain column in a predefined figure
    def plotModelInFigureColumn(self, Jall, colNum, title, axsArg = None):
        
        if axsArg == None:
            axs = self.axs
        else:
            axs = axsArg
            
        nx = axs.shape[1]
        
        plotInColumn = colNum
        subPlotNum   = plotInColumn + 1  

        nJ = len(Jall)
        #----------------flux-------------------------
        pyl.subplot(4, nx, subPlotNum)
        pyl.hold(True)
        xPlot = Jall
        yPlot = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis = transition['fluxcgs'] 
            yPlot[i] = yThis
            
        pyl.semilogy(xPlot, yPlot)
        pyl.axis([np.min(Jall), np.max(Jall),1e-10, 1])        
        
        #---------------Tex and T_R--------------------------
        subPlotNum += nx
        pyl.subplot(4, nx, subPlotNum)
        pyl.hold(True)

        xPlot = Jall
        yPlot1 = np.ndarray(nJ, dtype=np.float64)
        yPlot2 = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
    
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis1 = transition['Tex']
            yThis2 = transition['T_R'] 
    
            yPlot1[i] = yThis1
            yPlot2[i] = yThis2
        pyl.semilogy(xPlot, yPlot1, 'b')
        pyl.semilogy(xPlot, yPlot2, 'r')
        pyl.axis([np.min(Jall), np.max(Jall), 1, 10000])

        subPlotNum += nx 
        pyl.subplot(4, nx, subPlotNum)
        pyl.hold(True)

        xPlot = Jall
        yPlot = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
    
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis = transition['tau'] 
    
            yPlot[i] = yThis
        pyl.plot(xPlot, yPlot)
        pyl.axis([np.min(Jall), np.max(Jall), -1, np.max(yPlot)])


        subPlotNum += nx 
        pyl.subplot(4, nx, subPlotNum)
        pyl.hold(True)

        xPlot = Jall
        yPlot1 = np.ndarray(nJ, dtype=np.float64)
        yPlot2 = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
    
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis1 = transition['pop_up']
            yThis2 = transition['pop_down'] 
    
            yPlot1[i] = yThis1
            yPlot2[i] = yThis2
        pyl.semilogy(xPlot, yPlot1, 'b')
        pyl.semilogy(xPlot, yPlot2, 'r')
        pyl.axis([np.min(Jall), np.max(Jall), 1e-10, 1])

        # plotting the ylabels of the axes
        pyl.subplot(4, nx, 0*nx + 1);  pyl.ylabel('Flux')
        pyl.subplot(4, nx, 1*nx + 1);  pyl.ylabel('Tex, Trot')
        pyl.subplot(4, nx, 2*nx + 1);  pyl.ylabel('tau')
        pyl.subplot(4, nx, 3*nx + 1);  pyl.ylabel('pop dens')
        # plotting the xlabels
        for i in np.arange(nx+1):
            pyl.subplot(4, nx, 3*nx + i);  pyl.xlabel('J')
        # plotting the title
        pyl.subplot(4, nx, colNum + 1); pyl.title(title)

        # cleaning the tick labels and plotting the labels
        # removing the x label from the top rows
        for axsRow in axs[0:-1]:
            for ax in axsRow:
                for tick in ax.axes.xaxis.get_major_ticks():
                    tick.label1On = False

        # removing the y label from columns rows next to the first one
        for axsCol in axs[:, 1:nx]:
            for ax in axsCol:
                for tick in ax.axes.yaxis.get_major_ticks():
                    tick.label1On = False

        for ax in np.hstack( (axs[4-1], axs[: ,0]) ):
            xticks = ax.axes.xaxis.get_major_ticks()
            yticks = ax.axes.yaxis.get_major_ticks()
            for tick in [ xticks[0], xticks[-1], yticks[0], yticks[-1] ]:
                tick.label1On = False
