import numpy as np
import pylab as pyl
import subprocess
from string import *
#import re

class radex( ):
    """
    this class interfaces between radex. i.e run radex through python.
    
    after calling run(), this sets the flag determining whether the run
    was successfull or not : 
        status =     => run was successful
        status = 2    => run was successfull, but with warnings
        status = 3    => run was successfull, but with warnings
    
    methods :
        __init__
        setInFileParm(parm, value)
        
        setInFile(inFile)
        getInFileParm(parm)
        getInFile()
        setStatus( True|False )
        status = getStatus()
        genInputFileContentAsStr()
        getRawOutput()
        run( checkInput = False (default) | True )
        parseOutput()  
        tansitionInfo = getTransition( upper )
        warnings = getWarnings()
        nIter = getNIter()
        fig, axs = setupPlot(nx) # sets up the object to plot the radex output
    """

    def __init__(self, execPath):
        
        self.execPath     = execPath # absolute path to the radex executable
        self.moldataFiles =  { 'CO'   : 'co.dat'       ,
                               '13CO' : '13co.dat'     ,
                               'HCO+' : 'hco+@xpol.dat',
                               'HCN'  : 'hcn.dat'      ,
                               'HNC'  : 'hnc.dat'      ,
                               'CS'   : 'cs@xpol.dat'  ,
                               'CN'   : 'cn.dat'       }
        self.inFile       = None     # dictionary holding all the input parameters
        self.nCollPart    = None     # number of collision partners
        self.proccess     = None     # the subprocess.Popen object
        self.rawOutput    = None     # the output of the run
        self.lineInfo     = None     # line information from the output
        self.warnings     = None     # the warnings dumped by radex
        self.nIter        = None     # number of iterations 
        self.outputHdr    = None     # the header of the output, should be 
                                    # consistent with inFile
        self.transitions = None
        self.FLAGS       = {'DEFAULT'  : 0x000,   # default status, should be called before each run
                            'ERROR'    : 0x001,   # failed
                            'PARMSOK'  : 0x002,   # parameters are ok
                            'RUNOK'    : 0x004,   # radex did all the iterations (but not neccesarly converged)
                            'SUCCESS'  : 0x008,   # succeeded with no major warning (output should make sense)      
                            'WARNING'  : 0x010,   # success with warnings, warning are
                                                  # assigned to self.warnings
                            'ITERWARN' : 0x020}   # number of iterations warning
                            
                             
        self.status      = self.FLAGS['DEFAULT']

        #plotting attributes
        self.fig = None
        self.axs = None
        self.nx  = None
        self.ny  = 4
        
    def setInFileParm(self, parm, value):
        self.inFile[parm] = value
    def getInFileParm(self, parm):
        return self.inFile[parm]

    def setInFile(self, inFile):
        self.inFile = inFile
    def getInFile(self):
        return self.inFile
    def getStatus(self):
        return self.status
    def setStatus(self, status):
        self.status = status
    def setDefaultStatus(self):
        self.status = self.FLAGS['DEFAULT']
    
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

    def checkParameters(self):
        inFile = self.inFile

        for item in inFile:
            if inFile[item] == None:
                strng = 'Error : missing input parameter %s ' % (item)
                raise NameError(strng)
            
        # checking for correct range for the kinetic temperature
        if inFile['tKin'] < 0.1 or inFile['tKin'] > 1e4 :
            self.status |= self.FLAGS['ERROR']
            
        # checking for correct range for the densities of the collion partners
        for i, collPartner in enumerate(inFile['collisionPartners']):
            nDens = inFile['nDensCollisionPartners'][i]
            if nDens < 1e-4 or nDens > 1e12 :
                self.status |= self.FLAGS['ERROR']
                
        # checking for correct range for the species column density
        if inFile['molnDens'] < 1e5 or inFile['molnDens'] > 1e25:
            self.status |= self.FLAGS['ERROR']
              
        if self.status & self.FLAGS['ERROR']:
            return self.FLAGS['ERROR']
        else:
            self.status |= self.FLAGS['PARMSOK']
            return self.FLAGS['PARMSOK']

    # run radex with the set input 
    def run(self, checkInput = None ):
                        
        if checkInput == True:
            self.checkParameters()
        else:
            self.status |= self.FLAGS['PARMSOK']    
            
        if self.status &  self.FLAGS['PARMSOK']:
        
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
            self.parseOutput()
            self.status |= self.FLAGS['RUNOK']
        
            if self.nIter == 10000:  # ;;; get this from radex.inc
                self.status |= self.FLAGS['WARNING']
                self.status |= self.FLAGS['ITERWARN']
            else:
                self.status |= self.FLAGS['SUCCESS']
            
        return self.status
        
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
            
            #print self.rawOutput
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

    # set the axes and figure objects
    def setAxes(self, fig, axs):
        self.axs = axs
        self.fig = fig
        
    # sets up the object to plot the radex output ( nx rows and four coluns)
    def setupPlot(self, nx = None, fig = None, axs = None):
        
        self.nx = nx
        
        if fig == None and axs == None:
            self.makeAxes(nx)
        else:
            self.setAxes(fig, axs)
            
        return (self.fig, self.axs)
            
    
    def makeAxes(self, nx = None):

        if nx == None:
            nx = 1

        # axs[0,0] is the one in the top left corner
        fig, axs = pyl.subplots(4, nx, sharex = False, sharey = False, figsize=(8,8) )
        pyl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.0, hspace=0.0)
        self.setAxes(fig, axs)
        
    # plots the output of a certain model in a certain column in a predefined figure
    def plotModelInFigureColumn(self, Jall=None, inAxes=None, title=None):
                                
        nJ = len(Jall)
        #----------------flux-------------------------
        axes = inAxes[0]
        axes.lines = []
        
        xPlot = Jall
        yPlot = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis = transition['fluxcgs'] 
            
            yPlot[i] = yThis
            
        axes.semilogy(xPlot, yPlot, 'b')
        axes.axis([np.min(Jall), np.max(Jall), 1e-10, 1e-3])

        #---------------Tex and T_R--------------------------
        axes = inAxes[1]
        axes.lines = []        

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
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.semilogy(xPlot, yPlot2, 'r')
        axes.axis([np.min(Jall), np.max(Jall), 1, 10000])

        #------------------optical depth----------------------------------
        axes = inAxes[2]
        axes.lines = []        

        xPlot = Jall
        yPlot = np.ndarray(nJ, dtype=np.float64)
        for i in np.arange(nJ):
    
            Jthis = Jall[i]
    
            xPlot[i] = Jthis
    
            transition = self.getTransition( Jthis )
            yThis = transition['tau'] 
    
            yPlot[i] = yThis
        axes.plot(xPlot, yPlot, 'b')
        axes.axis([np.min(Jall), np.max(Jall), -1, np.max(yPlot)])

        #------------------population densities-----------------------------
        axes = inAxes[3]
        axes.lines = []        

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
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.semilogy(xPlot, yPlot2, 'r')
        axes.axis([np.min(Jall), np.max(Jall), 1e-10, 1])
    
    # set the appropriate labels of all the axes
    def setLabels(self):
        
        axs = self.axs
        def removeAll_xLabels(ax):
            for tick in ax.axes.xaxis.get_major_ticks():
                tick.label1On = False
        def removeAll_yLabels(ax):
            for tick in ax.axes.yaxis.get_major_ticks():
                tick.label1On = False
                        
        if self.nx == 1:
            axsLeft = axs.tolist()
            axsBotm = [axs[3]]
        else:
            axsLeft = axs[:,0].tolist()
            axsBotm = axs[3,:].tolist()

            # removing all y labels
            for axRow in axs[0:-1]:
                for ax in axRow[1:]:
                    removeAll_xLabels(ax)
                    removeAll_yLabels(ax)
        
        # removing the xlabeles of from the left axes            
        for ax in axsLeft[0:-1]:
            removeAll_xLabels(ax)
        # removing the xlabeles of from the left axes            
        for ax in axsBotm[1:]:
            removeAll_yLabels(ax)
                    
        # plotting the ylabels of the axes            
        axsLeft[0].set_ylabel('Flux')
        axsLeft[1].set_ylabel('Tex, Trot')
        axsLeft[2].set_ylabel('tau')
        axsLeft[3].set_ylabel('pop dens')
        # plotting the xlabels
        for ax in axsBotm:
            ax.set_xlabel('J')

        # removing the firt and last labels from the bottom and left axes
        for ax in axsBotm + axsLeft:
            xticks = ax.axes.xaxis.get_major_ticks()
            yticks = ax.axes.yaxis.get_major_ticks()
            for tick in [ xticks[0], xticks[-1], yticks[0], yticks[-1] ]:
                tick.label1On = False

        #removing all labels from axes with no lines (only for multiple models)
        if self.nx > 1:
            for ax in axsBotm:
                if len(ax.lines) == 0:
                    removeAll_xLabels(ax)
                    ax.set_xlabel('')
    
    # plotting a single model
    def plotModel(self):
        self.setupPlot(nx = 1)
        self.plotModelInFigureColumn(Jall=np.arange(10)+1, inAxes = self.axs, title='')
        self.setLabels()   
        pyl.show()
        