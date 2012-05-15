## @file radex.py
#  implementation of the radex interface in python

import numpy as np
import pylab as pyl
import subprocess

## @package Radex
# Radex   
# @author: Mher V. Kazandjian
# @date  : 2012-Mar-05
# @version: xxx
# @section Radex This is a python module through which RADEX can be called from python.
# Plotting features are also included to visualize the output. Minor modifications
# for the RADEX source are required to make it easier to the python module to parse
# the output. Otherwise, the numerics are left intact. See 
# www.strw.leidenuniv.nl/~moldata/radex.html 
# for more details on RADEX.
# @section Examples 
# A simple example on using this package/class is test_radex.py. A less trivial example
# is radexView.py

## Radex class. 
#  A wrapper class which runs radex parses its output into objects. see test_radex.py
#  and radexView.py for examples on using this class. 
#  @todo : docment the order in which things must be called
class radex( ):
    
    ## Initializes some of the instance variables needed to run radex.
    #  @param execPath (String) The path to the radex executable
    #  @param molDataDir (String) The Path to the directory containing the transition data for the species
    def __init__(self, execPath, molDataDir):
    
        ## string : Path to the radex executable 
        self.execPath     = execPath 
        ## string : Path to the directory containing the transition data for the species
        self.molDataDir   = molDataDir
        ## dict : dictonary for the species files, the keys of this dict are as follows 
        #  (see radex.py for an example) :\n
        #  SPECIE_STRING : FILENAME
        #  @todo implement a method to generate this dict automatically from all the files
        #  in the directory molDataFiles\n
        #  also implement a dict for the collision partners like 'e-' : 'e', but it might not
        #  be necessary since RADEX takes both 'e' and 'e-'
        self.molDataFiles =  { 'CO'   : 'co.dat'       ,
                               '13CO' : '13co.dat'     ,
                               'HCO+' : 'hco+@xpol.dat',
                               'HCN'  : 'hcn.dat'      ,
                               'HNC'  : 'hnc.dat'      ,
                               'CS'   : 'cs@xpol.dat'  ,
                               'C'    : 'catom.dat'    ,
                               'C+'   : 'c+.dat'       ,
                               'O'    : 'oatom.dat'    ,
                               'CN'   : 'cn.dat'       }
        ## dict : dictionary holding all the input parameters. It is used to construct the input
        #  parameter file that is piped to the radex executable. It should be of the form :\n
        #  \code
        #  inFile = { 'specStr'                : 'CO'        ,
        #             'freqRange'              : [0, 50000]  ,
        #             'tKin'                   : 10.0        ,
        #             'collisionPartners'      : ['H2']      ,
        #             'nDensCollisionPartners' : [1e3]       ,
        #             'tBack'                  : 2.73        ,
        #             'molnDens'               : 1e14        ,
        #             'lineWidth'              : 1.0         ,
        #             'runAnother'             : 1           }
        #  \endcode
        #  All of the items in this dict correspond to the same parameters in radex except
        #  'collisionPartners' and 'nDensCollisionPartners', 'molDataDir' and 'specStr'. 
        #  'collisionPartners' is a list of strings of species such as 'H2', or 'H' :\n
        #  \code 'collisionPartners' : ['H2'] \endcode
        #  or        
        #  \code 'collisionPartners' : ['H2','H+','e-'] \endcode
        #  'nDensCollisionPartners' should be a list of the same length holding the 
        #  number density of each collision partner respectively. For example :
        #  \code 'nDensCollisionPartners' : [1e3] \endcode
        #  or        
        #  \code 'nDensCollisionPartners' : [1e3, 1e1, 1e-2] \endcode
        #  the parameter molData in the input parameter file is constructed by appending
        #  the value of 'specStr' to self.#molDataDir 
        self.inFile       = None     
        ## integer : number of collision partners. This is the length of the list self.inFile['collisionPartners']
        self.nCollPart    = None     
        ## subproccess.popen object : Object used to communicate with the executable
        self.proccess     = None     
        ## string : the output of the run which would be dumped by RADEX when ran standalone 
        self.rawOutput    = None 
        ## list of strings : A list of strings containing the dumped by radex, if any. If there are
        #  any warnings, the appropriate flag is set in self.#FLAGS
        self.warnings     = None
        ## integer : number of iterations used when the run finishes.
        self.nIter        = None
        ## string : the header of the raw output. This is useful for inspecting wheather the 
        #  input parameters were constructed properly.
        self.outputHdr    = None  
        ## dict list : A list containing all the information of the computed lines.
        #  each item in the list is a dictionary with the following keys 
        #  (they are extracted from self.rawOutput) : 
        #    \verbatim
        #    'upper'    : string (the upper level of the transition) 
        #    'lower'    : string (the lower level of the transition) 
        #    'E_up'     : numpy.float64 ( Enery of the upper level ) 
        #    'Tex'      : numpy.float64 ( computed excitation temperatur )
        #    'tau'      : numpy.float64 ( computed optical depth )
        #    'T_R'      : numpy.float64 ( computed rotational temperature )
        #    'pop_up'   : numpy.float64 ( computed population density in the upper level )
        #    'pop_down' : numpy.float64 ( computed population density in the lower level )
        #    'fluxKkms' : numpy.float64 ( computed flux in \f$K.km^{-1}.s^{-1}\f$ )
        #    'fluxcgs'  : numpy.float64 ( computed flux in cgs )
        #    \endverbatim 
        self.transitions = None
        ## dict : Flags set when running radex which can be used to examing the output. The flags are :
        # \verbatim
        #  'DEFAULT'  : Default status, should be called before each run 
        #  'ERROR'    : failed 
        #  'PARMSOK'  : parameters are ok 
        #  'RUNOK'    : radex did all the iterations (but not neccesarly converged)
        #  'SUCCESS'  : succeeded with no major warning (output should make sense)
        #  'WARNING'  : success with warnings, warning are assigned to self.warnings 
        #  'ITERWARN' : number of iterations warning
        # \endverbatim See the source for the value of each flag
        self.FLAGS       = {'DEFAULT': 0x000, 'ERROR'  : 0x001, 'PARMSOK' : 0x002, 'RUNOK': 0x004,   
                            'SUCCESS': 0x008, 'WARNING': 0x010, 'ITERWARN': 0x020}   
                            
        ## intger : the default stauts set from self.FLAGS. The flags are set in a bitwise fashion.
        #  To check if a flag is set, it can be done in the following way :
        #  \code self.getStatus() & self.FLAGS['FLAG'] \endcode
        #  if the flag is set, it returns a number greater than 0, if not 
        #  it returns a zero. (for more details on the values of the flags, see radex.py)
        self.status      = self.FLAGS['DEFAULT']
        # plotting attributes
        ## matplotlib.figure object : when set, the output is plotted in this figure
        self.fig = None
        ## nump.ndarry matplotlib.axes object : run output are plotted in these axes objects
        # for a single model self.axs has the shape (4,). When nx is larger than 1, the shape
        # of self.#axs is  (4,nx). In other words, self.#axs is the object returned by
        # self.fig, self.axs = pylab.subplots(4, nx )
        self.axs = None
        ## integer : number of models to run and plot
        self.nx  = None
        ## integer : number of horizontal division in the figure (default)
        self.ny  = 4
        
    ## sets keys in the self.inFile dict
    #  @param parm  (String) A key to be modified or added in the self.inFile dict
    #  @param value (abitrary) The value to be set to the key  
    def setInFileParm(self, parm, value):
        self.inFile[parm] = value
    ## returns a parameter from self.inFile given the key in the dict
    #  @param parm (str) : the string used to extract the vlaue from the dict
    def getInFileParm(self, parm):
        return self.inFile[parm]
    def setInFile(self, inFile):
        self.inFile = inFile
    def getInFile(self):
        return self.inFile
    def getStatus(self):
        return self.status
    def getRawOutput(self):
        return self.rawOutput
    def setStatus(self, status):
        self.status = status
    ## sets the default status
    def setDefaultStatus(self):
        self.status = self.FLAGS['DEFAULT']
    ## This method can be used to extract individual transition information from the self.#transitions 
    #  list. For example, for CO, getTransition(0) would return the transition info for the 1-0 
    #  transition
    #  @param idx The index of the transition in the transition list. must be between 0 and len(self.#transitions) 
    #  @return (dict:self.#transitions item)
    def getTransition(self, idx):
        return self.transitions[idx]
    
    def getWarnings(self):
        return self.warnings
    def getNIter(self):
        return self.nIter

    ## set the axes and figure objects
    def setAxes(self, fig, axs):
        self.axs = axs
        self.fig = fig

    ## removes colliders from the dictionary of the self.#inFile if their density 
    #  is less than the minimum accepted value by RADEX (for now, this ie 1e-3 cm^{-3} )
    def filterColliders(self):
        # removing the colliders which have abundances less than the one range 
        # that radex accepts
        for (i, nDense) in enumerate(self.inFile['nDensCollisionPartners']):
            if nDense <= 1e-3:
                print 'poped ', self.inFile['nDensCollisionPartners'][i], self.inFile['collisionPartners'][i]
                self.inFile['nDensCollisionPartners'].pop(i)
                self.inFile['collisionPartners'].pop(i)
                 
                
    ## generates the parameter file contents from self.inFile, which can be passed to the 
    #   radex executable, as a string. 
    #   @return: (str)
    def genInputFileContentAsStr(self):
        self.nCollPart = len(self.inFile['collisionPartners'])
        
        strng = ''

        strng  += '%s/%s\n' % (self.molDataDir, self.molDataFiles[ self.inFile['specStr'] ])
        strng  += 'foo\n'
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
        #print '------------\n%s\n-----------------\n' % strng
        return strng

    ## chech whether the contents of self.inFile are within the ranges where
    #  RADEX can work.  
    #  @raise exception: NameException
    #  @attention This method sets the value of self.status. The flag 'PARMSOK' is set
    #  if the paramteres are ok, 'ERROR' is set otherwise
    def checkParameters(self):
        inFile = self.inFile
        
        for item in inFile:
            if inFile[item] == None:
                strng = 'Error : missing input parameter %s ' % (item)
                raise NameError(strng)
            
        # checking for correct range for the kinetic temperature
        if inFile['tKin'] <= 0.1 or inFile['tKin'] >= 1e4 :
            self.status |= self.FLAGS['ERROR']
            
        # checking for correct range for the densities of the collion partners
        for i, collPartner in enumerate(inFile['collisionPartners']):
            nDens = inFile['nDensCollisionPartners'][i]
            if nDens < 1e-3 or nDens > 1e12 :
                self.status |= self.FLAGS['ERROR']
                
        # checking for correct range for the species column density
        if inFile['molnDens'] < 1e5 or inFile['molnDens'] > 1e25:
            self.status |= self.FLAGS['ERROR']
              
        if self.status & self.FLAGS['ERROR']:
            return self.FLAGS['ERROR']
        else:
            self.status |= self.FLAGS['PARMSOK']
            return self.FLAGS['PARMSOK']

    ## run the radex executable.
    #  @param checkInput (bool): By default this is False. In this case, the input
    #  parameters are not checked and the flag 'PARMSOK' (see #FLAGS) is set to self.#status. 
    #  Otherwise, set this to True to force a paremter input check. 
    #  @param verbose (bool)   : By default this is False. Set it to True to
    #  to write the raw output to stdout 
    #  @return (int) #status \n upon a successful run, 'RUNOK' and 'SUCCESS' flags are set. If the number of iterations
    #  excceeds 10,000 'RUNOK','WARNING' and 'ITERWARN' flags are set. if the 'PARMSOK' is true
    #  self.#rawOutput and self.#transitions are set, otherwise they remaine None
    #  @todo: extract other warnings also, not just the one due to the max iterations
    #  @warning: when running the same radex instance multiple times, make sure to set the 
    #  status to the default before calling self.#run using self.#setDefaultStatus()
    def run(self, checkInput = None, verbose = None ):
        
        if checkInput == True:
            self.checkParameters()
        else:
            self.status |= self.FLAGS['PARMSOK']    

        if self.status &  self.FLAGS['PARMSOK']:
        
            radexInput = self.genInputFileContentAsStr()
            if verbose == True:
                print '---------------------------------------------'
                print radexInput
                print '---------------------------------------------'
            
            self.proccess = subprocess.Popen(self.execPath           , 
                                             stdin=subprocess.PIPE   ,  
                                             stdout=subprocess.PIPE  ,  
                                             stderr=subprocess.PIPE  )
            radexOutput = self.proccess.communicate(input=radexInput)[0]
            
            self.rawOutput = radexOutput
            self.parseOutput()
            self.status |= self.FLAGS['RUNOK']
        
            if self.nIter == 10000:  # ;;; get this from radex.inc
                self.status |= self.FLAGS['WARNING']
                self.status |= self.FLAGS['ITERWARN']
            else:
                self.status |= self.FLAGS['SUCCESS']
            
            if verbose == True:
                print '-------------raw output of Radex stdout------------------'
                print radexOutput
                print '---------------------------------------------------------------------'
                
        return self.status
        
    ## Once radex exectues and dumps transition information, this method is used
    #  to extract the line data.
    #  @return None\n The instance variable self.#transitions is set
    def parseOutput(self):
        #print self.rawOutput
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
            info['upper'   ] = lineSplt[0] 
            info['lower'   ] = lineSplt[1] 
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
        
    # -------------------------- Plotting Methods ---------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    
    ## sets up the object to plot the radex output 
    #  @param nx (integer) the number of vertical columns to setup (one model per column)
    #  @param fig ( matplotlib.figure object ) Set this object to use a figure which
    #  is pre-defined, otherwise, a new figure object is created. If this is not set
    #  the axs keyword is igonored and a self.#axs object is created.  
    #  @param axs ( numpy.ndarray of matplotlib.axes object ) The axes where the models will
    #  be plotted. Should have dimensions nx x 4
    #  @return  (self.#fig, self.#axs)\n in this method, also self.nx is set to nx
    def setupPlot(self, nx = None, fig = None, axs = None):
        
        self.nx = nx
        
        if fig == None and axs == None:
            self.makeAxes(nx)
        else:
            self.setAxes(fig, axs)
            
        return (self.fig, self.axs)
            
    ## creates and assings the self.#fig and self.#axs attributes
    #  @param nx (see #setupPlot)  
    def makeAxes(self, nx = None):

        if nx == None:
            nx = 1

        # axs[0,0] is the one in the top left corner
        fig, axs = pyl.subplots(4, nx, sharex = False, sharey = False, figsize=(8,8) )
        pyl.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.95, wspace=0.0, hspace=0.0)
        self.setAxes(fig, axs)
        
    ## plots the output of a certain model in a certain column in a predefined figure
    #  @param allTrans (integer list or numpy.int32) The indicies of the transitions in self.#transitions whose
    #  data should be plotted, ex : numpy. 
    #  @param inAxes ( nump.ndarray matplorlib.axes ) The axes colomn in which the model info will
    #  be plotted. Fluxes are plotted in inAxes[0], T_ex, T_rot are plotted in inAxes[1], 
    #  optical depth in inAxes[2] and population densiities in inAxes[3]
    #  @param title (string) The title to be written at the top of the axes column  
    def plotModelInFigureColumn(self, allTrans = None, inAxes = None, title = None):
        
        if allTrans == None:
            allTrans =  np.arange( len(self.transitions) )
            
        nTrans = len(allTrans)
        #----------------flux-------------------------
        axes = inAxes[0]
        axes.lines = []
        xticksStrs = ()
        
        xPlot = allTrans
        yPlot = np.ndarray(nTrans, dtype=np.float64)
        for i in np.arange(nTrans):
    
            thisTrans = allTrans[i]
            xPlot[i] = thisTrans
    
            transition = self.getTransition( thisTrans )
            yThis = transition['fluxcgs'] 
            
            yPlot[i] = yThis
            
            # construncting the latex string of the transitions
            upperStr = transition['upper'].split('_')
            if len(upperStr) == 2:
                upperStr = '%s_{%s}'% (upperStr[0], upperStr[1])
            else:
                upperStr = upperStr[0]
            lowerStr = transition['lower'].split('_')
            if len(lowerStr) == 2:
                lowerStr = '%s_{%s}'% (lowerStr[0], lowerStr[1])
            else:
                lowerStr = lowerStr[0]

            
            thisTickStr = '$%s-%s$' %( upperStr, lowerStr )
            xticksStrs = xticksStrs + (thisTickStr ,)
            
        axes.semilogy(xPlot, yPlot, 'b')
        axes.axis([np.min(allTrans), np.max(allTrans), 1e-10, 1e-3])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        

        #---------------Tex and T_R--------------------------
        axes = inAxes[1]
        axes.lines = []        

        xPlot = allTrans
        yPlot1 = np.ndarray(nTrans, dtype=np.float64)
        yPlot2 = np.ndarray(nTrans, dtype=np.float64)
        for i in np.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.getTransition( thisTrans )
            yThis1 = transition['Tex']
            yThis2 = transition['T_R'] 
    
            yPlot1[i] = yThis1
            yPlot2[i] = yThis2
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.semilogy(xPlot, yPlot2, 'r')
        axes.axis([np.min(allTrans), np.max(allTrans), 1, 10000])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        
        #------------------optical depth----------------------------------
        axes = inAxes[2]
        axes.lines = []        

        xPlot = allTrans
        yPlot = np.ndarray(nTrans, dtype=np.float64)
        for i in np.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.getTransition( thisTrans )
            yThis = transition['tau'] 
    
            yPlot[i] = yThis
        axes.plot(xPlot, yPlot, 'b')
        axes.axis([np.min(allTrans), np.max(allTrans), -1, np.max(yPlot)])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        
        #------------------population densities-----------------------------
        axes = inAxes[3]
        axes.lines = []        

        xPlot = allTrans
        yPlot1 = np.ndarray(nTrans, dtype=np.float64)
        yPlot2 = np.ndarray(nTrans, dtype=np.float64)
        for i in np.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.getTransition( thisTrans )
            yThis1 = transition['pop_up']
            yThis2 = transition['pop_down'] 
    
            yPlot1[i] = yThis1
            yPlot2[i] = yThis2
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.semilogy(xPlot, yPlot2, 'r')
        axes.axis([np.min(allTrans), np.max(allTrans), 1e-10, 1])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
    
    ## set the appropriate labels of all the axes
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
            ax.set_xlabel('Trans') #;;; use the up-down string (rotated 90 deg) to plot the trans labele

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
    
    ## plotting a single model
    #  @todo allow for a choice of the number of transitions to be plotted
    def plotModel(self):
        self.setupPlot(nx = 1)
        self.plotModelInFigureColumn( allTrans = None, inAxes = self.axs, title='')
        self.setLabels()   
        pyl.show()
        