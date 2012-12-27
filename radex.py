import numpy as np
import pylab as pyl
import subprocess
import logging, sys

class radex( ):
    """A wrapper class which runs radex parses its output into objects. see test_radex.py
    and radexView.py for examples on using this class. Upon a successull executation, all the
    transitions are stored in self.transitions.
     
    :todo: docment the order in which things must be called.
    
    :warning: make sure that the lower and upper string names for the transitions are not longer than
     10 character, or change the length accordinglt in self.generateTransitionDtype().
    
    :test: test_radex.py, radexView.py
    """

    def __init__(self, execPath, molDataDir):
        """Initializes some of the instance variables needed to run radex."""
        
        self.execPath     = execPath   #: string : The path to the radex executable 
        self.molDataDir   = molDataDir #: string : Path to the directory containing the transition data for the species
        self.molDataFiles =  { 'CO'   : 'co.dat'       ,
                               '13CO' : '13co.dat'     ,
                               'HCO+' : 'hco+@xpol.dat',
                               'HCN'  : 'hcn.dat'      ,
                               'HNC'  : 'hnc.dat'      ,
                               'CS'   : 'cs@xpol.dat'  ,
                               'C'    : 'catom.dat'    ,
                               'C+'   : 'c+.dat'       ,
                               'O'    : 'oatom.dat'    ,
                               'SiO'  : 'sio.dat'      ,
                               'CN'   : 'cn.dat'       }
        """ dict : dictonary for the species files, the keys of this dict are as follows 
          (see radex.py for an example) :
          
          SPECIE_STRING : FILENAME

        .. literalinclude:: radex.py
           :lines: 24-34
           :linenos:

        .. todo:: implement a method to generate this dict automatically from all the files
                  in the directory :data:`molDataFiles`.

        .. todo:: also implement a dict for the collision partners like 'e-' : 'e', but it might not
                  be necessary since RADEX takes both 'e' and 'e-'
        """
        self.inFile = None     
        """
        dict : dictionary holding all the input parameters. It is used to construct the input parameter file that is piped to the radex executable. It should be of the form.
          
          .. code-block:: python
             :linenos:

             inFile = { 'specStr'                : 'CO'        ,  
                        'freqRange'              : [0, 50000]  ,
                        'tKin'                   : 10.0        ,
                        'collisionPartners'      : ['H2']      ,
                        'nDensCollisionPartners' : [1e3]       ,
                        'tBack'                  : 2.73        ,
                        'molnDens'               : 1e14        , # column density of the species 'specStr'
                        'lineWidth'              : 1.0         ,
                        'runAnother'             : 1           }
             
             
        All of the items in this dict correspond to the same parameters in radex except:
        'collisionPartners' and 'nDensCollisionPartners', 'molDataDir' and 'specStr'. 
        'collisionPartners' is a list of strings of species such as 'H2', or 'H' :

        .. code-block:: python

           'collisionPartners' : ['H2'] 

        or

        .. code-block:: python

           'collisionPartners' : ['H2','H+','e-']

        'nDensCollisionPartners' should be a list of the same length holding the 
        number density of each collision partner respectively. For example :
        
        .. code-block:: python
          
           'nDensCollisionPartners' : [1e3]

        or
           
        .. code-block:: python
        
           'nDensCollisionPartners' : [1e3, 1e1, 1e-2]
           
        the parameter molData in the input parameter file is constructed by appending
        the value of 'specStr' to self.molDataDir 
        """
        self.nCollPart    = None
        """integer : number of collision partners. This is the length of the list self.inFile['collisionPartners']"""
     
        self.proccess     = None     
        """subproccess.popen object : Object used to communicate with the executable"""
         
        self.rawOutput    = None
        """string : the output of the run which would be dumped by RADEX when ran standalone"""
         
        self.warnings     = None
        """list of strings : A list of strings containing the dumped by radex, if any. If there are
        any warnings, the appropriate flag is set in self.#FLAGS
        """

        self.nIter        = None
        """integer : number of iterations used when the run finishes."""

        self.outputHdr    = None  
        """string : the header of the raw output. This is useful for inspecting wheather the 
        input parameters were constructed properly.
        """
         
        self.nTransitions = None
        """the number of transitions"""
        
        self.transitions = None
        """an ndarray of dtype self.transitionFormat holding all the info of the transitions. The dtype
           has the following keys :
        
        .. code-block:: python
            
                 'upper'    : string          # the upper level of the transition     
                 'lower'    : string          # the lower level of the transition 
                 'E_up'     : numpy.float64   # Enery of the upper level 
                 'Tex'      : numpy.float64   # computed excitation temperatur
                 'tau'      : numpy.float64   # computed optical depth
                 'T_R'      : numpy.float64   # computed rotational temperature
                 'pop_up'   : numpy.float64   # computed population density in the upper level
                 'pop_down' : numpy.float64   # computed population density in the lower level
                 'fluxKkms' : numpy.float64   # computed flux in $K.km^{-1}.s^{-1}$
                 'fluxcgs'  : numpy.float64   # computed flux in cgs
        """
        
        self.transitionDtype = self.generateTransitionDtype()
        """holds the numpy.dtype of the transitions. See self.generateTransitionDtype().
        """
        
        self.FLAGS       = {'DEFAULT': 0x000, 'ERROR'  : 0x001, 'PARMSOK' : 0x002, 'RUNOK': 0x004,   
                            'SUCCESS': 0x008, 'WARNING': 0x010, 'ITERWARN': 0x020}   
        """dict : Flags set when running radex which can be used to examing the output. The flags are :
         
           .. code-block:: python
           
             'DEFAULT'  : Default status, should be called before each run 
             'ERROR'    : failed 
             'PARMSOK'  : parameters are ok 
             'RUNOK'    : radex did all the iterations (but not NECCESARILY converged)
             'SUCCESS'  : succeeded with no major warning (output should make sense)
             'WARNING'  : success with warnings, warning are assigned to self.warnings.  
             'ITERWARN' : number of iterations warning
             
        See the source for the value of each flag
        """ 
                            
        self.status = self.FLAGS['DEFAULT']
        """intger : the default stauts set from self.FLAGS. The flags are set in a bitwise fashion.
        To check if a flag is set, it can be done in the following way :
        
        .. code-block:: python
        
           stt = self.getStatus() & self.FLAGS['FLAG']
           
        if the flag is set stt would be a number greater than 0, zero otherwise
        (for more details on the values of the flags, see radex.py)
        """
        
        # plotting attributes
        self.fig = None
        """matplotlib.figure object : when set, the output is plotted in this figure"""
        self.axs = None
        """nump.ndarry matplotlib.axes object : run output are plotted in these axes objects
        for a single model self.axs has the shape (4,). When nx is larger than 1, the shape
        of self.#axs is  (4,nx). In other words, self.#axs is the object returned by
        self.fig, self.axs = pylab.subplots(4, nx )
        """
        self.nx  = None
        """integer : number of models to run and plot"""
        self.ny  = 4
        """integer : number of horizontal division in the figure (default)"""
        self.logger = None
        
    def setInFileParm(self, parm, value):
        """This method sets values to the parameters to be passed to radex.

        :param parm: A key to be modified or added in the self.inFile dict.
        :type name: str.
        :param value: The value to be set to the key.
        :type value: arbitrary [depends on the key to be modified, see :data:`inFile` documentation]
        """
        self.inFile[parm] = value
    def getInFileParm(self, parm):
        """returns a parameter from self.inFile given the key in the dict as the argument 'parm'
           which is a string used to extract the vlaue from the dict.
        """
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
    def set_logger(self, logger):
        self.logger = logger
        
    def setDefaultStatus(self):
        """sets the default status and self.transitions, self.rawOutput and self.nTransitions to None"""
        self.status = self.FLAGS['DEFAULT']
        self.transitions = None
        self.nTransitions = None
        self.rawOutput = None
        
    def getTransition(self, idx):
        """This method can be used to extract individual transition information from the self.#transitions 
        list. For example, for CO, getTransition(0) would return the transition info for the 1-0 
        transition
        
        :param int32 idx: The index of the transition in the transition list. must be between 0 and len(self.#transitions) 
        :return: dict :data:`transitions` item
        """
        return self.transitions[idx]
    
    def getWarnings(self):
        return self.warnings
    def getNIter(self):
        return self.nIter
    def printSetFlags(self):
        """prints in a pretty way all the flags (whether they are set or not) in self.status"""
        print 'FLAG_NAME     SET'
        print '------------------'
        for i, key in enumerate(self.FLAGS):
            print '%-9s  %r' % (key, (self.getStatus() & self.FLAGS[key]) >= 1)
        
        
    def setAxes(self, fig, axs):
        """set the axes and figure objects"""

        self.axs = axs
        self.fig = fig

    def filterColliders(self):
        """removes colliders from the dictionary of the self.#inFile if their density 
           is less than the minimum accepted value by RADEX (for now, this ie 1e-3 cm^{-3} )
        """
        # removing the colliders which have abundances less than the one range 
        # that radex accepts
        for (i, nDense) in enumerate(self.inFile['nDensCollisionPartners']):
            if nDense <= 1e-3:
                print 'poped ', self.inFile['nDensCollisionPartners'][i], self.inFile['collisionPartners'][i]
                self.inFile['nDensCollisionPartners'].pop(i)
                self.inFile['collisionPartners'].pop(i)
                 
    def genInputFileContentAsStr(self):
        """generates the parameter file contents from self.inFile, which can be passed to the radex executable, as a string. 
        
        :return: (str)
        """
        
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

    def checkParameters(self):
        """chech whether the contents of self.inFile are within the ranges where RADEX can work.
        
        :raise:  exception NameException
        :attention: This method sets the value of self.status. The flag 'PARMSOK' is set if the paramteres are ok, 'ERROR' is set otherwise.
        """
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

    def run(self, checkInput = None, verbose = None ):
        """run the radex executable.
        
        :param bool checkInput: By default this is False. In this case, the input 
          parameters are not checked and the flag 'PARMSOK' (see #FLAGS) is set to :data:`status`. Otherwise, set this to True to force a paremter input check.
        :param bool verbose: By default this is False. Set it to True to to write 
          the raw output to stdout.
        :return: (int) :data:`status`. Upon a successful run, 'RUNOK' and 'SUCCESS' 
          flags are set. If the number of iterations excceeds 10,000 'RUNOK',
          'WARNING' and 'ITERWARN' flags are set. if the 'PARMSOK' is true 
          :data:`rawOutput` and :data:`transitions` are set, otherwise they 
          remaine None.
        :todo: extract other warnings also, not just the one due to the max 
          iterations.
        :warning: when running the same radex instance multiple times, make sure 
          to set the status to the default before calling :data:`run` using 
          :data:`setDefaultStatus()`.
        """
        
        self.warnings = []

        if checkInput == True:
            self.checkParameters()
        else:
            self.status |= self.FLAGS['PARMSOK']    

        if self.status &  self.FLAGS['PARMSOK']:
        
            radexInput = self.genInputFileContentAsStr()
            
            if verbose == True:
                print '---------------------------------------------------------------------------------------'
                print '----------------radex input parameters passed to the executable------------------------'
                print '---------------------------------------------------------------------------------------'
                print radexInput
                print '---------------------------------------------------------------------------------------'
            
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
        else:
            if verbose == True:
                print 'radex.py : parameters NOT ok. radex.status = ', self.status 
            
        return self.status
        
    def generateTransitionDtype(self):
        """generates the trasition dtype which will be used in assigning the transition info in self.transitionsNdArry.
        """
        fmt = [
               ('upper'   , np.str_, 10),   
               ('lower'   , np.str_, 10),   
               ('E_up'    , np.float64),
               ('Tex'     , np.float64), 
               ('tau'     , np.float64), 
               ('T_R'     , np.float64), 
               ('pop_up'  , np.float64), 
               ('pop_down', np.float64), 
               ('fluxKkms', np.float64), 
               ('fluxcgs' , np.float64), 
              ] # 2 * 10b + 8 * 8b = 84b per transition 
       
        return np.dtype(fmt)
    
    def flagSet(self, flag):
        """returns True 'flag' is set, false otherwise"""
        return self.status & self.FLAGS[flag]
        
    def parseOutput(self):
        """Once radex exectues and dumps transition information, this method is used
          to extract the line data.
         
        :return: None. The instance variable self.transitions is set.
        :note: in parsing the output, the transitions info is between the lines 
         containing the units:
          
             "(K)    (GHz) ..."
             
         and
         
             "Another calculation"
              
         thisway, we can get the number of transitions (the number of new lines, 
         the work on parsing the data without the need to append anythign to a 
         list..just preallocate the dtype and fill in the values.
        """
        
        try:
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
                info['upper'   ] = lineSplt[0].strip()
                info['lower'   ] = lineSplt[1].strip()
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
                
            #copying the content of self.transitions into self.transitionsNdarrya
            ##:note: do this at one go without storing things first in self.transitiosn
            #******************************************************
            self.nTransitions = len(transitions)
            transitionsNdarray = np.ndarray((self.nTransitions), dtype = self.transitionDtype)
            for i, trans in enumerate(transitions):
                for key in trans.keys():
                    transitionsNdarray[i][key] = trans[key]
            self.transitions = transitionsNdarray
            #******************************************************
            
            
        except (ValueError, IndexError):
            print self.rawOutput
            errStr = 'parsing the output failed'
            raise NameError(errStr)
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

        #plotting the intensities    
        axes.semilogy(xPlot, yPlot, 'b')
        axes.axis([np.min(allTrans), np.max(allTrans), 1e-10, 1e-1])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        ##np.savetxt('/home/mher/ism/tmp/intensities.out', np.array([xPlot, yPlot]).T, '%e')
        ##print '-----------> saved the file /home/mher/intensities.out'
        
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
        ##np.savetxt('/home/mher/ism/tmp/tau.out', np.array([xPlot, yPlot]).T, '%f')
        ##print '-----------> saved the file /home/mher/tau.out'
        
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

    def clearCurves(self):
        """Clears the curves in an axes column. This is usually used when radex fails
        and the data in a certain column need to be removed. Here it is assumes there
        is one column."""

        for ax in self.axs:            
            ax.lines = []
        
        
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

            # removing all x and y labels from axes to the left
            # of th zeroth column and above the bottom row
            for axRow in axs[0:-1]:
                for ax in axRow[1:]:
                    removeAll_xLabels(ax)
                    removeAll_yLabels(ax)
        
        #removing the x labeles of from axes on the zeroth column
        #excluding the one in the bottom one (lower-left corner)            
        for ax in axsLeft[0:-1]:
            removeAll_xLabels(ax)
        
        #removing the y labeles of from axes on the bottom row, to the right of the 
        #one on the lower-left corner
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
        if self.nx > 1:
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
        
    def clear(self):
        self.nCollPart    = None
        self.rawOutput    = None
        self.warnings     = []
        self.nIter        = None
        self.outputHdr    = None  
        self.nTransitions = None
        self.transitions = None
        self.status = None

    def setupLogger(self):
        """sets up the logger which will prepend info about the printed stuff. Assignes a value to self.logger."""
        # setting up the logger                                                                                                                                                                                                              
        # create logger                                                                                                                                                                                                                      
        self.logger = logging.getLogger('simple_example')
        self.logger.setLevel(logging.DEBUG)

        # create console handler and set level to debug                                                                                                                                                                                      
        ch = logging.StreamHandler( sys.stdout )  # setting the stream to stdout                                                                                                                                                             
        ch.setLevel(logging.DEBUG)

        # create formatter                                                                                                                                                                                                                   
        formatter = logging.Formatter('[%(asctime)s %(funcName)s() %(filename)s:%(lineno)s] %(message)s') # this was the original in the example                                                                              

        # add formatter to ch                                                                                                                                                                                                                
        ch.setFormatter(formatter)

        # add ch to logger                                                                                                                                                                                                                   
        self.logger.addHandler(ch)

        return self.logger
