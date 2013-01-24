import os, glob, sys, pickle
import numpy as np
import logging
import re

from time  import *
from ismUtils   import *
from radex      import *
from scipy      import interpolate
import chemicalNetwork
from scipy.interpolate import interp1d

class reader():
    """This is a class which provides an interface to the LAMDA database. 
    The files in the LAMBDA database are parsed into a dictioary of strings
    (for the names of the species) and numpy dtypes for the numerical values
    such as the rates coefficients.  The parsing is done based on the fields
    which specify the content of the following data. For example, first we
    look for the string 'MOLECULE' and extract its name form the following line
    since that is how the files are constructed.  It is not such a straight 
    forward task since sonetimes feilds which contain information also contain
    descriptive texts. For example the feild holding the name of the H20
    molecule also has :
    
           p-H2O spectroscopy from JPL ...
           
    The input units of the energy levels in the levels section is assumed to
    be in cm^-1 (which is the case in all  the files of the current version 
    of the database).
    
    It is assumed that the energy units read in the transitions (both radiative
    and collisional are in K.
    
    Creating an instance in the following way :
    
    .. code-block:: python
    
        reader = molData.reader(dirPath = '/home/mher/lambdaPath/')
        
    reads all the files with a .dat extention, parses them and stores the 
    parsed data in the attribute :data:`speciesInfo`.
    
    The method :data:`get_specie` can be used get the info for a cartain specie.
    Each files is parsed into a dictionary with the following keys:
    
    :warning: the quantum number reading is ignored for now.
    :test: Below are a few file which can be used as templates to read and use
       the data:
       
        - :file:`test_parsing_datafile.py`
        - :file:`test_leidenLambda.py`

    .. code-block:: python
            
              specDict = ['path']       #full path of the file
                         ['specStr']    #the string of the specie
                         ['info']       #the raw content of the string name feild (useful for extra checks)
                         ['weight']     #atomic weight
                         ['nlevels']    #number of levels 
                         ['levels'][0]  #an ndarray of dtypes holding info for each level
                                   [1]  #see the method *cleanSpecDict* for details of each key in the dtyp
                                    .   #the keys of the dtype are : 'n', 'E', 'g', 'qn'  
                                    . 
                                   [nlevels-1] 
                         ['nTransRad']         # number of radiative transitions
                         ['transRad'][0]       # an ndarray of dtypes holding the info of all the
                                     [1]       # radiative transitions. see the method *cleanSpecDict*
                                      .        # for details of each key in the dtyp 
                                      .        # the keys of the dtype are : 'n', 'u', 'l', 'A', 'nu', 'E'
                                     [nTransRad-1] 
                         ['transColl']['nPartners'] # number of collisional partners
                                      ['partnersList'][0] # a list of the collision partners (the raw content
                                                      [1] # of the collsion partner list, same as 'info'
                                                       .  # for each below)
                                                      [nPartners-1]  
                                      ['partner0']['nTrans'] # number of collisional transitions
                                           .      ['info']   # the raw content of the collion partner feild (useful for extra checks)
                                           .      ['nTemps'] # number of temperatures in the collisions table 
                                           .      ['temps']  # the temperatures values for which the rate coefficient are tabulated
                                           .      ['table'][0] #an ndarray holding the info of the collisional transitions with this  
                                           .               [1] #partner. The keys of the dtype are : 'n', 'u', 'l', 'rc'
                                           .                .  #'rc' is an interpolation function which given a temperature returns
                                           .                .  #the rate coefficient.
                                           .               [nTrans-1]   
                                      ['partner1']['nTrans']
                                           .      ['info']
                                           .      ['nTemps']
                                           .      ['temps']
                                           .      ['table'][0]
                                           .               [1] 
                                           .                .
                                           .                .
                                           .               [nTrans-1]
                                           .          
                                           .          
                                      ['partner(nPartners-1)']

    An example for getting a rate coefficient. Once the reader has been invoked, we
    can use the :data:`get_specie` method to extract a certain specie info. Lets say 
    we want the rate coefficient for the para-H2O molecule. 
    
    .. code-block:: python
            
        pNH3 = reader.get_specie(specStr = 'NH3', inInfo = 'p-')
        
    this retuns the species infor for para NH3. For H2O for example, it is better
    to provide the inPath parameter to decide which information to return
    based on the filename holding the data, since there are may o- and p- files
    in the current database.  
        
    :todo: add a method to return the index of a level given the quantum numbers
    :todo: add another method which feches the transition info give the upper and lower levels (all transitions -rad-col..etc..)
    :todo: add a method which dumps the graphvis stuff for the transitions (radiative, coll...)
    :todo: add a method which clears memory and check for memory leaks
    :todo: see what it is about the differences in the computed and tabulated deltaEs and energy levels
    """
    def __init__(self, *args, **kwargs):
        
        self.dataFiles = None
        """A list containing all the names of the files in :data:`dirPath`"""

        self.ignoreList = ('co@neufeld-old.dat')
        """A list of strings holding the names of the files to be ignored. Those files
        are old ones which are not combatible with the current parser of this class."""

        self.speciesInfo = ()
        """A dictionary holding the species string and the files containing data about them"""
                
        if 'dirPath' in kwargs:
            self.set_dirPath(kwargs['dirPath'])
            self.dataFiles = self.findFiles(**kwargs)
            self.collectAllFilesInfo(**kwargs) 
        else:
            self.dirPath = None
            """The path of the directory holding all the information about all the species"""
        
    def loadSpeciesData(self, specie = None):
        """Returns a numpy dtype holding all the information about the specie in the LMBDA database"""
        pass

    def collectAllFilesInfo(self, **kwargs):
        """Loops over all the data files in self.dirPath and collects the basic information
         
          - name of the species
          - file names containing info about that species
          
        The return values is a dict 'species' : [file1, file2,]   
        """ 
        
        for fPath in self.dataFiles:
            
            #ignoring proccessing file that can not be processed by default
            if fPath.split('/')[-1] in self.ignoreList:
                continue

            print fPath
            #if 'test.dat' not in fPath:  #;;;rmove this later
            #    continue                  #;;;remove this later
            #input('press 1 to process the next file %s' % fPath)
            specDict = self.parseDataFile(fPath)
            if 'specie' in kwargs:
                if kwargs['specie'] == specDict['specStr']:
                    specDict = self.cleanSpecDict(specDict)
            self.speciesInfo += (specDict,)
            
    def cleanSpecDict(self, specDict):
        """takes the output item of :data:`parseDataFile` and defines dtypes which hold 
           the info about the levels and the transitions. This method does the following
           modification to specDict:
           
             - specDict['levels'] : is replaced by a numpy dtype : 'n', 'E', 'g', 'qn'. 
             - specDict['transRad'] : is replaced by a numpy dtype : 'n', 'u', 'l', 'A', 'nu', 'E'. 
             - specDict['transColl']['partnerX'] : each collision partner X will have a new dict key ['trans']
               which is a dtype with keys 'n', 'u', 'l', 'rc'. 'rc' is an interpolation function which take 
               temperature as an argument and returns the rate coefficient for that temperature.
             - in the above : 'n' is a transition index.
             - ['levels'], ['transRad'], ['transColl']['partnerX']['dtype'] are ndarrays one entry for each level/transition.
        """ 
        
        #-----------------------------
        # extracting the energy levels
        #-----------------------------
        #            level index       energy (ev)       stat weight      level j             level k        leve inv            
        dt = np.dtype(
              [ ('n', np.int32), ('E', np.float64), ('g', np.float64), ('j', np.int32), ('k', np.int32), ('i', np.int32)]
             ) 
        levels = np.ndarray( specDict['nlevels'], dtype = dt)
        for idx, levelStr in enumerate(specDict['levels']):
            levelStrParsed = re.split('\s+|_', levelStr.strip())
            levels[idx]['n'] = np.int32(levelStrParsed[0]) - 1 
            levels[idx]['E'] = np.float64(levelStrParsed[1]) #*hPlank*cLight/kBoltz # energies in K
            levels[idx]['g'] = np.float64(levelStrParsed[2])
            """
            levels[idx]['j'] = np.int32(levelStrParsed[3])
            levels[idx]['k'] = np.int32(levelStrParsed[4])
            levels[idx]['i'] = np.int32(levelStrParsed[5])
            """
        #replacing the dict item with the numpy dtype array
        specDict.pop('levels')       
        specDict['levels'] = levels

        #-------------------------------------------------
        # extracting the radiative transitions information
        #-------------------------------------------------
        #                transition idx     upper           lower               frequency          energy      
        dt = np.dtype(
                      [ ('n', np.int32), ('u', np.int32), ('l', np.int32), ('A', np.float64), ('nu', np.float64), ('E', np.float64)]
                     ) 
        transRad = np.ndarray( specDict['nTransRad'], dtype = dt)
        for idx, transStr in enumerate(specDict['transRad']):
            transRadStrParsed = re.split('\s+', transStr.strip())
        
            transRad[idx]['n']  = np.int32(transRadStrParsed[0]) - 1 
            transRad[idx]['u']  = np.int32(transRadStrParsed[1]) - 1 
            transRad[idx]['l']  = np.int32(transRadStrParsed[2]) - 1
            transRad[idx]['A']  = np.float64(transRadStrParsed[3])
            transRad[idx]['nu'] = np.float64(transRadStrParsed[4])
            transRad[idx]['E']  = np.float64(transRadStrParsed[5])
            
        specDict.pop('transRad')       
        specDict['transRad'] = transRad

        #--------------------------------------------------------------------------
        #checking the read frequencies versus the computed ones from the transitions
        hPlank = 6.63e-27    # erg.s 
        cLight = 29979245800.0 # cm.s^-1
        kBoltz = 1.38e-16    # erg.K^-1 
        ev2erg = 1.602e-12   #erg

        totalRelErr = 0.0
        minRelErr = np.inf
        maxRelErr = 0.0
        for trans in transRad:
            lw =  trans['l']; up = trans['u'];
            # energy difference (converting from cm^-1 to K)
            dE = (levels[up]['E'] - levels[lw]['E'])*hPlank*cLight/kBoltz
            # transition frequency in GHz 
            nu = dE*kBoltz/hPlank * 1e-9
            relErr = (1.0 - trans['nu']/nu)
            minRelErr = np.minimum(abs(relErr), abs(minRelErr))
            maxRelErr = np.maximum(abs(relErr), abs(maxRelErr))
            totalRelErr += relErr
            #print relErr
        print 'avertage relative error = %e (min,max)=[%e,%e]' % (totalRelErr/specDict['nTransRad'], minRelErr, maxRelErr) 
        #--------------------------------------------------------------------------

        #-----------------------------------------------------------------------------
        # extracting the collisional transitions information (for the zeroth partner)
        #-----------------------------------------------------------------------------
        #                transition idx     upper           lower               frequency          energy      
        dt = np.dtype(
                      [ ('n', np.int32), ('u', np.int32), ('l', np.int32), ('rc', object) ]
                     )

        transColl = np.ndarray(specDict['transColl']['partner0']['nTrans'], dtype = dt)
        tKin =  specDict['transColl']['partner0']['temps']
        for idx, transCollStr in enumerate(specDict['transColl']['partner0']['table']):
            transCollStrParsed = re.split('\s+', transCollStr.strip())
            transColl[idx]['n']  = np.int32(transCollStrParsed[0]) - 1 
            transColl[idx]['u']  = np.int32(transCollStrParsed[1]) - 1
            transColl[idx]['l']  = np.int32(transCollStrParsed[2]) - 1
            # getting the interpolation function
            rc = np.array(transCollStrParsed[3:])
            
            if len(rc) == 1:
                f_rc = lambda x: rc if x == tKin else -1.0
            else:
                f_rc = interp1d(tKin, rc)
            transColl[idx]['rc'] = f_rc
            
        specDict['transColl']['partner0']['trans'] = transColl  

        return specDict
    
    def parseDataFile(self, fPath):
        """Parses the data file into its components and returns them as a nested dict object.
           The keys of the dict are :
           
           :return: (dict) specDict 
 
           .. code-block:: python
            
              specDict = ['path']
                         ['specStr']
                         ['specInfo']
                         ['weight']
                         ['nlevels']
                         ['levels']
                         ['nTransRad']
                         ['transRad']
                         ['transColl']['nPartners']
                                      ['partnersList']
                                      ['partner0']['nTrans']
                                                  ['info']   #todo: not implemented
                                                  ['nTemps']
                                                  ['temps']
                                                  ['table']
                                      ['partner1']['nTrans'] 
                                                  ['info']   #todo: not implemented
                                                  ['nTemps']
                                                  ['temps']
                                                  ['table']
                                                     .
                                                     .
                                                     .
                                                    etc
        """
        
        strng = open(fPath, 'r').read()
        
        #the dictonary object which will hold all the info about the specie
        specDict = {}
        specDict['path'] = fPath
        
        #regular expression which gets everything between a '!'
        #and a new line followed by a '!'
        tokens = re.findall(r"[\s]*!([\s\S]*?)\n\s*(?=[!]|$)", strng)
        #for t in tokens:
        #    print (t.split('\n'))[0:1]
        
        #token : the molecule
        #####################
        tok = tokens.pop(0)
        specDict['specStr'], specDict['info'] = self.get_specie_name(tok)
        #token : the molecular weight
        tok = tokens.pop(0)
        if 'MOLECULAR WEIGHT' in tok.upper() or 'MASS' in tok:
            specDict['weight'] = np.float64((tok.split('\n'))[-1])
        else:
            raise ValueError('expected the MOLECULE WEIGHT section of the file, got *%s* instead'  % tok)
        #token : n levels
        #####################
        tok = tokens.pop(0)
        if 'NUMBER OF ENERGY LEVELS' in tok.upper():
            tok = (tok.split('\n'))[-1]
            tok = re.findall(r"\w+", tok) #splitting to words 
            tok = tok[0] # picking the first (assuming it is the number we need)
            specDict['nlevels'] = np.int32(tok)
        else:
            raise ValueError('expected the the NUMBER OF ENERGY LEVELS section file, got *%s* instead' % tok)
        #token : the levels
        ####################
        tok = (tokens.pop(0)).split('\n')
        hdr  = tok[0].upper().strip()
        data = tok[1:]
        
        #checking the units of the energy levels
        hdrEnenrgy = hdr.split('+')[1]
        if 'CM' not in hdrEnenrgy:
            raise ValueError('Energy unit is not in CM^-1')
        else:
            specDict['energyUnit'] = 'cm^-1'
            
        if 'LEVEL' in hdr and 'ENERG' in hdr and 'WEIGHT' in hdr: 
            specDict['levels'] = data
        else:
            raise ValueError('expected the the LEVELs section file, got *%s* instead' % tok)
        #token : the number of radiative transitions
        ############################################
        tok = tokens.pop(0)
        if 'NUMBER OF RADIATIVE TRANSITIONS' in tok.upper():
            specDict['nTransRad'] = np.int32((tok.split('\n'))[-1])
        else:
            raise ValueError('expected the the NUMBER OF RADIATIVE TRANSITIONS section file, got *%s* instead' % tok)
        #token : the radiative transitions
        ##################################
        tok = (tokens.pop(0)).split('\n')
        hdr  = tok[0].upper()
        data = tok[1:]
        if 'TRANS' in hdr and 'U' in hdr and 'L' in hdr and 'FREQ' in hdr:
            specDict['transRad'] = data
        else:
            raise ValueError('expected the TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K) section file, got *%s* instead' % hdr)
        # collisional transitions
        ##########################
        specDict['transColl'] = {}
        specDict['transColl']['partnersList'] = []
        #token : the number of collisional partners
        tok = tokens.pop(0)
        if 'NUMBER OF COLL PARTNERS' in tok.upper() or 'NUMBER OF COLLISION PARTNERS' in tok.upper():
            specDict['transColl']['nPartners'] = np.int32((tok.split('\n'))[-1])
        else:
            raise ValueError('expected the NUMBER OF COLL PARTNERS section file, got *%s* instead' % tok)

        #parsing the collisional information for each partner
        for i in np.arange(specDict['transColl']['nPartners']):
            #token : the collisional partner
            tok = tokens.pop(0)
            """
            p1 = fPath.split('/')[-1]
            p2 = tok.split('\n')[1]
            p2 = re.sub(r'!.*','',p2)
            p2 = re.sub(r'from.*','',p2)
            p2 = re.sub(r',.*','',p2)
            p2 = re.sub(r'scaled.*','',p2)
            self.dummy += ((p1, p2,),)
            """
            if 'COLLISION' in tok.upper() or 'PARTNER' in tok.upper():
                specDict['transColl']['partnersList'].append((tok.split('\n'))[-1])
            else:
                raise ValueError('expected the COLLISIONS BETWEEN section file, got *%s* instead' % tok)
            partnerLabel = 'partner%s' % i
            specDict['transColl'][partnerLabel] = {}
            
            specDict['transColl'][partnerLabel]['info'] = tok
            
            #token : number of collisional transitions
            tok = tokens.pop(0)
            if 'NUMBER OF' in tok.upper() and 'COLL' in tok.upper() and 'TRANS' in tok.upper():
                tok = (tok.split('\n'))[-1]
                tok = re.findall(r"\w+", tok) #splitting to words     
                tok = tok[0] # picking the first (assuming it is the number we need)
                specDict['transColl'][partnerLabel]['nTrans'] = np.int32(tok)
            else:
                raise ValueError('expected the NUMBER OF COLL TRANS section file, got *%s* instead' % tok)
            #token : number of temperatures at which the rates are tabulated
            tok = tokens.pop(0)
            if 'NUMBER OF' in tok.upper() and 'TEMP' in tok.upper():
                specDict['transColl'][partnerLabel]['nTemps'] = np.int32((tok.split('\n'))[-1])
            else:
                raise ValueError('expected the NUMBER OF COLL TEMPS section file, got *%s* instead' % tok)
            #token : the temperature at which the rates are tabulated
            tok = tokens.pop(0)
            if 'TEMP' in tok.upper():
                temps = re.split('\n', tok, maxsplit=1)[1] # discaring the header
                temps = re.sub('\n', ' ', temps)           # replacing the newline (if any) with a whitespace char
                specDict['transColl'][partnerLabel]['temps'] = np.float64(temps.split())
            else:
                raise ValueError('expected the COLL TEMPS section file, got *%s* instead' % tok)
            #token : the collision rates data table
            tok = tokens.pop(0)
            if 'TRANS' in tok or 'RATE' in tok or 'cm^3 s^-1' in tok:
                specDict['transColl'][partnerLabel]['table'] = (tok.split('\n'))[1:]
            else:
                raise ValueError('expected the TRANS + UP + LOW + COLLRATES section file, got *%s* instead' % tok)

        return specDict 
         
    def findFiles(self, **kwargs):
        """returns a list of all the files ending with a .dat in :data:`dirPath`"""
        # getting the names of the meshes in that dir
        files = []
        for infile in glob.glob( os.path.join(self.dirPath, '*.dat') ):
            files.append(infile)
        return files
    
    def get_specie_name(self, tok):
        """Extracts the name of the specie from the token. In the current version
        of the database, the specie names also contain some text sometime like references
        and details about the specie. The specie name is in most cases the first item on
        the line after !MOLECULE, except for H20 where whenever H20 is found in the string
        the specie name is taken to be the second word in the string.
        
        :return: (specieName, info). Where info is some description in the string feild if any.
        :TODO: add a method which converts the string name to one which is compatible with the 
            specie class.
        """
        
        if 'MOLECULE' in tok.upper():
            fieldName, content = tok.split('\n')
            splt = content.strip().split(' ')
            if 'H2O' not in content.upper(): #only H2O has extra stuff in the name of the string feild
                specStr = splt[0]
            else:
                for word in splt:
                    if 'H2O' in word.upper():
                        specStr = word
                        break
                    else:
                        specStr = ''
            info = content
        else:
            raise ValueError('expected the MOLECULE section of the file, got *%s* instead' % tok)
        
        #removing the p-, o-, E- from the species names
        specStr = specStr.replace('o-','').replace('p-','').replace('E-','')
        return (specStr, info)
    
    def get_specie(self, specStr = None, inPath = None, inInfo = None):
        """This method can be used to return the data for a certain specie one
         the data has been read.  None of the keywords are mandatory. But
         at least one should be present. The checks are done using the 'AND'
         logic.  
         
         :param string specStr: the specie to be returned (not that if multiple 
           files contain info about the same specie, multiple species info
           will be returned. 
         :param string inFname: If this string is in the filename, the specie
           will be consider (assuming the remaining keyword checks restuls True).  
         :param string inInfo: if this string is in the MOLECULE feild, the check
           results as True.               
        """
        
        returnList = ()
        for spec in self.speciesInfo:
            
            consider = True
            
            if specStr in spec['specStr']:
                pass
            else:
                consider *= False
                
            if inPath != None:
                if inPath in spec['path'].split('/')[-1]:
                    pass
                else:
                    consider *= False
                    
            if inInfo != None:
                if inInfo in spec['info']:
                    pass
                else:
                    consider *= False
            
            if consider == True:                     
                returnList += (spec,)
    
        return returnList
    
    def set_dirPath(self, path):
        self.dirPath = path
    def get_dirPath(self):
        return self.dirPath
    