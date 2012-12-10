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

### - add a method to return the index of a level given the quantum numbers
### - add another method which feches the transition info give the upper and lower levels (all transitions -rad-col..etc..)

class reader():
    """This is a class which provides an interface to the LAMDA database. 
    """
    def __init__(self, *args, **kwargs):
        
        self.dataFiles = None
        """A list containing all the names of the files in self.dirPath"""

        self.speciesInfo = ()
        """A dictionary holding the species string and the files containing data about them"""
                
        if 'dirPath' in kwargs:
            self.set_dirPath(kwargs['dirPath'])
            self.dataFiles = self.findFiles()
            self.collectAllFilesInfo() 
        else:
            self.dirPath = None
            """The path of the directory holding all the information about all the species"""

    def loadSpeciesData(self, specie = None):
        """Returns a numpy dtype holding all the information about the specie in the LMBDA database"""
        pass

    def collectAllFilesInfo(self):
        """Loops over all the data files in self.dirPath and collects the basic information 
              - name of the species
              - file names containing info about that species
           The return values is a dict 'species' : [file1, file2,]   
        """ 
        
        for fPath in self.dataFiles:
            if 'p-nh3.dat' not in fPath:  #;;;rmove this later
                continue                  #;;;remove this later
            #if 'test.dat' not in fPath:  #;;;rmove this later
            #    continue                  #;;;remove this later
            #input('press 1 to process the next file %s' % fPath)
            specDict = self.parseDataFile(fPath)
            specDict = self.cleanSpecDict(specDict)
            self.speciesInfo += (specDict,)
            
    def cleanSpecDict(self, specDict):
        """takes the output item of self.parseDataFile() and defines dtypes which
           hold the info about the levels and the transitions. 
           
           this method does the following modification to specDict:
             -specDict['levels']
               is replaced by a numpy dtype : 'n', 'E', 'g', 'j', 'k', 'i' 
             -specDict['transRad']
               is replaced by a numpy dtype : 'n', 'u', 'l', 'A', 'nu', 'E' 
             -specDict['transColl']['partnerX']
               each collision partner X will have a new dict key ['trans']
               which is a dtype with keys 'n', 'u', 'l', 'rc'
               'rc' is an interpolation function which take temperature as an argument
               and returns the rate coefficient for that temperature.
             -in the above : 'n' is a transition index.
             -['levels'], ['transRad'], ['transColl']['partnerX']['dtype'] are ndarrays 
              one entry for each level/transition.    
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
            levels[idx]['j'] = np.int32(levelStrParsed[3])
            levels[idx]['k'] = np.int32(levelStrParsed[4])
            levels[idx]['i'] = np.int32(levelStrParsed[5])

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
            f_rc = interp1d(tKin, rc)
            transColl[idx]['rc'] = f_rc
            """
            if idx >= 0:    
                pyl.plot(tKin, rc, 'r+')
                TkinNew = np.arange(15,300,1)
                rcNew = transColl[idx]['rc'](TkinNew)
                pyl.plot(TkinNew, rcNew,'b-')
                pyl.show()
            """
        specDict['transColl']['partner0']['trans'] = transColl  
        
        return specDict
    
    def parseDataFile(self, fPath):
        """parses the data file into its components and returns them as a nested dict object.
            the keys of the dict are :
            
            specDict = ['path']
                       ['specStr']
                       ['weight']
                       ['nlevels']
                       ['levels']
                       ['nTransRad']
                       ['transRad']
                       ['transColl']['nPartners']
                                    ['partnersList']
                                    ['partner0']['nTrans']
                                                ['nTemps']
                                                ['temps']
                                                ['table']
                                    ['partner1']['nTrans']
                                                ['nTemps']
                                                ['temps']
                                                ['table']
                                                .
                                                .
                                                .
                                               etc
            returns specDict
        """
        
        #fPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles/ch3oh_merged_corrected.dat'
        print fPath
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
        tok = tokens.pop(0)
        if 'MOLECULE' in tok.upper():
            specDict['specStr'] = (tok.split('\n'))[-1]
        else:
            raise ValueError('expected the MOLECULE section of the file, got *%s* instead' % tok)
        #token : the molecular weight
        tok = tokens.pop(0)
        if 'MOLECULAR WEIGHT' in tok.upper() or 'MASS' in tok:
            specDict['weight'] = np.float64((tok.split('\n'))[-1])
        else:
            raise ValueError('expected the MOLECULE WEIGHT section of the file, got *%s* instead'  % tok)
        #token : n levels
        tok = tokens.pop(0)
        if 'NUMBER OF ENERGY LEVELS' in tok.upper():
            tok = (tok.split('\n'))[-1]
            tok = re.findall(r"\w+", tok) #splitting to words 
            tok = tok[0] # picking the first (assuming it is the number we need)
            specDict['nlevels'] = np.int32(tok)
        else:
            raise ValueError('expected the the NUMBER OF ENERGY LEVELS section file, got *%s* instead' % tok)
        #token : the levels
        tok = (tokens.pop(0)).split('\n')
        hdr  = tok[0].upper()
        data = tok[1:]
        if 'LEVEL' in hdr and 'ENERG' in hdr and 'WEIGHT' in hdr: #;;;not comparing the whole string as in the data files since it differers in some files
            specDict['levels'] = data
        else:
            raise ValueError('expected the the LEVELs section file, got *%s* instead' % tok)
        #token : the number of radiative transitions
        tok = tokens.pop(0)
        if 'NUMBER OF RADIATIVE TRANSITIONS' in tok.upper():
            specDict['nTransRad'] = np.int32((tok.split('\n'))[-1])
        else:
            raise ValueError('expected the the NUMBER OF RADIATIVE TRANSITIONS section file, got *%s* instead' % tok)
        #token : the radiative transitions
        tok = (tokens.pop(0)).split('\n')
        hdr  = tok[0].upper()
        data = tok[1:]
        if 'TRANS' in hdr and 'U' in hdr and 'L' in hdr and 'FREQ' in hdr:
            specDict['transRad'] = data
        else:
            raise ValueError('expected the TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K) section file, got *%s* instead' % hdr)
        ##########################
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
            if 'COLLISION' in tok.upper() or 'PARTNER' in tok.upper():
                specDict['transColl']['partnersList'].append((tok.split('\n'))[-1])
            else:
                raise ValueError('expected the COLLISIONS BETWEEN section file, got *%s* instead' % tok)
            partnerLabel = 'partner%s' % i
            specDict['transColl'][partnerLabel] = {}
            
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
            #token : the temperature at whcih the rate are tabulated
            tok = tokens.pop(0)
            if 'TEMP' in tok.upper():
                tok = (tok.split('\n'))[-1]
                specDict['transColl'][partnerLabel]['temps'] = np.float64(tok.split())
            else:
                raise ValueError('expected the COLL TEMPS section file, got *%s* instead' % tok)
            #token : the collision rates data table
            tok = tokens.pop(0)
            if 'TRANS' in tok or 'RATE' in tok or 'cm^3 s^-1' in tok:
                specDict['transColl'][partnerLabel]['table'] = (tok.split('\n'))[1:]
            else:
                raise ValueError('expected the TRANS + UP + LOW + COLLRATES section file, got *%s* instead' % tok)

        return specDict 
         
    def findFiles(self):
        """returns a list of all the files ending with a .dat in self.dirPath"""
        # getting the names of the meshes in that dir
        files = []
        for infile in glob.glob( os.path.join(self.dirPath, '*.dat') ):
            files.append(infile)
        return files
        
    def set_dirPath(self, path):
        self.dirPath = path
    def get_dirPath(self):
        return self.dirPath