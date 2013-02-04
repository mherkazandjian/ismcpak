import os, glob, fnmatch, sys, pickle
import numpy as np
import pylab as pyl
import logging
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d import Axes3D

from time  import *
from mesh  import *
from mylib.utils.ndmesh import ndmesh
from mylib.utils.misc   import fetchNestedDtypeValue
from ismUtils   import *
from radex      import *
from scipy      import interpolate
import chemicalNetwork

class meshArxv():
    """ this class generates and manipulates archives of PDR meshes.

     :param bool readDb: if this is set, the database is read during initialization.
      
     :param bool set_defaults: if this is passed as a keyword, self.set_default attributes and self.set_grid_axes_quantity_values are called upon intialization.
       
     :param bool no_init_from_run_parms: if this is passed as a keyword, nothing is initialized from used_parms.pkl.                

     FILES AND THEIR FORMATS: by default, the prefix name of the individual mesh files is assumed to be
      mesh.dat-id-xxxxxx
     
     The mesh database files are assumed to have the same prefix, for example, if the database file provided is foo, then this routine will assume (or write).
       
        - foo.info, foo.db
         
     for now the database can be stored into a single file and not split into
     multiple files.  
     
     the database for the meshes is constructed of two binary files:
     
         - foo.info : holding all the information about the meshes, locations and parameters
         - foo.db   : holds the data for all the meshes
       
     the info file has the following structure :
        
       - version number    : int32 int32 int32
       - nMeshes           : int32 
       - meshes info array : meshInfoArrayFormat (see also the method)
         ( nMeshes, 1)          

          the entries of the dtype are
           
            - mesh number     ( int64 ) :
            - data file index ( int64 ) : i.e in which .db file the mesh is located
            - offset ( int64 ) : the offset in bytes from the beginning of the file where the mesh is located
            - nSteps ( int64 ) : the number of slabs in the mesh
            - nSpecs ( int64 ) : the number of species in the mesh                                   .
                        
      the .db files have the following structure:
      
         - mesh_1 ( mesh dtype )
         - checkNum_1 ( int64 )
         - mesh_2 ( mesh dtype )
         - checkNum_2 ( int64 )

      For a mesh number i, the checkNum_i should be the same as the i^th entry 
      in the info array offset...i.e chechNum = infoAll[i]['info'][2] 
      
      features:
      
          - cliking on the panel of the temperature, sets the Av and gasT to the chemical
            network of the self.mshTmp object (which is the same as self.chemNet). Once
            those are set, the reactions can be examined, for example :
              
            .. code-block:: python
              
                ids = net.filter_reactions(withType='CP', 
                                           show=True,  
                                           fmt='id type rxn abg cst rate', 
                                           sort = True)
          
      :TODO: modify such that it can handle models with one slab only (the surface slab).
      
      see test_writeReadArxv.py for an example
    """
    def __init__(self, *args, **kwargs):
        
        # attributes
        if 'dirPath' in kwargs:
            self.dirPath = kwargs['dirPath']
        else:
            self.dirPath = None
            """path to the directory which contains all the database info and data"""
            
        self.nMeshes    = None
        """ np.long64 : number of meshes in the database (this is now as 1x1 array, convert it just into an np.int32 scalar"""
        
        self.ver        = None
        """ np.int32 array [3] version number of the database data"""
        
        self.meshes     = None
        """A list of all the meshes of 'mesh' dtypes (see mesh.py).
        
           For a mesh of at index 'i' in
            
               self.meshes[i]
               
           The corresponding info in the header are accessed as follows :
            
             - self.infoAll[i]['parms'][0]) which should be the same as self.meshes[i]['hdr']['G0'] 
             - self.infoAll[i]['parms'][1]) which should be the same as self.meshes[i]['hdr']['nGas']
             - self.infoAll[i]['parms'][2]) which should be the same as self.meshes[i]['hdr']['gammaMech']
        """
        
        self.infoAll    = None 
        """A numpy ndarray of dtype arxvHdrDtype (see below) which contains the info
           (headers) about all the meshes each entry in this array contains two things,
           information about the mesh and the parameters of the mesh. For example
           the elements x in self.infoAll[x] has the following contents : 
           self.infoAll[x]['info'][0]   mesh number\n 
           self.infoAll[x]['info'][1]   0 ( for now)\n
           self.infoAll[x]['info'][2]   offset from the start of file\n
           self.infoAll[x]['info'][3]   number of steps in the mesh\n
           self.infoAll[x]['info'][4]   number of species\n
         
           self.infoAll[x]['parms'][0]  G0\n
           self.infoAll[x]['parms'][1]  nGas\n
           self.infoAll[x]['parms'][2]  gammaMech\n
           self.infoAll[x]['parms'][3]  0.0, NOT USED\n
           self.infoAll[x]['parms'][4]  0.0  NOT USED\n
        """

        #---------------------------radex data storage attributes-----------------------------------
        self.verRadex = None
        """a 3 element np.int32 array holding the version of the radex database currently in use."""
        
        self.currentRadexDb = None
        """A string which holds the species name of the current radex database being used."""
        
        self.infoAllRadex = None
        """A numpy ndarray of length self.nMeshes of dtype returned by self.arxvRadexHdrFormat which
           holds the number of transitions computed for the mesh it corresponds to in self.meshes and
           self.infoAll. A value of zero means that there were was valid output from radex. The contents
           of each elements in the array are (describing the dtype entries):
        
           .. code-block:: python
            
             self.infoAllRadex[x]['info'][0]   mesh number 
             self.infoAllRadex[x]['info'][1]   number of transitions (0 => no radex data for this mesh)
             self.infoAllRadex[x]['info'][2]   offset from the start of file (computed and set only when the databse will be written, 0 otherwise)
             self.infoAllRadex[x]['info'][3]   warning code # holds radex.status after radex runs for the model
             self.infoAllRadex[x]['info'][4]   error code
         
             self.infoAllRadex[x]['parms'][0]  G0  (same as the entries in self.infoAll)
             self.infoAllRadex[x]['parms'][1]  nGas
             self.infoAllRadex[x]['parms'][2]  gammaMech
             self.infoAllRadex[x]['parms'][3]  0.0, NOT USED
             self.infoAllRadex[x]['parms'][4]  0.0  NOT USED


           :note: the transitions are stored even if there are warnings when running radex. So 
            take good care when analyzing the data. The info of this attribute is stored into 
            infoAllRadex.db.info.(specStr). Each specie will have its own .db.info.... file
            and a corresponding .db.(specStr) file which will hold all the data in self.meshesRadex.
        """
                    
        self.meshesRadex = None
        """A list holding all the info dumped by radex for each model. Each entry in the list is an ndarray
         of dtype :data:`radex.radex.transitionDtype`. The keys in each transition dtype are listed in 
         :data:`radex.radex.transitions` . 
        """
        
        self.radexDbs = {}
        """A dictionary holding the radex databases read so far. Each entry in the
        dictionary is a dictionary of two keys : 'infoAll' and 'meshes'. For example, 
        this could hold the databases for CO, HCN and HNC. Those can be mannualy
        to self.meshesRadex and self.infoAllRadex via :
        
        .. code-block:: python
        
            #setting the CO radex database
            self.meshesRadex  = self.radexDbs['CO']['meshes']
            self.infoAllRadex = self.radexDbs['CO']['infoAll']
            
            #setting the HCN radex database
            self.meshesRadex  = self.radexDbs['HCN']['meshes']
            self.infoAllRadex = self.radexDbs['HCN']['infoAll']
        """
        #-------------------------------------------------------------------------------------
        
        self.grid_qx = None
        """A list pointing to the quantity in the mesh.data dtype to be used for the x-axsi"""
        self.grid_x = None
        """the values for all the meshes corresponding to qx"""
        self.grid_qy = None
        """same as qx but for the y axis"""
        self.grid_y = None
        """same as qx but for the y axis"""
        self.grid_qz = None
        """same as qx but for the z axis"""
        self.grid_z = None
        """same as qx but for the z axis"""
        self.grid_z_unique = None
        """holds the unique values of grid_z"""

        self.nDbFiles   = None
        """ numpy.int32 number of database files"""
        
        self.chemNet    = None
        """ object of class type chemicalNetwork holds info about the chemical network used"""
        
        if 'plotRanges' in kwargs:
            self.ranges = kwargs['plotRanges']
        else:
            self.ranges = None
        """The ranges which will be used for displaying the x,y,z dimesntions for grid_x,yz"""
        
        # variables used for plotting
        self.pltGmSec    = None #: the value of the section in Gmech selected
        self.grds        = None #: 2x2 ndmesh array object
        self.pltRadex    = None
        self.resPltGrids = None
        
        #attributes used for storing the interpolation functions from which 
        #the 2D grids to be displayed
        self.grdInterp_f           = None
        self.abunGridInterp_f      = None
        self.colDensGridInterp_f   = None
        self.intensityGridInterp_f = None
        
        self.parms = kwargs
        """all the keywords used in setting up the database and running it."""
        
        """resolution of the interpolated grids"""
        
        if 'metallicity' in kwargs:
            self.set_metallicity( kwargs['metallicity'] )
        else:
            self.metallicity = None

        self.radexObj   = None
        #--------------------------GEN AND PLOT RADEX CURVES---------------------------------
        if 'radex' in self.parms:
            if self.parms['radex']['use']:
                #making the instance            
                self.radexObj = radex(self.parms['radex']['path'], self.parms['radex']['molDataDirPath'])
                #setting some default radex paraemeters
                inFile = {'specStr'                : self.parms['radex']['specStr'],
                          'outPath'                : 'foo' ,
                          'freqRange'              : self.parms['radex']['freqRange'],
                          'tKin'                   : None  ,  # depending on the model
                          'collisionPartners'      : self.parms['radex']['collisionPartners'],
                          'nDensCollisionPartners' : None  ,  # depending on the model
                          'tBack'                  : self.parms['radex']['tBack'],  # fixed
                          'molnDens'               : None  ,  # depending on the model
                          'lineWidth'              : self.parms['radex']['lineWidth'],  # fixed
                          'runAnother'             : 1     }
                self.radexObj.setInFile( inFile )

        self.mshTmp = None
        """a mesh object which is used to store single mesh data just for plotting purposes"""
        
        self.logger = self.setupLogger()
        """The logging object which will be used to output stuff to stdout"""
        
        self.gui = None
        self.parms_used_in_PDR_models = None
        """a dictionary holding the parameters used in modeling the PDRs.  This is read
        from self.dirPath/used_parms.pkl""" 
        
        #reading the used_parms.pkl file (if found) and setting the vaues read from it
        if 'no_init_from_run_parms' not in kwargs:
            self.read_used_parms_used_in_PDR_models()

        if 'readDb' in kwargs:
            kwargs.pop('readDb')
            self.readDb( check = True)

        if 'set_defaults' in kwargs:
            kwargs.pop('set_defaults')
            self.set_default_attributes()
            self.set_grid_axes_quantity_values()

    def read_used_parms_used_in_PDR_models(self):
        fpath_parms_used = os.path.join(self.dirPath,'used_parms.pkl')
        if os.path.isfile(fpath_parms_used):
            self.logger.debug('setting metallicity and chemical network from : %s' % fpath_parms_used)
            #loading the parameters
            parms = pickle.load(open(fpath_parms_used))
            #setting the attributes
            self.metallicity = parms['metallicity']
            self.setupChemistry(parms['chemistry']) 
            self.parms_used_in_PDR_models = parms 
        else:
            self.logger.debug('no used_parms.pkl file found : %s' % fpath_parms_used) 
    
    def construct(self, meshNamePrefix = None, writeDb = None ):
        """Construc the database anad write the .db files. If the meshNamePrefix 
            is not supplied all the files in :data:`dirPath` are assumed to be 
            data files and all of them are put in the database. The files with 
            matching name prefix are looked for resusively in the 'meshes' directory.
        """
        
        if meshNamePrefix == None:
            meshNamePrefix = ''
        
        sourceDir = os.path.join(self.dirPath, 'meshes')
        lookFor   = meshNamePrefix + '*'
        
        # getting the names of the meshes in that dir (recursively)
        files = []
        for root, dirnames, filenames in os.walk(sourceDir):
            for filename in fnmatch.filter(filenames, lookFor):
                files.append(os.path.join(root, filename))
        
        # setting variable for the 
        self.nMeshes = np.zeros( 1, dtype = np.int32 )
        self.nMeshes[0] = len(files) 
        self.logger.debug('found %s meshes' % (self.nMeshes))

        # defining the array holding the info about all the meshes
        # and their location in the database ..etc..
        arxvHdrDtype = np.dtype( self.arxvHdrFormat() )
        infoAll = np.ndarray( self.nMeshes, dtype = arxvHdrDtype )
        self.meshes = []
        
        # reading the individual mesh files and writing the files into .db file(s)
        for i in np.arange( self.nMeshes ):
    
            fName = files[i]
    
            # reading a mesh
            m = mesh(fName)
            mData = m.getData()

            self.meshes.append(mData)
        
            # filling the corresponding field in the header array
            infoAll[i]['info'][0] = i # mesh number 
            infoAll[i]['info'][1] = 0 # data file index
            ##infoAll[i]['info'][2] is filled in self.writeDb()
            infoAll[i]['info'][3] = mData['hdr']['nSteps']  
            infoAll[i]['info'][4] = mData['hdr']['nSpecs']
            infoAll[i]['info'][5] = mData['hdr']['version']
     
            infoAll[i]['parms'][0] = mData['hdr']['G0']
            infoAll[i]['parms'][1] = mData['hdr']['nGas']
            infoAll[i]['parms'][2] = mData['hdr']['gammaMech']
            infoAll[i]['parms'][3] = np.float64(0.0)
            infoAll[i]['parms'][4] = np.float64(0.0)
            
            del m
            
            if (i % np.int32(self.nMeshes / 100)) == 0:
                self.logger.debug('proccessed %d meshes...' % i) 
            
        self.infoAll = infoAll
        
        if writeDb != None:
            self.writeDb()

    def writeDb(self, dirName = None):
        """ writes the mesh.db and mesh.db.ino into the dir dirName. By default the files are
        written to the directorey self.dirName, unless a diffrent path is specified via
        the keyword dirName.
        
        :param string dirName: The directory where to write the database files.
        """
        
        if dirName == None:
            dirName = self.dirPath
        
        # writing the mesh to the database file
        dbDataFObj = file(dirName + 'meshes.db', 'wb')
        
        for i in np.arange( self.nMeshes ):
            self.meshes[i].tofile( dbDataFObj )
            self.infoAll[i]['info'][2] = dbDataFObj.tell()  # offset from the start of file
            np.array( self.infoAll[i]['info'][2] ).tofile(dbDataFObj)
        dbDataFObj.close()
        
        # writing the db info into a file
        dbInfoFObj = file(dirName + 'meshes.db.info', 'wb')
        self.ver.tofile( dbInfoFObj)
        self.nMeshes.tofile( dbInfoFObj )
        self.infoAll.tofile( dbInfoFObj )
        dbInfoFObj.close()
        
        self.logger.debug('wrote successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))        
        
    def readDb(self, check = None):
        """ reads the database and assigns the appropritate attributes (document)
            
            :param check: if this is set (to any value) the self.checkIntegrity() is called.
        """ 
        
        dbInfoFObj = file(self.dirPath + 'meshes.db.info', 'rb')
        dbDataFObj = file(self.dirPath + 'meshes.db'     , 'rb')

        self.ver = np.fromfile( dbInfoFObj, dtype = (np.int32, 3), count = 1)
        self.ver = self.ver[0]
    
        self.nMeshes = np.fromfile( dbInfoFObj, dtype = np.int32, count = 1)
        self.nMeshes = self.nMeshes[0]

        arxvHdrDtype = np.dtype( self.arxvHdrFormat() )
        self.infoAll = np.fromfile( dbInfoFObj, dtype = arxvHdrDtype, count = self.nMeshes)

        # reading the meshes in database into a list 
        self.meshes = []
        
        for i in np.arange(self.nMeshes):
            
            mDummy = mesh()
            
            nSteps  = self.infoAll[i]['info'][3]
            nSpecs  = self.infoAll[i]['info'][4]
            meshVer = self.infoAll[i]['info'][5]
            
            thisMeshDtype = mDummy.constructMeshDtype(nSpecs, nSteps, meshVer)
            thisMeshData = np.fromfile( dbDataFObj, dtype = thisMeshDtype, count = 1)
            
            self.meshes.append( thisMeshData[0] )                
            checkNum = np.fromfile( dbDataFObj, dtype = np.int64, count = 1)
            
            if checkNum != self.infoAll[i]['info'][2]:
                strng  = 'Error : checkpoint numbers do not match : database may be corrupt.\n'
                strng += 'read = %d, expected = %d' % (checkNum, self.infoAll[i]['info'][2])
                strng += ' when reading mesh %d' % i                
                raise NameError(strng)
                    
        self.logger.debug('read successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))
        
        if check:
            self.checkIntegrity()
            
    def writeDbRadex(self, dirName = None):
        """ writes the radex meshesRadex.db and meshesRadex.db.info files.  By default the files are
        written to the directorey self.dirName. , unless a diffrent path is specified via
        the keyword dirName. For each mesh, a block consists of : [log10_G0, log10_nGas, log10_gammaMech,],
        transitionData, filePosition....etc...
        if transitionData is none, then transitionData is not written in the block, but the rest of the
        block is written.
        
        
        :param string dirName: The directory where to write the database files.
        """
    
        if dirName == None:
            dirName = self.dirPath
        
        # writing the mesh to the database file
        dbDataFObj = file(dirName + 'meshesRadex.db.' + self.parms['radex']['specStr'], 'wb')
        
        for i in np.arange( self.nMeshes ):
            
            #writing the mesh parameters before writing the transition data (extra data integrirty check)
            thisMeshParms = np.array([self.meshes[i]['hdr']['G0'],
                                      self.meshes[i]['hdr']['nGas'],
                                      self.meshes[i]['hdr']['gammaMech']], dtype = np.float64)
            thisMeshParms.tofile( dbDataFObj )

            #writing the transition info to the transiotions db file
            if self.meshesRadex[i] != None:
                #writing the transition data
                self.meshesRadex[i].tofile( dbDataFObj )
                #filling the necessary values in infoAllRadex
                self.infoAllRadex[i]['info'][2] = dbDataFObj.tell()  # offset from the start of file
            else:
                self.infoAllRadex[i]['info'][2] = -1

            #writing the offset so far from the begining of the .db file
            np.array( self.infoAllRadex[i]['info'][2] ).tofile(dbDataFObj)

        dbDataFObj.close()

        # writing the db info into a file
        dbInfoFObj = file(self.parms['dirPath'] + 'meshesRadex.db.info.' + self.parms['radex']['specStr'], 'wb')
        self.verRadex.tofile( dbInfoFObj)
        self.nMeshes.tofile( dbInfoFObj )
        self.infoAllRadex.tofile( dbInfoFObj )
        dbInfoFObj.close()
        self.logger.debug('wrote successfully the radex database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))

    def readDbRadex(self, specStr, check = None):
        """ reads the database suffixed by specStr (i.e meshesRadex.db.(specStr)and assigns the
         appropritate attributes (document) and assigns the read data to self.meshesRadex and
         self.infoAllRadex. 
            
            :param bool check: if this is set (to any value) the self.checkIntegrity() is called.
            :param string specStr: the string of the specie whose database is to be read.
            :warning: the keyword specStr is not functional yet.
            :note: before calling this mehtod, an instance of the radex class should be     
              assigned to self.radexObj. 
        """ 
        
        #reading the database only if it is not already read. If it is read, it would
        #have a key in radexDbs
        if specStr not in self.radexDbs.keys():
            dbInfoFObj = file(self.dirPath + 'meshesRadex.db.info.' + specStr, 'rb')
            dbDataFObj = file(self.dirPath + 'meshesRadex.db.' + specStr, 'rb')
    
            self.verRadex = np.fromfile( dbInfoFObj, dtype = (np.int32, 3), count = 1)
            self.verRadex = self.verRadex[0]
    
            nMeshes = np.fromfile( dbInfoFObj, dtype = np.int32, count = 1)
            if self.nMeshes != nMeshes:
                raise ValueError('number of meshes in radex database is different from that of the meshes database')
            
            arxvRadexHdrDtype = np.dtype( self.arxvRadexHdrFormat() )
            infoAllRadex = np.fromfile( dbInfoFObj, dtype = arxvRadexHdrDtype, count = self.nMeshes)
            dbInfoFObj.close()
            
            # reading the meshes in database into a list 
            meshesRadex = []
            radexTransitionsDtype = self.radexObj.generateTransitionDtype()
    
            for i in np.arange(self.nMeshes):
    
                meshParms = np.fromfile( dbDataFObj, dtype = np.float64, count = 3)
    
                #overwriting the read parms into the infoAllRadex['parms']
                infoAllRadex[i]['parms'][0] = meshParms[0] #G0
                infoAllRadex[i]['parms'][1] = meshParms[1] #nGas
                infoAllRadex[i]['parms'][2] = meshParms[2] #gammaMech
                             
                # the number of transitions save in the database for this mesh
                nTransitions = infoAllRadex[i]['info'][1]
                
                if nTransitions > 0:
                    #reading the transition data for this mesh 
                    thisMeshRadexData = np.fromfile( dbDataFObj, dtype = radexTransitionsDtype, count = nTransitions)
                else:
                    thisMeshRadexData = None
                    
                #reading the check flag, which is the position in bytes in the database
                checkNum = np.fromfile( dbDataFObj, dtype = np.int64, count = 1)
    
                if checkNum != infoAllRadex[i]['info'][2]:
                    strng  = 'Error : While reading mesh %d \n' % i
                    strng += '           Error : checkpoint numbers do not match : radex database may be corrupt.\n'
                    strng += '           Error : read = %d, expected = %d' % (checkNum, infoAllRadex[i]['info'][2])
                    raise NameError(strng)
    
                #self.meshes.append( thisMeshRadexData[0] )                
                meshesRadex.append( thisMeshRadexData )                
            self.logger.debug('read successfully radex database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))
    
            self.infoAllRadex = infoAllRadex
            self.meshesRadex = meshesRadex
            
            if check:
                self.checkIntegrityRadex()            
        
            #storing the read Db into a dict
            self.currentRadexDb = specStr   
            radexDb = { 'meshes' : self.meshesRadex, 'infoAll' : self.infoAllRadex}
            self.radexDbs[specStr] = radexDb
            self.logger.debug('added %s radex database to memory' % specStr)
        else:
            self.logger.debug('%s radex database already read...skipping' % specStr)

    def readDbsRadex(self, species = None):
        """Loads radex database files from disk to a dictionary. The last database loaded is
        the one being used.
        
           :param species: A string or a list of strings holding the names of the species files to be loaded. 
        """
  
        if species == None:
            raise ValueError('keyword species is None, must be a list of strings.')

        if hasattr(species, '__iter__') == False:
            species = [species]      
            
        for specStr in species:
            self.readDbRadex(specStr, check = True)

    def use_radexDb(self, specStr):
        """A method which sets a radex database from the databases available in self.radexDbs."""
        self.meshesRadex    = self.radexDbs[specStr]['meshes']
        self.infoAllRadex   = self.radexDbs[specStr]['infoAll']
        self.currentRadexDb = specStr
        self.logger.debug('swithced to %s radex database' % specStr)
        
    def mergeDbs(self, newDbRunDirPath, outDirPath = None):
        """ merges two databases into one and write the resulting db and its info file
            :param string newDbRunDirPath: the dir in which the new db to be added to the current db is located.
        """
        
        arxv2 = meshArxv()
        arxv2.readDb( newDbRunDirPath )
        arxv2.checkIntegrity()

        if self.ver[0] != arxv2.ver[0] and self.ver[1] != arxv2.ver[1] and self.ver[2] != arxv2.ver[2]:
            raise ValueError('can not merge databases with different versions.')
        
        self.nMeshes[0] = self.nMeshes + arxv2.nMeshes
        
        # appending the meshes from the new arxv to the existing one
        for i in np.arange( arxv2.nMeshes ):
            self.meshes.append( arxv2.meshes.pop(0) )
        
        # merging the info arrays (need to update the numbers though
        self.infoAll = np.concatenate((self.infoAll, arxv2.infoAll),axis=1)
 
        # merging the info arrays
        for i in np.arange( self.nMeshes ):
            self.infoAll[i]['info'][0] = i # mesh number             
        
        if outDirPath != None:
            self.writeDb(outDirPath)
        
    def arxvHdrFormat(self):
        """ the format for the archive header from which the dtype will be constructed
         specifying the info and paramters for each line in the header
         in version 1 : the 'info' is defined as : ('info' , np.int64  , 5)
         in version 2 : the 'info' is defined as : ('info' , np.int64  , 6)
         """        
        
        # defining the version of the database
        self.ver = np.zeros([3], dtype = np.int32)
        self.ver[0] = 0
        self.ver[1] = 0
        self.ver[2] = 2

        return [ 
                  ('info' , np.int64  , 6),
                  ('parms', np.float64, 5)
               ]

    def arxvRadexHdrFormat(self):
        """ the format for the radex archive header from which the dtype will be constructed
         specifying the info and paramters for each line in the header
         in version 1 : the 'info' is defined as : ('info' , np.int64  , 5)
         """        
        
        # defining the version of the database
        self.verRadex = np.zeros([3], dtype = np.int32)
        self.verRadex[0] = 0
        self.verRadex[1] = 0
        self.verRadex[2] = 1

        return [ 
                  ('info' , np.int64  , 6),
                  ('parms', np.float64, 5)
               ]
        
    def checkIntegrity(self):
        """check archive integrity, this is a basic check where the information about the
             meshes in self.infoAll['info] is compared to those in self.data for each mesh.
             
             :NOTE: that this is a basic check.
             
             :TODO: see how to implement a checksum to check for the integrity of the data.
        """

        diff = 0.0
        for i in np.arange(self.nMeshes):
            x1 = np.array([np.log10(self.infoAll[i]['parms'][0]),
                           np.log10(self.infoAll[i]['parms'][1]),
                           np.log10(self.infoAll[i]['parms'][2])])

            x2 = np.array([np.log10(self.meshes[i]['hdr']['G0']),
                           np.log10(self.meshes[i]['hdr']['nGas']),
                           np.log10(self.meshes[i]['hdr']['gammaMech'])])

            xdv = x2 - x1
            diff +=  np.sqrt( np.dot(xdv,xdv) )
            
        if diff == 0.0:
            self.logger.debug('archive integrity test passed')
        else:
            strng = 'archive integrity test failed. database file may be corrupt' 
            raise NameError(strng)

    def checkIntegrityRadex(self):
        """checks the radex db archive integrity by comparing the read parameters for each mesh
        (3 float 64's before the transiotion data) which are stored into self.infoAllRadex['parms']
        once they are read from the meshesRadex.db.specStr and compares them to the mesh 
        parameters in self.infoAll['parms'].
             
        :TODO: see how to implement a checksum to check for the integrity of the data.
        """

        diff = 0.0
        for i in np.arange(self.nMeshes):
            x1 = np.array([np.log10(self.infoAllRadex[i]['parms'][0]),
                           np.log10(self.infoAllRadex[i]['parms'][1]),
                           np.log10(self.infoAllRadex[i]['parms'][2])])
            
            x2 = np.array([np.log10(self.infoAll[i]['parms'][0]),
                           np.log10(self.infoAll[i]['parms'][1]),
                           np.log10(self.infoAll[i]['parms'][2])])
            xdv = x2 - x1
            diff +=  np.sqrt( np.dot(xdv,xdv) )
            
        if diff == 0.0:
            self.logger.debug('archive integrity test passed')
        else:
            strng = 'archive integrity test failed. database file may be corrupt' 
            raise NameError(strng)
    
    def setupChemistry(self, parms = None):
        """sets up the attributes related to the chemistry assuming self.parms['chemistry'] is set.
           sets the attribute self.chemNet.
           
           :param parms: a dict holding the paths and definitions needed to setup the 
            chemical network.  An example of this dict :
            
            .. code-block:: python
            
                parms = {
                        'rxnFile'       : '/home/mher/rate99Fixed.inp',   # path of the reaction file
                        'specNumFile'   : '/home/mher/species.inp',       # the species in the rxnFile
                        'underAbunFile' : '/home/mher/underabundant.inp', # species to be excluded from the rxn file
                        'removeManual'  : ['13CH3'],                      # species to be removed explicitly not in the underabundant file
                        'baseSpecies'   : 'baseSpeciesDefault',           # name of the module holding the base species
                        'umistVer'      : 'umist99',                      # version of the umist database
                         }
        """
        
        if parms == None:
            parms =  self.parms['chemistry']
            
        #importing the module which holds the definitions of the base species
        baseSpecies = __import__(parms['baseSpecies'])
        baseSpecs = baseSpecies.baseSpecies()
        
        # settin up the orignial netowrk
        net = chemicalNetwork.chemicalNetwork(parms['rxnFile'], 
                                              baseSpecs,
                                              UMISTVER = parms['umistVer'])
        
        # reading the species to be removed from a file
        net.remove_species( underAbunFile = parms['underAbunFile'] )
        net.remove_species( species = parms['removeManual'] )
        # reading the species number and their corresponding indies and abundances from ascii files
        net.assign_numbers_to_species(fileName = parms['specNumFile'])
        
        self.setChemicalNetwork(net) # assiginig the chemical network to the archive
            
    def getQuantityFromAllMeshes(self, quantity, slabIdx = None, arrIdx = None):
        """ gets the quantity from all the meshes and returns it as a numpy array. 
            the quantity is mandatory, but no the slabIdx.
            
            :param int arrIdx: in case the quantity pointed to is an array, arrIdx would be the 
              index in the array we would like to retrieve.
        """
        
        values = np.zeros(self.nMeshes, dtype = np.float64)
        
        for i in np.arange(self.nMeshes):
            q = fetchNestedDtypeValue(self.meshes[i], quantity )
        
            if slabIdx != None and arrIdx != None:
                v = q[arrIdx][slabIdx]
            elif slabIdx != None and arrIdx == None:
                v = q[slabIdx]
            else:
                v = q

            values[i] = v  
                
        return values
        
    def construct3DInterpolationFunction(self, quantity = None, slabIdx = None, arrIdx = None, log10 = None, values = None, data = None, *args, **kwargs):
        """ returns a 3D interpolation function (interpolate.LinearNDInterpolator) which
             can be used to compute values determined by quantity give (nGas, G0, Gmech)
             (in log_10). The value to be interpolated upon is determined by the parameter
             quantity and the slab index determined by :data:`slabIdx`. For example, to construct an interpolation table for say surface
             temperature, this method can be invoked as :
             
             .. code-block:: python

                f = arvx.construct3DInterpolationFunction( quantity = ['state', 'gasT'], slabIdx = 0 )
             
             See :data:`computeInterpolated2DGrid` for an example on how to use this interpolation 
             function. Upon returning, this method sets the values of the attributes
             grid_x,y,z.
                
             :param list quantity: this is a list of strings which point to the data
                 to be extracted from the dtype :data:`mesh`. (see example above)
                 
             :param int32 slabIdx: the index of the slab from which the value will be extracted.
                 
             :param int32 arrIdx: The first index in the 2D array (ONLY necessary if the quantity pointed
               to is a 2D array, like ['state', 'abun']). In this case slabIdx is used as the seconds index. 
                
             :param bool log10: keyword which when passes as True, will generate the log10 of the quantity

             :param numpy.ndarray values: a numpy 1D array whose size if the same as the number of meshes which
               when passed as an argument, will be used as the value to be interpolated. In this case,
               quantity, slabIdx, arrIdx are ignored. 
             
             :param numpy.ndarray data: 
             
             :todo: modify this to construct the table over a selected range of the 3D 
                parameter space instead of using all the meshes. (optional)
                
             :warning: this fails if all the entries in one of the dimensions have the exact
               same value. in that case use :data:`construct2DInterpolationFunction`
            
        """
        
        # checking if the grid_x,y,z are set, if not, set them        
        self.set_grid_axes_quantity_values(**kwargs)
        self.set_attributes(**kwargs)
            
        if values == None:
            # the quantitiy we are intiesrested in showing
            values = self.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx, arrIdx = arrIdx)
        else:
            values = values

        if log10 != None and log10 == True:
            values[:] = np.log10(values[:])

        if data == None:
            data = np.array([self.grid_x, self.grid_y, self.grid_z]).T  #3D coordinates
        else:
            data = data
        
        ti = time()
        f = interpolate.LinearNDInterpolator(data, values) # getting the interpolation function     
        tf = time()
        self.logger.debug('constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti))

        return f

    def construct2DInterpolationFunction(self, quantity = None, slabIdx = None, *args, **kwargs):
        """ same as :data:`construct3DInterpolationFunction` but does it for n and G0
        
            :warning: make sure that the G_mech for all the models is the same. Other
              things will not make sense, bec the models would have diffrent G_mech
            :warning: make sure that the interpolated points are at the centroids of the grid points
        """
        
        x = np.log10( self.getQuantityFromAllMeshes( self.grid_qx ) )
        y = np.log10( self.getQuantityFromAllMeshes( self.grid_qy ) )
        
        values = self.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx)

        data = np.array([x, y]).T  #2D
        
        ti = time()
        f = interpolate.LinearNDInterpolator(data, values) # getting the interpolation function     
        tf = time()
        self.logger.debug('constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti))

        return f 

    def computeInterpolated2DGrid(self, ranges = None, res = None, zSec = None, 
                                  fInterp = None, *args, **kwargs):
        """ returns a 2D array ( a numpy ndarray ) of size res[0] and res[1] (check this if it is not the reverse) which holds
             the interpolated vlaues of 'quantity' over the domain determined by ranges for a slab whose index is slabIdx for
             a mechanical heating zSec (in log10). ( x is the horizontal direction of the mesh, and y is the vertical).
             
             An example would be (assuming the archive is constructed already):
             
             .. code-block:: python

                f   = arxv.construct3DInterpolationFunction(quantity = ['state', 'gasT'], 
                                                            slabIdx  = 0)
                grd = arxv.computeInterpolated2DGrid(slabIdx  = 0, 
                                                     ranges   = [[0,6],[0,6]],
                                                     res      = [100,100], 
                                                     zSec     = -30, 
                                                     fInterp  = f)


                import matplotlib.pyplot as plt
                plt.imshow(grd, extent=(0,1,0,1), origin='lower')
                plt.show()

             :param list ranges: a 2d list holding the ranges of the grid [[xmin, xmax],[ymin, ymax]]

             :param list res: the resolution in each dimension. [xResolution, yResolution]
             
             :param float zSec: the section in mechanical heating to be interpolated at.
             
             :param scipy.interpolate.LinearNDInterpolator fInterp: the interpolation function returned by
                 :data:`construct2DInterpolationFunction` or it 3D counterpart.        
        """
        
        # defining the points in the uniform 2D grid                                                                                                                                                                                                 
        grid_x, grid_y = np.mgrid[ranges[0][0]:ranges[0][1]:complex(0,res[0]),
                                  ranges[1][0]:ranges[1][1]:complex(0,res[1])]
        nPts = np.product(grid_x.shape)
        xNew = grid_x.reshape(nPts)
        yNew = grid_y.reshape(nPts)
        zNew = xNew.copy()
        
        # setting the values of mechanical heating to interpolate on
        zNew[:] = zSec      # all the grid points have the same mechanical heating
        
        dataNew = np.array( [xNew, yNew, zNew] ).T
        
        ti = time()
        tNew = fInterp(dataNew)
        tf = time()
        self.logger.debug('interpolated %d points in %f seconds at a rate of %e pts/sec' % (nPts, tf-ti, nPts / (tf-ti)))
        tNew = np.reshape(tNew, grid_x.shape)
        
        return tNew, grid_x, grid_y
        
    def showGrid(self, quantity = None, slabIdx = None, ranges = None, res = None, zSec = None, fInterp = None, *args, **kwargs):
        """This method displays a plot (useful for standalone use) of a grid of the quantity pointed by quantity in a new
           window. It makes use of :data:`computeInterpolated2DGrid` and :data:`construct3DInterpolationFunction`.
           For the documentation of the parameters see  :data:`computeInterpolated2DGrid`
           
           :param string quantity: (see construct3DInterpolationFunction documentation)
           :param interpolate.LinearNDInterpolator fInterp: an interpolation function function returned by
             :data:`construct3DInterpolationFunction` or it 2D analogue. If this is set, the quantity 
             keyword is meaningless.
           :param list res: (see self.computeInterpolated2DGrid)
           
           :return: aaa
        """

        self.set_attributes(**kwargs)
        self.set_default_attributes()
        
        # getting the grid
        if fInterp == None:
            f = self.construct3DInterpolationFunction(quantity = quantity, slabIdx  = slabIdx, *args, **kwargs)
        else:
            f = fInterp
        grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                        res      = res, 
                                                        zSec     = zSec, 
                                                        fInterp  = f, *args, **kwargs)
             
        # taking the transpose for plotting purposes
        grd = grd.T
        
        # plotting it
        if 'figure' not in kwargs:
            pyl.figure()
            pyl.subplot(111)
            
        im = pyl.imshow(grd, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        
        # adding levels and labels
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )
        #print grd
        CS = pyl.contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        pyl.clabel(CS, fmt = '%.1f' )
        
        pyl.colorbar(im, shrink = 0.8, extend = 'both')
        
        return ()

    def computeAndSetInterpolationFunctions(self, *args, **kwargs):
        """compute the 3D interpolation functions which will be used to compute the 2D grids.
            sets the attributes : self.grdInterp_f, self.abunGridInterp_f, self.colDensGridInterp_f
        """
        
        # computing the surface temperature interpolation function which will be used to make the high res maps
        #######################################################################################################
        if self.parms['gridsInfo']['00']['show']:
            self.grdInterp_f = self.construct3DInterpolationFunction(quantity = self.parms['gridsInfo']['00']['quantity'], 
                                                                     slabIdx  = self.parms['gridsInfo']['00']['slabIdx'],
                                                                     log10    = True,
                                                                     *args, **kwargs)
        
        # computing the abundance interpolation function which will be used to make the high res maps
        #############################################################################################
        if self.parms['gridsInfo']['01']['show']:
            self.abunGridInterp_f = self.construct3DInterpolationFunction(quantity = self.parms['gridsInfo']['01']['quantity'], 
                                                                          slabIdx  = self.parms['gridsInfo']['01']['slabIdx'],
                                                                          arrIdx   = self.chemNet.species[self.parms['gridsInfo']['01']['specStr']].num,
                                                                          log10    = True,
                                                                          *args, **kwargs)

        # computing the column density interpolation function which will be used to make the high res maps
        ##################################################################################################
        if self.parms['gridsInfo']['10']['show']:
            #getting the column densities for all the models
            m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
            values  = np.ndarray(self.nMeshes, dtype = np.float64)
             
            for i in np.arange(self.nMeshes):
                m.setData( self.meshes[i])
                specStr = self.parms['gridsInfo']['10']['specStr']
                colDens = m.getColumnDensity(specsStrs = [specStr], maxAv = self.parms['gridsInfo']['10']['maxAv'] )
                values[i] = colDens[0]
         
            self.colDensGridInterp_f = self.construct3DInterpolationFunction(values = values,
                                                                             log10  = True,
                                                                             *args, **kwargs)

        # computing the line intensity interpolation function which will be used to make the high res maps
        ##################################################################################################
        if self.parms['gridsInfo']['11']['show']:
             
            radiativeCodeType = self.parms['gridsInfo']['11']['type']
            
            #----------------------------emissions using radex---------------------------------------
            if radiativeCodeType == 'radex':

                values, grid_x, grid_y, grid_z = [], [], [], []

                transitionIdx = self.parms['gridsInfo']['11']['transitionIndx']
                quantity      = self.parms['gridsInfo']['11']['quantity']
            
                #looping over all the radex info of all the mesehes and getting the
                #transition data to be displayed
                for i in np.arange(self.nMeshes):
                    
                    if self.meshesRadex[i] != None:
                            
                            x, y, z = self.grid_x[i], self.grid_y[i], self.grid_z[i]
                            v = self.meshesRadex[i][transitionIdx][quantity]
                            
                            # appending the data and the value to a list to be used
                            # in making the interpolation function
                            
                            values.append( v )
                            grid_x.append( x )
                            grid_y.append( y )
                            grid_z.append( z )
                            
                # getting the data in the shape that is accepted by the interpolation construction        
                data   = np.array([grid_x, grid_y, grid_z], dtype = np.float64).T 
                values = np.array( values, dtype = np.float64)
                self.intensityGridInterp_f = self.construct3DInterpolationFunction(data   = data,
                                                                                   values = values,
                                                                                   log10  = True,
                                                                                   *args, **kwargs)

            #----------------------------emissions from the PDR slabs---------------------------------------
            if radiativeCodeType == 'pdr':

                #getting the column densities for all the models
                values  = np.ndarray(self.nMeshes, dtype = np.float64)
                #a temporary object used to calculate stuff from a pdr mesh object
                m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
                
                for i, data in enumerate(self.meshes):
                    
                    m.setData( data )
                    
                    quantity = self.parms['gridsInfo']['11']['quantity']
                    v = (1.0/(2.0*np.pi))*m.compute_integrated_quantity(quantity)

                    if v < 0:
                        self.logger.warning('negative intensitiy detected - setting it to 0')
                        v = np.float64(0)
                        #print self.grid_x[i], self.grid_y[i], self.grid_z[i], v 

                    values[i] = v 
        
                self.intensityGridInterp_f = self.construct3DInterpolationFunction(values = values,
                                                                                 log10  = True,
                                                                                 *args, **kwargs)
                    
    def show_pdr_mesh_quantity_grid(self, ranges = None, res = None, *args, **kwargs):
        """shows the surface temperature grid
            
           :todo: plot every other labeled contour as a tick in the colorbar
           
           :warning: there would be a memeory leak if things are plotted over and over
              the images, labels..etc..must be replaced and not new instance created
              and overplotted...take care of that later... 
        """
        panel = self.gui['maps2d']['00']
        
        grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                        res      = res,  
                                                        zSec     = self.pltGmSec, 
                                                        fInterp  = self.grdInterp_f, *args, **kwargs)

        grd = grd.T
        im00 = panel['axes'].imshow(grd, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )

        panel['contour'] = panel['axes'].contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        panel['axes'].clabel(panel['contour'], levels, fmt = '%.1f' )
        
        cbar = pyl.colorbar(im00, cax = panel['axesCbar'], ax = panel['axes'], 
                            orientation = 'horizontal', ticks = levels[::2], format = '%.1f')
        
        
    def showAbundancesGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the abundances grid"""
        panel = self.gui['maps2d']['01']
            

        grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                         res      = res,  
                                                         zSec     = self.pltGmSec, 
                                                         fInterp  = self.abunGridInterp_f, *args, **kwargs)
        
        grd = grd.T
        im01 = panel['axes'].imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )

        panel['contour'] = panel['axes'].contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        panel['axes'].clabel(panel['contour'],levels, fmt = '%.1f' )
        
        pyl.colorbar(im01, cax = panel['axesCbar'], ax = panel['axes'],
                     orientation = 'horizontal', ticks = levels[::2], format = '%.1f')

    def showColumnDensityGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the abundances grid"""
        panel = self.gui['maps2d']['10']

        grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                         res      = res,  
                                                         zSec     = self.pltGmSec, 
                                                         fInterp  = self.colDensGridInterp_f, *args, **kwargs)
        
        grd = grd.T
        im10 = panel['axes'].imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )
        
        panel['contour'] = panel['axes'].contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        panel['axes'].clabel(panel['contour'],levels, fmt = '%.1f' )
        
        pyl.colorbar(im10, cax = panel['axesCbar'], ax = panel['axes'],
                     orientation = 'horizontal', ticks = levels[::2], format = '%.1f')


    def showLineIntensityGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the line intensity grid"""            
        panel = self.gui['maps2d']['11']

        grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                         res      = res,  
                                                         zSec     = self.pltGmSec, 
                                                         fInterp  = self.intensityGridInterp_f, *args, **kwargs)
                                                        
        grd = grd.T
        im11 = panel['axes'].imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')

        mn = np.nanmin( grd[ np.where( np.fabs(grd) != np.inf ) ] )
        mx = np.nanmax( grd[ np.where( np.fabs(grd) != np.inf ) ] )
        nlevels = 10
        dl = (mx - mn)/nlevels            
        levels = np.arange( mn, mx, dl )
        
        # computing and displaying the contour levels to be plotted
        if self.parms['gridsInfo']['11']['showContours'] == True:
            panel['contour'] = panel['axes'].contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
            panel['axes'].clabel(panel['contour'],levels, fmt = '%.1f' )
        
        pyl.colorbar(im11, cax = panel['axesCbar'], ax = panel['axes'], orientation = 'horizontal', ticks = levels[::2], format = '%.1f')
        
        if self.parms['gridsInfo']['11']['type'] == 'radex':
            panel['axesCbar'].set_title('Intensity %s-%s' % (self.parms['gridsInfo']['11']['specStr'],self.parms['gridsInfo']['11']['transitionIndx'])  )
        elif self.parms['gridsInfo']['11']['type'] == 'pdr':
            panel['axesCbar'].set_title('Intensity %s-%s' % (self.parms['gridsInfo']['11']['quantity'][1], self.parms['gridsInfo']['11']['quantity'][3]) )
        
    def constructRadexDatabase(self, writeDb = None):
        """runs radex on all the models in self.meshes, and computes the line info 
           according to the parameters in self.parms['radex']. Once done computing,
           it stores all the generated info into self.infoAllRadex and :data:`meshesRadex`.
        
           :param bool writeDb: if this is set to true, self.infoAllRadex is writtent to a 
            file self.dirPath/infoAllRadex.db.info.'specStr' and self.meshesRadex is written
            to self.dirPath/meshesRadex.db.'specStr'
        """

        # defining the array holding the info about all the computed radex transitions
        # for each model in self.meshes and their location in the radex database ..etc..
        arxvRadexHdrDtype = np.dtype( self.arxvRadexHdrFormat() )
        infoAllRadex = np.ndarray( self.nMeshes, dtype = arxvRadexHdrDtype )
        #initiliizing to zero
        infoAllRadex[:] = 0
        meshesRadex = []
        
        radexObj = self.radexObj
        # mesh object which will be used as a utility to compute the stuff needed 
        # by radex to run, as in m.getRadexParameters() to get tKin, nDensCollisionPartners, molnDens
        m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
        
        specStr  = self.parms['radex']['specStr']
        xH2_Min  = self.parms['radex']['xH2_Min']
        
        for i in np.arange(self.nMeshes):
            
            self.logger.debug('pdr mesh index =  %d : ' % i)
            
            m.setData( self.meshes[i] )
            (gasTRadex, nColls, colDensThisSpec,) = m.getRadexParameters(speciesStr = specStr,
                                                                         threshold  = xH2_Min)
            # getting the collider densities in the same order of the supplied input spcie string list 
            nDensColls = [ nColls[collSpecStr] for collSpecStr in self.parms['radex']['collisionPartners'] ]
            collsStr   = list(self.parms['radex']['collisionPartners'])

            radexObj.setInFileParm('tKin', gasTRadex)
            radexObj.setInFileParm('collisionPartners', collsStr )
            radexObj.setInFileParm('nDensCollisionPartners', nDensColls )
            radexObj.setInFileParm('molnDens', colDensThisSpec)

            #remove colliders which are underabundant (below radex limits)
            radexObj.filterColliders()

            #all colliders have densities outsie the range accepted by radex
            #filling the data as empty and skipping this radex computation
            if len(radexObj.inFile['collisionPartners']) == 0:
                infoAllRadex[i]['info'][1] = 0 #no trainsitions
                meshesRadex.append(None)
                radexObj.printSetFlags()
                continue
            
            #----------------------------------running radex---------------------------------
            radexOutputMakesSense = False
            
            #parameters used for attempting to get radex to converge and vlaues make sense
            inFileOrig = radexObj.inFile.copy() # original parameters
            popDensTol = self.parms['radex']['popDensSumTol']  # ideally sum pop dens lower should be 1, for values less than tol
                                                              # the parameters are changed by 'changeFac' and another attempt is done
                                                              # until the sum of pop dens is close to 1 up to tol. 
            changeFac  = self.parms['radex']['changeFracTrial']  # the relative amount by which the input parameters to radex are changed
            nMaxTrial  = self.parms['radex']['nMaxTrial']     # number of time random changes are done before giving up
            
            nTried     = 0
            while not radexOutputMakesSense:
                self.logger.debug('radex trial %d' % nTried)
                radexObj.setDefaultStatus()
                status = radexObj.run( checkInput = True, verbose = False )

                if self.parms['radex']['checkOutputIntegrity']:
                    # checking the integrity of the solution in case it converged
                    if radexObj.getStatus() & radexObj.FLAGS['SUCCESS']:
                        
                        totalPopDensLower = np.sum(radexObj.transitions[:]['pop_down'])
                        if np.fabs(self.parms['radex']['popDensSumExpected'] - totalPopDensLower) > popDensTol:
                            
                            self.logger.warn('====> pop dens does NOT make sense, trial %d\n====> total pop dense lower = %20.10f' % (nTried, totalPopDensLower) )
    
                            #change in input parameters by a tiny amount 0.1% and run again
                            #(this is a numerical bug in radex which fails sometimes)
                            radexObj.inFile['tKin']     = inFileOrig['tKin']*(1.0 + np.random.rand()*changeFac)
                            radexObj.inFile['molnDens'] = inFileOrig['molnDens']*(1.0 + np.random.rand()*changeFac)
                            for j, dens in enumerate(radexObj.inFile['nDensCollisionPartners']):
                                radexObj.inFile['nDensCollisionPartners'][j] = inFileOrig['nDensCollisionPartners'][j]*(1.0 + np.random.rand()*changeFac)
                        else:
                            radexOutputMakesSense = True
                            print '====>', 'pop dense make sense, trial %d' % nTried
                    else: # radex did not succeed
                        break  # here radexOutputMakesSense woulbe be False
                    
                    nTried += 1
                    
                    if nTried >= nMaxTrial: #no meaningful output after nTried trials
                        break # here radexOutputMakesSense woulbe be False
                        
                else: #no need to check for validity of the radex pop densities
                    radexOutputMakesSense = True
                    break
            #------------------------------done running radex----------------------------------------
            
            #setting the basic radex run info into the attribute
            #-----------------------------------------------------------
            #saving the number of transitions into radex info attribute
            infoAllRadex[i]['info'][0] = i
            if radexObj.flagSet('ERROR') or (radexOutputMakesSense is False): # radex failed, no transition data available
                infoAllRadex[i]['info'][1] = 0 #no trainsitions
                meshesRadex.append(None)
                radexObj.printSetFlags()
            else: # radex succeeded, wirint the transition data
                infoAllRadex[i]['info'][1] = radexObj.nTransitions
                meshesRadex.append(radexObj.transitions)

            #infoAllRadex[i]['info'][2] = NOT SET HERE, IT IS SET WHEN WRITING THE DB TO A FILE
            
            #setting the status into the attribute
            infoAllRadex[i]['info'][3] = status
            #appending the transition data of this mesh to the attribute 
            #-----------------finished saving the info------------------
            
            
            if radexObj.flagSet('SUCCESS'):
                self.logger.debug('radexGrid : converged with no warnings')
            else:
                self.logger.debug('radexGrid : converged with warnings')
                self.logger.debug('------------------------------------')
                print radexObj.getWarnings()
                self.logger.debug('------------------------------------')
                continue

            self.logger.debug('---------------------------------------------------------')
            

        #copying the mesh parameters self.infoAll[:]['parms] to self.infoAllRadex[:]['parms']
        for i in np.arange( self.nMeshes ):
            infoAllRadex[i]['parms'][:] = self.infoAll[i]['parms']

        self.infoAllRadex = infoAllRadex
        self.meshesRadex = meshesRadex
        
        if writeDb == True:
            self.writeDbRadex()
        
    
    def save_radex_grids(self, relativeDirPath = None, basename = None, transitionInds = None, 
                         quantity = None, fileFormat = None, *args, **kwargs):
        """method that saves the radex grids (emission grids) for now into files.
        
           :param string fileFormat: When this is set to **'numpytxt'**, the grid 
            is dumped written to an ascii file using numpy.savetxt, one row on 
            each line.  When this parameter is set to **'3column'**, the 
            coordinates of the centroid of the grid cell and the value 
            (x_i, y_i, v_i) are written on each line.  
        """
        
        ranges          = self.parms['plotRanges'] 
        res             = self.resPltGrids 

        #self.set_default_attributes()
        #self.set_grid_axes_quantity_values()

        #dumping the parameters of meshutils used in generating the data
        timeStr = ctime().replace(' ','-')
        parmsOutputFile = self.dirPath + relativeDirPath + 'parms.out' #-' + timeStr 
        parmsOut = open(parmsOutputFile, 'wb')
        pickle.dump(self.parms, parmsOut)
        parmsOut.close()
        
        filesInfo = ()
        
        for zSec in self.grid_z_unique:
            
            for transitionIdx in transitionInds:
            
                nValidRadexInSec = 0 # number of meshes wiht radex data in this section
    
                ####################################################################################################################
                #getting the inerpolation function for this section in z
                ####################################################################################################################
                values, grid_x, grid_y, grid_z = [], [], [], []
                
                #looping over all the radex info of all the mesehes and getting the
                #transition data to be displayed
                for i in np.arange(self.nMeshes):
                    
                    if self.meshesRadex[i] != None:
                        
                        x, y, z = self.grid_x[i], self.grid_y[i], self.grid_z[i]
                        v = self.meshesRadex[i][transitionIdx][quantity]
                        
                        # appending the data and the value to a list to be used
                        # in making the interpolation function
                        
                        values.append( v )
                        grid_x.append( x )
                        grid_y.append( y )
                        grid_z.append( z )
                        
                        # getting the transition string
                        transStrng = '%s-%s' % (self.meshesRadex[i][transitionIdx]['upper'],self.meshesRadex[i][transitionIdx]['lower'])
                        nValidRadexInSec += 1
                
                if nValidRadexInSec > 0:
                    # getting the data in the shape that is accepted by the interpolation construction        
                    data   = np.array([grid_x, grid_y, grid_z], dtype = np.float64).T
                    values = np.array( values, dtype = np.float64)
                    intensityGridInterp_f = self.construct3DInterpolationFunction(data   = data,
                                                                                  values = values,
                                                                                  log10  = True,
                                                                                  *args, **kwargs)
        
                    ####################################################################################################################
                    #interpolating the grids
                    ####################################################################################################################
                    
        
                    grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                                     res      = res,  
                                                                     zSec     = zSec, 
                                                                     fInterp  = intensityGridInterp_f,
                                                                     *args, **kwargs)
                                                                
                    grd = grd.T
    
                    filename  = self.dirPath + relativeDirPath + basename
                    filename += '-' + quantity + '-' + transStrng
                    filename += '-log10zSec=%f' % zSec
                    filename += '.txt'
                    
                    #------------setting the format of the data to be written to a text file--------
                    if fileFormat == 'numpytxt': #writing the grid as a 2D table
                        fileData = grd
                    elif fileFormat == '3columns':
                        n = grd.size
                        # converting the grids into a 3 column format 
                        fileData = np.array([grdx.reshape(n), grdy.reshape(n), grd.reshape(n)]).T
                    #--------------------------------------------------------------------------------
                    
                    np.savetxt(filename, fileData, fmt = '%1.4e')
                    self.logger.debug('wrote the file %s' % filename)
                    
                    thisFileInfo = {'zSec'       : zSec,
                                    'transition' : transStrng,
                                    'filename'   : filename,
                                    }
                    filesInfo += (thisFileInfo,)
                else:
                    self.logger.warn('no valid radex data found in this section')
                    
        #dumping information about the files into a pickle file
        fileInfoPath = self.dirPath + relativeDirPath + 'filesInfo.out' #-' + timeStr
        fileInfoFobj = open(fileInfoPath, 'wb')
        pickle.dump(filesInfo, fileInfoFobj)
        fileInfoFobj.close()

    def save_PDR_emission_grids(self, relativeDirPath = None, basename = None, transitions = None, quantity = None, *args, **kwargs):
        """method that saves the emission grids computed from the PDR models (emission grids) for now into files"""
        
        ranges = self.parms['plotRanges'] 
        res    = self.resPltGrids 
        
        #dumping the parameters of meshutils used in generating the data
        timeStr = ctime().replace(' ','-')
        parmsOutputFile = self.dirPath + relativeDirPath + 'parms.out' 
        parmsOut = open(parmsOutputFile, 'wb')
        pickle.dump(self.parms, parmsOut)
        parmsOut.close()
        
        filesInfo = ()

        specStr = self.parms['gridsInfo']['11']['quantity'][1]
        quantity_PDR_em = ['fineStructureCoolingComponents', specStr, quantity, '']
        
        #a temporary object used to calculate stuff from a pdr mesh object
        m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)

        for zSec in self.grid_z_unique:
            
            for transStrng in transitions:
            
                ####################################################################################################################
                #getting the inerpolation function for this section in z
                ####################################################################################################################

                values  = np.ndarray(self.nMeshes, dtype = np.float64)
                quantity_PDR_em[3] = transStrng 
                
                #looping over all the radex info of all the mesehes and getting the
                #transition data to be displayed
                for i, data in enumerate(self.meshes):
                    
                    m.setData(data)
                    dxSlabs = m.compute_dx()
                    
                    q = fetchNestedDtypeValue(data, quantity_PDR_em )
                    v = np.sum( q*dxSlabs ) / (2.0 * np.pi)

                    if v < 0:
                        self.logger.warning('negative intensitiy detected - setting it to 0')
                        v = np.float64(0)
                        #print self.grid_x[i], self.grid_y[i], self.grid_z[i], v 

                    values[i] = v
                #--------done extracting the data from the meshes------- 
                
                # getting the data in the shape that is accepted by the interpolation construction        
                data   = np.array([self.grid_x, self.grid_y, self.grid_z], dtype = np.float64).T
                intensityGridInterp_f = self.construct3DInterpolationFunction(data   = data,
                                                                              log10  = True,
                                                                              values = values,
                                                                              *args, **kwargs)
        
                ####################################################################################################################
                #interpolating the grids
                ####################################################################################################################
        
                grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                                 res      = res,  
                                                                 zSec     = zSec, 
                                                                 fInterp  = intensityGridInterp_f,
                                                                 *args, **kwargs)
                                                                
                grd = grd.T

                filename  = self.dirPath + relativeDirPath + basename
                if quantity == 'rate': # just to make the filenames for emission the same as radex
                    filename += '-' + 'fluxcgs' + '-' + transStrng
                else:
                    filename += '-' + quantity + '-' + transStrng
                filename += '-log10zSec=%f' % zSec
                filename += '.txt'
                    
                np.savetxt(filename, grd, fmt = '%1.4e')
                self.logger.debug('wrote the file %s' % filename)
                    
                thisFileInfo = {'zSec'       : zSec,
                                'transition' : transStrng,
                                'filename'   : filename,
                               }
                filesInfo += (thisFileInfo,)
            #--------end loop over transitions-------------
        #-----end loop over zSections------------
        
        #dumping information about the files into a pickle file
        fileInfoPath = self.dirPath + relativeDirPath + 'filesInfo.out'# + timeStr
        fileInfoFobj = open(fileInfoPath, 'wb')
        pickle.dump(filesInfo, fileInfoFobj)
        fileInfoFobj.close()

    def save_PDR_quantity_grids(self, relativeDirPath = None, basename = None, quantity = None, slabIdx = None, *args, **kwargs):
        """method that saves a grid computed from the PDR models for now into files
        keywords : save, grid, pdr, quantity, slab
        """
        
        ranges = self.parms['plotRanges'] 
        res    = self.resPltGrids 
        
        #dumping the parameters of meshutils used in generating the data
        timeStr = ctime().replace(' ','-')
        parmsOutputFile = self.dirPath + relativeDirPath + 'parms.out' 
        parmsOut = open(parmsOutputFile, 'wb')
        pickle.dump(self.parms, parmsOut)
        parmsOut.close()
        
        filesInfo = ()
        
        for zSec in self.grid_z_unique:
            
            ####################################################################################################################
            #getting the inerpolation function for this section in z
            ####################################################################################################################
            grdInterp_f = self.construct3DInterpolationFunction(quantity = quantity, 
                                                                slabIdx  = slabIdx,
                                                                log10    = True,
                                                                *args, **kwargs)
            #--------done extracting the data from the meshes------- 
            
            ####################################################################################################################
            #interpolating the grids
            ####################################################################################################################
    
            grd, grdx, grdy = self.computeInterpolated2DGrid(ranges   = ranges,
                                                             res      = res,  
                                                             zSec     = zSec, 
                                                             fInterp  = grdInterp_f,
                                                             *args, **kwargs)
                                                            
            grd = grd.T

            filename  = self.dirPath + relativeDirPath + basename
            filename += '-log10zSec=%f' % zSec
            filename += '-%d-' % slabIdx
            filename += '.txt'
                
            np.savetxt(filename, grd, fmt = '%1.4e')
            self.logger.debug('wrote the file %s' % filename)
                
            thisFileInfo = {'zSec'       : zSec,
                            'slabIdx'    : slabIdx, 
                            'filename'   : filename,
                           }
            filesInfo += (thisFileInfo,)
        #-----end loop over zSections------------
        
        #dumping information about the files into a pickle file
        fileInfoPath = self.dirPath + relativeDirPath + 'filesInfo.out'# + timeStr
        fileInfoFobj = open(fileInfoPath, 'wb')
        pickle.dump(filesInfo, fileInfoFobj)
        fileInfoFobj.close()
                    
    #################################################################################
    def setupGui(self):
        """the method which setup the layout of the gui by defining a dict holding
        all the object and variables related to it. Returns the dict object and make
        use of self.parms"""
        
        gui = {}
        gui['figure']  = pyl.figure( figsize = (14,14) )
        gui['widgets'] = {} 
        gui['maps2d']  = {} # the grids for a certain section
        gui['ax2d']    = {} # the mesh points in a certain section
        gui['ax3d']    = {} # all the mehs points in 3D
        gui['AV']      = {} # plots as a function of Av for a certail slab
        
        # the grid plot in n,G0 showing the points where models are present in the section
        ax2d = {}
        #-----------------------------------------------------------------------------------
        ax2d['axes']  = gui['figure'].add_axes( [0.05, 0.7, 0.12, 0.12] )
        ax2d['pts1'], = ax2d['axes'].plot( [0], [0], color = 'b', marker = 'o', linestyle = '')
        ax2d['pts2'], = ax2d['axes'].plot( [1], [1], color = 'r', marker = 'o', linestyle = '')
        ax2d['pts3'], = ax2d['axes'].plot( [1], [1], color = 'w', marker = 'o', linestyle = '')
        ax2d['axes'].set_title('$\log_{10} n_{gas} = $ %4.2f\n$\log_{10} G_0 =$ %4.2f\n$\log_{10} \Gamma_{mech} = $  %5.2f\n' % (0, 0, 0) )
        ax2d['axes'].set_xlim(self.parms['plotRanges'][0])
        ax2d['axes'].set_ylim(self.parms['plotRanges'][1])
        ax2d['axes'].set_xlabel('$log_{10} n_{gas}$')
        ax2d['axes'].set_ylabel('$log_{10} G_0$')
        #-----------------------------------------------------------------------------------        
        gui['ax2d'] = ax2d

        # the axes to plot in the 3D grid (default axes are n,G0,gmech) showing the points
        # where models are present in the database
        ax3d = {}
        #-----------------------------------------------------------------------------------
        ax3d['axes'] = gui['figure'].add_subplot(111, projection='3d')
        ax3d['axes'].set_position((0.2, 0.68, 0.3, 0.3))
        #-----------------------------------------------------------------------------------        
        gui['ax3d'] = ax3d
                        
        # defining the axes which will be used to select the section in Z to show
        zSecSelector = {}
        #----------------------------------------------------------------------
        zSecSelector['axes'] = gui['figure'].add_axes( [0.05, 0.6, 0.4, 0.02] )
        zSecSelector['axes'].set_title('z section selector')
        zSecSelector['axes'].set_xlim( [self.ranges[2][0], self.ranges[2][1] ] )
        zSecSelector['axes'].set_ylim( [0, 1] )
        zSecSelector['axes'].get_yaxis().set_ticks([])
        zSecSelector['point'], = zSecSelector['axes'].plot([],[], 'r+', markersize = 20)
        zSecSelector['point'].set_xdata( self.pltGmSec )
        zSecSelector['point'].set_ydata( 0.5 )
        zSecSelector['pointsUnique'], = zSecSelector['axes'].plot([],[], 'yo')
        #----------------------------------------------------------------------
        gui['widgets']['zSecSelector'] = zSecSelector

        # defining the axes which will be used to select the transition index to be displayed in the emission panel
        transitionSelector = {}
        #----------------------------------------------------------------------
        transitionSelector['axes'] = gui['figure'].add_axes( [0.05, 0.53, 0.2, 0.02] )
        transitionSelector['axes'].set_title('transition selector')
        transitionSelector['axes'].set_xlim( [0, 40] )
        transitionSelector['axes'].set_ylim( [0, 1] )
        transitionSelector['axes'].get_yaxis().set_ticks([])
        transitionSelector['point'], = transitionSelector['axes'].plot([],[], 'r+', markersize = 20)
        if self.parms['gridsInfo']['11']['type'] == 'radex':
            transitionSelector['point'].set_xdata( self.parms['gridsInfo']['11']['transitionIndx'] )
            transitionSelector['point'].set_ydata( 0.5 )
        #----------------------------------------------------------------------
        gui['widgets']['transitionSelector'] = transitionSelector
        
        # defining the grid where the maps will be plotted (maps2d)
        #----------------------------------------------------------
        left  = 0.55
        bott  = 0.45
        sz    = 0.20
        vSpace = 0.06
        hSpace = 0.02
        
        # 0,0
        maps2d_00 = {}
        maps2d_00['axes'] = gui['figure'].add_axes([left, bott + sz + vSpace, sz, sz])
        maps2d_00['axesCbar'] = gui['figure'].add_axes([left, bott + sz + sz + vSpace + 0.017, 0.2, 0.01 ])
        maps2d_00['axesCbar'].set_title( '%s-%s' % tuple(self.parms['gridsInfo']['00']['quantity']) )
        maps2d_00['contour'] = None
        maps2d_00['axes'].set_ylabel( '$log_{10} G_0$' )
        for tick in maps2d_00['axes'].xaxis.get_major_ticks():
            tick.label1On = False
        maps2d_00['axes'].set_xlim(self.ranges[0][0], self.ranges[0][1])
        maps2d_00['axes'].set_ylim(self.ranges[1][0], self.ranges[1][1])
        gui['maps2d']['00'] = maps2d_00 
        
        # 0,1
        maps2d_01 = {}
        maps2d_01['axes'] = gui['figure'].add_axes([left + sz + hSpace, bott + sz + vSpace, sz, sz])
        maps2d_01['axesCbar'] = gui['figure'].add_axes( [left + sz + hSpace, bott + sz + sz + vSpace + 0.017, 0.2, 0.01] )
        maps2d_01['axesCbar'].set_title('$n_X/(n_H + 2 n_{H_2}$')
        maps2d_01['contour'] = None
        for tick in maps2d_01['axes'].xaxis.get_major_ticks():
            tick.label1On = False
        for tick in maps2d_01['axes'].yaxis.get_major_ticks():
            tick.label1On = False
        maps2d_00['axes'].set_xlim(self.ranges[0][0], self.ranges[0][1])
        maps2d_00['axes'].set_ylim(self.ranges[1][0], self.ranges[1][1])
        gui['maps2d']['01'] = maps2d_01 
        
        # 1,0
        maps2d_10 = {}
        maps2d_10['axes'] = gui['figure'].add_axes([left, bott, sz, sz])
        maps2d_10['axesCbar'] = gui['figure'].add_axes([left, bott + sz + 0.017, 0.2, 0.01 ])
        maps2d_10['axesCbar'].set_title('N(X)')
        maps2d_10['contour'] = None
        maps2d_10['axes'].set_xlabel( '$log_{10} n_{gas}$' )
        maps2d_10['axes'].set_ylabel( '$log_{10} G_0$' )
        maps2d_10['axes'].set_xlim(self.ranges[0][0], self.ranges[0][1])
        maps2d_10['axes'].set_ylim(self.ranges[1][0], self.ranges[1][1])
        (maps2d_10['axes'].yaxis.get_major_ticks())[-1].label1On = False
        (maps2d_10['axes'].xaxis.get_major_ticks())[-1].label1On = False
        gui['maps2d']['10'] = maps2d_10 
        
        # 1,1
        maps2d_11 = {}
        maps2d_11['axes'] = gui['figure'].add_axes([left + sz + hSpace, bott, sz, sz])
        maps2d_11['axesCbar'] = gui['figure'].add_axes( [left + sz + hSpace, bott + sz + 0.017, 0.2, 0.01] )
        if self.parms['gridsInfo']['11']['type'] == 'radex':
            maps2d_11['axesCbar'].set_title('Intensity %s-%s' % (self.parms['gridsInfo']['11']['specStr'],self.parms['gridsInfo']['11']['transitionIndx'])  )
        elif self.parms['gridsInfo']['11']['type'] == 'pdr':
            maps2d_11['axesCbar'].set_title('Intensity %s-%s' % (self.parms['gridsInfo']['11']['quantity'][1], self.parms['gridsInfo']['11']['quantity'][3]) )
        maps2d_11['contour'] = None
        maps2d_11['axes'].set_xlabel( '$log_{10} n_{gas}$' )
        for tick in maps2d_11['axes'].yaxis.get_major_ticks():
            tick.label1On = False
        maps2d_11['axes'].set_xlim(self.ranges[0][0], self.ranges[0][1])
        maps2d_11['axes'].set_ylim(self.ranges[1][0], self.ranges[1][1])

        gui['maps2d']['11'] = maps2d_11 
        
        
        # modifying ticks on the colorbars   
        for cbarAxs in [gui['maps2d']['00']['axesCbar'],
                        gui['maps2d']['01']['axesCbar'],
                        gui['maps2d']['10']['axesCbar'],
                        gui['maps2d']['11']['axesCbar']]:
            for tick in cbarAxs.xaxis.get_major_ticks():
                tick.label1On = True
                tick.label2On = False
            for tick in cbarAxs.yaxis.get_major_ticks():
                tick.label1On = False
                tick.label2On = False

        # the axes to plot in the mesh stuff as a function of AV.
        AV = {}
        #-----------------------------------------------------------------------------------
        left  = 0.07
        bott  = 0.05
        sz    = 0.20
        vSpace = 0.01
        hSpace = 0.01
        AV['00'] = {'axes' : gui['figure'].add_axes([left              , bott + sz + vSpace, sz, sz])}
        AV['01'] = {'axes' : gui['figure'].add_axes([left + sz + hSpace, bott + sz + vSpace, sz, sz])}
        AV['10'] = {'axes' : gui['figure'].add_axes([left              , bott              , sz, sz])}
        AV['11'] = {'axes' : gui['figure'].add_axes([left + sz + hSpace, bott              , sz, sz])}
        #-----------------------------------------------------------------------------------        
        gui['AV'] = AV


        # the axes stuff where the radex output wil be plotted
        radex = {}
        #-----------------------------------------------------------------------------------
        left  = 0.65
        bott  = 0.07
        sz    = 0.08
        vSpace = 0.0
        radex['0'] = {'axes' : gui['figure'].add_axes([left, bott + 3*sz + vSpace , 3*sz, sz])}
        radex['1'] = {'axes' : gui['figure'].add_axes([left, bott + 2*sz + vSpace , 3*sz, sz])}
        radex['2'] = {'axes' : gui['figure'].add_axes([left, bott + 1*sz + vSpace , 3*sz, sz])}
        radex['3'] = {'axes' : gui['figure'].add_axes([left, bott + 0*sz          , 3*sz, sz])}
        radex['title1'] = pyl.figtext(0.65, 0.4, '')
        radex['title2'] = pyl.figtext(0.9, 0.2 , '')
        #-----------------------------------------------------------------------------------        
        gui['radex'] = radex
        
        #-------------------------------------------------
        return gui
    
        
    def plotGrids(self, *args, **kwargs):
        """Main method for exploring the meshes in the database.
        
        :todo: change resGrids to a [res_x, res_y] insteads of it being just a scalar.
        """
        
        self.set_attributes(**kwargs)
        self.set_default_attributes()
        
        resGrids = self.parms['gridsRes']
        #setting the default zSec
        self.pltGmSec   = np.min(self.grid_z)
        self.logger.debug('set the current plotting section to %f, the minimum of the self.grid_z' % self.pltGmSec)
        self.radexParms = self.parms['radex']
        
        #setting up the gui attribute
        self.gui = self.setupGui()

        if self.parms['showGrids']:
            # getting the interpolation function which will be used to display the 2D grids
            self.computeAndSetInterpolationFunctions(*args, **kwargs)        
    
            # defining and intialising the ndmesh objects which will be used
            # for computing the grid properties and then displayed                                
                                    #xaxis    yaxis 
            self.grds = [ [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [self.ranges[0], self.ranges[1]], fill = 0.0 ), 
                           ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [self.ranges[0], self.ranges[1]], fill = 0.0 )], 
                          [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [self.ranges[0], self.ranges[1]], fill = 0.0 ),
                           ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [self.ranges[0], self.ranges[1]], fill = 0.0 )] ]  
    
            self.resPltGrids = [resGrids, resGrids] 
         
        # setting up a dummy mesh object to use its plotting functionalities
        self.mshTmp = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
        self.mshTmp.setFigureObjects(figObj = self.gui['figure'], 
                                     axObj  = np.array([
                                                        [self.gui['AV']['00']['axes'],self.gui['AV']['01']['axes']],
                                                        [self.gui['AV']['10']['axes'],self.gui['AV']['11']['axes']]
                                                        ])
                                     )
        self.mshTmp.setupFigures(plot_Av_range = self.parms['meshPltAvRng'])
        self.mshTmp.set_metallicity(self.metallicity)

        # plotting all the points in the archive in 3D
        self.plot_3D_grid_points(figure = self.gui['figure'], 
                                 axes   = self.gui['ax3d']['axes'], 
                                 ranges = self.ranges)
        
        # setting up the axes to plot the radex output for each selected model
        self.pltRadex = np.array([
                                  self.gui['radex']['0']['axes'],
                                  self.gui['radex']['1']['axes'],
                                  self.gui['radex']['2']['axes'],
                                  self.gui['radex']['3']['axes']
                                  ]
                                 )
        
        # getting and plotting the unique sections in z
        self.grid_z_unique = self.getUnique_grid_z_sections(1e-13)
        self.gui['widgets']['zSecSelector']['pointsUnique'].set_xdata(self.grid_z_unique)
        self.gui['widgets']['zSecSelector']['pointsUnique'].set_ydata(self.grid_z_unique*0.0 + 0.5)

        # attaching mouse click event to fig 1
        cid = self.gui['figure'].canvas.mpl_connect('button_press_event', self.on_b1_down)
        
        self.plotThisSec()
        
        # displaying
        pyl.draw()
        self.logger.debug('browing data....')

    def plotThisSec(self, *args, **kwargs):
        """updates the plotes once the z section value is changed
        """
        self.logger.debug('updating the grids for zSec = %f' % self.pltGmSec)

        #plotting the meshes n and G0 whose data is available for
        #this section
        indsThisSec = np.nonzero( np.fabs(self.grid_z - self.pltGmSec) < 1e-13 )[0]            
        self.gui['ax2d']['pts1'].set_xdata( self.grid_x[indsThisSec] )
        self.gui['ax2d']['pts1'].set_ydata( self.grid_y[indsThisSec] )
        
        #overplotting the radex points available for this section
        if self.parms['gridsInfo']['11']['show']:
            cond1 = np.fabs(self.grid_z - self.pltGmSec) < 1e-13
            for i in np.arange(cond1.size): 
                c0 = (self.infoAllRadex[i]['info'][2] == -1)
                cond1[i] *= c0
            indsThisSec = np.nonzero( cond1 )[0]
            self.gui['ax2d']['pts3'].set_xdata( self.grid_x[indsThisSec] )
            self.gui['ax2d']['pts3'].set_ydata( self.grid_y[indsThisSec] )
        
        # plotting the grids
        #-------------------
        if self.parms['showGrids']:
            # temperature grid (top left grid)
            if self.parms['gridsInfo']['00']['show']:
                self.show_pdr_mesh_quantity_grid(ranges = self.parms['plotRanges'], res = self.resPltGrids, *args, **kwargs)
    
            # abundances (top left grid)
            if self.parms['gridsInfo']['01']['show']:
                self.showAbundancesGrid(ranges = self.parms['plotRanges'], res = self.resPltGrids, *args, **kwargs)
    
            # column densities (bottom left grid)
            if self.parms['gridsInfo']['10']['show']:
                self.showColumnDensityGrid(ranges = self.parms['plotRanges'], res = self.resPltGrids, *args, **kwargs)
    
            # line intensity (bottom right grid)
            if self.parms['gridsInfo']['11']['show']:
                self.showLineIntensityGrid(ranges = self.parms['plotRanges'], res = self.resPltGrids, *args, **kwargs)
        
    def on_b1_down(self, event):
        """method called on the event of a mouse button down in self.fig. See self.setupGui()
           for the layout of the window.
        """
       
        ti = time()

        # get the x and y coords, flip y from top to bottom
        xd, yd = event.xdata, event.ydata
        
        if event.button == 1:

            clickedInAxes = False
            
            # getting the value of the section in z to display
            if event.inaxes is self.gui['widgets']['zSecSelector']['axes']:
                clickedInAxes = True
                # setting the section closest to the data available
                inds = np.argmin(np.fabs( self.grid_z - xd ) )
                self.pltGmSec = self.grid_z[inds]
                self.gui['widgets']['zSecSelector']['point'].set_xdata( self.pltGmSec )
                # updating the title
                newTitle = 'z section selector : %f (%e)' % (self.pltGmSec, 10**self.pltGmSec)
                self.gui['widgets']['zSecSelector']['axes'].set_title(newTitle)
                                
                #for all the axes of the 2d maps
                if self.parms['showGrids']:
                    for panel in self.gui['maps2d'].values():
                        #deleting the countour lines (if they are found)
                        if 'collections' in dir(panel['contour']):
                            for c in panel['contour'].collections:
                                paths = c.get_paths()
                                del paths[:]
                            #deleting the labels
                            for txt in panel['contour'].labelTexts:
                                txt.set_text('')
                        #deleting the images
                        del panel['axes'].images[:]
                    
                self.plotThisSec() #;;; rename this to update 2D grids
                pyl.draw()

            # getting the value of the transition to display in the emission panel
            if event.inaxes is self.gui['widgets']['transitionSelector']['axes']:
                clickedInAxes = True
                # setting the section closest to the data available
                self.parms['gridsInfo']['11']['transitionIndx'] = np.round(xd)
                self.gui['widgets']['transitionSelector']['point'].set_xdata( self.parms['gridsInfo']['11']['transitionIndx'] )
                # updating the title
                newTitle = 'transition selector, index = %d' % self.parms['gridsInfo']['11']['transitionIndx']
                self.gui['widgets']['transitionSelector']['axes'].set_title(newTitle)

                #for all the axes of the 2d maps
                for panel in self.gui['maps2d'].values():
                    #deleting the countour lines (if they are found)
                    if 'collections' in dir(panel['contour']):
                        for c in panel['contour'].collections:
                            paths = c.get_paths()
                            del paths[:]
                        #deleting the labels
                        for txt in panel['contour'].labelTexts:
                            txt.set_text('')
                    #deleting the images
                    del panel['axes'].images[:]

                self.computeAndSetInterpolationFunctions()
                self.plotThisSec() #;;; rename this to update 2D grids
                pyl.draw()

            # getting the coordinates inthe grid to display as a single mesh as a function of Av
            # and maybe display the radex calculations                
            if event.inaxes is self.gui['ax2d']['axes']:
                clickedInAxes = True
                self.logger.debug('####################################button 1 clicked###########################################')
                indMin = self.get_mesh_index(x = xd, y = yd, z = self.pltGmSec)
                self.mshTmp.setData( self.meshes[indMin] )
                
                #updating the title of ['ax2d']['axes']
                strng = '$\log_{10} n_{gas} = $ %4.2f\n$\log_{10} G_0 =$ %4.2f\n$\log_{10} \Gamma_{mech} = $  %5.2f\n' %  (np.log10(self.mshTmp.data['hdr']['nGas']), np.log10(self.mshTmp.data['hdr']['G0']), np.log10(self.mshTmp.data['hdr']['gammaMech']))
                strng += '$\Gamma_{mech}$ = %.2e\n' % self.mshTmp.data['hdr']['gammaMech']   
                self.gui['ax2d']['axes'].set_title(strng)
                
                self.gui['ax2d']['pts2'].set_xdata( self.grid_x[indMin] )
                self.gui['ax2d']['pts2'].set_ydata( self.grid_y[indMin] )
                self.gui['ax2d']['pts2'].set_color('r')
                                    
                self.mshTmp.plot() #plotting the PDR mesh curves in the panels
                                    
                if self.parms['radex']['use']:
                    self.computeAndSetRadexCurves(meshObj = self.mshTmp)
                
                #;;;remove later (to plot the emissions along the grid)
                """
                pyl.figure(figsize=(6,6))
                pyl.plot(self.mshTmp.data['state']['Av'],
                         self.mshTmp.data['fineStructureCoolingComponents']['C']['rate']['1-0'],'r')
                pyl.plot(self.mshTmp.data['state']['Av'],
                         self.mshTmp.data['fineStructureCoolingComponents']['C']['rate']['2-1'],'b')
                """
                #;;;end remove later
                pyl.draw()

            # setting the Av of the position clicked on the plots of the current mesh
            # to the chemical netowrk object.
            if event.inaxes is self.gui['AV']['00']['axes'] or event.inaxes is self.gui['AV']['01']['axes'] or event.inaxes is self.gui['AV']['10']['axes'] or event.inaxes is self.gui['AV']['11']['axes']:
                
                clickedInAxes = True
                Av_clicked = xd

                self.logger.debug('####################################button 1 clicked###########################################')
                self.logger.debug('Av clicked = %.4f' % Av_clicked)
                
                #getting the index along the mesh which is closest in Av to the place
                #where the mouse was clicked, the closes values will be used in the
                #chemical netowrk of the mesh
                mshAv = self.mshTmp.data['state']['Av']
                indAv = np.argmin( np.fabs(mshAv - Av_clicked) )
                
                #the values to be assigned to the chemical network
                Av_use = self.mshTmp.data['state']['Av'][indAv]
                gasT_use = self.mshTmp.data['state']['gasT'][indAv]
                abun_use = self.mshTmp.data['state']['abun'][:,indAv]

                self.logger.debug('Av   used = %.4f' % Av_use)
                self.logger.debug('Tgas used = %.4f' % gasT_use)
                
                #setting the parameters of the gas and the environment at that slab
                self.mshTmp.chemNet.set_environment_state(gasT_use, 
                                                          self.parms_used_in_PDR_models['zeta'], 
                                                          Av_use,
                                                          self.parms_used_in_PDR_models['albedo'],
                                                          self.mshTmp.data['hdr']['nGas'],
                                                          self.mshTmp.data['hdr']['G0'])
                
                #plotting the vertical lines on the gui indicating the positions
                #in the slab used for the chemistry
                self.mshTmp.plot_v_lines_used_in_chemnet()
                pyl.draw()
                
                #set the abundances at that slab
                self.mshTmp.chemNet.set_abundances(fromArray = abun_use) 
                self.mshTmp.chemNet.compute_rxn_constants()
                self.mshTmp.chemNet.compute_rxn_rates()

                self.logger.debug('set the environment variable to the chemical netowrk of the mesh.')
            
            if clickedInAxes == True:
                tf = time()
                self.logger.debug('time elapased after click : %f\n' % (tf - ti))
        
    def computeAndSetRadexCurves(self, meshObj = None):
        """This is a utilty method (make it a private method), for populating the radex axes
          with the plots for the specie passed in the global parameter self.parms['radex']['specStr'].
          It takes as a paremeter the :data:`mesh` which will be used for doing the radex computations.
        """
        
        msh = meshObj
        #assigning the axes in which the radex curves will be plotted
        self.radexObj.setupPlot(nx = 1, fig = self.gui['figure'], axs = self.pltRadex)
        self.radexObj.set_logger(self.logger)

        (gasTRadex, nColls, colDensThisSpec,) = msh.getRadexParameters(speciesStr = self.parms['radex']['specStr'], 
                                                                       threshold  = self.parms['radex']['xH2_Min'])

        # getting the collider densities in the same order of the supplied input spcie string list 
        nDensColls = [ nColls[collSpecStr] for collSpecStr in self.parms['radex']['collisionPartners'] ]
        collsStr   = list(self.parms['radex']['collisionPartners'])
        #print 'input coll species', self.radexParms['collisionPartners'] 
        #print 'nColls after putting them in the right order = ', nDensColls
        
        self.radexObj.setInFileParm('specStr', self.parms['radex']['specStr'])
        self.radexObj.setInFileParm('tKin', gasTRadex)
        self.radexObj.setInFileParm('collisionPartners', collsStr )
        self.radexObj.setInFileParm('nDensCollisionPartners', nDensColls )
        self.radexObj.setInFileParm('molnDens', colDensThisSpec)
        
        self.radexObj.filterColliders()
            
        if len(self.radexObj.inFile['collisionPartners']) == 0:
            self.logger.debug('not enough colliders')
        else:
            #updating the ['radex']['title1'] and ['radex']['title2']
            #title1
            strng = 'radex LVG data'
            self.gui['radex']['title1'].set_text(strng)
            #title2
            strng = 'gasT\n%f\nN(specie)\n%e\n' % (gasTRadex, colDensThisSpec)
            for i, specStr in enumerate(self.radexObj.inFile['collisionPartners']):
                strng += '%s\n%e\n' % (specStr, self.radexObj.inFile['nDensCollisionPartners'][i])                
            self.gui['radex']['title2'].set_text(strng)

            # writing the parameters to a file (for debugging purposes)
            # this file can be used to re-run radex as standalone
            fName = self.parms['radex']['path'] + '-debug.inp'
            fObj = open(fName, 'w')
            fObj.write(self.radexObj.genInputFileContentAsStr() )
            self.logger.debug('input radex file written to %s' % fName)
            
            #----------------------------------running radex---------------------------------
            radexOutputMakesSense = False
            
            #parameters used for attempting to get radex to converge and vlaues make sense
            inFileOrig = self.radexObj.inFile.copy() # original parameters
            popDensTol = self.parms['radex']['popDensSumTol']  # ideally sum pop dens lower should be 1, for values less than tol
                                                              # the parameters are changed by 'changeFac' and another attempt is done
                                                              # until the sum of pop dens is close to 1 up to tol. 
            changeFac  = self.parms['radex']['changeFracTrial']  # the relative amount by which the input parameters to radex are changed
            nMaxTrial  = self.parms['radex']['nMaxTrial']     # number of time random changes are done before giving up

            nTried     = 0
            while not radexOutputMakesSense:
                self.radexObj.setDefaultStatus()
                self.radexObj.run( checkInput = True, verbose = True )

                if self.parms['radex']['checkOutputIntegrity']:
                    # checking the integrity of the solution in case it converged
                    if self.radexObj.getStatus() & self.radexObj.FLAGS['SUCCESS']:
                                            
                        totalPopDensLower = np.sum(self.radexObj.transitions[:]['pop_down'])
                        if np.fabs(self.parms['radex']['popDensSumExpected'] - totalPopDensLower) > popDensTol:
                            
                            self.logger.warn('====> pop dens does NOT make sense, trial %d\n====> total pop dense lower = %20.10f' % (nTried, totalPopDensLower) )
    
                            #change in input parameters by a tiny amount 0.1% and run again
                            #(this is a numerical bug in radex which fails sometimes)
                            self.radexObj.inFile['tKin']     = inFileOrig['tKin']*(1.0 + np.random.rand()*changeFac)
                            self.radexObj.inFile['molnDens'] = inFileOrig['molnDens']*(1.0 + np.random.rand()*changeFac)
                            for i, dens in enumerate(self.radexObj.inFile['nDensCollisionPartners']):
                                self.radexObj.inFile['nDensCollisionPartners'][i] = inFileOrig['nDensCollisionPartners'][i]*(1.0 + np.random.rand()*changeFac)
                        else:
                            radexOutputMakesSense = True
                            print '====>', 'pop dense make sense, trial %d' % nTried
                    else: # radex did not succeed
                        self.logger.debug('radex failed')
                        self.logger.debug('---------------------------------------------------------')
                        self.logger.debug('------------------------radex raw output-----------------')
                        self.logger.debug('---------------------------------------------------------')
                        print self.radexObj.rawOutput
                        self.logger.debug('---------------------------------------------------------')
                        self.radexObj.clearCurves()
                        break  # here radexOutputMakesSense woulbe be False
                    
                    nTried += 1
                    
                    if nTried >= nMaxTrial: #no meaningful output after nTried trials
                        break # here radexOutputMakesSense woulbe be False
                        
                else: #no need to check for validity of the radex pop densities
                    radexOutputMakesSense = True
                    break
            #------------------------------done running radex----------------------------------------
        
            # plotting the data
            # ----------------
            if self.radexObj.getStatus() & self.radexObj.FLAGS['SUCCESS']:
                
                self.radexObj.plotModelInFigureColumn(allTrans = np.arange(self.parms['radex']['maxDisplayTranistion']),
                                                      inAxes = self.pltRadex, 
                                                      title='')
                self.radexObj.setLabels()
            else:
                self.logger.debug('--------------warnings---------------')
                for warning in self.radexObj.warnings:
                    print warning
                self.logger.debug('--------------warnings---------------')
                    

        
    def clear(self):
        """clears all the bufferes allocated in the instance"""
        
        del self.infoAll
        del self.nDbFiles
        del self.chemNet
        del self.pltGmSec
        del self.radexParms
        del self.grds  
        del self.pltRadex
        pyl.close()

    def plotCurvesFromMeshes( self, inds = None, qx = None, qy = None, qz = None, *args, **kwargs ):
        """Plots curves from the meshes in the database on top of each other in the
          same window. This method inherits all keword argument from pyl.plot()
          
          :param list inds: an optional list of the indicies of the meshes in the database. If
           this is not passed, all the meshes in the database are considered.
          :param list qx: a list of strings pointing to the quantinity in the mesh dtype to
            be used as the absissa of the points. 
          :param list qy: same as qx but will be used as the ordinates.  
          :param list qz: a thrid quantity used to identity each curve (used in the legend)
          :warning: check and test this methods and document it better.
        """
        
        fig = pyl.figure()  # make the figure
        ax  = fig.add_axes([0.1,0.1, 0.6, 0.6])    # make new axes
        ax.set_xlabel(qx[-1])
        ax.set_ylabel(qy[-1])
        
        plts  = ()
        names = ()
        
        for m in self.meshes:
            x = fetchNestedDtypeValue(m, qx)
            y = fetchNestedDtypeValue(m, qy)
            plt, = ax.plot(x, y)
            
            if qz != None:
                z = fetchNestedDtypeValue(m, qz)
                plts  = plts  + (plt,)
                names = names + (z,)
        
        #print plts
        fig.legend(plts, names)
        pyl.show()
    
    def plot_3D_grid_points(self, **kwargs):
        """plots in 3D the parameters of the meshes in the database. By default
           the x,y,z coordinates are the log10 of nGas, G0, and gammaMech
           
           :param matplotlib.pyplot.figure figure: if this keyword is passed, a new figure is not 
               created and the axes are created in the figure object passed.
           :param  mpl_toolkits.mplot3d.Axes3D axes: if this keyword is passed, things are plotted
               in these axes, otherwise new ones are created in the figure.
           :return: (figure, axes, plt) 
        """

        if self.grid_x == None or self.grid_y == None or self.grid_z == None:
            raise TypeError('data of one or more coordinates not set!!.')
         
        if 'figure' not in kwargs:
            fig = pyl.figure()
        else:
            fig = kwargs['figure'] 
            
        # setting the x,y,z coordinate variables
        if self.grid_x == None or self.grid_x == None or self.grid_x == None:
            self.set_grid_axes_quantity_values()            
        
        if 'axes' not in kwargs:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = kwargs['axes']
        
        zPlt = self.grid_z
        
        plt = ax.plot(self.grid_x, self.grid_y, zPlt, 'o')
        ax.set_xlim( kwargs['ranges'][0][0], kwargs['ranges'][0][1] )
        ax.set_ylim( kwargs['ranges'][1][0], kwargs['ranges'][1][1] )
        ax.set_zlim( kwargs['ranges'][2][0], kwargs['ranges'][2][1] )
        ax.set_xlabel(self.grid_qx[-1])
        ax.set_ylabel(self.grid_qy[-1])
        ax.set_zlabel(self.grid_qz[-1])

        return (fig, ax, plt)
        
    ###############################setter and getters##############################
    def setChemicalNetwork(self, chemNet):
        self.chemNet = chemNet 
    def set_grid_qx(self, qx):
        self.grid_qx = qx
    def set_grid_qy(self, qy):
        self.grid_qy = qy
    def set_grid_qz(self, qz):
        self.grid_qz = qz
    def set_metallicity(self, metallicity):
        self.metallicity = metallicity

    def set_attributes(self, **kwargs):
        """set values of attributes from provided keywords
        """
        
        if 'grid_qx' in kwargs:
            self.set_grid_qx( kwargs['grid_qx'] )
        if 'grid_qy' in kwargs:
            self.set_grid_qy( kwargs['grid_qy'] )
        if 'grid_qz' in kwargs:
            self.set_grid_qz( kwargs['grid_qz'] )
        
    def set_default_attributes(self):
        """sets the default values of the attributes if they are not set
        """

        if self.grid_qx == None:
            self.set_grid_qx( ['hdr', 'nGas'] )
        if self.grid_qy == None:
            self.set_grid_qy( ['hdr', 'G0'] )
        if self.grid_qz == None:
            self.set_grid_qz( ['hdr', 'gammaMech'] )

    def set_grid_axes_quantity_values(self, *args, **kwargs):
        """assigns the values of self.grid_x,y,z from self.grid_qx,qy,qz. The values
         are extracted from all the meshes in the database.
        
        :param bool relativeGmech: if this is set, the z quantity data is set as the
          ratio of mechanicalHeating of the mesehs to the surface heating when the
          mechanical heating is zero. If this is set, the keyword 'referenceDatabase'
          sould be provided. It is assumed that the reference database has exact same
          parameters as the current one. It is assumed that the minimum gmech in 
          reference database is low enough to be assumed to be zero. An interpolation
          function is constructed which picks the lowest gMech of the reference
          database. This keyword cuases self.grid_qz to be set to 
          'gMech/gSurface(gMech=0)'

        :param string referenceDatabasePath: a string containing the path of the 
             reference database.
        
        :param float64 min_gMech: If this is present, the gMech of the reference models
          used would be this values instead of the minimum.
         
        :note: by default, the log of the quantities from the meshes are used. Only
          when the relative keyword is present, the value of the heating ratios is 
          set, and NOT the log of the heating ratios. 
        """
        
        if 'relativeGmech' in kwargs and 'grid_qz' in kwargs:
            raise ValueError('Either relativeGmech or grid_qz can be passed, not both') 
        
        #making sure self.qx,y,z are set
        self.set_attributes(**kwargs)
        self.set_default_attributes()

        if self.grid_x == None:
            self.grid_x = np.log10(self.getQuantityFromAllMeshes(self.grid_qx))
        if self.grid_y == None:
            self.grid_y = np.log10(self.getQuantityFromAllMeshes(self.grid_qy))
        
        if self.grid_z == None:
            
            # check whether we need to set the z quantity as the ration of heating rates
            try: 
                setRelativeGmech = kwargs['relativeGmech']
            except:
                setRelativeGmech = False                
                        
            # if relativeGmech is present AND it is True
            if setRelativeGmech == True:
            
                #setting the path of the reference database
                if 'referenceDatabasePath' in kwargs:
                    referenceDbDirPath = kwargs['referenceDatabasePath']
                elif 'runDirPath2' in self.parms_used_in_PDR_models:
                    referenceDbDirPath =  self.parms_used_in_PDR_models['runDirPath2']
                else:
                    raise ValueError("""missing the path of the reference database. Either "referenceDatabasePath"
                                        should be passed as a keyword or self.parms_used_in_PDR_models['runDirPath2']
                                        should be set.""")

                #computing (for the mesh points in the current database) the ratio of 
                #the mechanical heating to the surface heating when gmech = 0
                
                self.set_grid_qz(['','gMech/gSurface(gMech=0)'])
                
                # reading the reference archive
                self.logger.debug('setting up the reference archive')
                t0 = time()
                arxvRef = meshArxv(dirPath = referenceDbDirPath, readDb = True )
                self.logger.debug('time reading %f' % (time() - t0))
                
                #finding the minimum gMech of the reference database, and constructing
                #an interpolation function which returns the surface heating corresponding
                #to that mechanical heating rate.
                if 'min_gMech' in self.parms:
                    min_gMech = self.parms['min_gMech']
                else: 
                    min_gMech = np.min( np.array( 
                                                [mesh['hdr']['gammaMech'] for mesh in arxvRef.meshes] 
                                                )
                                       )
                
                gMechZero = self.grid_x.copy()
                #setting the values of gMech where the surface heating will be computed
                gMechZero[:] = np.log10(min_gMech)
                
                #getting the surface heating that the models in the current would have
                #from the reference database which has gmech = 0
                arxvRef.set_grid_qx( self.grid_qx )
                arxvRef.set_grid_qy( self.grid_qy )
                f = arxvRef.construct3DInterpolationFunction(quantity = ['therm', 'heating'], 
                                                             slabIdx  = 0, 
                                                             log10 = True,
                                                             grid_qz = ['hdr','gammaMech'])
                dataNew   = np.array( [self.grid_x, self.grid_y, gMechZero] ).T
                
                # the surface heating that the models in the current would have if the 
                # mechanical heating were zero
                gammaSurf = f(dataNew) 

                z = np.log10(self.getQuantityFromAllMeshes( ['hdr', 'gammaMech']) )
                self.grid_z = np.log10(10.0**z / 10.0**gammaSurf)
            else:
                # just use gMech as the 3rd axis
                self.grid_z = np.log10( self.getQuantityFromAllMeshes(self.grid_qz) )
    
    def setupLogger(self):
        """sets up the logger which will prepend info about the printed stuff. Assignes a value to self.logger."""
        # setting up the logger                                                                                                                                                                                                              
        # create logger                                                                                                                                                                                                                      
        logger = logging.getLogger('simple_example')
        if not len(logger.handlers):
            logger.setLevel(logging.DEBUG)

            # create console handler and set level to debug                                                                                                                                                                                      
            ch = logging.StreamHandler( sys.stdout )  # setting the stream to stdout                                                                                                                                                             
            ch.setLevel(logging.DEBUG)
    
            # create formatter                                                                                                                                                                                                                   
            formatter = logging.Formatter('[%(asctime)s %(funcName)s() %(filename)s:%(lineno)s] %(message)s') # this was the original in the example                                                                              
    
            # add formatter to ch                                                                                                                                                                                                                
            ch.setFormatter(formatter)
    
            # add ch to logger                                                                                                                                                                                                                   
            logger.addHandler(ch)

        self.logger = logger
        return self.logger
    
    def writeRadexDbAscii(self):
        """a utlity method which can be used to output stuff to an ascii file. 
        """
        
        fileName = '/home/mher/ism/marissa-0.csv'
        
        fObj = open(fileName, 'wb')
        for i in np.arange(self.nMeshes):
                
            if self.meshesRadex[i] != None:
                    
                x = self.grid_x[i]
                y = self.grid_y[i]
                z = self.grid_z[i]
                
                strng = '%f,%f,%f' % (x,y,z)
                for trans in self.meshesRadex[i]:
                    strng += ',%s-%s,%e' % (trans['upper'],trans['lower'],trans['fluxcgs'])
                strng += '\n'
                if z < -29.5:
                    print strng
                    fObj.write(strng)
                #v = self.meshesRadex[i][transitionIdx][quantity]
                

        fObj.close()
        #writing some of the output of the DB in to an ascii file
        
    def getUnique_grid_z_sections(self, relDiff):
        """returns the unique section in z up to a relative difference of relDiff"""
        
        z = self.grid_z
        
        sections = []

        while z.size != 0:
            mn = np.min(z)
            #print 'min = ', mn, 10.0**mn
            if mn == 0.0:
                inds = np.where( z != mn)
            else:
                inds = np.where( abs(1.0 - z/mn) > relDiff)
            z = z[inds]
            sections.append(mn)
    
        return np.array(sections, dtype=np.float64)
    
    def do_something_for_all_meshes(self):
        """A template method which loops over all meshes and gets/computes something 
        from each mesh.
        """
        
        fName = '/home/mher/ism/marissa/CO-ladder.out'
        fObj = open(fName, 'w')
        
        #a temporary object used to calculate stuff from a pdr mesh object
        m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
                        
        #--------------------------looping over all the meshes------------------
        for i, pdrData in enumerate(self.meshes):
            
            x, y, z = self.grid_x[i], self.grid_y[i], self.grid_z[i] # model parms 
            pdrMeshData = pdrData                                    # pdf mesh data
            rdxMeshData = self.meshesRadex[i]                        # radex data
            
            m.setData( pdrMeshData )
            dxSlabs = m.compute_dx()

            if z == -10:
                print 'mesh x,y,z = %.2f, %.2f, %.2f' % (x, y, z),
                
                """
                #----------------dumping the total cooling--------------------
                quantity = ['therm','heating']
                v = (1.0/(2.0*np.pi))*m.compute_integrated_quantity(quantity)
                fObj.write('%d %.2f %.2f %.5e\n' % (i, x, y, v))
                #----------------done dumping the total cooling---------------
                """
                
                #----------------dumping the radex CO lines--------------------
                fObj.write('%d %.2f %.2f ' % (i, x, y))

                if rdxMeshData != None:
                    for trans in self.meshesRadex[i]:
                        fObj.write('%.5e %.5e ' % (trans['pop_up'], trans['fluxcgs']))
                else:
                    fObj.write('nan')
                fObj.write('\n')
                #----------------done dumping the total cooling---------------
                
                
                
        #--------------------------done looping over the meshes---------------       
                
        fObj.close()
        print 'wrote the file %s'
                
    def write_pdr_code_guess_database(self, path):
        """writes the guess database that is used by the pdr code. See 
        the struct 'database' in oneSided/src/vars.h.
        
        Steps involved in preparing the data before calling this method:
        
          - collect a bunch of mesh files in simulation folder
          - construct the database by running : constructReadArchive.py
          - then call this method and save the file to be used later 
          
        """
        from amuse.community.pdr.pyUtils.guess_db import guess_db
        
        pdrGuessDb = guess_db( path = path )

        # setting the x,y,z coordinate variables
        if self.grid_x == None or self.grid_x == None or self.grid_x == None:
            self.set_grid_axes_quantity_values()            

        #getting the data from all the meshes
        #------------------------------------
        logDens   = self.grid_x
        logG0     = self.grid_y
        logGammaM = self.grid_z
        Teq       = self.getQuantityFromAllMeshes(quantity = ['state','gasT'], slabIdx = 0)
        #collecting the abundance (taking the transpose to make it compatible with 
        #the shape acceptable by the 'guess' class
        abun      = [ m['state']['abun'][:,0] for m in self.meshes],
        abun      = np.array(abun, dtype = np.float64)[0].T
        #------------------------------------
        
        #print pdrGuessDb.get_path()
        pdrGuessDb.set_db_data(logDens, logG0, logGammaM, Teq, abun)
        pdrGuessDb.write()    
    
    def get_mesh_index(self, x = None, y = None, z = None):
        """Returns the index of the mesh in the x,y,z parmaerter coordinate lests self.grid_x, self.grid_y, self.grid_z
        which are closest to the input x,y,z."""
        
        l2Distance  = np.sqrt( (y - self.grid_y)**2 + (x - self.grid_x)**2 + (z - self.grid_z)**2 )
        indMin = l2Distance.argmin()
    
        return indMin
    
    def get_mesh_data(self, x = None, y = None, z = None):
        """Returns the data (of dtype :data:`mesh`.data) of the pdr mesh in the x,y,z parmaerter 
        coordinate lests self.grid_x, self.grid_y, self.grid_z which are closest to the input x,y,z."""
        
        ind = self.get_mesh_index(x, y, z)
        return self.meshes[ind] 
    
    def grid_interpolator(self, n_gas, fuv_heating, g_mech,
                          quantity = None, slabIdx = None):
        """Returns interpolated vlaues for input quantities. The input
        values are gas density, fuv heating and mechanical heating as
        numpy arrays with amuse units.  These are converted to the 
        same units as the PDR meshes.  An interpolation function is 
        constructed and the interpolation is done for a quantity pointing
        to the dtype of data of a mesh."""
        
        from amuse.units import units

        #----------------------------------------------------------------
        # converting the quantities to units compatible with PDR models
        # getting the quantities as numpy arrays
        #----------------------------------------------------------------
        #converting gas density to cgs
        meanmwt=1.3|units.amu
        n_gas_cgs = (n_gas/meanmwt).value_in( units.cm **-3 )
        
        #converting the G0 to correct units
        G0 = 6.54*fuv_heating.value_in( units.none )
        
        #converting the mechanical heating to the correct unit
        #from m2 / s3 (Inti Units) to  erg / (cm^3 s)
        g_mech = g_mech.as_quantity_in( units.erg / (units.g * units.s ) )
        # converting the mechanical heating rate from per unit mass (of H gas)     
        # to per unit volume (see notesISM.odt)
        g_mech = g_mech.value_in( g_mech.unit ) # gMech in erg / (g s )
        g_mech = g_mech * 1.6474e-24 * n_gas_cgs # gMech in erg / (cm^3 s) 

        #talking the log of the input quantities at which the interpolated
        #quantities will be returned
        lx = np.log10(n_gas_cgs)
        ly = np.log10(G0)
        lz = np.log10(g_mech)
        #-----------------------done converting units----------------------
        
        #-----------------------generate the interpolation function--------
        fInterp = self.construct3DInterpolationFunction(quantity, slabIdx)

        #picking the SPH particles which are within the parametrer space
        #of the PDR grids 
        lx = lx.clip(self.grid_x.min(), self.grid_x.max())
        ly = ly.clip(self.grid_y.min(), self.grid_y.max())
        lz = lz.clip(self.grid_z.min(), self.grid_z.max())
        
        dataNew = np.array( [lx, ly, lz] ).T
        qInterp = fInterp(dataNew)

        return qInterp