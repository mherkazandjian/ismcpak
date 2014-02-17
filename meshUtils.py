import os, glob, fnmatch, sys, pickle
import numpy as np
from numpy import array, vstack, save, load
import pylab as pyl
import logging
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d import Axes3D

from time  import *
from mesh  import *
from mylib.utils.ndmesh import ndmesh
from mylib.utils.misc   import fetchNestedDtypeValue, scale
from ismUtils   import *
from radex      import *
from scipy      import interpolate
import chemicalNetwork
import lineDict
from despotic import cloud
from multiprocessing import Process, Pipe
from mylib.utils.interpolation import interpolator_sectioned

class meshArxv(object):
    """ this class generates and manipulates archives of PDR meshes. By defaults self.set_default attributes and self.set_grid_axes_quantity_values are called upon intialization.
    
     :param bool readDb: if this is set, the database is read during initialization.
      
     :param bool no_init_from_run_parms: if this is passed as a keyword, nothing is initialized from used_parms.pkl.                

     .. code-block:: python
     
         ## basic usage
         arxvPDR = meshUtils.meshArxv(dirPath = '/home/mher/ism/runs/oneSidedDB-z-1.0')
         
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
        self.logger = None
        
        if 'dirPath' in kwargs:
            self.dirPath = kwargs['dirPath']
        else:
            self.dirPath = None
            """path to the directory which contains all the database info and data"""
            
        self.nMeshes    = None
        """ np.int32 : number of meshes in the database (this is now as 1x1 array, convert it just into an np.int32 scalar"""
        
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
        """A dict which holds the name of the species self.currenRadexDb['specStr'] whose radex 
        db is in the attributes self.meshesRadex and self.infoAllRadexspecies name of the current 
        radex database being used and the Av depth upto which the data was computed
        self.currenRadexDb['Av']""" 

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


           .. note:: the transitions are stored even if there are warnings when running radex. So 
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
        dictionary is a dictionary holding the Dbs of the species for different Avs.
        For example :
    
        .. code-block:: python
        
            Dbs_Av10 = self.radexDbs['10.00']
               
        in this case Dbs_Av10 is a dict holding the databases for the species at Av=10

        .. warning:: the resolution in storing Av is %.2f.
        
        The keys of Dbs_Av10 are species string. For example:

        .. code-block:: python
        
            Db_CO_Av10 = self.radexDbs['10.00']['CO']
        
        Db_CO_sAv10 is a dict of two keys : 'infoAll' and 'meshes'. As a full example 
        for accessing the Dbs for CO and HCN for Av=10 and Av=30 could be set manually 
        to self.meshesRadex and self.infoAllRadex via:
        
        .. code-block:: python
        
            #setting the CO radex database for Av=10
            self.meshesRadex = self.radexDbs['10.00']['CO']['meshes']
            self.infoAllRadex = self.radexDbs['10.00']['CO']['infoAll']
            #this is equiveletly done using
            self.use_radexDb(Av=10, specStr=CO)
            
            #setting the CO radex database for Av=30
            self.meshesRadex = self.radexDbs['30.00']['CO']['meshes']
            self.infoAllRadex = self.radexDbs['30.00']['CO']['infoAll']
            #this is equiveletly done using
            self.use_radexDb(Av=30, specStr=CO)
    
            #similarly for other species and Avs
        
        This could be achived also using the method 
        """
        self.radexObj = None
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
        self.grid_x_unique = None
        """holds the unique values of grid_x"""
        self.grid_y_unique = None
        """holds the unique values of grid_y"""
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
            self.set_metallicity(kwargs['metallicity'])
        else:
            self.metallicity = None
        
        self.gui = None
        self.parms_used_in_PDR_models = None
        """a dictionary holding the parameters used in modeling the PDRs.  This is read
        from self.dirPath/used_parms.pkl""" 

        #sections in gm within which saperate interpolation functions will be built.
        # .. note:: make sure there are enough points in the database within each section
        self.intervals_z = [ 
                         [-50.0, -35.0],  
                         [-35.0, -34.0], [-34.0, -33.0], [-33.0, -32.0], [-32.0, -31.0], [-31.0, -30.0],    
                         [-30.0, -29.0], [-29.0, -28.0], [-28.0, -27.0], [-27.0, -26.0], [-26.0, -25.0],    
                         [-25.0, -24.0], [-24.0, -23.0], [-23.0, -22.0], [-22.0, -21.0], [-21.0, -20.0],    
                         #[-20.0, -13.0],
                       ]
        self.ghost_z = 1e-6 #1.1
    
        #sections in Av within which saperate interpolation functions will be built.
        # .. note:: make sure there are enough points in the database within each section
        self.intervals_t = [  [0.01, 1.0],  [1.0, 2.0],   [2.0, 3.0],
                         [3.0, 4.0],   [4.0, 5.0],   [5.0, 6.0],   [6.0, 7.0],   [7.0, 8.0],   [8.0, 9.0], [9.0, 10.0],
                         [10.0, 12.0], [12.0, 14.0], [14.0, 16.0], [16.0, 18.0], [18.0, 20.0], 
                         [20.0, 22.0], [22.0, 24.0], [24.0, 26.0], [26.0, 28.0], #, [28.0, 30.0],
                       ] 
    
        self.ghost_t = 1e-6 #1.1
        
        if 'no_init' not in kwargs:

            #The logging object which will be used to output stuff to stdout
            self.logger = self.setupLogger()
                        
            #reading the used_parms.pkl file (if found) and setting the vaues read from it
            if 'no_init_from_run_parms' not in kwargs:
                self.read_used_parms_used_in_PDR_models()
    
            #setting up the auxiliary attributes if the parms are available
            #self.mshTmp and self.radexObj
            if 'radex' in self.parms and  self.parms['radex']['use'] == True:
                if self.parms['radex']['use']:
                    #making the instance            
                    self.setup_default_radex_instance(self.parms['radex'])
    
            self.mshTmp = None
            """a mesh object which is used to store single mesh data just for plotting purposes"""
            if self.chemNet != None and self.metallicity != None:
                self.mshTmp = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
    
            if 'readDb' in kwargs and kwargs['readDb'] == True:
                kwargs.pop('readDb')
                self.readDb(check=True)
    
            #reading|computing radex emission databases
            #-------------------------------------------
            if 'radex' in self.parms:
                if self.parms['radex']['use']:
                    if 'gridsInfo' in self.parms and 'show' in self.parms['gridsInfo']['11'] and self.parms['gridsInfo']['11']['show'] == True:
                        #----
                        if self.parms['gridsInfo']['11']['type'] == 'radex':
                            if self.parms['radex']['compute']:
                                self.constructRadexDatabase(writeDb = self.parms['radex']['writeDb'])
                            else:
                                self.readDbsRadex(Av = self.parms['gridsInfo']['11']['Av_max'], 
                                                  species = self.parms['gridsInfo']['11']['specStr']
                                                  )
                        #----
                        else:
                            self.logger.debug('pdr emissions will be computed on the spot')
                    #----
                    if self.parms['radex']['loadAllDbs']:
                        #reading all the databases
                        self.readDbsRadex(allDbs = True)
                        #set the current Db to the one specieified in the input parms
                        if 'gridsInfo' in self.parms and self.parms['gridsInfo']['11']['type'] == 'radex': 
                            self.use_radexDb(self.parms['gridsInfo']['11']['Av_max'], 
                                             self.parms['gridsInfo']['11']['specStr'])
                #----
            #-------------------------------------------
    
            #set some default attributes (the quantities of the grids) if meshes have been read
            if self.meshes != None:
                
                if 'grid_qx' in self.parms and 'grid_qy' in self.parms and 'grid_qz' in self.parms:
                    #specifying the grid quantities explicityly
                    self.set_attributes(**kwargs)
                    self.set_grid_axes_quantity_values()
                else:
                    if 'relativeGmech' in self.parms:
                        # setting the x,y,z quantities to be used for ploting
                        self.set_grid_axes_quantity_values(relativeGmech = self.parms['relativeGmech'])
                    else:
                        self.set_default_attributes()
                        self.set_grid_axes_quantity_values()
            else:
                self.logger.debug('no meshes read, not setting defualy attributes for the grids')
        #
        #finished intializing and setting up some attributes
        
    def read_used_parms_used_in_PDR_models(self):
        fpath_parms_used = os.path.join(self.dirPath, 'used_parms.pkl')
        if os.path.isfile(fpath_parms_used):
            self.logger.debug('setting metallicity and chemical network from : %s' % fpath_parms_used)
            #loading the parameters
            parms = pickle.load(open(fpath_parms_used))
            #setting the attributes
            self.metallicity = parms['metallicity']
            self.setupChemistry(parms['chemistry']) 
            self.parms_used_in_PDR_models = parms
            
            if 'min_gMech' not in self.parms and 'min_gMech' in parms:
                self.parms['min_gMech'] = parms['min_gMech']
                
            if 'relativeGmech' not in self.parms and 'relativeGmech' in parms:
                self.parms['relativeGmech'] = parms['relativeGmech']
                                 
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
     
            #setting the variables which make the run unique
            #-----------------------------------------------
            infoAll[i]['parms'][0] = mData['hdr']['G0']
            infoAll[i]['parms'][1] = mData['hdr']['nGas']
            infoAll[i]['parms'][2] = mData['hdr']['gammaMech']
            infoAll[i]['parms'][3] = np.float64(0.0)
            infoAll[i]['parms'][4] = np.float64(0.0) #not used
            #------------------------------------------------
            
            
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

    def make_copy(self, dirName = None, x = None, y = None, z = None):
        """makes a copy of all the databases with some selections specified by the arrays x,y,z 
           picking a meshes whose coordinates are 'exactly' in x,y,z. This is done for all the Avs
           available for the Radex meshes.
           
           example : 
           
               arxv = meshUtils.meshArxv(**parms)
               arxv.make_copy(dirName = '/home/mher/foo/',
                              x = arxv.grid_x_unique[::2], 
                              y = arxv.grid_y_unique[::2], 
                              z = arxv.grid_z_unique[::2],
                             ) 
        """
        
        self.writeDb_selection(dirName, x, y, z)

        for key_Av in self.radexDbs:
            
            for key_specStr in self.radexDbs[key_Av]:
                
                self.use_radexDb(numpy.float64(key_Av), key_specStr)
                
                self.writeDbRadex_selection(dirName, x, y, z)

        
    def writeDb_selection(self, dirName = None, x = None, y = None, z = None):
        """ writes the mesh.db and mesh.db.ino into the dir dirName.  Only meshes whose
            self.grid_x,y,z[i] are in x,y,z are written.
        
        :param string dirName: The directory where to write the database files.
        """
        
        if dirName == None:
            raise ValueError('keyword dirName not set...please provide it')
        if dirName.replace('/','') == self.parms['dirPath'].replace('/',''):
            raise ValueError("keyword dirName should be different from self.params['dirPath']")
        
        # writing the mesh to the database file
        dbDataFObj = file(dirName + 'meshes.db', 'wb')
        
        #a temporary copy of self.infoAll which will store the info that will be written
        #to meshes.db.info in the new dirName 
        infoAll_new = self.infoAll.copy()
        
        counter = 0 #to keep track of the number of meshes with the specified coords
        
        #looping over all the coordinates of the grids...
        for i in np.arange( self.nMeshes ):
            
            #writing only meshes whose coordinates are in the specified ones
            if self.grid_x[i] in x and self.grid_y[i] in y and self.grid_z[i] in z:
                
                self.meshes[i].tofile( dbDataFObj )
                
                infoAll_new[counter] = self.infoAll[i]
                
                infoAll_new[counter]['info'][0] = counter # mesh number 

                infoAll_new[counter]['info'][2] = dbDataFObj.tell()  # offset from the start of file
                
                np.array( infoAll_new[counter]['info'][2] ).tofile(dbDataFObj)
                
                counter += 1
                
        dbDataFObj.close()
        
        nMeshes_new = numpy.int32(counter)
        
        # writing the db info into a file
        dbInfoFObj = file(dirName + 'meshes.db.info', 'wb')
        self.ver.tofile( dbInfoFObj)
        nMeshes_new.tofile( dbInfoFObj )
        infoAll_new[0:counter].tofile( dbInfoFObj )        
        dbInfoFObj.close()
        self.logger.debug('wrote successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))        

    def writeDbRadex_selection(self, dirName = None, x = None, y = None, z = None):
        """ writes the meshesRadex.db and meshRadex.db.info into the dir dirName.  Only meshes whose
            self.grid_x,y,z[i] are in x,y,z are written. 
        
        .. todo:: add an example of a safe way to use this.
        
        :param string dirName: The directory where to write the database files.
        """
        
        if dirName == None:
            raise ValueError('keyword dirName not set...please provide it')
        if dirName.replace('/','') == self.parms['dirPath'].replace('/',''):
            raise ValueError("keyword dirName should be different from self.params['dirPath']")

        Av_end, specStr = self.currentRadexDb['Av'], self.currentRadexDb['specStr']
        
        dbFname, dbInfoFname = self.getRadexDbPath(Av_end, specStr, dirpath = dirName) 

        #creating the Av_XX.XX dir is neccessary
        dirPath_this_db = dbFname.replace(dbFname.split('/')[-1],'')
        if os.path.exists(dirPath_this_db) == False:
            os.mkdir(dirPath_this_db)
            self.logger.debug('created directory : %s' % dirPath_this_db)
        
        # writing the mesh to the database file
        dbDataFObj = file(dbFname, 'wb')
        
        #a temporary copy of self.infoAll which will store the info that will be written
        #to meshes.db.info in the new dirName 
        infoAllRadex_new = self.infoAllRadex.copy()
        
        counter = 0 #to keep track of the number of meshes with the specified coords
        
        #looping over all the coordinates of the grids...
        for i in np.arange( self.nMeshes ):
            
            #writing only meshes whose coordinates are in the specified ones
            if self.grid_x[i] in x and self.grid_y[i] in y and self.grid_z[i] in z:

                #writing the mesh parameters before writing the transition data (extra data integrirty check)
                thisMeshParms = np.array([self.meshes[i]['hdr']['G0'],
                                          self.meshes[i]['hdr']['nGas'],
                                          self.meshes[i]['hdr']['gammaMech']], dtype = np.float64)
                thisMeshParms.tofile( dbDataFObj )

                infoAllRadex_new[counter] = self.infoAllRadex[i]
                infoAllRadex_new[counter]['info'][0] = counter # mesh number
                
                #writing the transition info to the transiotions db file
                if self.meshesRadex[i] != None:
                    #writing the transition data
                    self.meshesRadex[i].tofile( dbDataFObj )
                    #filling the necessary values in infoAllRadex
                    infoAllRadex_new[counter]['info'][2] = dbDataFObj.tell()  # offset from the start of file
                else:
                    infoAllRadex_new[counter]['info'][2] = -1                
                    
                #writing the offset so far from the begining of the .db file
                np.array( infoAllRadex_new[counter]['info'][2] ).tofile(dbDataFObj)
                
                counter += 1
                
        dbDataFObj.close()
        
        nMeshesRadex_new = numpy.int32(counter)
        
        # writing the db info into a file
        dbInfoFObj = file(dbInfoFname, 'wb')
        self.ver.tofile( dbInfoFObj)
        nMeshesRadex_new.tofile( dbInfoFObj )
        infoAllRadex_new[0:counter].tofile( dbInfoFObj )
        dbInfoFObj.close()
        self.logger.debug('wrote successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))
        
    def readDb(self, check = None):
        """ reads the database and assigns the appropritate attributes (document)
            
            :param check: if this is set (to any value) the self.checkIntegrity() is called.
        """ 
        
        dbInfoFObj = file(self.dirPath + 'meshes.db.info', 'rb')
        dbDataFObj = file(self.dirPath + 'meshes.db'     , 'rb')

        self.ver = np.fromfile( dbInfoFObj, dtype = (np.int32, 3), count = 1)[0]        
        self.nMeshes = np.fromfile( dbInfoFObj, dtype = np.int32, count = 1)[0]

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
    
    def getRadexDbPath(self, Av, specStr, dirpath = None):
        """Give the Av and the species string, returns the paths string pointing to the file (see self.radexDbs).
        The path to the .db and to the .info.db files
        """

        if dirpath == None:
            dirName = self.parms['dirPath']
        else:
            dirName = dirpath
            
        radexDbsPath_thisAv = dirName + 'radexDbs/Av_%.2f/' % Av
        
        dbFname = radexDbsPath_thisAv + 'meshesRadex.db.%s' % specStr
        dbInfoFname = radexDbsPath_thisAv + 'meshesRadex.db.info.%s' % specStr
        
        return (dbFname, dbInfoFname)

 
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

        Av_end, specStr = self.currentRadexDb['Av'], self.currentRadexDb['specStr']
        
        dbFname, dbInfoFname = self.getRadexDbPath(Av_end, specStr) 

        #creating the Av_XX.XX dir is neccessary
        dirPath_this_db = dbFname.replace(dbFname.split('/')[-1],'')
        if os.path.exists(dirPath_this_db) == False:
            os.mkdir(dirPath_this_db)
            self.logger.debug('created directory : %s' % dirPath_this_db)
             
        #---------------------writing the mesh to the database file------------------
        dbDataFObj = file(dbFname, 'wb')
        
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
        #----------------------done writing the meshes data---------------------------------
        
        #-------------------------writing the db info into a file---------------------------
        dbInfoFObj = file(dbInfoFname, 'wb')
        self.verRadex.tofile( dbInfoFObj)
        numpy.array(self.nMeshes, dtype=numpy.int32).tofile( dbInfoFObj )
        self.infoAllRadex.tofile( dbInfoFObj )
        dbInfoFObj.close()
        self.logger.debug('wrote successfully the radex database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name))        
        
    def radex_db_has_been_read(self, Av, specStr):
        """defining a function which checks if the specified Db has been
           already read. If found returns True, Flase otherwise
        """

        if '%.2f' % Av in self.radexDbs.keys():
            if specStr in self.radexDbs['%.2f' % Av].keys():
                return True
            else:
                return False #specie with this Av is not present
        else:
            return False  #no db with this Av is present

    def readDbRadex(self, Av, specStr, check = None):
        """Reads the database suffixed by specStr (i.e meshesRadex.db.(specStr)and assigns the
           appropritate attributes (document) and assigns the read data to self.meshesRadex and
           self.infoAllRadex. 
            
           :param float Av: this is the Av for the corresponding database at which the emissions were computed 
           :param bool check: if this is set (to any value) the self.checkIntegrity() is called.
           :param string specStr: the string of the specie whose database is to be read.
           
           .. warning:: the keyword specStr is not functional yet.
           .. note:: before calling this mehtod, an instance of the radex class should be     
               assigned to self.radexObj. 
        """ 
                
        #reading the database only if it is not already read. If it is read, it would
        #have a key in radexDbs
        if self.radex_db_has_been_read(Av, specStr) == False:
            
            #paths of the files of this db
            dbFname, dbInfoFname = self.getRadexDbPath(Av, specStr)
            
            #file object of this db
            dbDataFObj, dbInfoFObj = file(dbFname, 'rb'), file(dbInfoFname, 'rb')
    
            #reading the db info
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
            radexTransitionsDtype = radex(None,None).generateTransitionDtype()
    
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
                            
            #setting state attributes related to the current radexDb
            self.append_radex_db_related_attributes(meshesRadex, infoAllRadex, Av, specStr)

            if check:
                self.checkIntegrityRadex()

        else:
            self.logger.debug('radex database for %s and Av = %f already read...skipping' % (specStr, Av))

    def readDbsRadex(self, Av = None, species = None, allDbs = None, in_Av_rng = None):
        """Loads radex database files from disk to a dictionary. The last database loaded is
        the one being used.
        
        The structure of the dir where the Dbs are present should be as follows:
        
            - dirPath/radexDbs/AvXX.XX/meshRadex.db.CO
            - dirPath/radexDbs/AvXX.XX/meshRadex.db.infoCO
            
              and so on for each specie

            - dirPath/radexDbs/AvYY.YY/meshRadex.db.CO
            - dirPath/radexDbs/AvYY.YY/meshRadex.db.infoCO

        where XX.XX and YY.YY are the Av's for which the Dbs have been computed. For example :
         
            - dirPath/radexDbs/Av10.00/meshRadex.db.CO
            - dirPath/radexDbs/Av10.0/meshRadex.db.infoCO
        
        
        :param species: A string or a list of strings holding the names of the species files to be loaded. 
        :param Av: A float or a list floats holding the Avs to be restored (if they are present). It
         the floats are formatted into %.2f when they are included in the database path for each path. For
         example if Av = 12.5433 and the species = 'CO', this is formatted to :
           
             dirPath/radexDbs/Av12.54/meshRadex.db.CO
        :param bool all: if this is set, the directory radexDbs is read recursively and
         all the radex Dbs in it are loaded into self.radexDbs. If this is passed, the
         parameters Av and species are ignored 
        """
  
        if allDbs == None:
            #loading the specified radex Dbs
            if species == None and Av == None and  in_Av_rng == None:
                raise ValueError('one of the keywords "species", "Av" or "in_Av_rng" must be provided.')
            
            if hasattr(species, '__iter__') == False:
                species = [species]
            
            if Av != None and in_Av_rng != None:
                raise ValueError('either Av or in_Av_rng must be provided, not both')
            elif Av != None and hasattr(Av, '__iter__') == False:
                Av = [Av]
            elif in_Av_rng != None:
                Avs_on_disk = self.query_available_radex_dbs()
                inds_in_rng = numpy.where((Avs_on_disk >= in_Av_rng[0])*(Avs_on_disk <= in_Av_rng[1]))[0]
                
                if inds_in_rng.size != 0:
                    Av = Avs_on_disk[inds_in_rng]
                else:
                    raise ValueError('no radex Dbs were found in the Av range [%f, %f]' % (in_Av_rng[0],in_Av_rng[1]))
            
            for AvThis in Av:   
                for specStr in species:
                    self.readDbRadex(AvThis, specStr, check = True)
        else:
            #loading all the Dbs in radexDs/

            dirPath = self.parms['dirPath'] + 'radexDbs'
            pAv = subprocess.Popen(['ls', dirPath], stdout=subprocess.PIPE)
            Avlsout = pAv.stdout.read()
            
            print 'directores which might contain Dbs : ', Avlsout
            
            #getting the Avs in the radexDbs dir
            Avs = []
            for Avline in Avlsout.split('\n'):
                Avline = Avline.replace('Av_','').strip()
                if Avline != '':
                    Av = float(Avline)
                    Avs.append(Av)
                    
                    print 'looking for Dbs in the dir whose Av is %.2f' % Av
                    
                    #getting the species names
                    subDirName = dirPath + '/Av_%.2f' % Av
                    print subDirName
                    pSpcs = subprocess.Popen(['ls', subDirName], stdout=subprocess.PIPE)
                    speclsout = pSpcs.stdout.read()
                    
                    print 'Db files found ', speclsout
                    
                    for specline in speclsout.split('\n'):
                        if 'info' in specline:
                            continue
                        else:
                            if specline != '':
                                specStr = specline.replace('meshesRadex.db.', '').strip()
                                
                                self.readDbRadex(Av, specStr, check = True)#read the Db

    def query_available_radex_dbs(self):
        """method which looks for the available radex dbs, and for each Av returns a dict of the 
        available dbs in the dirpath"""
    
        sourceDir = os.path.join(self.dirPath, 'radexDbs')
        lookFor   = 'meshesRadex.db.*'

        matches = []
        for root, dirnames, filenames in os.walk(sourceDir):
            for filename in fnmatch.filter(filenames, lookFor):
                matches.append(os.path.join(root, filename))
        
        #getting the Avs of the matched databases from the dirs in which the matches were found
        if len(matches) != 0:
            
            Avs = []
            for match in matches:
                match = match.split('/')[-2]  #the name of the dir
                AvStr = match.replace('Av_','')
                Avs.append(AvStr)
            
            Avs = numpy.unique(numpy.float64(Avs))
            
            return Avs
        
        else:
            return None
        
                                        
    def use_radexDb(self, Av=None, specStr=None, silent=None, load_if_not_in_memory=None):
        """A method which sets a radex database from the databases available in self.radexDbs."""

        if self.radex_db_has_been_read(Av, specStr) == False:
            self.readDbRadex(Av, specStr)
        
        if self.radex_db_has_been_read(Av, specStr) == True:

            if self.currentRadexDb['Av'] != Av or self.currentRadexDb['specStr'] != specStr:
                #setting the new Db only if it is different from the specified one 
                self.meshesRadex    = self.radexDbs['%.2f' % Av][specStr]['meshes']
                self.infoAllRadex   = self.radexDbs['%.2f' % Av][specStr]['infoAll']
                self.currentRadexDb = {'Av':Av, 'specStr':specStr}
                if silent != None and silent == False:
                    self.logger.debug('swithced to %s radex database at Av = %f' % (specStr, Av))
        else:
            raise ValueError('radex Db for Av %.2f for the specie %s has not been read into memory' % (Av, specStr))
        
    def mergeDbs(self, newDbRunDirPath=None, new_arxv=None, outDirPath = None):
        """ merges two databases into one and write the resulting db and its info file
            :param string newDbRunDirPath: the dir in which the new db to be added to the current db is located.
            :param object new_arxv: The meshArxv object to be added to self
        """
        
        if newDbRunDirPath!=None and new_arxv==None:
            arxv2 = meshArxv()
            arxv2.readDb( newDbRunDirPath )
            arxv2.checkIntegrity()
        else:
            if new_arxv != None:
                arxv2=new_arxv 
            
        if self.ver[0] != arxv2.ver[0] and self.ver[1] != arxv2.ver[1] and self.ver[2] != arxv2.ver[2]:
            raise ValueError('can not merge databases with different versions.')
        
        self.nMeshes += arxv2.nMeshes
        
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
            
            diff_this = np.sqrt( np.dot(xdv,xdv) )
             
            if diff_this != 0.0:
                raise NameError(
                                """"archive integrity test failed. database file may be corrupt.
                                    parms in infoAll[i]     : %e %e %e
                                    parms in mesh[i]['hdr'] : %e %e %e
                                """ % (
                                       self.infoAll[i]['parms'][0], self.infoAll[i]['parms'][1], self.infoAll[i]['parms'][2],
                                       self.meshes[i]['hdr']['G0'], self.meshes[i]['hdr']['nGas'], self.meshes[i]['hdr']['gammaMech']
                                      )
                                ) 

            diff += diff_this 
            
        if diff == 0.0:
            self.logger.debug('archive integrity test passed')

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

            diff_this = np.sqrt( np.dot(xdv,xdv) )
             
            if diff_this != 0.0:
                raise NameError("""archive integrity test failed. database file may be corrupt.
                                   parms in infoAllRadex[i] : %e %e %e
                                   parms in infoAll[i]      : %e %e %e
                                """ % (
                                      self.infoAllRadex[i]['parms'][0], self.infoAllRadex[i]['parms'][1], self.infoAllRadex[i]['parms'][2],
                                      self.infoAll[i]['parms'][0], self.infoAll[i]['parms'][1], self.infoAll[i]['parms'][2]
                                      )
                                )
                                
            diff += diff_this 
            
        if diff == 0.0:
            self.logger.debug('archive integrity test passed')

    
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

        print '------------------------------------------------------------------------------'
        ###----------------------------setup the network as in the pdr code---------------
        
        #scratching out reactions  (check if 3707 is removed in the original PDR code)
        net.pop_reactions( [799, 3706, 3707, 4076, 4412, 4413] )
        net.pop_reactions( net.filter_reactions(withType='TR') )
        
        #changing temperature ranges of some of the reactions
        net.setattr_rxn(5   , 'Tu', 700.0 ); net.setattr_rxn(6   , 'Tl', 700.0);
        net.setattr_rxn(240 , 'Tu', 4150.0); net.setattr_rxn(241 , 'Tl', 4150.0);
        net.setattr_rxn(240 , 'Tu', 4150.0); net.setattr_rxn(241 , 'Tl', 4150.0);
        
        #setting the lower bounds of the temperature of some reactions
        # [[ID1, Tlb1], [ID2, Tlb2]...[]]
        ID_Tlb = [ [63 , 10.0 ], [64, 10.0  ], [80, 300.0 ],  [131, 100.0], 
                   [196, 300.0], [210, 200.0], [234, 200.0],  [236, 200.0],
                   [273, 200.0], [275, 300.0], [298, 200.0],  [337, 250.0],
                   [343, 200.0], [348, 500.0], [354, 80.0 ],  [360, 200.0],
                   [363, 300.0], [366, 200.0], [374, 25.0 ],  [379, 150.0],
                   [387, 300.0], [393, 13.0 ], [400, 300.0],  [414, 200.0],
                   [425, 200.0], [432, 200.0], [4278,300.0]
                 ]
        
        for ID, Tlb in ID_Tlb:
            print 'setting a lower bound to the T to be used in computing the rxn cst'
            net.setattr_rxn(ID, 'Tlb', Tlb)
            print '       ',
            net.print_reactions(ID, fmt='type rxn trng Tlb')
        
        #in looking at chemcial_balance.c we see that all the reactions which have no
        #Tl and Tu set, the rates are valid for all temperatures, so we set those to
        #the min and max of UMIST Tl=10 and Tu=41000
        for rxn in net.reactions:
            if rxn.Tl == -1 and rxn.Tu == -1:
                rxn.Tl = 10.0
                rxn.Tu = 41000.0
                print 'set trange manually ', 
                rxn.display(fmt='rxn trng')
        
        #setting the functions which compute the reactions constants
        print 'setting the functions which will compute the reaction constants'
        net.set_rxn_cst_functions()
        
        #checking if all the reactions have the functions set
        found = 0
        for rxn in net.reactions:
            if rxn.compute_rxn_cst_func == None:
                found = 0
                rxn.display(fmt='id rxn')
                
        if found > 0:
            raise ValueError('one or more reaction does not have a rxn cst function set...see above..')
        else:
            print 'reaction constant functions for all the reactions are set'
            
        #now that the ranges of some of the reactions are set/modified, we can merge them
        #and set the functions which compute the constants give a temperature 
        net.merge_identical_reactions()
        net.set_all_rxn_bounds()
        ###------------------done setting up the network as in the pdr code---------------
        
        self.setChemicalNetwork(net) # assiginig the chemical network to the archive
            
    def getQuantityFromAllMeshes(self, quantity, slabIdx = None, arrIdx = None):
        """ gets the quantity from all the meshes and returns it as a numpy array. 
            the quantity is mandatory, but no the slabIdx.
            
            :param int arrIdx: in case the quantity pointed to is an array, arrIdx would be the 
              index in the array we would like to retrieve.
        """

        values = np.zeros(self.nMeshes, dtype = np.float64)
        
        if 'from_meshes_info' in quantity:
            #getting the values from the info corresponding to each mesh in self.infoMeshes. In this case 
            # quantity should be as follows ['from_meshes_info', 'info|parms', index]
            
            q_name, q_indx = quantity[1], quantity[2] 
            
            for i in np.arange(self.nMeshes):

                q = fetchNestedDtypeValue(self.infoAll[i], [q_name] )
            
                values[i] = q[q_indx]  
        
        else:
            #getting the values from the meshes data self.meshes
             
            for i in np.arange(self.nMeshes):
                
                q = fetchNestedDtypeValue(self.meshes[i], quantity )
            
                if slabIdx != None and arrIdx != None:
                    v = q[arrIdx][slabIdx]
                elif slabIdx != None and arrIdx == None:
                    v = q[slabIdx]
                else:
                    v = q
    
                values[i] = v  
        #
        
        return values

    def apply_function_to_all_meshes(self, func, func_kw = None):
        """Applies func to all the meshes in self.meshes. func_kw are passed to func. The results are 
        returned as a list which is the same length as self.nMeshes"""
        
        if func_kw == None: func_kw = {}
            
        dataRet = []

        for mesh in self.meshes:            
            self.mshTmp.setData(mesh)
            v = func(self.mshTmp, **func_kw)
            dataRet.append(v)
        return dataRet

    def apply_function_to_all_radex_meshes(self, func, func_kw = None):
        """Applies func to all the radex meshes in self.meshesRadex. func_kw are passed to func. The 
        results are returned as a list which is the same length as self.nMeshes"""

        if func_kw == None: func_kw = {}
        
        dataRet = []
        
        for mesh in self.meshesRadex:
            v = func(mesh, **func_kw)
            dataRet.append(v)
        return dataRet
        
        
    def construct3DInterpolationFunction(self, quantity = None, slabIdx = None, arrIdx = None, log10 = None, values = None, 
                                         data = None, interpolator = None, remove_nan_inf = None, 
                                         f_interp_dim = None, *args, **kwargs):
        """Returns a 3D interpolation function (interpolate.LinearNDInterpolator) which
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
           :param numpy.ndarray data: the coordinates of the mesh points x,y,z to be used 
           :param string interpolator: 'nearest'|'linear'|'both' uses scipy.interpolate.NearestNDInterpolator or 
             scipy.interpolate.LinearNDInterpolator respectively if provided.  By default, linear interpolation
             is used. 
             
           .. todo:: modify this to construct the table over a selected range of the 3D 
               parameter space instead of using all the meshes. (optional)
                
           .. warning:: this fails if all the entries in one of the dimensions have the exact same value. in that case use :data:`construct2DInterpolationFunction`
        """
        
        # checking if the grid_x,y,z are set, if not, set them        
        self.set_grid_axes_quantity_values(**kwargs)
        self.set_attributes(**kwargs)
        
        if f_interp_dim == None:
            f_interp_dim = '3D'
        
        if values == None:
            # the quantitiy we are intiesrested in showing
            values = self.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx, arrIdx = arrIdx)
        else:
            values = values

        if log10 != None and log10 == True:
            values[:] = np.log10(values[:])

        if data == None:
            if f_interp_dim == '3D':
                data = np.array([self.grid_x, self.grid_y, self.grid_z]).T  #3D coordinates
            elif f_interp_dim == '2D':
                data = np.array([self.grid_x, self.grid_y]).T  #2D coordinates                                
            else:
                raise ValueError('dim of interpolation function is unknow.')
        else:
            data = data
        
        #make sure the data points and the values are of compatible shapes
        if data.shape[0] != values.size:
            raise ValueError('data and values are of incompatible shape')
            
        if remove_nan_inf != None and remove_nan_inf == True:
            inds_valid = numpy.where(numpy.isfinite(values))[0]
            
            if inds_valid.size == 0:
                raise ValueError('no valid points were found')
            
            values = values[inds_valid]
            data = data[inds_valid, :]
             
        # getting the interpolation function
        ti = time()
        if interpolator == None or interpolator == 'linear':        
            f = interpolate.LinearNDInterpolator(data, values) 
        elif interpolator == 'nearest':
            f = interpolate.NearestNDInterpolator(data, values) # getting the interpolation function
        tf = time()
        self.logger.debug('constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti))

        return f

    def computeInterpolated2DGrid(self, ranges = None, res = None, zSec = None, 
                                  fInterp = None, *args, **kwargs):
        """Returns a 2D array ( a numpy ndarray ) of size res[0] and res[1] (check this if it is not the reverse) which holds
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
        
        if zSec != None:
            # setting the values of mechanical heating to interpolate on
            zNew[:] = zSec      # all the grid points have the same mechanical heating
            dataNew = np.array( [xNew, yNew, zNew] ).T            
        else:
            dataNew = np.array( [xNew, yNew] ).T

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

    def construct_4D_interpolation_function_sectioned(self, ):
        pass 
    
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
            quantity = self.parms['gridsInfo']['11']['quantity']

            #----------------------------emissions using radex---------------------------------------
            if radiativeCodeType == 'radex':

                values, grid_x, grid_y, grid_z = [], [], [], []

                transitionIdx = self.parms['gridsInfo']['11']['transitionIndx']
            
                #looping over all the radex info of all the mesehes and getting the
                #transition data to be displayed
                for i in np.arange(self.nMeshes):
                    
                    if self.meshesRadex[i] != None:
                            
                            x, y, z = self.grid_x[i], self.grid_y[i], self.grid_z[i]
                            v = self.meshesRadex[i][transitionIdx][quantity]
                            
                            # appending the data and the value to a list to be used
                            # in making the interpolation function
                            
                            values.append(v)
                            grid_x.append(x)
                            grid_y.append(y)
                            grid_z.append(z)
                            
                # getting the data in the shape that is accepted by the interpolation construction        
                data   = np.array([grid_x, grid_y, grid_z], 'f8').T 
                values = np.array(values, 'f8')
                
                self.intensityGridInterp_f = self.construct3DInterpolationFunction(data   = data,
                                                                                   values = values,
                                                                                   log10  = True,
                                                                                   remove_nan_inf = True,
                                                                                   *args, **kwargs)

            #----------------------------emissions from the PDR slabs---------------------------------------
            if radiativeCodeType == 'pdr':

                Av_max = self.parms['gridsInfo']['11']['Av_max']

                #getting the column densities for all the models
                values  = np.ndarray(self.nMeshes, dtype = np.float64)
                #a temporary object used to calculate stuff from a pdr mesh object
                m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
                
                for i, data in enumerate(self.meshes):
                    
                    m.setData(data)
                    v = m.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])

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
            
           .. todo:: plot every other labeled contour as a tick in the colorbar
           
           .. warning:: there would be a memeory leak if things are plotted over and over
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
        
        self.grds[0][0][:] = grd        
        
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
        
        self.grds[0][1][:] = grd        

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

        self.grds[1][0][:] = grd        

    def showLineIntensityGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the line intensity grid"""            
        panel = self.gui['maps2d']['11']
        pltParms = self.parms['gridsInfo']['11']
        
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
        if pltParms['showContours'] == True:
            panel['contour'] = panel['axes'].contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
            panel['axes'].clabel(panel['contour'],levels, fmt = '%.1f' )
        
        pyl.colorbar(im11, cax = panel['axesCbar'], ax = panel['axes'], orientation = 'horizontal', ticks = levels[::2], format = '%.1f')

        #formatting the string of the title of the emission panel
        titleStr  = 'Intensity '
        
        if pltParms['type'] == 'radex':
            titleStr += '%s-%s' % (self.currentRadexDb['specStr'], pltParms['transitionIndx']) 
        elif pltParms['type'] == 'pdr':
            titleStr += '%s-%s' % (pltParms['quantity'][1], pltParms['quantity'][3])

        titleStr += ' Av = [0,%2.f]' % self.parms['gridsInfo']['11']['Av_max']
        panel['axesCbar'].set_title(titleStr)

        self.grds[1][1][:] = grd
        
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
        
        #list which will hold the transitions for the models (intiialized to None)
        meshesRadex = [None for i in np.arange(self.nMeshes)] 

        radex_parms = self.parms['radex']
        
        Av_range = radex_parms['Av_range']
        Av_end   = Av_range[1]
        specStr  = radex_parms['specStr']
        
        #
        #
        def put_radex_transitions_in_meshesRadex_array(mesh_idx):
            
            i = mesh_idx

            #utility object used to get info from the PDR and feeding them to Radex to get emissions 
            radex_obj_utility = self.radexObj.copy()              
            pdr_mesh_obj_utility = self.mshTmp.copy()
            
            pdr_mesh_obj_utility.setData( self.meshes[i] )

            has_lines, radex_parm_from_pdr_mesh = self.run_radex_on_pdr_model(pdr_mesh_obj = pdr_mesh_obj_utility, 
                                                                              radex_obj = radex_obj_utility,
                                                                              radex_parms = radex_parms,    
                                                                              Av_range = Av_range)
            
            if has_lines == None:
                infoAllRadex[i]['info'][1] = 0 #no trainsitions
                meshesRadex[i] = None
                radex_obj_utility.printSetFlags()
            else:                
                #setting the basic radex run info into the attribute
                #-----------------------------------------------------------
                #saving the number of transitions into radex info attribute
                infoAllRadex[i]['info'][0] = i
                if radex_obj_utility.flag_is_set('ERROR') or (has_lines is False): # radex failed, no transition data available
                    infoAllRadex[i]['info'][1] = 0 #no trainsitions
                    meshesRadex[i] = None
                    radex_obj_utility.printSetFlags()
                else: # radex succeeded, wirint the transition data
                    infoAllRadex[i]['info'][1] = radex_obj_utility.nTransitions
                    meshesRadex[i] = radex_obj_utility.transitions
    
                #infoAllRadex[i]['info'][2] = NOT SET HERE, IT IS SET WHEN WRITING THE DB TO A FILE
                
                #setting the status into the attribute
                infoAllRadex[i]['info'][3] = radex_obj_utility.status
                #appending the transition data of this mesh to the attribute 
                #-----------------finished saving the info------------------
                
                
                if radex_obj_utility.flag_is_set('SUCCESS'):
                    self.logger.debug('radexGrid : converged with no warnings')
                else:
                    self.logger.debug('radexGrid : converged with warnings')
                    self.logger.debug('------------------------------------')
                    print radex_obj_utility.getWarnings()
                    self.logger.debug('------------------------------------')
    
                self.logger.debug('---------------------------------------------------------')
                
            del radex_obj_utility               
            del pdr_mesh_obj_utility
        #
        #
        
        t0 = time()

        #running radex on all the pdr meshes and storing them into meshesRadex
        if self.parms['nThreads'] == 1:
            for i in np.arange(self.nMeshes):
                self.logger.debug('pdr mesh index =  %d : ' % i)
                put_radex_transitions_in_meshesRadex_array(i)
        else:
        #running radex instances in parallel on all the pdr meshes and storing them into meshesRadex
        
            nThreads = self.parms['nThreads']
            nEachThread = self.nMeshes/nThreads
            
            #
            #
            def run_bunch(tidx, conn):

                start_idx = tidx*nEachThread
                                
                #setting the loads                
                if tidx != nThreads - 1:
                    end_idx = (tidx+1)*nEachThread 
                else: 
                    end_idx = self.nMeshes 

                #runing the bunch of meshes for each thread
                for mesh_idx in numpy.arange(start_idx, end_idx):
                    put_radex_transitions_in_meshesRadex_array(mesh_idx)

                #sending the array infoAllRadex                
                conn.send([start_idx, end_idx])
                conn.send(infoAllRadex[start_idx:end_idx])
                #sending the transitions info
                for mesh_idx in numpy.arange(start_idx, end_idx):
                    conn.send([meshesRadex[mesh_idx]])
            #
            #
                        
            #making the pipes for communicating the results, pipes[i][0] = parent conn, pipes[i][1] child conn
            pipes = [Pipe() for tidx in range(nThreads)]
            #making the proccesses
            ps = [Process(target=run_bunch, args=(tidx, pipes[tidx][1])) for tidx in range(nThreads) ]
            #running  
            for p in ps: p.start()
            
            #receiving the data (note : if data grater 32MB might raise value errors)
            for tidx in range(nThreads):
                #recieving the indicies of the meshes and the meshInfoRadex array for thread tidx
                start_idx, end_idx = pipes[tidx][0].recv()
                infoAllRadex[start_idx:end_idx] = pipes[tidx][0].recv()
                #recieving the transitions info
                for mesh_idx in numpy.arange(start_idx, end_idx):
                    meshesRadex[mesh_idx] = pipes[tidx][0].recv()[0]
                
            #joining the threads
            for p in ps: p.join()
            
        #done getting the transitions from all the PDR models
        
        self.logger.debug('Time running in parallel = %.2f seconds ' % (time()-t0))
        
        #copying the mesh parameters self.infoAll[:]['parms] to self.infoAllRadex[:]['parms'] (info are set in the loop above)
        for i in np.arange( self.nMeshes ):
            infoAllRadex[i]['parms'][:] = self.infoAll[i]['parms']

        #setting state attributes related to the current radexDb
        self.append_radex_db_related_attributes(meshesRadex, infoAllRadex, Av_end, specStr)
        
        if writeDb == True:
            self.writeDbRadex()
    
    def append_radex_db_related_attributes(self, meshesRadex, infoAllRadex, Av_end, specStr):
        """This method organizes and append a newly computed/read radex Db to the attributes
        and in memory. 
        """
        
        self.infoAllRadex = infoAllRadex
        self.meshesRadex = meshesRadex

        #storing the read Db into a nested dict
        self.currentRadexDb = {'Av' : Av_end, 'specStr': specStr}   
        radexDb = { 'meshes' : self.meshesRadex, 'infoAll' : self.infoAllRadex}
            
        #if the Av is not in self.radexDbs, make the empty dict object which will hold this Db and its info
        if '%.2f' % Av_end not in self.radexDbs:
            self.radexDbs['%.2f' % Av_end] = {}
                
        self.radexDbs['%.2f' % Av_end][specStr] = radexDb
        self.logger.debug('added %s radex database to memory' % specStr)
        
    def check_radex_databases_for_bogus_meshes(self):
        '''This method checks for radex meshes (self.meshesRadex) which could hold incorrect results. For example, 
        it checks for jumps in the population densities'''
                    
        inds_bogus = []
        
        for i in numpy.arange(self.nMeshes):
            
            if self.meshesRadex[i] != None:
                
                self.radexObj.transitions = self.meshesRadex[i] 
 
                check = True
                check &= self.radexObj.check_for_pop_density_continuity()
                check &= self.radexObj.check_for_pop_density_positivity() 
                check &= self.radexObj.check_for_pop_density_continuity()
                check &= self.radexObj.check_for_pop_density_all_non_zero()
                
                if check :
                    continue
                else:
                    #print i, 'possibly incorrect radex run'
                    #print self.grid_x[i], self.grid_y[i], 10.0**self.grid_z[i]
                    inds_bogus.append(i)
                    
        return inds_bogus
    
    def fix_bogus_radex_meshes(self):
        '''Checks for radex models which might be bogus, and re-runs radex with stricter convergence 
        conditions and re-updates the databases.  This method will look for all the radex
        meshes for a specific specie for all the Av loaded into memory.
        How to use this method? for example in the gui we see some jumps in the maps for the species CO for Av=4.0
        a) load the radex databases for CO using: 
            arxv.readDbs(allDbs=True)  #or arxv.readDbRadex(Av=4.0, specStr='CO')
        b) set the 'species' parameter below to the one we need to check, species='CO' in this case
        c) arxv.fix_bogus_radex_meshes()
        d) if all goes well, running arxv.fix_bogus_radex_meshes() shud not complain
        e) as a test, first do these steps with self.writeDbRadex() commented
        '''
        
        #species = ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS']
        species = ['CN']        
        
        for specCheck in species:
            
            self.parms['radex']['specStr']=specCheck
            self.setup_default_radex_instance(self.parms['radex'])
            
            for AvStrDb in self.radexDbs.keys():
                for specStrDb in self.radexDbs[AvStrDb]:
                    
                    if specCheck != specStrDb:
                        continue
                        
                    #loading a database
                    self.use_radexDb(Av=float(AvStrDb), specStr=specCheck, silent=False)
    
                    #finding the meshes whihc might have anomaleous data (if any)                   
                    inds = self.check_radex_databases_for_bogus_meshes()
                   
                    if len(inds) == 0:
                        continue
                   
                    self.logger.debug('indicies of anomaleous radex meshes')
                    print '                 ', inds
                   
                    #re-running radex on the false meshes
                    for ind in inds:
                       
                        print '               ',
                        print 'x, y, 10**z = ', self.grid_x[ind], self.grid_y[ind], 10.0**self.grid_z[ind]
                       
                        self.mshTmp.setData(self.meshes[ind])
                        self.compute_and_set_radex_curves(pdr_mesh_obj = self.mshTmp, compute_only=True)
                        self.meshesRadex[ind][:] = self.radexObj.transitions[:]
                   
                    #writing the updated database
                    self.writeDbRadex()
                    
    #################################################################################
    def setup_gui(self):
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
        transitionSelector['axes'] = gui['figure'].add_axes( [0.05, 0.55, 0.2, 0.02] )
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
        maps2d_10['axesCbar'].set_title('N(%s)' % self.parms['gridsInfo']['10']['specStr'])
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
        left  = 0.10
        bott  = 0.035
        sz    = 0.20
        vSpace = 0.01
        hSpace = 0.01
        AV['00'] = {'axes' : gui['figure'].add_axes([left              , bott + sz + vSpace, sz, sz])}
        AV['01'] = {'axes' : gui['figure'].add_axes([left + sz + hSpace, bott + sz + vSpace, sz, sz])}
        AV['10'] = {'axes' : gui['figure'].add_axes([left              , bott              , sz, sz])}
        AV['11'] = {'axes' : gui['figure'].add_axes([left + sz + hSpace, bott              , sz, sz])}
        #-----------------------------------------------------------------------------------        
        gui['AV'] = AV


        # the axes stuff where the emission lines will be plotted
        em_lines = {}
        #-----------------------------------------------------------------------------------
        left  = 0.65
        bott  = 0.07
        sz    = 0.08
        vSpace = 0.0
        em_lines['0'] = {'axes' : gui['figure'].add_axes([left, bott + 3*sz + vSpace , 3*sz, sz])}
        em_lines['1'] = {'axes' : gui['figure'].add_axes([left, bott + 2*sz + vSpace , 3*sz, sz])}
        em_lines['2'] = {'axes' : gui['figure'].add_axes([left, bott + 1*sz + vSpace , 3*sz, sz])}
        em_lines['3'] = {'axes' : gui['figure'].add_axes([left, bott + 0*sz          , 3*sz, sz])}
        em_lines['title1'] = pyl.figtext(0.65, 0.4, '')
        em_lines['title2'] = pyl.figtext(0.9, 0.2 , '')
        #-----------------------------------------------------------------------------------        
        gui['em_lines'] = em_lines
        
        #-------------------------------------------------
        return gui
    
    def update_gui_em_lines_titles(self, solver_instance, Av_range):
        '''
        Sets the title of the line emission panels in the gui. 
        
        :param solver_instance: The instanance of the code used to get the emissions (radex or despotic)
        :param Av_range: The range in Av used in the in PDR mesh to compute the emissions
        '''
        
        if isinstance(solver_instance, radex):
            #updating the display strings on the gui related to Radex
            #updating the ['radex']['title1'] and ['radex']['title2']
            #title1
            
            radexObj = self.radexObj
            
            strng = 'radex LVG data for species %s, Av = [%.2f, %.2f]' %  (radexObj.inFile['specStr'], Av_range[0], Av_range[1])
            self.gui['em_lines']['title1'].set_text(strng)
            #title2
            strng = 'gasT\n%f\nN(specie)\n%e\n' % (radexObj.inFile['tKin'], radexObj.inFile['molnDens'])
            for i, specStr in enumerate(radexObj.inFile['collisionPartners']):
                strng += '%s\n%e\n' % (specStr, radexObj.inFile['nDensCollisionPartners'][i])                
            self.gui['em_lines']['title2'].set_text(strng)
        
    def plotGrids(self, *args, **kwargs):
        """Main method for exploring the meshes in the database.
        
           .. todo:: change resGrids to a [res_x, res_y] insteads of it being just a scalar.
        """
        
        self.set_attributes(**kwargs)
        self.set_default_attributes()
        
        resGrids = self.parms['gridsRes']
        #setting the default zSec
        self.pltGmSec   = np.min(self.grid_z)
        self.logger.debug('set the current plotting section to %f, the minimum of the self.grid_z' % self.pltGmSec)
        self.radexParms = self.parms['radex']
        
        #setting up the gui attribute
        self.gui = self.setup_gui()

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
                                  self.gui['em_lines']['0']['axes'],
                                  self.gui['em_lines']['1']['axes'],
                                  self.gui['em_lines']['2']['axes'],
                                  self.gui['em_lines']['3']['axes']
                                  ]
                                 )
        
        # getting and plotting the unique sections in z
        self.set_unique_grid_sections(1e-13) 
        self.gui['widgets']['zSecSelector']['pointsUnique'].set_xdata(self.grid_z_unique)
        self.gui['widgets']['zSecSelector']['pointsUnique'].set_ydata(self.grid_z_unique*0.0 + 0.5)

        # attaching mouse click event to fig 1
        self.gui['figure'].canvas.mpl_connect('button_press_event', self.on_button_down)
        
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
        if self.parms['gridsInfo']['11']['type'] == 'radex' and self.parms['gridsInfo']['11']['show']:
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
    
    def on_button_down(self, event):
        """method called on the event of a mouse button down in self.fig. See self.setup_gui()
           for the layout of the window.
        """
       
        ti = time()

        # get the x and y coords, flip y from top to bottom
        xd, yd = event.xdata, event.ydata
        
        if event.button == 1:

            clickedInAxes = False
            
            # getting the value of the section in z to display
            #---------------------------------------------------------------
            if self.gui['widgets']['zSecSelector']['axes'].contains(event)[0]:
                clickedInAxes = True
                # setting the section closest to the data available
                inds = np.argmin(np.fabs( self.grid_z - xd ) )
                self.pltGmSec = self.grid_z[inds]
                self.gui['widgets']['zSecSelector']['point'].set_xdata(self.pltGmSec)
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
                    
                self.plotThisSec() #TODO:: rename this to update 2D grids
                pyl.draw()

            # updating the grid to the selected transition cliked
            #----------------------------------------------------------------------
            if self.gui['widgets']['transitionSelector']['axes'].contains(event)[0]:
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

            # selecting a mesh from points: getting the coordinates inthe grid to display 
            # as a single mesh as a function of Av and maybe display the radex calculations
            #------------------------------------------------------------------------------------
            if self.gui['ax2d']['axes'].contains(event)[0] or self.gui['maps2d']['00']['axes'].contains(event)[0] or self.gui['maps2d']['01']['axes'].contains(event)[0] or self.gui['maps2d']['10']['axes'].contains(event)[0] or self.gui['maps2d']['11']['axes'].contains(event)[0]:
                clickedInAxes = True
                self.logger.debug('####################################button 1 clicked###########################################')
                indMin = self.get_mesh_index(x = xd, y = yd, z = self.pltGmSec)
                self.mshTmp.setData(self.meshes[indMin])
                
                #updating the title of ['ax2d']['axes']
                strng = '$\log_{10} n_{gas} = $ %4.2f\n$\log_{10} G_0 =$ %4.2f\n$\log_{10} \Gamma_{mech} = $  %5.2f\n' %  (np.log10(self.mshTmp.data['hdr']['nGas']), np.log10(self.mshTmp.data['hdr']['G0']), np.log10(self.mshTmp.data['hdr']['gammaMech']))
                strng += '$\Gamma_{mech}$ = %.2e\n' % self.mshTmp.data['hdr']['gammaMech']   
                
                self.gui['ax2d']['axes'].set_title(strng)
                
                self.gui['ax2d']['pts2'].set_xdata( self.grid_x[indMin] )
                self.gui['ax2d']['pts2'].set_ydata( self.grid_y[indMin] )
                self.gui['ax2d']['pts2'].set_color('r')

                self.mshTmp.plot() #plotting the PDR mesh curves in the panels (vs Av)
                
                if self.parms['radex']['use']:
                    
                    if self.parms['radex']['plot_model_from_Db'] == False:
                        self.compute_and_set_radex_curves(pdr_mesh_obj = self.mshTmp)
                    else:
                        if self.meshesRadex[indMin] != None:
                            self.set_radex_curves_from_db(indMin)
                        else:
                            self.logger.debug('No radex data available for this PDR model clicked on.')
                    
                pyl.draw()

            # setting the Av of the position clicked on the plots of the current mesh
            # to the chemical netowrk object.
            #----------------------------------------------------------------------------

            if self.gui['AV']['00']['axes'].contains(event)[0] or self.gui['AV']['01']['axes'].contains(event)[0] or self.gui['AV']['10']['axes'].contains(event)[0] or self.gui['AV']['11']['axes'].contains(event)[0]:
                
                clickedInAxes = True
                Av_clicked = xd
                chemNet = self.mshTmp.chemNet
                
                self.logger.debug('####################################button 1 clicked###########################################')
                self.logger.debug('Av clicked = %.4f' % Av_clicked)

                #computing the reaction rates at the closes Av available in the pdr mesh
                self.compute_rxn_rates_at_Av(Av_clicked)
                #the used Av in the pdr mesh
                Av_use = self.mshTmp.chemNet.Av

                #plotting the vertical lines on the gui indicating the positions
                #in the slab used for the chemistry
                self.mshTmp.plot_v_lines_used_in_chemnet()
                
                specStr = self.parms['gridsInfo']['10']['specStr']

                self.logger.debug('dominat reactions destroyin %s' % specStr)
                IDs_filtered = chemNet.filter_reactions(withReacts=specStr, 
                                                        fmt='id type rxn trng cst rate', 
                                                        show=10, 
                                                        sort=True)
                self.logger.debug('dominat reactions producing %s' % specStr)
                IDs_filtered = chemNet.filter_reactions(withProds=specStr, 
                                                        fmt='id type rxn trng cst rate', 
                                                        show=10, 
                                                        sort=True)
                
                self.logger.debug('set the environment variable to the chemical netowrk of the mesh.')
            
            
                if self.parms['radex']['use']:

                    self.compute_and_set_radex_curves(meshObj = self.mshTmp,
                                                      Av_range = [0, Av_use])

                pyl.draw()

            if clickedInAxes == True:
                tf = time()
                self.logger.debug('time elapased after click : %f\n' % (tf - ti))

        if event.button == 3: #upon a right click
            print 'button 3 pressed'

            self.plot_integrated_emissions()


    def replace_radexObj_transitions_into_radexDb(self, no_write=None, no_refresh=None):
        '''This is an auxiliary method which replaces the transitions data in self.radexObj 
        (corresponding to a mesh with coordinates self.mshTmp into the array position corresponding 
        to self.meshesRadex.  The z quantity used is the plotting gmech section
        
        If the optional keyword write is passed, the Db is written to the appropriate file
        corresponding to it Av and specie.
        
        Howto use : 
          1- assigm a the mesh data to arxv.mshTmp by clicking somewhere on a grid
          2- call arxv.replace_radexObj_transitions_into_radexDb(write=True)
          3- call arxv.computeAndSetInterpolationFunctions()
          4- press on the z selectior of the current section to refresh the maps
        '''
        
        x, y = np.log10(self.mshTmp.data['hdr']['nGas']), np.log10(self.mshTmp.data['hdr']['G0'])
        z = self.pltGmSec #using as 'z' the plotting gmech
        
        ind = self.get_mesh_index(x = x, y = y, z = z)
        
        self.meshesRadex[ind][:] = self.radexObj.transitions[:]
        
        if no_write != None and no_write == False:
            pass
        else:
            self.writeDbRadex()

        if no_refresh != None and no_refresh == False:
            pass
        else:
            self.computeAndSetInterpolationFunctions()
    
    def run_radex_on_pdr_model(self, pdr_mesh_obj = None, radex_obj = None, radex_parms = None, Av_range = None):
        '''
        Runs radex on a PDR model and return the output in self.radexObj.
        
        :param pdr_mesh_obj: The mesh object from which the parameters needed by radex will be extracted and used.
        :param radex_parms: The paramteres needed by radex (the parameters should be a dict of the form self.parms['radex'])
        :param Av_range: The range in Av of the pdr model to be considered.
        '''
        
        write_debug_info = False
        
        #getthing the stuff radex needs from the PDR mesh
        radex_parm_from_pdr_mesh = pdr_mesh_obj.get_radex_parameters(speciesStr = radex_parms['specStr'], 
                                                                     threshold  = radex_parms['xH2_Min'],
                                                                     Av_range   = Av_range)
        
        (gasTRadex, nColls, colDensThisSpec, Av_range_used,) = radex_parm_from_pdr_mesh
         
        #determining what gas density to use for the collider for H2
        if radex_parms['use_pdr_gas_den_H2']:
            #using nGas/2 as the density of H2 (allowed when the only collider is H2)
            nDensColls = [pdr_mesh_obj.data['hdr']['nGas']/2.0]
            collsStr   = list(radex_parms['collisionPartners'])
            
            if len(radex_parms['collisionPartners']) > 1:
                raise ValueError("""collisionsPartners can be only H2 when setting the density 
                of the collider H2 to the full gas density""")
        else:
            #using the weighted densities extracted from the PDR
            
            # getting the collider densities in the same order of the supplied input spcie string list 
            nDensColls = [ nColls[collSpecStr] for collSpecStr in radex_parms['collisionPartners'] ]
            collsStr   = list(radex_parms['collisionPartners'])
            #print 'input coll species', self.radexParms['collisionPartners'] 
            #print 'nColls after putting them in the right order = ', nDensColls

        #print '================='
        #print collsStr
        #print nDensColls

        radex_obj.setInFileParm('specStr', radex_parms['specStr'])
        radex_obj.setInFileParm('tKin', gasTRadex)
        radex_obj.setInFileParm('collisionPartners', collsStr)
        radex_obj.setInFileParm('nDensCollisionPartners', nDensColls)
        radex_obj.setInFileParm('molnDens', colDensThisSpec)
        
        #remove colliders which are underabundant (below radex limits)
        radex_obj.filterColliders()
        
        if len(radex_obj.inFile['collisionPartners']) == 0:
            self.logger.debug('not enough colliders')
            return None, None
        else:
            
            if write_debug_info:
                # writing the parameters to a file (for debugging purposes)
                # this file can be used to re-run radex as standalone
                fName = radex_parms['path'] + '-debug.inp'
                fObj = open(fName, 'w')
                fObj.write(radex_obj.genInputFileContentAsStr() )
                self.logger.debug('input radex file written to %s' % fName)
            
            if radex_parms['checkOutputIntegrity'] == False:
                radex_obj.setDefaultStatus()
                radex_obj.run(checkInput = True, verbose = radex_parms['verbose'])
                
                if radex_obj.flag_is_set('RUNOK'):
                    has_lines = True
                else:
                    return None, None
                
            else:
                #running radex (multiple times if necessary) for it to converge for this set of parms     
                status, has_lines = radex_obj.run_mutiple_trials(
                                                                 expected_sum_pop_dens = radex_parms['popDensSumExpected'],
                                                                 rel_pop_dens_tol = radex_parms['popDensSumTol'],   
                                                                 change_frac_trial = radex_parms['changeFracTrial'],    
                                                                 max_trials = radex_parms['nMaxTrial'],
                                                                 verbose = radex_parms['verbose'],
                                                                 strict = radex_parms['strict']
                                                                 )
        return has_lines, radex_parm_from_pdr_mesh 
    
    def compute_and_set_radex_curves(self, pdr_mesh_obj = None, radex_parms = None, Av_range = None, compute_only = None, radex_obj = None):
        """This is a utilty method (make it a private method), for populating the radex axes
          with the plots for the specie passed in the global parameter self.parms['radex']['specStr'].
          It takes as a paremeter the :data:`mesh` which will be used for doing the radex computations.
          
          if meshObj is not passed, self.mshTmp is used as the pdr mesh object.
          if Av_range is not passed, the maximum Av is used.
          if compute_only is passed, nothing is plotted, only the emissions are computed and set
          to self.radexObj
        """
        if pdr_mesh_obj == None: pdr_mesh_obj = self.mshTmp
        if radex_obj    == None: radex_obj = self.radexObj
        if radex_parms  == None: radex_parms = self.parms['radex']
        
        if compute_only == None:
            #setting the axes
            radex_obj.setupPlot(nx = 1, fig = self.gui['figure'], axs = self.pltRadex)
            radex_obj.set_logger(self.logger)
            #clearing them in case there was anything
            radex_obj.clearCurves()

        radexOutputMakesSense, radex_parm_from_pdr_mesh = self.run_radex_on_pdr_model(pdr_mesh_obj = pdr_mesh_obj,
                                                                                      radex_obj = radex_obj, 
                                                                                      radex_parms = radex_parms,    
                                                                                      Av_range = Av_range)

        #plotting only when the radex solution makes sense (pop dens dont add up to 
        #what they should be up to a certain tolerence)
        if radexOutputMakesSense:
            
            if compute_only == None:
                # plotting the data (even if it does not converge)                
                if radex_obj.flag_is_set('SUCCESS'):
                    radex_obj.plotModelInFigureColumn(
                                                      allTrans = np.arange(radex_parms['maxDisplayTranistion']),
                                                      inAxes = self.pltRadex,     
                                                      title = '',
                                                      em_unit=radex_parms['quantity']
                                                     )                             
                radex_obj.setLabels()
            else:
                radex_obj.print_warnings() #printing the warnings
        else:
            self.logger.debug('radex output doesnt make sesne')

        if compute_only == None:
            
            if radex_parm_from_pdr_mesh != None:
                Av_range_used = radex_parm_from_pdr_mesh[3]
            else:
                Av_range_used = [-1, -1]
                
            self.update_gui_em_lines_titles(radex_obj, Av_range_used)
        
    def set_radex_curves_from_db(self, ind, radex_obj = None):
        """This is a utilty method which plots the radex transitions from a the precomputed radex
          database. The Av of the current radex Db is used.
          
          if Av_range is not passed, currenlty used Av of the gui is used.
        """
        if radex_obj == None:
            
            if self.radexObj == None:
                self.setup_default_radex_instance(self.parms['radex'])
            else:   
                radex_obj = self.radexObj
        
        radex_parms = self.parms['radex']

        radex_obj.transitions = self.meshesRadex[ind] 
        
        radex_obj.plotModelInFigureColumn(allTrans = np.arange(radex_parms['maxDisplayTranistion']),
                                                     inAxes = self.pltRadex,     
                                                     title = '')                             
        #radex_obj.setLabels()
        strng = 'transitions from DB [index = %d], Av = %.2f of clicked mesh' % (ind, self.currentRadexDb['Av'])
        self.gui['em_lines']['title1'].set_text(strng)
        #Av_range_used = [0.0, self.currentRadexDb['Av']]
        #self.update_gui_em_lines_titles(radex_obj, Av_range_used)

    def compute_and_set_despotic_curves(self, meshObj = None, radex_parms = None, Av_range = None, compute_only = None):
        """This is a utilty method (make it a private method), for populating the radex axes
          with the plots for the specie passed in the global parameter self.parms['radex']['specStr'].
          It takes as a paremeter the :data:`mesh` which will be used for doing the radex computations.
          
          if meshObj is not passed, self.mshTmp is used as the pdr mesh object.
          if Av_range is not passed, the maximum Av is used.
          if compute_only is passed, nothing is plotted, only the emissions are computed and set
          to self.radexObj
          
          .. todo:: continue working on this
        """

        if meshObj == None:
            meshObj = self.mshTmp
          
        radexObj = self.radexObj
        
        radex_parms = self.parms['radex']
        
        #getting the parameters needed by despotic to compute the emissions
        weigted_parms = meshObj.get_weighted_averages(speciesStr = radex_parms['specStr'], 
                                                      threshold  = radex_parms['xH2_Min'],
                                                      Av_range   = Av_range)
        
        (TMean, nDenseColl, N_spec, nSpecMean, Av_range_used,) = weigted_parms
        
        nGas = meshObj.data['hdr']['nGas']
        specStr = radex_parms['specStr']
        
        mycloud = cloud()
        mycloud.nH =  nGas                                         #gas density
        mycloud.colDen = Av2NH(Av_range_used[1], self.metallicity) #cloud column density
        mycloud.sigmaNT = radex_parms['lineWidth']                 #non-thermal velocity despersion 
        mycloud.comp.xoH2 = 0.25                                   #ortho-H2 composition, xoH2 molecule per H nucleus
        mycloud.comp.xpH2 = 0.25                                   #para-H2 composition, xpH2 molecule per H nucleus
        mycloud.Tg = TMean                                         #cloud gas kinetic temperature
        mycloud.Td = 0.0                                           #cloud dust temperature
        
        mycloud.addEmitter(specStr, nSpecMean/nGas)  
        
        #computing the lines
        lines = mycloud.lineLum(specStr)
        
        #assinging the lines info to a radex instance (just for plotting purposes)
        radexObj.set_attributes_from_despotic(specStr, mycloud, lines)

        if compute_only == None:
            #setting the axes
            radexObj.setupPlot(nx = 1, fig = self.gui['figure'], axs = self.pltRadex)
            radexObj.set_logger(self.logger)
            #clearing them in case there was anything
            radexObj.clearCurves()
        
        """
        #determining what gas density to use for the collider for H2
        if radex_parms['use_pdr_gas_den_H2']:
            pass
            #using nGas/2 as the density of H2 (allowed when the only collider is H2)
            nDensColls = [meshObj.data['hdr']['nGas']/2.0]
            collsStr   = list(radex_parms['collisionPartners'])
            
            if len(radex_parms['collisionPartners']) > 1:
                raise ValueError(''''collisionsPartners can be only H2 when setting the density 
                of the collider H2 to the full gas density''')
        else:
            #using the weighted densities extracted from the PDR
            
            # getting the collider densities in the same order of the supplied input spcie string list 
            nDensColls = [ nColls[collSpecStr] for collSpecStr in radex_parms['collisionPartners'] ]
            collsStr   = list(radex_parms['collisionPartners'])
            #print 'input coll species', self.radexParms['collisionPartners'] 
            #print 'nColls after putting them in the right order = ', nDensColls
        """
        
        #print '================='
        #print collsStr
        #print nDensColls
        
        """
        radexObj.setInFileParm('specStr', radex_parms['specStr'])
        radexObj.setInFileParm('tKin', gasTRadex)
        radexObj.setInFileParm('collisionPartners', collsStr)
        radexObj.setInFileParm('nDensCollisionPartners', nDensColls)
        radexObj.setInFileParm('molnDens', colDensThisSpec)
        
        #remove colliders which are underabundant (below radex limits)
        radexObj.filterColliders()
        """
        
        if compute_only == None:
            #updating the display strings on the gui related to Radex
            #updating the ['radex']['title1'] and ['radex']['title2']
            #title1
            strng = 'radex LVG data for species %s, Av = [%.2f, %.2f]' %  (radex_parms['specStr'], Av_range_used[0], Av_range_used[1])
            self.gui['em_lines']['title1'].set_text(strng)
            #title2
            strng = 'gasT\n%f\nN(specie)\n%e\n' % (TMean, nSpecMean)
            #for i, specStr in enumerate(radexObj.inFile['collisionPartners']):
            #    strng += '%s\n%e\n' % (specStr, radexObj.inFile['nDensCollisionPartners'][i])                
            self.gui['em_lines']['title2'].set_text(strng)
            
        radexObj.plotModelInFigureColumn(allTrans = np.arange(radex_parms['maxDisplayTranistion']),
                                         inAxes = self.pltRadex,     
                                         title = '')                             
        radexObj.setLabels()            
    
        """
        if len(radexObj.inFile['collisionPartners']) == 0:
            self.logger.debug('not enough colliders')
        else:
            # writing the parameters to a file (for debugging purposes)
            # this file can be used to re-run radex as standalone
            fName = radex_parms['path'] + '-debug.inp'
            fObj = open(fName, 'w')
            fObj.write(radexObj.genInputFileContentAsStr() )
            self.logger.debug('input radex file written to %s' % fName)
            
            if radex_parms['checkOutputIntegrity'] == False:
                radexObj.setDefaultStatus()
                radexObj.run(checkInput = True, verbose = radex_parms['verbose'])
                radexOutputMakesSense = True
            else:
                #running radex (multiple times if necessary) for it to converge for this set of parms     
                status, radexOutputMakesSense = radexObj.run_mutiple_trials(expected_sum_pop_dens = radex_parms['popDensSumExpected'],
                                                                            rel_pop_dens_tol = radex_parms['popDensSumTol'], 
                                                                            change_frac_trial = radex_parms['changeFracTrial'],    
                                                                            max_trials = radex_parms['nMaxTrial'],
                                                                            verbose = radex_parms['verbose'])
            
            #plotting only when the radex solution makes sense (pop dens dont add up to 
            #what they should be up to a certain tolerence)
            if radexOutputMakesSense:
                
                if compute_only == None:
                    # plotting the data (even if it does not converge)
                    if radexObj.flag_is_set('SUCCESS'):
                        radexObj.plotModelInFigureColumn(allTrans = np.arange(radex_parms['maxDisplayTranistion']),
                                                         inAxes = self.pltRadex,     
                                                         title = '')                             
                    radexObj.setLabels()
                else:
                    radexObj.print_warnings() #printing the warnings
        """
        print 'bbbbbbbbbbbbbbbbbbbbb'

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

          .. warning:: check and test this methods and document it better.
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
        """Assigns the values of self.grid_x,y,z from self.grid_qx,qy,qz. The values
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
         
           .. note:: by default, the log of the quantities from the meshes are used. Only
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
                parms_ref_grid = {'radex': {'use': False}, }
                arxvRef = meshArxv(readDb = True, dirPath = referenceDbDirPath, **parms_ref_grid)
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
        
    def get_unique_grid_x_sections(self, relDiff):
        """returns the unique section in x up to a relative difference of relDiff"""
        return self.get_unique_grid_sections('x', relDiff)
    
    def get_unique_grid_y_sections(self, relDiff):
        """returns the unique section in y up to a relative difference of relDiff"""
        return self.get_unique_grid_sections('y', relDiff)
    
    def get_unique_grid_z_sections(self, relDiff):
        """returns the unique section in z up to a relative difference of relDiff"""
        return self.get_unique_grid_sections('z', relDiff)
    
    def set_unique_grid_sections(self, relDiff):
        """sets the attributes self.grid_x_unique, self.grid_y_unique, self.grid_z_unique"""
        self.grid_x_unique = self.get_unique_grid_x_sections(relDiff) 
        self.grid_y_unique = self.get_unique_grid_y_sections(relDiff) 
        self.grid_z_unique = self.get_unique_grid_z_sections(relDiff) 
        
    def get_unique_grid_sections(self, axis, relDiff):
        """returns the unique section in up to a relative difference of relDiff of one of the
        axes x,y or z specified by the axis argument ('x' for self.grid_x, 'y' for self.grid_y
        and 'z' for self.grid_z)."""
        
        if axis == 'x': v = self.grid_x
        if axis == 'y': v = self.grid_y
        if axis == 'z': v = self.grid_z
        
        sections = []

        while v.size != 0:
            mn = np.min(v)
            #print 'min = ', mn, 10.0**mn
            if mn == 0.0:
                inds = np.where( v != mn)
            else:
                inds = np.where( abs(1.0 - v/mn) > relDiff)
            v = v[inds]
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
                v = m.compute_integrated_quantity(quantity)
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
        
        pdrGuessDb = guess_db(path = path)

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
        self.mshTmp.data = self.meshes[ind] 
        return  self.mshTmp
    
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
    
    def setup_default_radex_instance(self, radex_parms):
        """sets up an instance of the object 'radex' to after using
           the parameters from radex_parms and setting 'tKin', 'nDensCollisionPartners'
           and 'molnDens' as None. Also sets the attribute self.radexObj
        """
        #making the instance            
        radexObj = radex(radex_parms['path'], radex_parms['molDataDirPath'], logger = self.logger)
        #setting some default radex paraemeters
        inFile = {'specStr'                : radex_parms['specStr'],
                  'outPath'                : 'foo' ,
                  'freqRange'              : radex_parms['freqRange'],
                  'tKin'                   : None  ,  # depending on the model
                  'collisionPartners'      : radex_parms['collisionPartners'],
                  'nDensCollisionPartners' : None  ,  # depending on the model
                  'tBack'                  : radex_parms['tBack'],  # fixed
                  'molnDens'               : None  ,  # depending on the model
                  'lineWidth'              : radex_parms['lineWidth'],  # fixed
                  'runAnother'             : 1     }
        radexObj.setInFile( inFile )
        
        self.radexObj = radexObj
        
    def plot_integrated_emissions(self):
        """plots quantities of a PDR mesh as a function of Av, also some
        quantites from the radexDbs (self.radexDbs)
        """
        xrng = [0.01, 30.0]  
        yrng = [1e-12, 1.0]
        xscale = 'log' #'linear' #'log'
        yscale = 'log' #'linear' #'log'
        
        fig = pylab.figure(figsize=(8,6))
        ax = fig.add_axes([0.15, 0.1, 0.6, 0.8])
        
        m = self.mshTmp

        lG0 = numpy.log10(m.data['hdr']['G0'])
        lnGas = numpy.log10(m.data['hdr']['nGas'])
        lGmech = numpy.log10(m.data['hdr']['gammaMech'])
        plotTitle =  '$\log_{10} n_{gas} = $ %4.2f  $\log_{10} G_0 =$ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f' %  (lnGas, lG0, lGmech)
        
        titles = []  #will store legend entries of the plots 
        plots = []   #will store the line object of each plot
        symbols = [] #the symbols to be used in plotting the lines
        #-------------------computing stuff from the PDR meshes--------------        
        quantities = []
    
        #-------------------------cooling components--------------------
        quantities.append(['therm','cooling'])
        titles.append(r'$\Lambda_{total}$'); symbols.append('-x')
        quantities.append(['cooling','metaStable'])
        titles.append(r'$\Lambda_{MS}$');  symbols.append('-')
        quantities.append(['cooling','fineStructure'])
        titles.append(r'$\Lambda_{FS}$');  symbols.append('-')
        quantities.append(['cooling','roVib'])
        titles.append(r'$\Lambda_{RV}$');  symbols.append('-')
        quantities.append(['cooling','lymanAlpha'])
        titles.append(r'$\Lambda_{LyA}$');  symbols.append('-')
        #--------------------------cooling components--------------------

        #-------------------------heating components--------------------
        quantities.append(['therm','heating'])
        titles.append(r'$\Gamma_{total}$'); symbols.append('-o')        
        quantities.append(['heating','photo'])
        titles.append(r'$\Gamma_{photo}$'); symbols.append('->')
        quantities.append(['heating','cr'])
        titles.append(r'$\Gamma_{CR}$'); symbols.append('-<')
        #--------------------------cooling components--------------------
        
        #-------------------------fine structure lines--------------------------
        '''
        quantities.append(['fineStructureCoolingComponents','C+','rate','1-0'])
        titles.append(r'CII');  symbols.append('-')
        quantities.append(['fineStructureCoolingComponents','C','rate','1-0'])
        titles.append(r'CI1-0');  symbols.append('-')
        quantities.append(['fineStructureCoolingComponents','C','rate','2-1'])
        titles.append(r'CI2-1');  symbols.append('-')
        quantities.append(['fineStructureCoolingComponents','O','rate','1-0'])
        titles.append(r'OI');  symbols.append('-')
        '''
        #-------------------------fine structure lines---------------------------
        
        Avs = m.data['state']['Av']
        
        for i, quantity in enumerate(quantities):
            
            vs = []
            for Av in Avs: 
                v1 = m.compute_integrated_quantity(quantity, Av_range = [0.0, Av])
                vs.append(v1)

            if 'fineStructureCoolingComponents' in quantity[0]:
                flux_total = m.compute_integrated_quantity(quantity, Av_range = [0.0, 30.0])
                flux_upto_Av = m.compute_integrated_quantity(quantity, Av_range = [0.0, 10.0])
                print 'contributions from %s(%s) 0  up to Av=10 is %.2f percent' % (quantity[1],quantity[3],100.0*flux_upto_Av/flux_total)
                            
            vs = numpy.array(vs)
            
            plt, = ax.semilogy(Avs, vs, symbols[i])

            plots.append(plt)
        #-------------------done computing stuff from the PDR meshes--------------        

        #-------------------extracting and plotting stuff from the radex Dbs------
        """
        specStrs    = ['CO' , '13CO'  , 'HCN' , 'HNC', 'CS']
        transitions = [0    ,   0   ,   0   ,   0  ,  0  ]
        titleStrs   = ['1-0',  '1-0', '1-0', '1-0' , '1-0']
        quantity    = 'fluxcgs'
                
        for i,specStr in enumerate(specStrs):
            Avs, vs = self.get_quantity_from_radex_meshes_vs_Av(mesh_indx,
                                                                specStr, 
                                                                quantity,
                                                                transitions[i])
            #print Avs, vs
            #plots.append(ax.semilogy(Avs, vs,'--o')[0])
            plots.append(ax.plot(Avs, vs,'--o')[0])
            titles.append(specStr + titleStrs[i])
        """
        #-------------done extracting and plotting stuff from the radex Dbs-------
        

        #pylab.legend(plots, titles)
        pylab.xlabel('Av', size='large')
        pylab.ylabel('cooling', size='large')
        pylab.title(plotTitle)
        pylab.xlim(xrng)
        pylab.ylim(yrng)
        pylab.xscale(xscale)
        pylab.yscale(yscale)
        
        pylab.legend(plots, titles, 
                    bbox_to_anchor = (1.3, 0.1, 0.1, 1))
        
                
        pylab.draw()
        pylab.show()
        self.logger.debug('exitting auxiliaray method')
        
    def get_quantity_from_radex_meshes_vs_Av(self, mesh_indx, specStr, quantity, transition_idx):
        """returns a quantity vs Av from a certain mesh indicated by 'mesh_indx'
        from readeDbs. A tuple is reutured (Avs, quantity) where Avs is the
        visual extinction at which the quantities returned. Both are numpy arrays.
        This was primarily written to get the flux as a function of Av from the 
        radex meshes. """
        
        #making a copy of the current Db to be reset to it upon exitting this method
        currRadexDb = self.currentRadexDb.copy() 

        Avs = []
        vs = [] 
        
        for AvKey in self.radexDbs.keys():
            Av = float(AvKey)
            #setting the radexDb to use
            self.use_radexDb(Av=Av, specStr=specStr)
            #getting the value we need and appending it
            if self.meshesRadex[mesh_indx] != None:
                vs.append(self.meshesRadex[mesh_indx][quantity][transition_idx])
                Avs.append(Av)
            print 'gMech radex = ', mesh_indx, numpy.log10(self.infoAllRadex[mesh_indx]['parms'][2])
        
        Avs = numpy.array(Avs)    
        vs = numpy.array(vs)

        self.use_radexDb(Av=currRadexDb['Av'], specStr=currRadexDb['specStr'])
        
        indsSorted = Avs.argsort()
        return (Avs[indsSorted], vs[indsSorted])  
    
    def compute_rxn_rates_at_Av(self, Av):
        """computes the reation rates at the closest Av of the available Av in the
        pdr mesh object self.mshTmp. The rates are set to self.mshTmp.chemNet 
        """
        
        #getting the index along the mesh which is closest in Av to the place
        #where the mouse was clicked, the closes values will be used in the
        #chemical netowrk of the mesh
        mshAv = self.mshTmp.data['state']['Av']
        indAv = np.argmin( np.fabs(mshAv - Av) )
        
        #the values to be assigned to the chemical network
        Av_use    = self.mshTmp.data['state']['Av'][indAv]
        gasT_use  = self.mshTmp.data['state']['gasT'][indAv]
        dustT_use = self.mshTmp.data['state']['dustT'][indAv]
        beta_CO   = self.mshTmp.data['selfSheilding']['CO'][indAv]
        beta_13CO = self.mshTmp.data['selfSheilding']['13CO'][indAv]
        beta_H2   = self.mshTmp.data['selfSheilding']['H2'][indAv]
        abun_use  = self.mshTmp.data['state']['abun'][:,indAv]

        #setting the parameters of the gas and the environment at that slab
        self.chemNet.set_environment_state(T           = gasT_use, 
                                           zeta        = self.parms_used_in_PDR_models['zeta'], 
                                           Av          = Av_use,
                                           albedo      = self.parms_used_in_PDR_models['albedo'],
                                           nDens       = self.mshTmp.data['hdr']['nGas'],
                                           G0          = self.mshTmp.data['hdr']['G0'],
                                           Tdust       = dustT_use,
                                           metallicity = self.metallicity,
                                           PHI_PAH     = 0.5,
                                           beta_CO     = beta_CO,
                                           beta_13CO   = beta_13CO,
                                           beta_H2     = beta_H2)
                        
        #cleaning previously computerd rxn rates and constants (if they were computed)
        self.chemNet.set_all_rxn_rates_and_cst_to_none()
         
        #set the abundances at that slab
        self.chemNet.set_abundances(fromArray = abun_use) 
        self.chemNet.compute_rxn_constants()
        self.chemNet.compute_rxn_rates() 
        
    def plot_rxns_rates_vs_Av(self, rxnIDs = None, Avs = None, ylim = None):
        """plots the reaction rates for the reactions whose IDs are passed in the list
        rxnIDs. The rates are computes for the values of Av passed in the list Avs.
        
        example :
         
            arxv.plot_rxns_rates_vs_Av([1559, 2100], numpy.arange(0,30,5), ylim=[1e-24, 1e-5])
        """
        
        if hasattr(rxnIDs, '__iter__') == False:
            raise TypeError('rxnIDs must be a list')
        if hasattr(Avs, '__iter__') == False:
            raise TypeError('rxnIDs must be a list')
        
        nRxns = len(rxnIDs)
        nAv   = len(Avs)
        
        #array which will store the rates of the reactions with rxnIDs
        rates = numpy.zeros((nRxns, nAv), dtype=numpy.float64)
        
        for j, Av in enumerate(Avs):
            
            #for each Av computing the rates
            self.compute_rxn_rates_at_Av(Av)
            
            #extractin the rates
            for i, rxnID in enumerate(rxnIDs):
                rate = self.chemNet.getattr_rxn(rxnID, 'rate')
                if rate == None:       
                    rates[i, j] = -1
                else:
                    rates[i, j] = rate

        pylab.figure(figsize=(5,5))
        plots = []
        titles = []
        
        #plotting the reaction rates
        for i, rxnID in enumerate(rxnIDs):
            rxnObj = self.chemNet.get_reactions(IDs = rxnID)[0]
            rxnStr = rxnObj.get_compact_rxn_string()
            
            plt, = pylab.semilogy(Avs, rates[i,:])
            plots.append(plt)
            titles.append(rxnStr)

        pylab.legend(plots, titles)
        
        if ylim != None:
            pylab.ylim(ylim)
        pylab.show()
    
    def plot_single_grid(self):
        
        #z_sec = -10.0
        z_sec = numpy.log10(1e-10)
        
        ## getting all the emission for a certain line at a certain Av
        v = self.get_emissions_from_databases(
                                              line={
                                                     'type'    : 'radex-lvg', 
                                                     'code'    : 'CO1-0',
                                                     'em_unit' : 'cgs',
                                                   }, 
                                              Av_use=10.0
                                             )
        
        ## the coordinates of the modesl in the DB (log10(n), log10(G0), log10(gmech))
        x, y, z = self.grid_x, self.grid_y, self.grid_z
        
        ## keeping only a section in z
        inds = numpy.where( (z >= (z_sec - 1e-6))*(z <= (z_sec + 1e-6)) )
        
        x_in, y_in, v_in = x[inds], y[inds], v[inds]
        
        ## constructing the interpolation function
        f_interp = interpolate.LinearNDInterpolator(np.array([x_in, y_in]).T, v_in)
        
        xu, yu = numpy.unique(x_in), numpy.unique(y_in)
        
        xx, yy = numpy.meshgrid(xu, yu, indexing='ij')
        
        data_i = np.array( [xx.flatten(), yy.flatten()] ).T
        
        vi = f_interp(data_i)
        v_xy = vi.reshape(xx.shape).T
        pylab.figure()
        pylab.imshow(v_xy, origin='lower', vmin=-15.0, vmax=-6.0)
        
        inds_nan = numpy.where(numpy.isfinite(vi) == False)[0]
        if inds_nan.size > 0:
            vi[inds_nan] = numpy.nanmin(v_in)
        
        
        ############
        #vi = vi.reshape(xi.shape)
        #f_cubic = interpolate.RectBivariateSpline(numpy.unique(xi), numpy.unique(yi), vi)
        
        #f_cubic = interpolate.SmoothBivariateSpline(xi.flatten(), yi.flatten(), vi.flatten(), w=None,
        #                                            kx=3, ky=3, s=0.05, eps=1e-2)
        
        f_cubic = interpolate.LSQBivariateSpline(xx.flatten(), yy.flatten(), vi, [0,1],[0,1], w=None,
                                                 kx=3, ky=3, eps=1e-3)
        ###########

        #xic, yic = numpy.linspace(0.0, 6.0, 100), numpy.linspace(0.0, 6.0, 100), 
        xic, yic = numpy.mgrid[0.0:6.0:100j, 0.0:6.0:100j] 
        
        #if using RectBivariateSpline
        vic = f_cubic.ev(xic.flatten(), yic.flatten())   
        vic = vic.reshape(xic.shape)
        pylab.figure()
        pylab.imshow(vic.T, origin='lower', vmin=-15.0, vmax=-6.0)
        
        #asdad
    
    def plot_grid_n_G0_standalone(self, maps2d_idx, title = None, levels = None):
        """Plots a grid in a new standalone figure instead of the main gui. One of the maps2d index
        (tuple) (for (0,0), or (0,1), (1,0), (1,1) ) ... sould be passed."""
        
        i, j = maps2d_idx
        fig = pylab.figure(figsize=(6,6))
        
        grid_info = self.gui['maps2d']['%s%s' % (i,j)]
        """

        maps2d_10['axesCbar'].set_title('N(%s)' % self.parms['gridsInfo']['10']['specStr'])
        maps2d_10['contour'] = None
        maps2d_10['axes'].set_xlabel( '$log_{10} n_{gas}$' )
        maps2d_10['axes'].set_ylabel( '$log_{10} G_0$' )
        maps2d_10['axes'].set_xlim(self.ranges[0][0], self.ranges[0][1])
        maps2d_10['axes'].set_ylim(self.ranges[1][0], self.ranges[1][1])
        (maps2d_10['axes'].yaxis.get_major_ticks())[-1].label1On = False
        (maps2d_10['axes'].xaxis.get_major_ticks())[-1].label1On = False
        gui['maps2d']['10'] = maps2d_10 
        """
    
        left = 0.14
        bottom = 0.12
        sz = 0.7
        vSpace = 0.07
        nlevels = 12
        
        ax = fig.add_axes([left, bottom, sz, sz])
        ax_cbar = fig.add_axes([left, bottom + sz + vSpace, sz, 0.02])

        grd_data = self.grds[i][j]
        ranges = (grid_info['axes'].get_xlim()[0], grid_info['axes'].get_xlim()[1], grid_info['axes'].get_ylim()[0], grid_info['axes'].get_ylim()[1])
        
        ax.set_xlim(ranges[0],ranges[1])
        ax.set_ylim(ranges[2],ranges[3])
        
        im = ax.imshow(grd_data, extent=ranges, origin='lower')
        
        if levels == None:
            dl = (np.nanmax(grd_data) - np.nanmin(grd_data))/nlevels
            levels = np.arange( np.nanmin(grd_data), np.nanmax(grd_data), dl )
        
        contours = ax.contour(grd_data, levels, extent=ranges, origin='lower', colors = 'black')

        
        ax.clabel(contours, levels, fmt = '%.1f' )
        
        cbar = pyl.colorbar(im, cax = ax_cbar, ax = ax, 
                            orientation = 'horizontal', ticks = levels[::2], format = '%.1f')
        
        ax.set_xlabel(r'log$_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='xx-large')
        ax.set_ylabel(r'log$_{10}$ G$_0$', size='xx-large')

        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.tick_params(axis='both', which='minor', labelsize=15)

        ax_cbar.tick_params(axis='both', which='major', labelsize=15)
        ax_cbar.tick_params(axis='both', which='minor', labelsize=15)

        vmax, vmin = np.nanmax(grd_data), np.nanmin(grd_data)  
        nTicks = 5
        cbar_levels = numpy.linspace(vmin, vmax, nTicks)
        cbar.set_ticks(cbar_levels)
        
        title = r'log$_{10}$ [ CO(1-0)' + r'/ erg cm$^{-3}$ s$^{-1}$ ]'                     
        #title = r'log$_{10}$ [ HCN(1-0)' + r'/ erg cm$^{-3}$ s$^{-1}$ ]'       
        ax_cbar.set_title(title, size = 'xx-large')
        
        pylab.show()
        
        
    def print_available_radex_dbs(self):
        """prints the keys and the nested keys of the attribute self.radexDbs"""
        
        for key_Av in self.radexDbs:
            print 'Av = %-5s :' % key_Av,
            for key_specStr in self.radexDbs[key_Av]:
                print '%-8s ' % key_specStr,
            print 
        
    def get_emissions_from_databases(self, line=None, Av_use=None):
        '''gets the emissions from the pdr database or radex models determined by the keyword line
        
        .. code-block:: python
        
            CO1_0 = arxv.get_emissions_from_databases(
                                                      line   = {'type':'pdr', 'code':'CO1-0'}
                                                      Av_use = 10.0,
                                                     )
        '''
    
        lineInfo = lineDict.lines[line['code']]
        
        if line['type'] == 'pdr' or line['type'] == 'radex-lvg':
            pass
        else:
            raise ValueError('unknow line type, please specify "pdr" or "radex-lvg"')
         
        #setting the pointer to the quantity which will be used to access the emission data 
        #depending on the line type requested
        #---------------------------------------------------------------------------------------------------
        if line['type'] == 'pdr':

            v = self.apply_function_to_all_meshes(
                                                  pdr_mesh_log_intensity, 
                                                  func_kw = {
                                                             'quantity': lineInfo['ismcpak'],
                                                             'up_to_Av': Av_use,
                                                            }
                                                  )
        #---------------------------------------------------------------------------------------------------        
        if line['type'] == 'radex-lvg':

            self.use_radexDb(Av=Av_use, specStr=lineInfo['specStr'], load_if_not_in_memory=True)

            func_kw = {}
            func_kw['transitionIdx'] = lineInfo['radexIdx']
            if 'em_unit' in line:
                func_kw['em_unit'] = line['em_unit'] 
            
            v = self.apply_function_to_all_radex_meshes(
                                                        radex_mesh_log_intensity, 
                                                        func_kw = func_kw,
                                                       )
        return  np.array(v, 'f8') 

    def get_emission_from_all_radex_dbs_for_Av_range(self, line=None, Avs=None, quantity=None, keep_nans=None):
        '''gets the emissions from the pdr database or radex models determined by the keyword line for
         the specified Avs in the Avs array. line should be an entry for a radex line in lineDict.lines
        
        :param Avs: array specifiying the Avs to use, or 'all' indicating the 
         use of all available Avs. 
        
        .. code-block:: python
        
            v, data = arxv.get_emission_from_all_databases_for_Av_range(
                                                                        line = '13CO1-0',
                                                                        Avs  = [1.0, 2.0, 5.0, 10.0],
                                                                        quantity = 'fluxKkms',  
                                                                       )
        The available Avs can be queried using :    
        
            self.print_available_radex_dbs()  # print available DB on disk for all species 
            self.query_available_radex_dbs()  # pring the available Avs for the current DB
                              
        see also get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range()
        '''
        
        specStr   = lineDict.lines[line]['specStr']
        radex_idx = lineDict.lines[line]['radexIdx']
        
        if type(Avs) == type(''):
            assert(Avs == 'all')
            Avs_use = self.query_available_radex_dbs()
        else:
            Avs_use = numpy.array(Avs)
        
        ## collecting the data from the database corresponding to all the data in Av avaiable
        for i, Av in enumerate(Avs_use):
            
            #using the database of the specific Av
            self.use_radexDb(specStr=specStr, Av=Av)
            
            #getting the data corresponding to this Av
            v = self.apply_function_to_all_radex_meshes(
                                                        radex_mesh_quantity, 
                                                        func_kw = {
                                                                   'transitionIdx':radex_idx,
                                                                   'quantity' : quantity
                                                                   }
                                                        )

            v = numpy.array(v)
                        
            if keep_nans != None and keep_nans == True: 
                xGrd, yGrd, zGrd = self.grid_x, self.grid_y, self.grid_z
            else:
                #keeping the points which are useful (finite ones)
                inds_valid = numpy.isfinite(v)
                v = v[inds_valid]
                xGrd, yGrd, zGrd = self.grid_x[inds_valid], self.grid_y[inds_valid], self.grid_z[inds_valid]

            AvGrd = numpy.ones(xGrd.shape, 'f')*Av
            data = numpy.array([xGrd, yGrd, zGrd, AvGrd], dtype = numpy.float64).T
                    
            #soting the x,y,z,t, and v into arrays to be used later to construct the interpoaltion function
            if i == 0:
                data_all_Av = data
                v_all_Av = v
            else:
                data_all_Av = numpy.vstack( (data_all_Av, data) )
                v_all_Av = numpy.hstack( (v_all_Av, v) )

        return v_all_Av, data_all_Av

    def get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range(self,
                                                                    sectioned=False,                                                                     
                                                                    **kwargs):
        '''Retuns an interpolation function of the emission from radex DBs.  This method takes the same
         parmaeters as get_emission_from_all_radex_dbs_for_Av_range and takes the same keywords.  The 
         returned interpolation function returns the emission info give logn, logG0, logGmech and Av.  
        
        .. code-block:: python
        
            F = arxv.get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range(
                                                                                line = '13CO1-0',
                                                                                Avs  = [1.0, 2.0, 5.0, 10.0],
                                                                                quantity = 'fluxKkms',  
                                                                                )

            # the quantity for a given set of parameters can be returned
            #        logn  logG0 logGM   Av  
            parms = [ 0.0, 0.0  , -50.0, 9.0]
            q = F(array(parms).reshape(1,4))

            #or for a vector of inputs using
            logn  = array([0.0  ,   1.0,   2.0,   3.0,   4.0])
            logG0 = array([3.0  ,   3.1,   3.2,   3.3,   3.4])
            logGM = array([-30.0, -29.0, -28.0, -27.0, -26.0])
            Avs   = array([1.0  ,   2.0,   3.0,   4.0,   5.0])
            qv = F.get(vstack((logn, logG0, logGM, Avs)).T)
            
            #constructing the interpolation might take about 10 seconds depending on the amount of 
            #points in the DBs, so the constructed function can be saved into a file and loaded later 
            #using
            save('/home/mher/foo.npy', F)
            #or
            numpy.savez_compressed('/home/mher/tmp/foo.npz', F)

            F = load('/home/mher/foo.npy')
             
        see also get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range
        '''
        ## getting the data
        v, data = self.get_emission_from_all_radex_dbs_for_Av_range(**kwargs)

        ## constructing the linear interpolation function        
        if sectioned == False:
            F = interpolate.LinearNDInterpolator(data, v)
            setattr(F, 'get', F) 
        else:
            F = interpolator_sectioned(data, v, ## not tested
                                       intervals_z=self.intervals_z,     
                                       ghost_z=self.ghost_z,
                                       intervals_t=self.intervals_t,
                                       ghost_t=self.ghost_t,
                                       scipy_interpolator=interpolate.LinearNDInterpolator)

        return F
     
    def get_4D_interp_quantity(self, info=None, save=None, **kwargs):
        ''' returns an interpolation function for different quantities defined by the info
        keyword.  
        
        :param dict info: a dict holding the neccessary info which is used to extract data from
         the meshes or the DBs to construct the interpolation function.
          
        .. code-block:: python
        
            # uses get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range to get the data
            info = {'source' : 'radex' }
             
        .. code-block:: python
        F = self.get_4D_interp_quantity(
                                        info={'source':'radex', 
                                        save=False,
                                        line='13CO1-0',
                                        Avs='all',
                                        quantity='fluxKkms'
                                        sectioned=False,
                                       )
        '''
        
        if info['source'] == 'radex':
            F = self.get_4D_interp_emission_func_from_all_radex_dbs_for_Av_range(**kwargs)
            

        ## saving the interpolation function to the disk            
        if save == True:

            fpath = self.get_interp_func_path(info, **kwargs)
                        
            ## if the object is of type scipy interpolator....
            '''
            numpy.savez_compressed(fpath, 'not_sectioned', F)
            print 'interpolation function saved to \n\t%s:' % fpath
            '''
            
            ## if the object is of type interpolator_sectioned
            if type(F) == interpolator_sectioned:
                F.save(fpath)
                        
        return F
       
    def load_interp_func(self, info=None, **kwargs):
        '''reads a the saved interpolation function and returns it
        
        .. code-block:: python

            arxvPDR.load_interp_func(info={'source':'radex'}, 
                                     line='13CO1-0', 
                                     quantity='fluxKkms')

        '''
        
        fpath = self.get_interp_func_path(info, **kwargs) + '.npz'
        
        fd = numpy.load(fpath)
        
        method = fd['method'].flatten()[0]
 
        if method == 'not_sectioned':
            F =  fd['data'].flatten()[0]
            
        if method == 'sectioned':
            F = interpolator_sectioned(None, None, load_from = fpath)
 
        return F
    
    def get_interp_func_path(self, info=None, **kwargs):
        '''sets up the path of the interpolation function based on the 
        quantities it provides'''
        
        if info['source'] == 'radex':
            
            fpath =  os.path.join(self.dirPath, 'radexDbs', 'interp', 
                                    '%s__%s' % (kwargs['line'],  
                                                kwargs['quantity'])
                                 )
        
         
        return fpath
    
    def get_emission_grid_from_databases(self, 
                                         line=None, Av_use=None, ranges=None, res=None, z_sec=None, dz=None,
                                         f_interp_dim=None, interp='linear', zoom=5, 
                                         get_section_data_using='select',
                                         **kwargs):
        '''returns an interpolated grid from the database for a pdr or radex models determined by the 
         keyword line.
         
        :param string get_section_data_using: 'select' or '3d-interpolation'. gets the values section in z either by 
         selecting poitns in the 3D grid and interpolating over that 2D surfrace. Or gets the values at the 
         section first by interpolating at that section using 3D interpolation. Once the data in that z section
         is available, a 2D interpolation is done over that data to get an interpolated values at that z section. 
         
        '''

        ## the coordinates of the models for ex: log10(n), log10(G0), log10(gmech)
        x, y, z = self.grid_x, self.grid_y, self.grid_z

        ## the emission extracted from the database for the specified line for a certain AV
        v = self.get_emissions_from_databases(line, Av_use)

        def get_coords_inds_in_sec():
            '''returns the indicies of x,y coordinates of the available data withing a section in z '''
            
            #selecting a section in z
            if dz != None:
                inds = numpy.where( (z >= (z_sec - dz))*(z <= (z_sec + dz)) )
            else:
                inds = numpy.where(z == z_sec)[0]
            
            return inds
        #
        
        def get_z_section_within_z_range():
            '''returns the x,y,v of valid points in a z section in databse within a z section''' 
        
            inds = get_coords_inds_in_sec() 
                
            #keeping only the points withing the z section
            x_in, y_in, v_in = x[inds], y[inds], v[inds]
            
            return x_in, y_in, v_in
        #

        def get_z_section_by_3D_interpolation():
            '''constructing the 3D interpolation function which will be used to interpolate at the 
               grid points of this section.  Only the points which do not return an nan are kept.
            '''
            data=np.array([x, y, z], 'f8').T
            fInterp = self.construct3DInterpolationFunction(data=data, values=v, remove_nan_inf=True)

            inds = get_coords_inds_in_sec() 

            ## grids in x and y in the z section with the same resolution as the original points
            xxGrdc, yyGrdc = numpy.meshgrid(numpy.unique(x[inds]), numpy.unique(y[inds]), indexing='ij')
            
            v_new = fInterp(
                           np.array(
                                    [
                                     xxGrdc.flatten(), 
                                     yyGrdc.flatten(), 
                                     numpy.ones(xxGrdc.size, 'f8')*z_sec
                                    ]
                                   ).T
                          )
            
            return xxGrdc.flatten(), yyGrdc.flatten(), v_new 
        #
        
        if f_interp_dim != None and f_interp_dim == '2D':
 
            ## getting a section of values in z by keeping only the points withing the z section
            if get_section_data_using == 'select':
                x_sec, y_sec, v_sec = get_z_section_within_z_range()
            elif get_section_data_using == '3d-interp':
                x_sec, y_sec, v_sec = get_z_section_by_3D_interpolation()
            else:
                raise ValueError("invalid z selection method, use 'select' or '3d-interp'")
            
            if interp in ['cubic', 'linear', 'nearest']:

                # getting the coordinates of the grid in each dime
                xx, yy = numpy.unique(x_sec), numpy.unique(y_sec)
                xxc = numpy.linspace(xx.min(), xx.max(), res[0])
                yyc = numpy.linspace(yy.min(), yy.max(), res[1])
                xxGrdc, yyGrdc = numpy.meshgrid(xxc, yyc, indexing='ij')

                # removing all nan's from the data
                inds_not_nan = numpy.isfinite(v_sec)
                x_sec_finite, y_sec_finite, v_sec_finite = x_sec[inds_not_nan], y_sec[inds_not_nan], v_sec[inds_not_nan] 

                '''
                x_sec_finite, y_sec_finite, v_sec_finite = x_new, y_new, v_new 
                '''
                
                # interpolating
                grd = interpolate.griddata(
                                           numpy.vstack((x_sec_finite, y_sec_finite)).T, v_sec_finite,
                                           (xxGrdc, yyGrdc), 
                                           method=interp, 
                                           fill_value=numpy.nan
                                           )
                
                grd_x, grd_y = xxGrdc, yyGrdc
            

            
        else:
            ## gettomg the data interpolated from a 3D interpolator
            data=np.array([self.grid_x, self.grid_y, self.grid_z], 'f8').T
            
            fInterp = self.construct3DInterpolationFunction(data=data, values=v, remove_nan_inf=True)  #interpolator='nearest',
            grd, grd_x, grd_y = self.computeInterpolated2DGrid(ranges=ranges, res=res, zSec=z_sec, fInterp = fInterp)            
    
        return grd

    def as_dict(self):
        '''returns a dict of a selection of the attributes. Mainly the ones which can be pickeled
         and serialized. This method was written so that the components of this object can be
         pushed to remote machines/hosts/processes using IPython's parallel capabilitie'''
        
        #the attributes to be returned as keys of a dict
        ret_attr =  ('dirPath', 'nMeshes', 'ver', 'meshes', 'infoAll', 'verRadex', 'currentRadexDb',)
        ret_attr += ('infoAllRadex', 'meshesRadex', 'radexDbs', 'grid_qx', 'grid_x', 'grid_qy',)
        ret_attr += ('grid_qy', 'grid_y', 'grid_qz', 'grid_z', 'grid_x_unique', 'grid_y_unique',)
        ret_attr += ('grid_z_unique', 'nDbFiles', 'chemNet', 'ranges', 'parms', 'metallicity',)
        ret_attr += ('parms_used_in_PDR_models',)

        #the dict to be returned
        dict_ret = dict()
        
        #setting the keys of the dict to be returned
        for attr in ret_attr:
            dict_ret[attr] = getattr(self, attr) 
             
        return dict_ret
    
    def make_copy_from_data(self, data_as_dict):
        '''returns a new :data:`meshUtils.meshArxv` object and sets the attributes from the data
         passed. The data should be a dict returned by the method :data:`meshUtils.as_dict`
         '''
        
        arxv_ret = meshArxv(no_init=True)
        
        for key in data_as_dict:
            setattr(arxv_ret, key, data_as_dict[key])
            
        return arxv_ret
    
######################################################################################################################
#                utility functions               utility functions               utility functions
######################################################################################################################
def radex_mesh_log_intensity(mesh_radex, **kwargs):
    '''returns the log10 intentisty emitted from a radex model applied on a pdr mesh with the Av of the current Db'''
    
    transitionIdx = kwargs['transitionIdx']
    if 'em_unit' in kwargs:        
        quantity = 'flux' + kwargs['em_unit']
    else:
        quantity = 'fluxcgs'         
        
    
    if mesh_radex == None:
        return numpy.nan
    else:
        return numpy.log10(mesh_radex[transitionIdx][quantity])

def pdr_mesh_log_intensity(meshObj, **kwargs):
    '''returns the log10 intentisty (per sr) emitting from a pdr mesh with a certain Av
    
    .. todo:: to be deprecated
    '''
    
    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = meshObj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    
    return numpy.log10(value)


def radex_mesh_quantity(mesh_radex, **kwargs):
    '''This is a utility function which returns a quantity from a radex mesh. The log of the flux is returned, for the rest
     of the quantites, the quantity as is is returned. If the mesh_radex is None, nan is returned. This function can be passed
     to  meshArxv.apply_function_to_all_radex_meshes(radex_mesh_quantity, func_kw = kwargs)
     where kwargs are transitionIdx = integer (the transition index) and quantity = string (the quantity in the gtype to be returned)  
    '''
    
    if mesh_radex == None:
        return numpy.nan
    else:
        transition_info = mesh_radex[kwargs['transitionIdx']] 
        quantity = transition_info[kwargs['quantity']]  #a radex quantity
        
        if kwargs['quantity'] == 'fluxcgs' or kwargs['quantity'] == 'fluxKkms':
            return numpy.log10(quantity) 
        else:
            return quantity    
        
def pdr_mesh_integrated_quantity(mesh_obj, **kwargs):
    '''This is a utility function which returns an intgerated quantity from a mesh up to a certain Av. This function can be passed
     to  meshArxv.apply_function_to_all_meshes(pdr_mesh_integrated_quantity, kwargs) where kwargs are the 'quantity' in a :data:`mesh.mesh`
     object and 'up_to_Av', which is the Av upto which the integration will be done. Also the keyword 'as_log10' can be passed
     so that the log10 of the quantity is retured (this keyword is optional).
     mesh_obj should be an object of type :data:`mesh.mesh`.
     
     This is not so useful in computing column densities since there is not one quantity that points to the abundance of 
     species, since they are stored in a 2D array. For that purpose use, pdr_mesh_column_density()
     
    .. todo:: this should be used instead of pdr_mesh_log_intensity. 
    '''

    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = mesh_obj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    
    if 'as_log10' in kwargs and kwargs['as_log10'] == True:
        return numpy.log10(value)
    else:
        return value
    
def pdr_mesh_column_density(mesh_obj, **kwargs):
    '''Same as pdr_mesh_integrated_quantity, but returns a column density of a certain species.
    '''

    specStr  = kwargs['specStr']
    up_to_Av = kwargs['up_to_Av']

    value = mesh_obj.getColumnDensity(specsStrs = [specStr], maxAv = up_to_Av)

    
    if len(value) == 1:
        value = value[0]
        
    if 'as_log10' in kwargs and kwargs['as_log10'] == True:
        return numpy.log10(value)
    else:
        return value
