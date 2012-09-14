import os, glob, sys
import numpy as np
import pylab as pyl
from matplotlib.widgets import Button
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d import Axes3D

from time     import *
from mesh     import *
from utils    import *
from ndmesh   import *
from ismUtils import *
from radex    import *
from utils    import fetchNestedDtypeValue
from scipy    import interpolate
import chemicalNetwork

class meshArxv():
    """ this class generates and manipulates archives of PDR meshes.
               
     FILES AND THEIR FORMATS:\n
     by default, the prefix name of the individual mesh files is assumed to be
     mesh.dat-id-xxxxxx
     
     the mesh database files are assumed to have the same prefix, for example, if the
     database file provided is foo, then this routine will assume (or write)  
         foo.info, foo.db
     for now the database can be stored into a single file and not split into
     multiple files.  
     
     the database for the meshes is constructed of two binary files:
         foo.info : holding all the information about the meshes, locations and parameters
         foo.db   : holds the data for all the meshes
       
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

      for a mesh number i, the checkNum_i should be the same as the i^th entry 
      in the info array offset...i.e chechNum = infoAll[i]['info'][2] 
        
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
        """ a list of all the meshes of 'mesh' dtypes (see mesh.py)
        
            for a mesh of at index 'i' in 
               self.meshes[i]
            the corresponding info in the header are accessed as follows :
            
            self.infoAll[i]['parms'][0]) which should be the same as self.meshes[i]['hdr']['G0'] 
            self.infoAll[i]['parms'][1]) which should be the same as self.meshes[i]['hdr']['nGas']
            self.infoAll[i]['parms'][2]) which should be the same as self.meshes[i]['hdr']['gammaMech']
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

        #---------------------------radex data storage attributes0000-----------------------------------
        self.verRadex = None
        """a 3 element np.int32 array holding the version of the radex database currently in use."""
        
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
        of dtype radex.transitionDtype. 
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

        self.nDbFiles   = None
        """ numpy.int32 number of database files"""
        
        self.chemNet    = None
        """ object of class type chemicalNetwork holds info about the chemical network used
        """
    
        # variables used for plotting
        self.pltGmSec = None #: the value of the section in Gmech selected
        self.pltGrds  = None 
        self.grds     = None #: 2x2 ndmesh array object
        self.grdsCbarAxs = None
        self.pltRadex    = None
        self.fig         = None
        self.grdPltPts1  = None
        self.grdPltPts2  = None
        self.grdPltTitle = None
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

    # read all the meshes files in the dir and construct the
    # database
    def construct(self, meshNamePrefix = None, writeDb = None ):
        """ construc the database anad write the .db files. If the meshNamePrefix is
              not supplied all the files in self.dirPath are assumed to be data files
              and all of them are put in the database."""
        
        if meshNamePrefix == None:
            meshNamePrefix = ''
        
        # getting the names of the meshes in that dir
        files = []
        for infile in glob.glob( os.path.join(self.dirPath, 'meshes/'+meshNamePrefix + '*') ):
            files.append(infile)
        
        # setting variable for the 
        self.nMeshes = np.zeros( 1, dtype = np.int32 )
        self.nMeshes[0] = len(files) 
        print 'found %s meshes' % (self.nMeshes)

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
                print 'proccessed %d meshes...' % i 
            
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
        
        print 'wrote successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)        
        
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
                    
        print 'read successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)
        
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

            np.array( self.infoAllRadex[i]['info'][2] ).tofile(dbDataFObj)

        dbDataFObj.close()
        
        # writing the db info into a file
        dbInfoFObj = file(self.parms['dirPath'] + 'meshesRadex.db.info.' + self.parms['radex']['specStr'], 'wb')
        self.verRadex.tofile( dbInfoFObj)
        self.nMeshes.tofile( dbInfoFObj )
        self.infoAllRadex.tofile( dbInfoFObj )
        dbInfoFObj.close()
        print 'wrote successfully the radex database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)

    def readDbRadex(self, specStr, check = None):
        """ reads the database sufficed by specStr (i.e meshesRadex.db.(specStr)and assigns the appropritate attributes (document)
            
            :param bool check: if this is set (to any value) the self.checkIntegrity() is called.
            :param string specStr: the string of the specie whose database is to be read.
            :warning: the keyword specStr is not functional yet.
            :note: before calling this mehtod, an instance of the radex class should be     
              assigned to self.radexObj. 
        """ 
        
        dbInfoFObj = file(self.dirPath + 'meshesRadex.db.info.' + self.parms['radex']['specStr'], 'rb')
        dbDataFObj = file(self.dirPath + 'meshesRadex.db.' + self.parms['radex']['specStr'], 'rb')

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
                strng  = 'Error : checkpoint numbers do not match : database may be corrupt.\n'
                strng += 'read = %d, expected = %d' % (checkNum, infoAllRadex[i]['info'][2])
                strng += ' when reading mesh %d' % i                
                raise NameError(strng)

            #self.meshes.append( thisMeshRadexData[0] )                
            meshesRadex.append( thisMeshRadexData )                
        print 'read successfully radex database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)

        self.infoAllRadex = infoAllRadex
        self.meshesRadex = meshesRadex
        
        if check:
            self.checkIntegrityRadex()
            
            
    
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
            print 'archive integrity test passed'
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
            print 'archive integrity test passed'
        else:
            strng = 'archive integrity test failed. database file may be corrupt' 
            raise NameError(strng)
    
    def setupChemistry(self):
        """sets up the attributes related to the chemistry assuming self.parms['chemistry'] is set.
           sets the attribute self.chemNet.
        """
        
        #importing the module which holds the definitions of the base species
        baseSpecies = __import__(self.parms['chemistry']['baseSpecies'])
        baseSpecs = baseSpecies.baseSpecies()
        
        # settin up the orignial netowrk
        net = chemicalNetwork.chemicalNetwork(self.parms['chemistry']['rxnFile'], 
                                              baseSpecs,
                                              UMISTVER = self.parms['chemistry']['umistVer'])
        # reading the species to be removed from a file
        net.removeSpecies( underAbunFile = self.parms['chemistry']['underAbunFile'] )
        net.removeSpecies( species = self.parms['chemistry']['removeManual'] )
        # reading the species number and their corresponding indies and abundances from ascii files
        net.assignNumbersToSpecies(fileName = self.parms['chemistry']['specNumFile'])
        
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
                 
             :param int32 slabIdx: the index of the slab from which the value will be
                 extracted
                 
             :param int32 arrIdx: The first index in the 2D array (in case the quantity points to
              a 2D array) to be extracted first (usually for ['state','abun']). And slabIdx is used
              as the seconds index. 
                
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
        print 'constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti)

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
        print 'constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti)

        return f 

    def compute3DInterpolatedData(self, quantity = None, slabIdx = None, x = None, y = None, z = None, 
                                  fInterp = None, *args, **kwargs):
        """this method is the same as computeInterpolated2DGrid but it returns interpolated values based on
           the input values, i.e """
        pass
                                 
    
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
             
             :param interpolate.LinearNDInterpolator fInterp: the interpolation function returned by
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
        print 'interpolated %d points in %f seconds at a rate of %e pts/sec' % (nPts, tf-ti, nPts / (tf-ti))
        tNew = np.reshape(tNew, grid_x.shape)
        
        return tNew
        
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
        grd = self.computeInterpolated2DGrid(ranges   = ranges,
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
        #pyl.set_xlable('log_10 nGas')
        #pyl.set_ylable('log_10 G0')
        
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
            #getting the column densities for all the models
            m = mesh(chemNet = self.chemNet, metallicity = self.metallicity)
            values  = np.ndarray(self.nMeshes, dtype = np.float64)
            grid_x = []
            grid_y = []
            grid_z = []
             
            for i in np.arange(self.nMeshes):
                if self.meshesRadex[i] != None:
                    transitionIdx = self.parms['gridsInfo']['11']['transitionIndx']
                    quantity = self.parms['gridsInfo']['11']['quantity']
                    print self.meshesRadex[i][transitionIdx][quantity]
                    values[i] = self.meshesRadex[i][transitionIdx][quantity]
                    asdadad
                    ####get the self.grid_x, grid_y, grid_z corresonding to the 
                    ####meshes which have radex data and pass them to construct3DInterpolationFunction
                    grid_x.append(self.grid_x[i])
                    grid_y.append(self.grid_y[i])
                    grid_z.append(self.grid_z[i])
                    data = np.array([grid_x, grid_y, grid_z]).T # getting the data in the shape that is 
                                                                # accepted by the interpolation construction
                else:
                    values[i] = np.nan
                    
            self.intensityGridInterp_f = self.construct3DInterpolationFunction(data = data,
                                                                               values = values,
                                                                             log10  = True,
                                                                             *args, **kwargs)

    def showSurfaceTemperatureGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the surface temperature grid
            
           :todo: plot every other labeled contour as a tick in the colorbar
           
           :warning: there would be a memeory leak if things are plotted over and over
              the images, labels..etc..must be replaced and not new instance created
              and overplotted...take care of that later... 
        """
        grd = self.computeInterpolated2DGrid(ranges   = ranges,
                                             res      = res,  
                                             zSec     = self.pltGmSec, 
                                             fInterp  = self.grdInterp_f, *args, **kwargs)
        
        grd = grd.T
        im00 = pyl.imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )
        
        CS = pyl.contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        pyl.clabel(CS,levels, fmt = '%.1f' )
        
        pyl.colorbar(im00, cax=self.grdsCbarAxs[0][0], ax=pyl.gca(), orientation = 'horizontal')

    def showAbundancesGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the abundances grid"""            

        grd = self.computeInterpolated2DGrid(ranges   = ranges,
                                             res      = res,  
                                             zSec     = self.pltGmSec, 
                                             fInterp  = self.abunGridInterp_f, *args, **kwargs)
        
        grd = grd.T
        im01 = pyl.imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )

        CS = pyl.contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        pyl.clabel(CS,levels, fmt = '%.1f' )
        
        pyl.colorbar(im01, cax=self.grdsCbarAxs[0][1], ax=pyl.gca(), orientation = 'horizontal')

    def showColumnDensityGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the abundances grid"""            

        grd = self.computeInterpolated2DGrid(ranges   = ranges,
                                             res      = res,  
                                             zSec     = self.pltGmSec, 
                                             fInterp  = self.colDensGridInterp_f, *args, **kwargs)
        
        grd = grd.T
        im10 = pyl.imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )
        
        CS = pyl.contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        pyl.clabel(CS,levels, fmt = '%.1f' )
        
        pyl.colorbar(im10, cax=self.grdsCbarAxs[1][0], ax=pyl.gca(), orientation = 'horizontal')


    def showLineIntensityGrid(self, ranges = None, res = None, *args, **kwargs):
        """shows the line intensity grid"""            

        grd = self.computeInterpolated2DGrid(ranges   = ranges,
                                             res      = res,  
                                             zSec     = self.pltGmSec, 
                                             fInterp  = self.intensityGridInterp_f, *args, **kwargs)
                                                        
        grd = grd.T

        im11 = pyl.imshow(grd,extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower')
        nlevels = 10
        dl = (np.nanmax(grd) - np.nanmin(grd))/nlevels
        levels = np.arange( np.nanmin(grd), np.nanmax(grd), dl )

        CS = pyl.contour(grd, levels, extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), origin='lower', colors = 'black')
        pyl.clabel(CS,levels, fmt = '%.1f' )
        
        pyl.colorbar(im11, cax=self.grdsCbarAxs[1][1], ax=pyl.gca(), orientation = 'horizontal')
        
    
    def constructRadexDatabase(self, writeDb = None):
        """runs radex on all the models in self.meshes, and computes the line info according
        to the parameters in self.parms['radex']. Once done computing, it stores
        all the generated info into self.infoAllRadex and self.meshesRadex.
        
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
            
            print 'pdr mesh index =  %d : ' % i
            
            m.setData( self.meshes[i] )
            (gasTRadex, nColls, colDensThisSpec,) = m.getRadexParameters(speciesStr = specStr,
                                                                         threshold  = xH2_Min)
            print gasTRadex, nColls, colDensThisSpec
            
            # getting the collider densities in the same order of the supplied input spcie string list 
            nDensColls = [ nColls[collSpecStr] for collSpecStr in self.parms['radex']['collisionPartners'] ]
            collsStr   = list(self.parms['radex']['collisionPartners'])
            print 'input coll species', self.parms['radex']['collisionPartners'] 
            print 'nColls after putting them in the right order = ', nDensColls

            print 'radexGrid : radex parms : ', gasTRadex, nDensColls, colDensThisSpec
            radexObj.setInFileParm('tKin', gasTRadex)
            radexObj.setInFileParm('collisionPartners', collsStr )
            radexObj.setInFileParm('nDensCollisionPartners', nDensColls )
            radexObj.setInFileParm('molnDens', colDensThisSpec)

            #remove colliders which are underabundant (below radex limits)
            radexObj.filterColliders()

            if len(radexObj.inFile['collisionPartners']) == 0:
                print 'all colliders have densities outsie the range accepted by radex'
                continue
            
            #setting radex to the defaul status before running it
            radexObj.setDefaultStatus()
            
            #run radex for this model
            status = radexObj.run( checkInput = True )

            #setting the basic radex run info into the attribute
            #-----------------------------------------------------------
            #saving the number of transitions into radex info attribute
            infoAllRadex[i]['info'][0] = i
            if not radexObj.flagSet('ERROR'):
                infoAllRadex[i]['info'][1] = radexObj.nTransitions
                meshesRadex.append(radexObj.transitions)
            else:
                infoAllRadex[i]['info'][1] = 0 #no trainsitions
                meshesRadex.append(None)
                radexObj.printSetFlags()

            #infoAllRadex[i]['info'][2] = NOT SET HERE, IT IS SET WHEN WRITING THE DB TO A FILE
            
            #setting the status into the attribute
            infoAllRadex[i]['info'][3] = status
            #appending the transition data of this mesh to the attribute 
            #-----------------finished saving the info------------------
            
            
            if radexObj.flagSet('SUCCESS'):
                print 'radexGrid : converged with no warnings'
            else:
                print 'radexGrid : converged with warnings'
                print '------------------------------------'
                print radexObj.getWarnings()
                print '------------------------------------'
                continue

            print '---------------------------------------------------------'
            

        #copying the mesh parameters self.infoAll[:]['parms] to self.infoAllRadex[:]['parms']
        for i in np.arange( self.nMeshes ):
            infoAllRadex[i]['parms'][:] = self.infoAll[i]['parms']

        self.infoAllRadex = infoAllRadex
        self.meshesRadex = meshesRadex
        
        if writeDb == True:
            self.writeDbRadex()
        

    # sets up a grid, and runs the LVG models over the grid
    # the computed grid of emission stuff is saved in the LVG 
    # method.
    def saveGridsToFiles( self, resGrids, lgammaMechSec, radexParms, ranges = None ):
        self.pltGmSec   = lgammaMechSec
        self.radexParms = radexParms

        if ranges == None:
            raise ValueError('missing value of the parameter ranges\n')
        
                                #xaxis    yaxis 
        self.grds = [ [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 ), 
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 )], 
                      [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 ),
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 )] ]  

        # setting up a dummy mesh object to use its plotting functionalities
        msh = mesh()
        msh.set_chemNet( self.chemNet )
        msh.set_metallicity( self.metallicity )
        
        lG0All   = np.log10(self.infoAll['parms'][:,0])
        lnGasAll = np.log10(self.infoAll['parms'][:,1])
        lgmAll   = np.log10(self.infoAll['parms'][:,2])
        
        indsThisSec = np.nonzero( np.fabs(lgmAll - self.pltGmSec) < 1e-6 )[0]
        #self.computeSurfaceTemperatureGrid( indsThisSec )
        #self.computeAbundanceAtSurfaceGrid( indsThisSec, self.radexParms['specStr'] )
        #self.computeColumnDensityGrid( indsThisSec, self.radexParms['specStr'] )
        if self.radexParms['use'] == True:
            self.computeLineEmissionLvgGrid(indsThisSec, self.radexParms['specStr'])

                    
    
    #################################################################################         
    def plotGrid(self, resGrids, lgammaMechSec, radex = None, ranges = None, *args, **kwargs):
        """Main method for exploring the meshes in the database.
        
        :todo: change resGrids to a [res_x, res_y] insteads of it being just a scalar.
        """
        
        self.set_attributes(**kwargs)
        self.set_default_attributes()

        self.pltGmSec   = lgammaMechSec
        self.radexParms = radex

        if ranges == None:
            raise ValueError('missing value of the parameter ranges\n')
        else:
            nMin  = ranges[0][0]; nMax  = ranges[0][1];
            G0Min = ranges[1][0]; G0Max = ranges[1][1];

        # definig plotting windows and setting the locations of subplots
        fig1, axs1, = pyl.subplots(3, 3, sharex=False, sharey=False, figsize=(14,14))
        self.fig = fig1
        
        #---------------------------------------------------------------------------------------------------
        # the grid plot in n,G0 showing the points where models are present in the
        # database
        # + - -
        # - - - 
        # - - -
        axsGrd = axs1[0,0];  axsGrd_n = 331;
        axsGrd.set_position((0.05, 0.8, 0.12, 0.12))
        pyl.subplot(axsGrd_n)
        pyl.hold(True)
        self.grdPltPts1, = pyl.plot( [0], [0], 'bo' )
        self.grdPltPts2, = pyl.plot( [1], [1], 'ro')
        self.grdPltTitle  = pyl.title('$\log_{10} n_{gas} = $ %4.2f\n$\log_{10} G_0 =$ %4.2f\n$\log_{10} \Gamma_{mech} = $  %5.2f\n' % (0, 0, 0) )
        pyl.xlim( xmin = -1, xmax = 7.0)
        pyl.ylim( ymin = -1, ymax = 7.0)
        pyl.xlabel('$log_{10} n_{gas}$')
        pyl.ylabel('$log_{10} G_0$')
        #---------------------------------------------------------------------------------------------------
        # the subplots where things are plotted as a function of Av       
        # - - -
        # + + -
        # + + - 
        left  = 0.07
        bott  = 0.05
        sz    = 0.20
        vSpace = 0.01
        hSpace = 0.01
        axsAvPlts = np.array([ [axs1[1,0],axs1[1,1]], [axs1[2,0],axs1[2,1]] ])        
        axsAvPlts[0,0].set_position((left              , bott + sz + vSpace, sz, sz))
        axsAvPlts[0,1].set_position((left + sz + hSpace, bott + sz + vSpace, sz, sz))
        axsAvPlts[1,0].set_position((left              , bott              , sz, sz))
        axsAvPlts[1,1].set_position((left + sz + hSpace, bott              , sz, sz))
        axsAvPlts_n = np.array( [[334,335], [337,338]])

        # the subplots where things are plotted in n,G0 plane (grid properties)               
        # - + +
        # - - +
        # - - +
        left  = 0.55
        bott  = 0.45
        sz    = 0.20
        vSpace = 0.035
        hSpace = 0.02
        # defining the new axes array
        axsGrds = np.array( [ [axs1[0,1], axs1[0,2] ], [axs1[1,2], axs1[2,2] ] ])
        axsGrds[0,0].set_position((left              , bott + sz + vSpace, sz, sz))
        axsGrds[0,1].set_position((left + sz + hSpace, bott + sz + vSpace, sz, sz))
        axsGrds[1,0].set_position((left              , bott              , sz, sz))
        axsGrds[1,1].set_position((left + sz + hSpace, bott              , sz, sz))
        axsGrds_n = np.array( [[332,333], [336,339]])
        # setting up the labels of the axes and the major ticks
        pyl.subplot( axsGrds_n[0,0] )
        pyl.ylabel( '$log_{10} G_0$' )
        for tick in pyl.gca().xaxis.get_major_ticks():
            tick.label1On = False
        pyl.xlim( xmin = nMin , xmax = nMax )
        pyl.ylim( ymin = G0Min, ymax = G0Max)
        
        pyl.subplot( axsGrds_n[0,1] )
        for tick in pyl.gca().xaxis.get_major_ticks():
            tick.label1On = False
        for tick in pyl.gca().yaxis.get_major_ticks():
            tick.label1On = False
        pyl.xlim( xmin = nMin , xmax = nMax )
        pyl.ylim( ymin = G0Min, ymax = G0Max)

        pyl.subplot( axsGrds_n[1,0] )
        pyl.xlabel( '$log_{10} n_{gas}$' )
        pyl.ylabel( '$log_{10} G_0$' )
        pyl.xlim( xmin = nMin , xmax = nMax )
        pyl.ylim( ymin = G0Min, ymax = G0Max)
        (pyl.gca().yaxis.get_major_ticks())[-1].label1On = False
        (pyl.gca().xaxis.get_major_ticks())[-1].label1On = False
        
        pyl.subplot( axsGrds_n[1,1] )
        pyl.xlabel( '$log_{10} n_{gas}$' )
        for tick in pyl.gca().yaxis.get_major_ticks():
            tick.label1On = False
        pyl.xlim( xmin = nMin , xmax = nMax )
        pyl.ylim( ymin = G0Min, ymax = G0Max)
        
        # getting the interpolation function which will be used to display the 2D grids
        self.computeAndSetInterpolationFunctions(*args, **kwargs)        
        
        # defining and intialising the ndmesh objects which will be used
        # for computing the grid properties and then displayed                                
                                #xaxis    yaxis 
        self.grds = [ [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [[nMin, nMax], [G0Min, G0Max]], fill = 0.0  ), 
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [[nMin, nMax], [G0Min, G0Max]], fill = 0.0 )], 
                      [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [[nMin, nMax], [G0Min, G0Max]], fill = 0.0 ),
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = [[nMin, nMax], [G0Min, G0Max]], fill = 0.0 )] ]  
                        
        # creating the axes for the colorbars
        self.grdsCbarAxs = [ [pyl.axes([left, bott + sz + sz + vSpace + 0.017, 0.2, 0.01 ]), pyl.axes( [left + sz + hSpace, bott + sz + sz + vSpace + 0.017, 0.2, 0.01] ) ],
                             [pyl.axes([left, bott + sz +               0.017, 0.2, 0.01 ]), pyl.axes( [left + sz + hSpace, bott + sz +               0.017, 0.2, 0.01] ) ] ]
                
        for cbarAxsSubList in self.grdsCbarAxs:
            for cbarAxs in cbarAxsSubList:
                
                for tick in cbarAxs.xaxis.get_major_ticks():
                    tick.label1On = True
                    tick.label2On = False
                for tick in cbarAxs.yaxis.get_major_ticks():
                    tick.label1On = False
                    tick.label2On = False
        
        # setting up the axes to plot the radex output for each selected model
        left  = 0.65
        bott  = 0.07
        sz    = 0.08
        vSpace = 0.0
        self.pltRadex = [pyl.axes([left, bott + 3*sz + vSpace , 3*sz, sz ]),
                         pyl.axes([left, bott + 2*sz + vSpace , 3*sz, sz ]),
                         pyl.axes([left, bott + 1*sz + vSpace , 3*sz, sz ]), 
                         pyl.axes([left, bott + 0*sz          , 3*sz, sz ])]
        self.pltRadex = np.array(self.pltRadex)                           
        
        # setting up a dummy mesh object to use its plotting functionalities
        msh = mesh()
        msh.setFigureObjects(fig1, axs1, axsAvPlts_n)
        msh.setupFigures()
        msh.set_chemNet( self.chemNet )
        msh.set_metallicity( self.metallicity )
        
        lG0All   = np.log10(self.infoAll['parms'][:,0])
        lnGasAll = np.log10(self.infoAll['parms'][:,1])
        lgmAll   = np.log10(self.infoAll['parms'][:,2])
        
        #text object which will show the section in gamma mech
        tt = fig1.text(0.55, 0.02, '$Log_{10}\Gamma_{mech}$ = %5.2f' % self.pltGmSec )
        
        # the axes to plot in the 3D grid (default axes are n,G0,gmech) 
        # showing the points where models are present in the database
        ax3d = fig1.add_subplot(111, projection='3d')
        ax3d.set_position((0.2, 0.68, 0.3, 0.3))
        self.plot_3D_grid_point(figure = fig1, axes = ax3d, ranges = ranges, log10z = self.parms['plotLog10zAxs'])

        self.resPltGrids = [resGrids, resGrids] 
         
        def plotThisSec():
            """updates the plotes once the z section value is changed
            """
            
            print 'plotThisSec(): updating the grids'
            
            tt.set_text('$log_{10}\Gamma_{mech}$ = %5.2f' % self.pltGmSec)            
            indsThisSec = np.nonzero( np.fabs(lgmAll - self.pltGmSec) < 1e-6 )[0]            
            self.grdPltPts1.set_xdata( lnGasAll[indsThisSec] )
            self.grdPltPts1.set_ydata( lG0All[indsThisSec]   )
            
            # plotting the grids
            #-------------------
            # temperature grid (top left grid)
            if self.parms['gridsInfo']['00']['show']:
                pyl.subplot( axsGrds_n[0, 0] )
                self.showSurfaceTemperatureGrid(ranges = ranges, res = self.resPltGrids, *args, **kwargs)
            
            # abundances (top left grid)
            if self.parms['gridsInfo']['01']['show']:
                pyl.subplot( axsGrds_n[0, 1] )
                self.showAbundancesGrid(ranges = ranges, res = self.resPltGrids, *args, **kwargs)

            # column densities (bottom left grid)
            if self.parms['gridsInfo']['10']['show']:
                pyl.subplot( axsGrds_n[1, 0] )
                self.showColumnDensityGrid(ranges = ranges, res = self.resPltGrids, *args, **kwargs)

            # line intensity (bottom right grid)
            if self.parms['gridsInfo']['11']['show']:
                pyl.subplot( axsGrds_n[1, 1] )
                self.showLineIntensityGrid(ranges = ranges, res = self.resPltGrids, *args, **kwargs)
            
                         
        # defining the buttons to control mechanical heating section        
        def nextSec(event):
            self.pltGmSec += 1.0
            plotThisSec()
            pyl.draw()
    
        def prevSec(event):
            self.pltGmSec -= 1.0
            plotThisSec()
            pyl.draw()        
        
        # defining the event when a point in a section is clicked
        def onB1Down(event):
            
            ti = time()

            # get the x and y coords, flip y from top to bottom
            b      = event.button 
            x, y   = event.x, event.y
            xd, yd = event.xdata, event.ydata
        
            if event.button==1:
                if event.inaxes is not None:
            
                    l2Distance  = np.sqrt( (yd - lG0All )**2 + (xd - lnGasAll)**2 + (self.pltGmSec - lgmAll)**2 )
                    rMin = min(l2Distance)
                    indMin = l2Distance.argmin()
                    msh.setData( self.meshes[indMin] )
                                        
                    self.grdPltTitle.set_text('$\log_{10} n_{gas} = $ %4.2f\n$\log_{10} G_0 =$ %4.2f\n$\log_{10} \Gamma_{mech} = $  %5.2f\n' % (np.log10(msh.data['hdr']['nGas']), np.log10(msh.data['hdr']['G0']), np.log10(msh.data['hdr']['gammaMech'])))
                    
                    self.grdPltPts2.set_xdata( lnGasAll[indMin] )
                    self.grdPltPts2.set_ydata( lG0All[indMin] )
                    self.grdPltPts2.set_color('r')
                                        
                    msh.plot()
                                        
                    if self.parms['radex']['use']:
                        self.computeAndSetRadexCurves(meshObj = msh)
                    
                    pyl.draw()
                    
            tf = time()
            print tf - ti
        #----------------------------end method onB1Down---------------------------------
        
        # attaching mouse click event to fig 1
        cid = fig1.canvas.mpl_connect('button_press_event', onB1Down)
        
        plotThisSec()

        # attaching the gamma mech button events
        axPrev = pyl.axes([0.70, 0.02, 0.02, 0.02])
        bPrev = Button(axPrev, '-')
        bPrev.on_clicked(prevSec)

        axNext = pyl.axes([0.73, 0.02, 0.02, 0.02])
        bNext  = Button(axNext, '+')
        bNext.on_clicked(nextSec)
        
        # displaying
        pyl.draw()
        print 'browing data....'

    def computeAndSetRadexCurves(self, meshObj = None):
        """This is a utilty method (make it a private method), for populating the radex axes
          with the plots for the specie passed in the global parameter self.parms['radex']['specStr'].
          It takes as a paremeter the :data:`mesh` which will be used for doing the radex computations.
        """
        
        msh = meshObj
        #assigning the axes in which the radex curves will be plotted
        self.radexObj.setupPlot(nx = 1, fig = self.fig, axs = self.pltRadex)
        
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
            print 'not enough colliders'
        else:
    
            print 'Radex input parms computed from the slab: ', gasTRadex, nColls, colDensThisSpec
                
            self.radexObj.setDefaultStatus()
            self.radexObj.run( checkInput = True, verbose = True )

            if self.radexObj.getStatus() & self.radexObj.FLAGS['SUCCESS']:
                
                self.radexObj.plotModelInFigureColumn(allTrans = None,
                                                      inAxes = self.pltRadex, 
                                                      title='')
                self.radexObj.setLabels()
            else:
                for warning in self.radexObj.warnings:
                    print 'meshUtils.py : ', warning  

        
    def clear(self):
        """clears all the bufferes allocated in the instance"""
        
        del self.infoAll
        del self.nDbFiles
        del self.chemNet
        del self.pltGmSec
        del self.radexParms
        del self.fig
        del self.grdPltPts1
        del self.grdPltPts2
        del self.grdPltTitle
        del self.grds  
        del self.grdsCbarAxs
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
        
        print plts
        fig.legend(plts, names)
        pyl.show()
    
    def plot_3D_grid_point(self, log10z = None, **kwargs):
        """plots in 3D the parameters of the meshes in the database. By default
           the x,y,z coordinates are the log10 of nGas, G0, and gammaMech
           
           :param matplotlib.pyplot.figure figure: if this keyword is passed, a new figure is not 
               created and the axes are created in the figure object passed.
           :param  mpl_toolkits.mplot3d.Axes3D axes: if this keyword is passed, things are plotted
               in these axes, otherwise new ones are created in the figure.
           :param bool log10z: when this is set and it is true, the log10 of the z quantity is plotted.  
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
        

        if log10z != None and log10z == True:
            zPlt = np.log10(self.grid_z)
        else:
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
          reference database is low enough to be assumed to be zero. This keyword
          cuases self.grid_qz to be set to 'gMech/gSurface(gMech=0)'

        :param string referenceDatabasePath: a string containing the path of the 
             reference database.
         
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
                if 'referenceDatabasePath' not in kwargs:
                    raise ValueError('missing the keyword argument "referenceDatabasePath"')
                else:
                    #computing (for the mesh points in the current database) the ratio of 
                    #the mechanical heating to the surface heating when gmech = 0
                    
                    self.set_grid_qz(['','gMech/gSurface(gMech=0)'])
                    
                    # reading the reference archive
                    print 'setting up the reference archive'
                    t0 = time()
                    arxvRef = meshArxv( dirPath = kwargs['referenceDatabasePath'], 
                                        metallicity = self.metallicity )
                    arxvRef.readDb( check = True )
                    print 'time reading %f' % (time() - t0)
                    arxvRef.setChemicalNetwork(self.chemNet) # assiginig the chemical network to the archive
    
                    gMechZero = self.grid_x.copy()
                    gMechZero[:] = np.log10(arxvRef.meshes[0]['hdr']['gammaMech']) # assumig all the meshes have this value as well
                    
                    arxvRef.set_grid_qx( self.grid_qx )
                    arxvRef.set_grid_qy( self.grid_qy )
                    f         = arxvRef.construct3DInterpolationFunction(quantity = ['therm', 'heating'], 
                                                                         slabIdx  = 0, 
                                                                         log10 = True,
                                                                         grid_qz = ['hdr','gammaMech'])
                    dataNew   = np.array( [self.grid_x, self.grid_y, gMechZero] ).T
                    
                    # the surface heating that the models in the current would have if the 
                    # mechanical heating were zero
                    gammaSurf = f(dataNew) 
            
                    z = np.log10(self.getQuantityFromAllMeshes( ['hdr', 'gammaMech']) )
                    self.grid_z = 10.0**z / 10.0**gammaSurf
            else:
                # just use gMech as the 3rd axis
                self.grid_z = np.log10( self.getQuantityFromAllMeshes(self.grid_qz) )
    
    """
    def computeSurfaceTemperatureGrid( self, res = None, ranges = None ):
        #generates color map of the surface temperature. if meshInds is not provided,
        #    all the meshes in the arxive are used. In that case, res and ranges shoudl
        #    be provided.
        #    
        #    :deprecated: remove this function...to slow..
        # checking stuff
        if meshInds == None:
            meshInds = np.arange(self.nMeshes)
        
        if self.grds == None: # for standalone use of this function
            if res == None:
                raise ValueError('res = None, i need to have the resolution of the grid first!!')
            if ranges == None:
                raise ValueError('ranges = None, i need to have the ranges!!!')
            grd = ndmesh( (res, res), dtype = np.float64, ranges = ranges, fill = 0.0) 
        else:
            grd = self.grds[0][0] # should have called self.plot first
            ranges = grd.ranges 
            grd.fill(0.0) 
        
        # doing the computations to fill the grid 
        tGasGrid = grd     
        spcs = self.chemNet.species
        nInCells = tGasGrid.copy()
        nx, ny = tGasGrid.shape
        
        # computing the surface temperature grid
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            #gasT = self.meshes[i]['state']['gasT']
            gasT = fetchNestedDtypeValue(self.meshes[i], ['state', 'gasT'] )
            zThis = gasT[0]
            
            indxInGrid = scale(xThis, 0, nx, ranges[0][0], ranges[0][1], integer = True) 
            indyInGrid = scale(yThis, 0, ny, ranges[1][0], ranges[1][1], integer = True) 
        
            tGasGrid[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1
        
        #print tGasGrid
        #print nInCells
        tGasGrid[:] = np.log10(tGasGrid / nInCells)
        
        return tGasGrid
    """