import os, glob, sys
import numpy as np
import pylab as pyl
from matplotlib.widgets import Button
import matplotlib.cm as cm

from time import *
from mesh import *
from utils import *
from ndmesh import *
from ismUtils import *
from radex import *
from utils import fetchNestedDtypeValue
from scipy import interpolate

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
    def __init__(self, *args, **kwrds):
        
        # instance variables
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
        """A numpy array of dtype arxvHdrDtype (see below) which contains the info
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
        
        if 'metallicity' in kwrds:
            self.set_metallicity( kwrds['metallicity'] )
        else:
            self.metallicity = None
            
        self.radexObj   = None
        self.radexParms = None
        
    # read all the meshes files in the dir and construct the
    # database
    def construct(self, dirName , meshNamePrefix = None, writeDb = None ):
        """ construc the database anad write the .db files. If the meshNamePrefix is
              not supplied all the files in dirName are assumed to be data files
              and all of them are put in the database."""
        
        if meshNamePrefix == None:
            meshNamePrefix = ''

        # getting the names of the meshes in that dir
        files = []
        for infile in glob.glob( os.path.join(dirName, 'meshes/'+meshNamePrefix + '*') ):
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
            self.writeDb(dirName)

    def writeDb(self, dirName):
        """ writes the db files into the dir dirName"""
                
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

        
    def readDb(self, dirName ):
        """ reads the database and assigns the appropritate attributes (document)""" 
        
        dbInfoFObj = file(dirName + 'meshes.db.info', 'rb')
        dbDataFObj = file(dirName + 'meshes.db'  , 'rb')

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
                
    def constructTemperatureGrid(self):
        return 1
    
    def constructModelsGrid(self, log_nGas=None, log_G0=None):
        
        if log_nGas == None:
            errStr = 'denisty range of the grid not set'
            raise NameError(errStr)
        if log_G0 == None:
            errStr = 'FUV G0 range of the grid not set'
            raise NameError(errStr)            

    def getQuantityFromAllMeshes(self, quantity, slabIdx = None):
        """ gets the quantity from all the meshes and returns it as a numpy array. 
            the quantity is mandatory, but no the slabIdx.
        """
        
        values = np.zeros(self.nMeshes, dtype = np.float64)
        for i in np.arange(self.nMeshes):
            q = fetchNestedDtypeValue(self.meshes[i], quantity )
            if slabIdx != None:
                values[i] = q[slabIdx]
            else:
                values[i] = q  
                
        return values
        
    def construct3DInterpolationFunction(self, quantity = None, slabIdx = None, log10 = None, *args, **kwargs):
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
                 
             :param bool log10: keyword which when passes as True, will generate the log10 of the quantity
                 
             :todo: modify this to construct the table over a selected range of the 3D 
                parameter space instead of using all the meshes. (optional)
                
             :warning: this fails if all the entries in one of the dimensions have the exact
               same value. in that case use :data:`construct2DInterpolationFunction`
            
        """
        
        x = np.log10( self.getQuantityFromAllMeshes(self.grid_qx) )
        y = np.log10( self.getQuantityFromAllMeshes(self.grid_qy) )
        z = np.log10( self.getQuantityFromAllMeshes(self.grid_qz) )
        
        values = self.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx)

        if log10 != None and log10 == True:
            values[:] = np.log10(values[:])
        
        data = np.array([x, y, z]).T  #3D
        #data = np.array([x, y]).T  #2D
        
        ti = time()
        f = interpolate.LinearNDInterpolator(data, values) # getting the interpolation function     
        tf = time()
        print 'constructed the interpolation function from %d points in %f seconds' % (len(values), tf-ti)

        self.grid_x, self.grid_y, self.grid_z = (x, y, z)
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
                                 
    
    def computeInterpolated2DGrid(self, quantity = None, slabIdx = None, ranges = None, res = None, zSec = None, 
                                  fInterp = None, *args, **kwargs):
        """ returns a 2D array ( a numpy ndarray ) of size res[0] and res[1] (check this if it is not the reverse) which holds
             the interpolated vlaues of 'quantity' over the domain determined by ranges for a slab whose index is slabIdx for
             a mechanical heating zSec (in log10). ( x is the horizontal direction of the mesh, and y is the vertical).
             
             An example would be (assuming the archive is constructed already):
             
             .. code-block:: python

                f   = arxv.construct3DInterpolationFunction(quantity = ['state', 'gasT'], 
                                                            slabIdx  = 0)
                grd = arxv.computeInterpolated2DGrid(quantity = ['state', 'gasT'], 
                                                     slabIdx  = 0, 
                                                     ranges   = [[0,6],[0,6]],
                                                     res      = [100,100], 
                                                     zSec     = -30, 
                                                     fInterp  = f)


                import matplotlib.pyplot as plt
                plt.imshow(grd, extent=(0,1,0,1), origin='lower')
                plt.show()

             
             :param list quantity: (see definition in :data:`construct3DInterpolationFunction`
                 
             :param int32 slabIdx: (see definition in :data:`construct3DInterpolationFunction`
             
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
    
    def getInterpFunctionGmechToSurfaceHeating(self, quantity, referenceDatabasePath = None, slabIdx = None, log10 = None, *args, **kwargs):
        """ This methods returns a 3D interpolation function (n,G0,r) where r is
            gMech/GammaSurface(Gmech=0) i.e the ratio of the mechanical heating used
            for a grid, devided by the surface heating when the mechanical heating is zero.
            So given n,G0,r the interpoaltion function computes the quantity at the locations
            we desire. 
            GammaSurface(gMech=0) is obtained from a reference database whose path is passed
            as a parameter.
            
            :WARNING: it is assumed that the reference database has the same chemical network
            and metallicity as the current one. it is also assumed that all the models in the 
            reference database have the same mechanical heating rate which is too small and
            can be neglected and treated as zero.
        """
        
        self.set_grid_qz('gMech/gSurface(gMech=0)')
        
        self.set_attributes(**kwargs)
        self.set_default_attributes()
        
        # reading the reference archive
        print 'setting up the reference archive'
        t0 = time()
        arxvRef = meshArxv( metallicity = self.metallicity )
        arxvRef.readDb( referenceDatabasePath )
        arxvRef.checkIntegrity()
        print 'time reading %f' % (time() - t0)
        arxvRef.setChemicalNetwork(self.chemNet) # assiginig the chemical network to the archive
        
        # the x and y of the current database which will be used as input to get gSurface(gMech=0)
        # from the referece database 
        x = np.log10(self.getQuantityFromAllMeshes( self.grid_qx ) )
        y = np.log10(self.getQuantityFromAllMeshes( self.grid_qy ) )
        gMechZero = x.copy()
        gMechZero[:] = np.log10(arxvRef.meshes[0]['hdr']['gammaMech']) # assumig all the meshes have this value as well
        
        arxvRef.set_grid_qx( self.grid_qx )
        arxvRef.set_grid_qy( self.grid_qy )
        f         = arxvRef.construct3DInterpolationFunction(quantity = ['therm', 'heating'], slabIdx  = 0, log10 = True)
        dataNew   = np.array( [x, y, gMechZero] ).T
        
        # the surface heating that the models in the current would have if the 
        # mechanical heating were zero
        gammaSurf = f(dataNew) 

        z = np.log10(self.getQuantityFromAllMeshes( ['hdr', 'gammaMech']) )
        r = 10.0**z / 10.0**gammaSurf

        # the quantity to be interpolated
        values = self.getQuantityFromAllMeshes( quantity, slabIdx = slabIdx)
        
        if log10 != None and log10 == True:
            values[:] = np.log10(values[:])
        
        # coordinates where the interpolation will be done
        data = np.array([x, y, r]).T
        
        f = interpolate.LinearNDInterpolator(data, values)      

        # changin values of attributes
        self.grid_x, self.grid_y, self.grid_z = (x, y, r)
        
        return f
    
    def showGrid(self, quantity = None, slabIdx = None, ranges = None, res = None, zSec = None, fInterp = None, *args, **kwargs):
        """This method displays a plot of a grid of the quantity pointed by quantity in a new
           window. It makes use of :data:`computeInterpolated2DGrid` and :data:`construct3DInterpolationFunction`.
           For the documentation of the parameters see  :data:`computeInterpolated2DGrid`
        """

        self.set_attributes(**kwargs)
        self.set_default_attributes()
        
        # getting the grid
        if fInterp == None:
            f = self.construct3DInterpolationFunction(quantity = quantity, slabIdx  = slabIdx, *args, **kwargs)
        else:
            f = fInterp
        grd = self.computeInterpolated2DGrid(quantity = quantity, 
                                             slabIdx  = slabIdx, 
                                             ranges   = ranges,
                                             res      = res, 
                                             zSec     = zSec, 
                                             fInterp  = f, *args, **kwargs)
        
        # checking/setting the quantities to be consided as the x and y axes of the grids 
        if self.grid_qx == None:
            self.grid_qx = ['hdr', 'nGas']
        if self.grid_qy == None:
            self.grid_qy = ['hdr', 'G0']
        
        # taking the transpose for plotting purposes
        grd = grd.T
        
        # plotting it
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
                
    def computeSurfaceTemperatureGrid( self, meshInds = None, res = None, ranges = None ):
        """generates color map of the surface temperature. if meshInds is not provided,
            all the meshes in the arxive are used. In that case, res and ranges shoudl
            be provided.
        """
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

    # produce the grid for heating or cooling determined with
    #    whichThermal = 'heating' | 'cooling'
    # and the specific process is determined by 
    #   whichProcsess = string
    def computeHeatingCoolingGrid( self, slabIndex = None, meshInds = None, whichThermal = None, whichProcess = None, res = None, ranges = None):
        """generates color map of the various heating and cooling processes for a certain slab.
        if meshInds is not provided, all the meshes in the arxiv are used. In that case, res and ranges shoudl
            be provided.
        """
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
            grd = self.grds[1][1] # should have called self.plot first
            ranges = grd.ranges 
            grd.fill(0.0) 
        
        # doing the computations to fill the grid 
        nInCells = grd.copy()
        nx, ny = grd.shape

        # computing the surface temperature grid
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            
            proc  = self.meshes[i][whichThermal][whichProcess]
            #proc1  = self.meshes[i]['therm']['heating']  
            #proc2  = self.meshes[i]['therm']['cooling']
            #proc   = proc1 / proc2
            
            zThis = proc[slabIndex]

            indxInGrid = scale(xThis, 0, nx, ranges[0][0], ranges[0][1], integer = True)
            indyInGrid = scale(yThis, 0, ny, ranges[1][0], ranges[1][1], integer = True) 
        
            grd[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1
        
        #print tGasGrid
        #print nInCells
        grd[:] = np.log10(grd / nInCells)
        
        return grd
        
    def computeAbundanceAtSurfaceGrid( self, meshInds, specStr ):
        abunGrid = self.grds[0][1] 
        nInCells = abunGrid.copy()
        abunGrid.fill(0.0)
        nInCells.fill(0.0)
        nx, ny = abunGrid.shape
        
        # computing the abundace of a specie
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            
            abunAllSpcs = self.meshes[i]['state']['abun']
            specIdx = self.chemNet.species[specStr].num
            slabIdx = 0
            zThis = abunAllSpcs[specIdx][slabIdx]  #<----------- 

            indxInGrid = scale(xThis, 0, nx, 0, 6.0, integer = True) 
            indyInGrid = scale(yThis, 0, ny, 0, 6.0, integer = True) 
        
            abunGrid[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1
        
        abunGrid[:] = np.log10(abunGrid / nInCells)

    def computeColumnDensityGrid( self, meshInds, specStr ):            
        colDensGrid = self.grds[1][0] 
        nInCells = colDensGrid.copy()
        colDensGrid.fill(0.0)
        nInCells.fill(0.0)
        nx, ny = colDensGrid.shape

        # computing the abundace of a specie
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            abunAllSpcs = self.meshes[i]['state']['abun']
            specIdx = self.chemNet.species[specStr].num
            Av = self.meshes[i]['state']['Av']
            nDensThisSpec = (10**xThis)*abunAllSpcs[ specIdx ][:]
            # setting the thickness of the last slab to the one before it
            dxSlabs          =  getSlabThicknessFromAv(Av, 10**xThis, self.metallicity)
            dxSlabsNew       =  np.ndarray( len(dxSlabs)+1, dtype = np.float64 )
            dxSlabsNew[0:-1] =  dxSlabs
            dxSlabsNew[-1]   =  dxSlabs[-1]
            dxSlabs          =  dxSlabsNew
            
            colDensThisSpec = np.sum( dxSlabs * nDensThisSpec )     
            zThis = colDensThisSpec  # <---------------
            #print zThis 

            indxInGrid = scale(xThis, 0, nx, 0, 6.0, integer = True) 
            indyInGrid = scale(yThis, 0, ny, 0, 6.0, integer = True) 
        
            colDensGrid[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1

        colDensGrid[:] = np.log10(colDensGrid / nInCells)

    def computeLineEmissionLvgGrid( self, meshInds, specStr):
        lineIntense = self.grds[1][1] 
        nInCells = lineIntense.copy()
        lineIntense.fill(0.0)
        nInCells.fill(0.0)
        nx, ny = lineIntense.shape
        
        radexObj = radex(self.radexParms['radexPath'], self.radexParms['molDataDirPath'])
        inFile = { 'specStr'                : self.radexParms['specStr'],
                   'freqRange'              : [0, 50000]      ,
                   'tKin'                   : None            ,
                   'collisionPartners'      : None            ,
                   'nDensCollisionPartners' : None            ,
                   'tBack'                  : 2.73            ,
                   'molnDens'               : None            ,
                   'lineWidth'              : 1.0             ,
                   'runAnother'             : 1               }
        radexObj.setInFile( inFile )

        transitionNum = self.radexParms['plotTransitionInGrid']
        every = 100
        nDone = 0
        upper = None
        lower = None
        # computing the abundace of a specie
        for i in meshInds[0::every]:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            
            thisMeshObj = mesh( None, self.chemNet, self.metallicity)
            thisMeshObj.setData( self.meshes[i] )
            
            print 'pdr mesh index =  %d : ' % i

            (gasTRadex, nColls, colDensThisSpec,) = thisMeshObj.getRadexParameters('H2',  # ;;; this parameter is redundant
                                                                                   self.radexParms['specStr'],  
                                                                                   self.radexParms['xH2_Min'])  # ;;; this parameter is redundant

            # getting the collider densities in the same order of the supplied input spcie string list 
            nDensColls = [ nColls[collSpecStr] for collSpecStr in self.radexParms['collisionPartners'] ]
            collsStr   = list(self.radexParms['collisionPartners'])
            #print 'input coll species', self.radexParms['collisionPartners'] 
            #print 'nColls after putting them in the right order = ', nDensColls

            print 'radexGrid : radex parms : ', gasTRadex, nDensColls, colDensThisSpec
            radexObj.setInFileParm('tKin', gasTRadex)
            radexObj.setInFileParm('collisionPartners', collsStr )
            radexObj.setInFileParm('nDensCollisionPartners', nDensColls )
            radexObj.setInFileParm('molnDens', colDensThisSpec)

            radexObj.filterColliders()
            
            if len(radexObj.inFile['collisionPartners']) == 0:
                print 'not enough colliders'
                continue
            
            radexObj.setDefaultStatus()
            status = radexObj.run( checkInput = True )

            if status & radexObj.FLAGS['SUCCESS']:
                print 'radexGrid : converged with no warnings'
            else:
                print 'radexGrid : converged with warnings'
                print '------------------------------------'
                print radexObj.getWarnings()
                print '------------------------------------'
                continue
                           
            transition = radexObj.getTransition( transitionNum )
            upperStr = transition['upper']
            lowerStr = transition['lower']
            zThis = transition['fluxcgs']
                            
            indxInGrid = scale(xThis, 0, nx, 0, 6.0, integer = True) 
            indyInGrid = scale(yThis, 0, ny, 0, 6.0, integer = True) 
        
            lineIntense[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1
            #print 'press a key to continue...'
            #print input()
            nDone += 1
            print 100.0*np.float64(nDone)/np.float64(len(meshInds)), '%'
            print '----------------------------------------------------' 
        
        lineIntense[:] = lineIntense / nInCells

        if nDone > 0:
            #-----writing the grid to a file -------------------------
            # ;;;; dump all the transitions
            fName = ('%s/%s-%s-%s-%.1f%d.npy') % ('/home/mher/ism/docs/paper02/lineData',
                                                  self.radexParms['specStr'],
                                                  upperStr,
                                                  lowerStr,
                                                  self.metallicity,
                                                  self.pltGmSec)
            print fName 
            np.save(fName, lineIntense )
            #------ done writing the data to a file -------------------
        
        lineIntense[:] = np.log10(lineIntense)
        #print lineIntense

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
        self.computeLineEmissionLvgGrid(indsThisSec, self.radexParms['specStr'])

                    
    
    #################################################################################         
    def plotGrid(self, resGrids, lgammaMechSec, radexParms, ranges = None):

        self.pltGmSec   = lgammaMechSec
        self.radexParms = radexParms

        if ranges == None:
            raise ValueError('missing value of the parameter ranges\n')
        else:
            nMin  = ranges[0][0]; nMax  = ranges[0][1];
            G0Min = ranges[1][0]; G0Max = ranges[1][1];

        # definig plotting windows and setting the locations of subplots
        fig1, axs1, = pyl.subplots(3, 3, sharex=False, sharey=False, figsize=(12,12))
        self.fig = fig1
        
        # the grid plot in n,G0 showing the points where models are present in the
        # database
        # + - -
        # - - - 
        # - - -
        axsGrd = axs1[0,0];  axsGrd_n = 331;
        axsGrd.set_position((0.1,0.7,0.25,0.25))
        pyl.subplot(axsGrd_n)
        pyl.hold(True)
        self.grdPltPts1, = pyl.plot( [0], [0], 'bo' )
        self.grdPltPts2, = pyl.plot( [1], [1], 'ro')
        self.grdPltTitle  = pyl.title('$\log_{10} n_{gas} = $ %4.2f $\log_{10}  = G_0$ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n' % (0, 0, 0) )
        pyl.xlim( xmin = -1, xmax = 7.0)
        pyl.ylim( ymin = -1, ymax = 7.0)
        pyl.xlabel('$log_{10} n_{gas}$')
        pyl.ylabel('$log_{10} G_0$')

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
        

        # plot a section in gamma mech
        def plotThisSec():
            
            tt.set_text('$log_{10}\Gamma_{mech}$ = %5.2f' % self.pltGmSec)            
            indsThisSec = np.nonzero( np.fabs(lgmAll - self.pltGmSec) < 1e-6 )[0]            
            self.grdPltPts1.set_xdata( lnGasAll[indsThisSec] )
            self.grdPltPts1.set_ydata( lG0All[indsThisSec]   )
            
            # plotting the grids
            # ---> temperature grid (top left grid)
            #--------------------------------------------------------------
            pyl.subplot( axsGrds_n[0, 0] )
            ###
            self.computeSurfaceTemperatureGrid( indsThisSec, ranges = ranges )
            ###
            im00 = self.grds[0][0].imshow(interpolation='nearest') 
            cbar00 = pyl.colorbar(im00, cax=self.grdsCbarAxs[0][0], ax=pyl.gca(), orientation = 'horizontal')
            cbar00.set_ticks([0.0, 1.0, 2.0, 3.0, 4.0])

            # some other diagnostic (top left grid)
            # ---> plotting abundances
            #--------------------------------------------------------------
            pyl.subplot( axsGrds_n[0, 1] )
            ###
            self.computeAbundanceAtSurfaceGrid( indsThisSec, self.radexParms['specStr'] )
            ###
            im01 = self.grds[0][1].imshow(interpolation='nearest')
            cbar01 = pyl.colorbar(im01, cax=self.grdsCbarAxs[0][1], ax=pyl.gca(), orientation = 'horizontal')
            cbar01.set_ticks([-4.0, -3.0, -2.0, -1.0, 0.0])

            # some other diagnostic (bottom left grid)
            # ---> plotting column densities
            #--------------------------------------------------------------
            pyl.subplot( axsGrds_n[1, 0] )
            ###
            self.computeColumnDensityGrid( indsThisSec, self.radexParms['specStr'] )
            ###
            im10 = self.grds[1][0].imshow(interpolation='nearest')
            cbar10 = pyl.colorbar(im10, cax=self.grdsCbarAxs[1][0], ax=pyl.gca(), orientation = 'horizontal')
            cbarTickValues =  [15, 16, 17, 18, 19, 20]
            contourValues  =  [15, 16, 17, 18, 18.5, 18.7, 18.9, 19, 19.05, 20] #CO
            cbar10.set_ticks( cbarTickValues )
            self.grds[1][0].plotContour( levels = contourValues )
                   
            # some other diagnostic (bottom right grid)
            # ---> plotting line intensities
            #--------------------------------------------------------------
            pyl.subplot( axsGrds_n[1, 1] )
            ###
            """
            self.computeLineEmissionLvgGrid(indsThisSec, self.radexParms['specStr'])
            #cbarTickValues =  [-12, -11, -10, -9, -8, -7, -6, -5, -4]
            cbarTickValues =  [-2, -1, 0, 1, 2]
            """
            ###    
            self.computeHeatingCoolingGrid(slabIndex = 0, meshInds = indsThisSec, whichThermal = 'heating', whichProcess = 'photo', ranges = ranges)
            cbarTickValues =  [-26, -24, -22, -20, -18, -16, -14]
            im11 = self.grds[1][1].imshow(interpolation='nearest')
            cbar11 = pyl.colorbar(im11, cax=self.grdsCbarAxs[1][1], ax=pyl.gca(), orientation = 'horizontal')            
            cbar11.set_ticks( cbarTickValues )
            self.grds[1][1].plotContour( levels = cbarTickValues, colors = 'black' )
            print 'done'
        
                         
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
                                        
                    self.grdPltTitle.set_text('$\log_{10} n_{gas} = $ %4.2f $\log_{10}  = G_0$ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n' % (np.log10(msh.data['hdr']['nGas']), np.log10(msh.data['hdr']['G0']), np.log10(msh.data['hdr']['gammaMech'])))
                    
                    self.grdPltPts2.set_xdata( lnGasAll[indMin] )
                    self.grdPltPts2.set_ydata( lG0All[indMin] )
                    self.grdPltPts2.set_color('r')
                                        
                    msh.plot()
                    
                    
                    #----------------------------------------------------
                    if self.radexObj == None:       
                        
                        self.radexObj = radex(self.radexParms['radexPath'], self.radexParms['molDataDirPath'])                 
                        inFile = {'specStr'                : self.radexParms['specStr'],
                                  'outPath'                : 'foo'       ,
                                  'freqRange'              : [0, 0]  ,
                                  'tKin'                   : None        ,
                                  'collisionPartners'      : None      ,
                                  'nDensCollisionPartners' : None      ,
                                  'tBack'                  : 2.73        ,
                                  'molnDens'               : None        ,
                                  'lineWidth'              : 1.0         ,
                                  'runAnother'             : 1           }
                        self.radexObj.setInFile( inFile )
                       
                        
                    self.radexObj.setupPlot(nx = 1, fig = self.fig, axs = self.pltRadex)
                    
                    (gasTRadex, nColls, colDensThisSpec,) = msh.getRadexParameters('H2',   # ;;; this parameter is redundant 
                                                                                   self.radexParms['specStr'], 
                                                                                   self.radexParms['xH2_Min'])  # ;;; this parameter is redundant
    
                    # getting the collider densities in the same order of the supplied input spcie string list 
                    nDensColls = [ nColls[collSpecStr] for collSpecStr in self.radexParms['collisionPartners'] ]
                    collsStr   = list(self.radexParms['collisionPartners'])
                    #print 'input coll species', self.radexParms['collisionPartners'] 
                    #print 'nColls after putting them in the right order = ', nDensColls

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
                        #----------------------------------------------------
                    
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

    # setters and getters
    def set_metallicity(self, metallicity):
        self.metallicity = metallicity
        
    ## clears all the bufferes allocated in the instance
    def clear(self):
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
    
    def plot_3D_grid_point(self, **kwargs):
        """plots in 3D the parameters of the meshes in the database. By default
           the x,y,z coordinates are the log10 of nGas, G0, and gammaMech
        """
        
        if self.grid_x == None or self.grid_y == None or self.grid_z == None:
            raise TypeError('data of one or more coordinates not set!!.')
         
        fig = pyl.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.grid_x, self.grid_y, self.grid_z, 'o')
        ax.set_xlim( kwargs['ranges'][0][0], kwargs['ranges'][0][1] )
        ax.set_ylim( kwargs['ranges'][1][0], kwargs['ranges'][1][1] )
        ax.set_zlim( kwargs['ranges'][2][0], kwargs['ranges'][2][1] )
        ax.set_xlabel(self.grid_qx[-1])
        ax.set_ylabel(self.grid_qy[-1])
        ax.set_zlabel(self.grid_qz[-1])

        
    ###############################setter and getters##############################
    def setChemicalNetwork(self, chemNet):
        self.chemNet = chemNet
        
    def set_grid_qx(self, qx):
        self.grid_qx = qx
    def set_grid_qy(self, qy):
        self.grid_qy = qy
    def set_grid_qz(self, qz):
        self.grid_qz = qz

    def set_attributes(self, **kwargs):
        """set values of attributes from provided keywords
        """
        
        if 'grid_qx' in kwargs:
            self.set_grid_qx( kwargs['grid_qx'] )
        if 'grid_qy' in kwargs:
            self.set_grid_qy( kwargs['grid_qy'] )
        
    def set_default_attributes(self):
        """sets the default values of the attributes if they are not set
        """
        
        if self.grid_qx == None:
            self.set_grid_qx( ['hdr', 'nGas'] )
        if self.grid_qy == None:
            self.set_grid_qy( ['hdr', 'G0'] )
        if self.grid_qz == None:
            self.set_grid_qz( ['hdr', 'gammaMech'] )
            
        