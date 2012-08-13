import os, glob, sys
import numpy as np
import pylab as pyl
from matplotlib.widgets import Button
from time import *
from mesh import *
from utils import *
from ndmesh import *
from ismUtils import *
from radex import *

class meshArxv():
    """ this class generates and manipulates archives of meshes.
               
     FILES AND THEIR FORMATS:
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
            - data file index ( int64 ) :
               i.e in which .db file the mesh is located
            - offset ( int64 ) :
               the offset in bytes from the beginning of the file where the mesh is located
            - nSteps ( int64 ):
               the number of slabs in the mesh
            - nSpecs ( int64 ) :
               the number of species in the mesh                                   .
                        
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
        """ nump.long64 : number of meshes in the database """

        self.ver        = None
        """ numpy int32 array [3] version number of the database data"""
        
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
        """ a numpy array of dtype arxvHdrDtype (see below) which contains the info
         (headers) about all the meshes each entry in this array contains two things,
          information about the mesh and the parameters of the mesh. For example
          the elements x in self.infoAll[x] has the following contents : 
          self.infoAll[x]['info'][0]   mesh number 
          self.infoAll[x]['info'][1]   0 ( for now)
          self.infoAll[x]['info'][2]   offset from the start of file
          self.infoAll[x]['info'][3]   number of steps in the mesh  
          self.infoAll[x]['info'][4]   number of species
         
          self.infoAll[x]['parms'][0]  G0
          self.infoAll[x]['parms'][1]  nGas
          self.infoAll[x]['parms'][2]  gammaMech
          self.infoAll[x]['parms'][3]  0.0, NOT USED
          self.infoAll[x]['parms'][4]  0.0  NOT USED
        """

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
        self.pltRadex = None
        self.fig      = None
        self.grdPltPts1 = None
        self.grdPltPts2 = None
        self.grdPltTitle = None
        
        if 'metallicity' in kwrds:
            self.set_metallicity( kwrds['metallicity'] )
        else:
            self.metallicity = None
            
        self.radexObj   = None
        self.radexParms = None
        
    # read all the meshes files in the dir and construct the
    # database
    def construct(self, dirName ):

        meshNamePrefix = 'mesh.dat-id-'        

        # defining the files objects for the database fiels        
        dbInfoFObj = file(dirName + 'meshes.db.info', 'wb')
        dbDataFObj = file(dirName + 'meshes.db'  , 'wb')
       
        # getting the names of the meshes in that dir
        files = []
        #for infile in glob.glob( os.path.join(dirName, meshNamePrefix + '*0000*') ):
        for infile in glob.glob( os.path.join(dirName, 'meshes/'+meshNamePrefix + '*') ):
            files.append(infile)
        #print files
        
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
            infoAll[i]['info'][2] = dbDataFObj.tell()  # offset from the start of file
            infoAll[i]['info'][3] = mData['hdr']['nSteps']  
            infoAll[i]['info'][4] = mData['hdr']['nSpecs']
            infoAll[i]['info'][5] = mData['hdr']['version']
     
            infoAll[i]['parms'][0] = mData['hdr']['G0']
            infoAll[i]['parms'][1] = mData['hdr']['nGas']
            infoAll[i]['parms'][2] = mData['hdr']['gammaMech']
            infoAll[i]['parms'][3] = np.float64(0.0)
            infoAll[i]['parms'][4] = np.float64(0.0)
             
            # writing the mesh to the database file
            mData.tofile( dbDataFObj )
            np.array( infoAll[i]['info'][2] ).tofile(dbDataFObj)
            del m
            
        self.infoAll = infoAll
        
        # writing basic info to the .info file
        self.ver.tofile( dbInfoFObj)
        self.nMeshes.tofile( dbInfoFObj )
        self.infoAll.tofile( dbInfoFObj )
        
        # closing the files
        dbDataFObj.close()
        dbInfoFObj.close()
        print 'wrote successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)

    def readDb(self, dirName ):
        
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
        
    # check archive integrity, this is a basic check where the information about the
    # meshes in self.infoAll['info] is compared to those in self.data for each mes
    def checkIntegrity(self):
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
        
    def setChemicalNetwork(self, chemNet):
        self.chemNet = chemNet
        
    def constructTemperatureGrid(self):
        return 1
    
    def constructModelsGrid(self, log_nGas=None, log_G0=None):
        
        if log_nGas == None:
            errStr = 'denisty range of the grid not set'
            raise NameError(errStr)
        if log_G0 == None:
            errStr = 'FUV G0 range of the grid not set'
            raise NameError(errStr)            

    def computeSurfaceTemperatureGrid( self, meshInds = None, res = None, ranges = None ):
        
        if not meshInds:
            meshInds = np.arange(self.nMeshes)
        
        if self.grds == None: # for standalone use of this function
            if res == None:
                raise ValueError('res = None, i need to have the resolution of the grid first!!')
            if ranges == None:
                raise ValueError('ranges = None, i need to have the ranges!!!')
            tGasGrid = ndmesh( (res, res), dtype = np.float64, ranges = ranges, fill = 0.0) 
        else:
            tGasGrid = self.grds[0][0] # should have called self.plot first 
                    
        spcs = self.chemNet.species
        nInCells = tGasGrid.copy()
        nInCells.fill(0.0)
        tGasGrid.fill(0.0)
        nx, ny = tGasGrid.shape
        
        # computing the surface temperature grid
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            gasT = self.meshes[i]['state']['gasT']
            
            zThis = gasT[0]

            indxInGrid = scale(xThis, 0, nx, 0, 6.0, integer = True) 
            indyInGrid = scale(yThis, 0, ny, 0, 6.0, integer = True) 
        
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
    def computeHeatingCoolingGrid( self, meshInds, slabIndex, whichThermal, whichProcess ):
        
        grid = self.grds[1][1]
        nInCells = grid.copy()
        nInCells.fill(0.0)
        grid.fill(0.0)
        nx, ny = grid.shape

        # computing the surface temperature grid
        for i in meshInds:
            xThis = np.log10(self.meshes[i]['hdr']['G0'])
            yThis = np.log10(self.meshes[i]['hdr']['nGas'])
            
            #proc  = self.meshes[i][whichThermal][whichProcess]
            proc1  = self.meshes[i]['therm']['heating']  
            proc2  = self.meshes[i]['therm']['cooling']
            proc   = proc1 / proc2
            
            print proc
            
            zThis = proc[slabIndex]

            indxInGrid = scale(xThis, 0, nx, 0, 6.0, integer = True) 
            indyInGrid = scale(yThis, 0, ny, 0, 6.0, integer = True) 
        
            grid[indyInGrid][indxInGrid] += zThis
            nInCells[indyInGrid][indxInGrid] += 1
        
        #print tGasGrid
        #print nInCells
        grid[:] = np.log10(grid / nInCells)
        
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

                                #xaxis    yaxis 
        self.grds = [ [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 ), 
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 )], 
                      [ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 ),
                       ndmesh( (resGrids, resGrids), dtype=np.float64, ranges = ranges, fill = 0.0 )] ]  

        if ranges == None:
            nMin  = 0.0; nMax  = 6.0;
            G0Min = 0.0; G0Max = 6.0;
        else:
            nMin  = ranges[0,0]; nMax  = ranges[0,1] 
            G0Min = ranges[1,0]; G0Max = ranges[1,1]
                        
        #for grdSubList in self.grds:
        #    for grd in grdSubList:
        #        # setting up the 2D meshes for the grids and initializing them
        #        grd.fill(0)
        #        grd.setup( [[nMin, nMax], [G0Min, G0Max]] )
        
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
    def plotGrid(self, resGrids, lgammaMechSec, radexParms):

        self.pltGmSec   = lgammaMechSec
        self.radexParms = radexParms

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
        nMin  = 0.0; nMax  = 6.0;
        G0Min = 0.0; G0Max = 6.0;
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
                
        #for grdSubList in self.grds:
        #    for grd in grdSubList:
        #        # setting up the 2D meshes for the grids and initializing them
        #        grd.fill(0)
        #        grd.setup(  )

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
            self.computeSurfaceTemperatureGrid( indsThisSec )
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
            self.computeHeatingCoolingGrid(indsThisSec, 0, 'heating', 'photo' )
            ###
            im11 = self.grds[1][1].imshow(interpolation='nearest')
            cbar11 = pyl.colorbar(im11, cax=self.grdsCbarAxs[1][1], ax=pyl.gca(), orientation = 'horizontal')            
            cbar11.set_ticks( cbarTickValues )
            self.grds[1][1].plotContour( levels = cbarTickValues )
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
