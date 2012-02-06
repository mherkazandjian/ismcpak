"""
 this class generates and manipulates archives of meshes
  
 methods :
   arxvHdrFormat()

 instance variables :
    self.infoAll  : a numpy array of dtype arxvHdrDtype (see below) which contains 
                    the info (headers) about all the meshes
                    each entry in this array contains two things, information about the
                    mesh and the parameters of the mesh. For example the elements
                    x in self.infoAll[x] has the following contents : 
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
                      
    self.meshes  : a list of all the meshes of 'mesh' dtypes 
   
 by default, the prefix name of the individual mesh files is assumed to be
 mesh.dat-id-xxxxxx
 
 the mesh database files are assumed to have the same prefix, for example, if the
 database file provided is foo, then this routine will assume (or write) 
    foo.info
    foo.db
 for now the database can be stored into a single file and not split into
 multiple files.  
 
  the database for the meshes is constructed of two binary files:
   foo.info : holding all the information about the meshes, locations and parameters
   foo.db   : holds the data for all the meshes
   
    the info file has the following structure :
    
    - version number    : int32 int32 int32
    - nMeshes           : int32 
    - meshes info array : meshInfoArrayFormat (see also the method)
       ( nMeshes, 1)          .
                              .
                              .  

            the entries of the dtype are 
                - mesh number     ( int64 ) :
                - data file index ( int64 ) :
                     i.e in which .db file the mesh is located
                - offset ( int64 ) :
                     the offset in bytes from the beginning of the file where the mesh
                     is located
                - nSteps ( int64 ):
                     the number of slabs in the mesh
                - nSpecs ( int64 ) :
                    the number of species in the mesh                                   .
                    
    the .db files have the following structure:
    
     - mesh_1 ( mesh dtype )
     - checkNum_1 ( int64 )
     - mesh_2 ( mesh dtype )
     - checkNum_2 ( int64 )
              .
              .
              .
    for a mesh number i, the checkNum_i should be the same as the i^th entry 
    in the info array offset...i.e chechNum = infoAll[i]['info'][2] 
    
    see test_writeReadArxv.py for an example
"""

import os
import glob
import numpy as np
from matplotlib.widgets import Button

from mesh import *

class meshArxv():
    
    def __init__(self):
        
        # instance variables
        self.nMeshes    = None
        self.ver        = None
        self.meshes     = None
        self.infoAll    = None
        self.nDbFiles   = None
        self.chemNet    = None
        
        # variables used for plotting
        self.pltGmSec = None
        
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
            
            nSteps = self.infoAll[i]['info'][3]
            nSpecs = self.infoAll[i]['info'][4]
        
            thisMeshDtype = mDummy.constructMeshDtype(nSpecs, nSteps)
            thisMeshData = np.fromfile( dbDataFObj, dtype = thisMeshDtype, count = 1)

            self.meshes.append( thisMeshData )                
            checkNum = np.fromfile( dbDataFObj, dtype = np.int64, count = 1)
            
            if checkNum != self.infoAll[i]['info'][2]:
                strng  = 'Error : checkpoint numbers do not match : database may be corrupt.\n'
                strng += 'read = %d, expected = %d' % (checkNum, self.infoAll[i]['info'][2])
                strng += ' when reading mesh %d' % i                
                raise NameError(strng)
                    
        print 'read successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)
        
    # the format for the archive header from which the dtype will be constructed
    # specifying the info and paramters for each line in the header        
    def arxvHdrFormat(self):
        
        # defining the version of the database
        self.ver = np.zeros([3], dtype = np.int32)
        self.ver[0] = 0
        self.ver[1] = 0
        self.ver[2] = 1

        return [ 
                  ('info' , np.int64  , 5),
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

            x2 = np.array([np.log10(self.meshes[i]['hdr']['G0'][0]),
                           np.log10(self.meshes[i]['hdr']['nGas'][0]),
                           np.log10(self.meshes[i]['hdr']['gammaMech'][0])])

            xdv = x2 - x1
            diff +=  np.sqrt( np.dot(xdv,xdv) )
            
        if diff == 0.0:
            print 'archive integrity test passed'
        else:
            strng = 'archive integrity test failed. database file may be corrupt' 
            raise NameError(strng)
        
    def setChemicalNetwork(self, chemNet):
        self.chemNet = chemNet
        
    
    def plotGrid(self):
        
        self.pltGmSec = -24.0
        
        fig1 = pyl.figure(1, figsize=(6, 6))
        axs1 = fig1.add_subplot(111)
        msh = mesh()

        lG0All   = np.log10(self.infoAll['parms'][:,0])
        lnGasAll = np.log10(self.infoAll['parms'][:,1])
        lgmAll   = np.log10(self.infoAll['parms'][:,2])
        
        
        def plotThisSec():
            pyl.figure(1)
            pyl.subplot(111)
            pyl.hold(False)
            indsThisSec = ( np.fabs(lgmAll - self.pltGmSec) < 1e-6 ).nonzero()
            pyl.plot( lnGasAll[indsThisSec], lG0All[indsThisSec] , 'bo' )
            pyl.xlim( xmin = -1, xmax = 7.0)
            pyl.ylim( ymin = -1, ymax = 7.0)
            pyl.hold(True)
            pyl.draw()

        plotThisSec()
        
        def nextSec(event):
            self.pltGmSec += 1.0
            print self.pltGmSec
            plotThisSec()
    
        def prevSec(event):
            self.pltGmSec -= 1.0
            print self.pltGmSec
            plotThisSec()

        axNext = pyl.axes([0.81, 0.02, 0.1, 0.05])
        bNext  = Button(axNext, 'Next')
        bNext.on_clicked(nextSec)
        axPrev = pyl.axes([0.61, 0.02, 0.1, 0.05])
        bPrev = Button(axPrev, 'Prev')
        bPrev.on_clicked(prevSec)
        
        fig2, axs2 = pyl.subplots(2, 2, sharex=True, sharey=False, figsize=(8,8))
        
        def onB1Down(event):
            # get the x and y coords, flip y from top to bottom
            b      = event.button 
            x, y   = event.x, event.y
            xd, yd = event.xdata, event.ydata
        
            if event.button==1:
                if event.inaxes is not None:
            
                    l2Distance  = np.sqrt( (yd - lG0All )**2 + (xd - lnGasAll)**2 + (self.pltGmSec - lgmAll)**2 )
                    rMin = min(l2Distance)
                    indMin = l2Distance.argmin()
                    msh.setData( self.meshes[indMin][0] )
             
                    pyl.figure(1)
                    pyl.subplot(111)
                    pyl.hold(False)
                    plotThisSec()
                    pyl.hold(True)
                    pyl.plot( lnGasAll[indMin], lG0All[indMin], 'ro')
                    pyl.title('$\log_{10} n_{gas} = $ %4.2f $\log_{10}  = G_0$ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n' % (np.log10(msh.data['hdr']['nGas']), np.log10(msh.data['hdr']['G0']), np.log10(msh.data['hdr']['gammaMech'])))
                    pyl.xlabel('$log_{10} n_{gas}$')
                    pyl.ylabel('$log_{10} n_{Gas}$')
                    pyl.draw()

                    
                    pyl.figure(2)            
                    msh.setFigure(fig2, axs2)
                    msh.plot(self.chemNet)
                    
            pyl.draw()

        cid = fig1.canvas.mpl_connect('button_press_event', onB1Down)

        pyl.figure(1)
        pyl.show()
        print 'browing data....'
