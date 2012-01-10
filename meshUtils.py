"""
 this class generates and manipulates archives of meshes
  
 methods :
   arxvHdrFormat()

   
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

from mesh import *

class meshArxv():
    
    def __init__(self):
        
        self.nMeshes    = None
        self.ver        = None
        self.meshes     = None
        self.infoAll    = None
        self.nDbFiles   = None
         
    # read all the meshes files in the dir and construct the
    # database
    def construct(self, dirName ):

        meshNamePrefix = 'mesh.dat-id-'        

        # defining the files objects for the database fiels        
        dbInfoFObj = file(dirName + 'meshes.db.info', 'wb')
        dbDataFObj = file(dirName + 'meshes.db'  , 'wb')

        
        # getting the names of the meshes in that dir
        files = []
#        for infile in glob.glob( os.path.join(dirName, meshNamePrefix + '*0000*') ):
        for infile in glob.glob( os.path.join(dirName, 'meshes/'+meshNamePrefix + '*') ):
            files.append(infile)
        #print files
        
        # setting variable for the 
        self.nMeshes = np.zeros( 1, dtype = np.int32 )
        self.nMeshes[0] = len(files)
        print 'found %s meshes' % (self.nMeshes)

        # defining the array holding the info about the meshes
        # and their location in the database ..etc..
        arxvHdrDtype = np.dtype( self.arxvHdrFormat() )
        infoAll = np.ndarray( self.nMeshes, dtype = arxvHdrDtype )

        # reading the individual mesh files and writing the files into .db file(s)
        for i in np.arange( self.nMeshes ):
    
            fName = files[i]
            #print fName
    
            # reading a mesh
            m = mesh(fName)
            mData = m.getData()

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
        #print self.ver
    
        self.nMeshes = np.fromfile( dbInfoFObj, dtype = np.int32, count = 1)
        self.nMeshes = self.nMeshes[0]

        arxvHdrDtype = np.dtype( self.arxvHdrFormat() )
        self.infoAll = np.fromfile( dbInfoFObj, dtype = arxvHdrDtype, count = self.nMeshes)

        # reading the meshes in database into a list 
        mDummy = mesh()
        self.meshes = []
        
        for i in np.arange(self.nMeshes):
                
            nSteps = self.infoAll[i]['info'][3]
            nSpecs = self.infoAll[i]['info'][4]
    
            thisMeshDtype = mDummy.constructMeshDtype(nSpecs, nSteps)
            thisMeshData = np.fromfile( dbDataFObj, dtype = thisMeshDtype, count = 1)
            self.meshes.append( thisMeshData )                
            checkNum = np.fromfile( dbDataFObj, dtype = np.int64, count = 1)
                
            if checkNum != self.infoAll[i]['info'][2]+1:
                str = 'Error : checkpoint numbers do not match : database may be corrupt.'                 
                raise NameError(str)
                    
        print 'read successfully database files : \n  %s\n  %s' % (dbInfoFObj.name, dbDataFObj.name)
        
    # the format for the archive header from which the dtype will be constructed        
    def arxvHdrFormat(self):
        
        self.ver = np.zeros([3], dtype = np.int32)
        self.ver[0] = 0
        self.ver[1] = 0
        self.ver[2] = 1

        return [ 
                  ('info' , np.int64  , 5),
                  ('parms', np.float64, 5)
               ]
