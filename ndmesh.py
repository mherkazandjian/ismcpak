import numpy as np
import pylab as pyl
from matplotlib.artist import setp

class ndmesh(np.ndarray):
    """This class provides utilities for manipulating mgrids such as
     plotting, filling data...etc.. FOR NOW the plotting stuff work 
     for 2D only it works only with 2D meshes.
     methods :
    
    NOTES : see http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    for proper subclassing of the class
      
    Inherits all from ndarray and the attributes defined here. Below is a 
    demonstration of the data structure of a 2D ndarray defined as:
     
    .. code-block:: python
        
        x = ndmesh((2,5), ranges=[[0,1], [0,1]] )  #the range values here can be arbitrary  
        
        #when printed it looks like this :   
        #two rows and 5 columns
        #   |-----|-----|-----|-----|-----|
        #   | 0,0 | 0,1 | 0,2 | 0,3 | 0,4 |  x[0,:] = x[0]
        #   |-----|-----|-----|-----|-----|  
        #   | 1,0 | 1,1 | 1,2 | 1,3 | 1,4 |  x[1,:] = x[1]
        #   |-----|-----|-----|-----|-----|
        #   x[:,0]      x[:,0]      x[:,0]
        #   
        #since this was implemented mainly for 2D meshes which might be needed to
        #be plotted, with x increasing from left to right and y increasing from bottom
        #to top. the method plotGrid() and plotCntrd() would produce
        # 
        #(0,n1-1)             (n1, n2)
        #   |--------|--------|
        #   |        |        |
        #   | x[0,4] | x[1,4] |
        #   |        |        |
        #   |--------|--------|
        #   |        |        |
        #   | x[0,3] | x[1,3] |
        #   |        |        |
        #   |--------|--------|
        #   |        |        |
        #   | x[0,2] | x[1,2] |
        #   |        |        |
        #   |--------|--------|
        #   |        |        |
        #   | x[0,1] | x[1,1] |
        #   |        |        |
        #   |--------|--------|
        #   |        |        |
        #   | x[0,0] | x[1,0] |
        #   |        |        |
        #   |--------|--------|
        #(0,0)  x[0]     x[1]  (n1-1,0)
        #
        #          Fig 1
                   
    ndmesh.imshow() plots the 2D color coded data in the same orientation as
    Fig 1. which flips the axes and the origin would be at the top which is what
    pylab.imshow(x) does. In getting the orientation the same as Fig 1, x.transpose
    is should be used in ndmesh.imshow in adition to having the axis from bottom 
    to top instead of the default used in pylab.imshow() which is from top to bottom. 

    An example : creating a 3D mesh.
    
    .. code-block:: python
    
        nx, ny, nz = 10, 10, 10
        xr, yr, zr = [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]
        msh3D = ndmesh((nx, ny, nz),  ranges = [xr, yr, zr])
        
    See test_ndmesh.py and test_ndmesh2.py for a demonstration
        
    author : mher 2012-03-11
    """
    
    def __new__(cls, *args, **kwrds):
        newArgs = ['ranges', 'fill']
        
        for arg in newArgs: # removing keyword arguments not compatible with ndarray
            if arg in kwrds:
                kwrds.pop(arg)
        return np.ndarray.__new__(cls, *args, **kwrds)
    
    def __init__(self, *args, **kwrds):
        
        self.ranges  = None #: numpy.ndarray : [ [min max] [min max] ... [min max] ]
        self.cntrd   = None #: nunmpy.mgrid  : ( nDim, n1, n2, ..., nn )
        self.spos    = None #: numpy.mgrid   : ( nDim, n1, n2, ..., nn )
        self.epos    = None #: numpy.mgrid   : ( nDim, n1, n2, ..., nn )
        self.dl      = None #: numpy.ndarray : [ (max-min)/n1, (max-min)/n2,...., (max-min)/nn ], increment size in each dimension
        self.l       = None #: numpy.ndarray : [  max-min, max-min,...., max-min ], ength in each dimenion
        self.nCells  = None #: numpy.float   : float64 : n1*n2*....*nn
        self.image   = None #: image  : object reference returned by pylab.imshow(), bject returned by self.imshow()
        self.contour = None #: contour: object reference returned by pylab.contour()

        #self.data   = self
        
        if 'ranges' in kwrds:
            self.setup( kwrds['ranges'])
        if 'fill' in kwrds:
            self.fill(kwrds['fill'])
            
    # setup the locations of the cells as mgid objects
    def setup(self, ranges):
        
        # setting the ranges
        self.ranges = np.array(ranges, dtype = np.float64)
        #print 'ranges = ', self.ranges
        
        if len(self.ranges) != len(self.shape):
            raise ValueError("""The input ranges do not correspond to the same dimensions
                                of the ndmesh to be created.""")
            
        shape = self.shape
        #print 'shape = ', shape
        #print 'dim   = ', self.ndim

        # setting the increment size in each dimension ( xMax - xMin )/dx
        self.dl = (self.ranges[:,1] - self.ranges[:,0])/np.array(self.shape, dtype = np.float64 )
        #print 'dl = ', self.dl 
        #print self.dl

        # defining the spos grid in terms of the indicies
        evalStr = "np.mgrid["        
        for i, rangeDim_i in list(enumerate(self.ranges)):
            evalStr += "%d:%d:1.0" % (0, self.shape[i])
            if i != self.ndim-1:
                evalStr += ','
        evalStr += ']'
        
        #the grid with the starting positions defined as the indicies (this ensures the desried shape of the grid)
        spos_norm = eval(evalStr)
        
        #the actual starting coordinates of the cells
        self.spos = ((spos_norm.T * self.dl + self.ranges[:,0]).T)
        
        # computing the epos grid
        self.epos = (self.spos.T + self.dl).T
        # computing the centroids
        self.cntrd = (self.spos + self.epos)/2.0
        # computing the total number of cells
        self.nCells = np.prod(self.shape)

    def plotGrid(self):
        """2D only : plots the grid cells boundaries on the current axis"""
        
        gs = np.resize(self.spos , (2, self.nCells))
        ge = np.resize(self.epos , (2, self.nCells))

        allLines = ()

        for i in np.arange( self.nCells ):
            poss = gs[:,i] # [xs, ys, zs...]
            pose = ge[:,i] # [xe, ye, ze...]
    
            # extracting the starting and ending coords as individual variables for
            # clarity
            xs = poss[0]; xe = pose[0]
            ys = poss[1]; ye = pose[1]
            # defining the path for each reactangle
            xPathRect = [ xs, xe, xe, xs, xs ]
            yPathRect = [ ys, ys, ye, ye, ys ]
            # plotting the rectangle
            lines = pyl.plot( xPathRect , yPathRect )
            allLines = allLines + (lines,)
    
            setp(allLines, linewidth=0.5, color='r')
            
        self.pltGrd = allLines
        
    def plotCntrd(self):
        """2D only : plots the centroid son the current axis"""

        gc = np.resize(self.cntrd, (2, self.nCells))
        self.pltCnt = pyl.plot( gc[0,:], gc[1,:], '.')

    # shows the color coded and scaled data of the object
    def imshow(self, *args, **kwrds):
        """2D only : plots the color coded image of the data, also inherits all the arguments 
         and keywords of pylab. imshow() """
        if self.image != None: # deleting the previous displayed image 
            del self.image

        self.image = pyl.imshow( self.transpose(), origin = 'bottom', extent=[ self.ranges[0,0], self.ranges[0,1], self.ranges[1,0], self.ranges[1,1] ], aspect='auto', *args, **kwrds)
        return self.image

    def imUpdate(self, *args, **kwrds):
        """2D only : updates the displayed image with the current data of self"""

        if self.image == None:
            errStr = 'no image is set, self.imshow was not called before'
            raise NameError(errStr)
        else:
            self.image.set_data(self.transpose() )

    def plotContour(self, *args, **kwrds):
        """2D only : draws the contour plot of the 2D data. shows the color coded and 
        scaled data of the object"""
        self.contour = pyl.contour( self.transpose(), extent=[ self.ranges[0,0], self.ranges[0,1], self.ranges[1,0], self.ranges[1,1] ], *args, **kwrds)
        return self.contour

    def getCntrd(self):
        """returns the centroids mgrid object"""
        return self.cntrd

    def __del__(self):
        """deleting the attributes
        .. todo:: implement a good cleaning.
        """
        pass
        #del self.ranges
        #del self.cntrd
        #del self.spos
        #del self.epos
        #del self.dl
        #del self.l
        #del self.nCells
        #del self.image
        #del self.contour