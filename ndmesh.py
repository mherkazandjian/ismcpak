import numpy as np
import pylab as pyl

class ndmesh( np.ndarray ):
    
    def __init__(self, *args, **kwrds):
        
        self.ranges = None  # ranges for each dimension ( a list [ min1, max1, min2, max2...] )
        self.ldim   = None  # ndarray, length of each dimenstion max1-min1, max2-min2
        self.dl     = None  # size of each bin in each dim
        self.nbins  = None  # number of bins in each coord
        
        self.cntrds = None  # centroids of each bin in each dim
        self.scrds  = None  # starting coords of each bin in each dim
        self.ecrds  = None  # ending coords of each bin in each dim
    
    # set the ranges in each dimesnsion
    def set_ranges(self, ranges):
        
        if len(ranges) != 2*self.ndim:
            errStr = "number of ranges does not agree with the dimensions of the array"
            raise ValueError(errStr)
        else:
            self.ranges = np.float64(ranges)
        
    def get_ranges(self):
        
        if self.ranges == None:
            errStr = "attribute not set"
            raise ValueError(errStr)
        else:
            return self.ranges
    
    def centroids(self, **kwrds):
        
        if 'ranges' in kwrds:
            self.set_ranges( kwrds['ranges'])
                    
        if self.ranges == None:
            errStr = "ranges not set, cannot compute centroids"
            raise ValueError(errStr)
        else:
            
            ndim  = self.ndim
            shape = np.float64(self.shape)
            
            self.ldim  = np.ndarray( ndim, dtype = np.float64 )
            self.dl    = np.ndarray( ndim, dtype = np.float64 )            
            

            self.cntrds = np.ndarray( shape, dtype = np.float64 )
            self.scrds  = np.ndarray( shape, dtype = np.float64 )
            self.ecrds  = np.ndarray( shape, dtype = np.float64 )
            
            # computing the length in each dim and the bin size and the number of bins                        
            for dimInd in np.arange(0, ndim, 1):
                
                mn = self.ranges[2*dimInd]
                mx = self.ranges[2*dimInd + 1]
                
                self.ldim[dimInd] = mx - mn
                self.dl[dimInd]   = (mx - mn)/shape[dimInd]
                
                # computing bins start,center,end positions
                ##### CONTINUE IMPLEMENTING THESE
                
                dimInd += 1
                
            
            #dl = np.
            #print self.shape
            #print self.ndim
        
            #asdasdads
            #self.centroids =  np.ndarray( self.shape, dtype = np.float64 )
            #print self.centroids
        
    def misc(self):
        print self.shape
    