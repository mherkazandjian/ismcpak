import time
import numpy
from scipy.interpolate import griddata as gdata

def sectioned_4D_interpolator(xyzt, v, scipy_interpolator):
    """Returns an interpolator function which makes use of sub-interpolation function 
    for ranges in the 3rd and 4th dinemnsion. xyzt is an array whose columns are : x, y, z, t.
    
    .. warning:: make sure the input intervals in z and t are within the bounds of the data.
    
    .. todo:: add a keyword which takes uniform sections in both as (nz, nt) and construct a grid
     of intervals automatically
    """
    
    #sections in gm within which saperate interpolation functions will be built.
    # .. note:: make sure there are enough points in the database within each section
    intervals_z = [ 
                     [-50.0, -35.0],  
                     [-35.0, -34.0], [-34.0, -33.0], [-33.0, -32.0], [-32.0, -31.0], [-31.0, -30.0],    
                     [-30.0, -29.0], [-29.0, -28.0], [-28.0, -27.0], [-27.0, -26.0], [-26.0, -25.0],    
                     [-25.0, -24.0], [-24.0, -23.0], [-23.0, -22.0], [-22.0, -21.0], [-21.0, -20.0],    
                     #[-20.0, -13.0],
                   ]
    ghost_z = 1e-6 #1.1

    #sections in Av within which saperate interpolation functions will be built.
    # .. note:: make sure there are enough points in the database within each section
    intervals_t = [  [0.01, 1.0],  [1.0, 2.0],   [2.0, 3.0],
                     [3.0, 4.0],   [4.0, 5.0],   [5.0, 6.0],   [6.0, 7.0],   [7.0, 8.0],   [8.0, 9.0], [9.0, 10.0],
                     [10.0, 12.0], [12.0, 14.0], [14.0, 16.0], [16.0, 18.0], [18.0, 20.0], 
                     [20.0, 22.0], [22.0, 24.0], [24.0, 26.0], [26.0, 28.0], #, [28.0, 30.0],
                   ] 

    ghost_t = 1e-6 #1.1
    
    #splitting data into sections specified by intervals_z and constructing an interpolation
    #function for each interval (the intevals include an extra dex in each bound of the specified
    #bound to make sure interpolation is done correctly...i.e the scipy interpolator doesnt return
    # nan's)
    
    z, t = xyzt[:, 2], xyzt[:, 3]

    #numpy 2D array whose dimesnions are the number of intervals in gm and Av (in each dimensions)
    #for instance all_fInterps[0,3] corresponds to an interpolation function within the range 
    #intervals_z[0] and intervals_t[3]    
    all_fInterps = numpy.ndarray((len(intervals_z), len(intervals_t)), dtype='object')
    
    #lookup points which lie in each section of gm and Av and get the points in there
    #and construct an interpolation function using those points
    for i, interval_z in enumerate(intervals_z):
        
        interval_z_low  = interval_z[0] - ghost_z
        interval_z_high = interval_z[1] + ghost_z        
        
        for j, interval_t in enumerate(intervals_t):

            interval_t_low  = interval_t[0] - ghost_t
            interval_t_high = interval_t[1] + ghost_t

            #finding the points in the DB in this section
            inds_in_this_interval = numpy.where(
                                               (z >= interval_z_low)*(z <= interval_z_high)*
                                               (t >= interval_t_low)*(t <= interval_t_high)
                                              )[0]
                               
            if inds_in_this_interval.size != 0:
                
                #if enough points are found, make the interpolation funcion
                xyzt_this_box = xyzt[inds_in_this_interval, :]
                v_in_this_box = v[inds_in_this_interval]
                
                z_min, z_max = numpy.min(xyzt_this_box[:,2]), numpy.max(xyzt_this_box[:,2])
                t_min, t_max = numpy.min(xyzt_this_box[:,3]), numpy.max(xyzt_this_box[:,3])

                print '(i,j) = (%02d, %02d) %06d DB points in box gm = [%.2f, %.2f ], Av = [%04.2f,%04.2f ] ' % (i, j, inds_in_this_interval.size, interval_z_low, interval_z_high, interval_t_low, interval_t_high)
                print '\t\t\t ranges_i = [%.2f %.2f], ranges_j = [%.2f %.2f]' % (z_min, z_max, t_min, t_max)
                t0 = time.time()
                fInterp_this_interval = scipy_interpolator(xyzt_this_box, v_in_this_box)
                print '\t\t\t interpolation function constructed in %.2e sec' % (time.time()-t0)
                all_fInterps[i,j] = fInterp_this_interval
            
            else:                
                print 'no points found in box gm = [%f, %f ], Av = [%f,%f ]' % (interval_z_low, interval_z_high, interval_t_low, interval_t_high)
        #
        
    def fInterp_from_section_interpolator(data):
        """devides data (values at which we would like to get interpolated values) into sections of 
        gm and Av (3rd and 4th colomn) and returns the interpolated values by passing them to the 
        sub-interpolation function in all_fInterps.
        """
        
        nPts = data.shape[0]
        
        #array which will hold the interpolated values from all the sections
        all_values = numpy.zeros(nPts, dtype=numpy.float64)
        
        data_z, data_t = data[:, 2], data[:, 3]
        inds_orig = numpy.arange(nPts) 
        
        #interpolating from the sectioned interpolation function. Lookup the points
        #which are in the _all_Avbounds of each interpolation function and interpolate the
        #values and set them in the array to be returned.
        for i, interval_z in enumerate(intervals_z):
            
            for j, interval_t in enumerate(intervals_t):
            
                inds_in_box = numpy.where( 
                                           (data_z >= interval_z[0])*
                                           (data_z <  interval_z[1])* 
                                           (data_t >= interval_t[0])*
                                           (data_t <  interval_t[1]) 
                                         )[0]
                                                                                
                if inds_in_box.size == 0:
                    continue
                
                data_in_this_interval = data[inds_in_box, :]

                print '(i,j) = (%02d, %02d) %06d points in box gm = [%.2f, %.2f ], Av = [%4.1f,%4.1f ] ' % (i, j, inds_in_box.size, interval_z[0], interval_z[1], interval_t[0], interval_t[1]),                
                t0 = time.time()
                v = all_fInterps[i,j](data_in_this_interval)
                #v = 0.0
                dt = time.time()-t0
            
                #setting the values in the array to be returned
                all_values[inds_orig[inds_in_box]] = v
                
                print 'interpolated in %5.2f sec, %05d points at %e points/sec' % (dt, inds_in_box.size, inds_in_box.size/dt)

        return all_values
    #
    
    return fInterp_from_section_interpolator

class interpolator_sectioned(object):
    """Returns an interpolator function which makes use of sub-interpolation function 
    for ranges in the 3rd and 4th dinemnsion. xyzt is an array whose columns are : x, y, z, t.
    
    .. note:: this only works for 4D interpoaltion functions i.e x,y,z,t and v.
     The first four are what refer to as the 4D variables and v are the values 
     which are interpolated.
    
    .. warning:: make sure the input intervals in z and t are within the bounds of the data.
    
    .. todo:: add a keyword which takes uniform sections in both as (nz, nt) and construct a grid
     of intervals automatically
    """
    
    def __init__(self, xyzt, v, intervals_z=None, ghost_z=None,
                 intervals_t=None, ghost_t=None,
                 scipy_interpolator=None, load_from=False, verbose=False):

        self.intervals_z = intervals_z
        self.ghost_z = ghost_z
        self.intervals_t = intervals_t
        self.ghost_t = ghost_t
        self.all_fInterps = None
        #-----------------------        
                 
        if type(load_from) == type(''):
            self.load(load_from)
        else:
            #splitting data into sections specified by intervals_z and constructing 
            #an interpolation function for each interval (the intevals include an 
            #extra dex in each bound of the specified bound to make sure interpolation 
            #is done correctly...i.e the scipy interpolator doesnt return # nan's)
            z, t = xyzt[:, 2], xyzt[:, 3]
        
            #numpy 2D array whose dimesnions are the number of intervals in gm and Av 
            #(in each dimensions) for instance all_fInterps[0,3] corresponds to an 
            #interpolation function within the range intervals_z[0] and intervals_t[3]    
            self.all_fInterps = numpy.ndarray((len(intervals_z), 
                                               len(intervals_t)), dtype='object')
            
            #lookup points which lie in each section of gm and Av and get the points 
            #in there and construct an interpolation function using those points
            for i, interval_z in enumerate(intervals_z):
                
                interval_z_low  = interval_z[0] - ghost_z
                interval_z_high = interval_z[1] + ghost_z        
                
                for j, interval_t in enumerate(intervals_t):
        
                    interval_t_low  = interval_t[0] - ghost_t
                    interval_t_high = interval_t[1] + ghost_t
        
                    #finding the points in the DB in this section
                    inds_in_this_interval = numpy.where(
                                                       (z >= interval_z_low)*(z <= interval_z_high)*
                                                       (t >= interval_t_low)*(t <= interval_t_high)
                                                      )[0]
                                       
                    if inds_in_this_interval.size != 0:
                        
                        #if enough points are found, make the interpolation funcion
                        xyzt_this_box = xyzt[inds_in_this_interval, :]
                        v_in_this_box = v[inds_in_this_interval]
                        
                        z_min, z_max = numpy.min(xyzt_this_box[:,2]), numpy.max(xyzt_this_box[:,2])
                        t_min, t_max = numpy.min(xyzt_this_box[:,3]), numpy.max(xyzt_this_box[:,3])
                        
                        if verbose == True:
                            print '(i,j) = (%02d, %02d) %06d DB points in box z = [%.2f, %.2f ], t = [%04.2f,%04.2f ] ' % (i, j, inds_in_this_interval.size, interval_z_low, interval_z_high, interval_t_low, interval_t_high)
                            print '\t\t\t ranges_i = [%.2f %.2f], ranges_j = [%.2f %.2f]' % (z_min, z_max, t_min, t_max)
                            
                        t0 = time.time()
                        fInterp_this_interval = scipy_interpolator(xyzt_this_box, v_in_this_box)
                        dt = time.time()-t0
                        
                        if verbose == True:
                            print '\t\t\t interpolation function constructed in %.2e sec' % dt
                             
                        self.all_fInterps[i,j] = fInterp_this_interval                
                    else:                
                        print 'no points found in box gm = [%f, %f ], Av = [%f,%f ]' % (interval_z_low, interval_z_high, interval_t_low, interval_t_high)
                #

    def get(self, data, verbose=False):
        """devides data (values at which we would like to get interpolated values) into sections of 
        gm and Av (3rd and 4th colomn) and returns the interpolated values by passing them to the 
        sub-interpolation function in all_fInterps.
        
        .. todo:: add a check which when the number of points is low, then an efficient method is 
         used in determining the box where the point is and the interpolation function in that
         box is used, or use some clustering technique, meshing, gridding the points and using 
         reverse indicies instead of looping over all the boxes and calling where every time
        """
        
        nPts = data.shape[0]
        
        #array which will hold the interpolated values from all the sections
        all_values = numpy.zeros(nPts, dtype=numpy.float64)
        
        data_z, data_t = data[:, 2], data[:, 3]
        inds_orig = numpy.arange(nPts) 
        
        #interpolating from the sectioned interpolation function. Lookup the 
        #points which are in the _all_Avbounds of each interpolation function 
        #and interpolate the values and set them in the array to be returned.
        for i, interval_z in enumerate(self.intervals_z):
            
            for j, interval_t in enumerate(self.intervals_t):
            
                inds_in_box = numpy.where(
                                           (data_z >= interval_z[0])*
                                           (data_z <  interval_z[1])* 
                                           (data_t >= interval_t[0])*
                                           (data_t <  interval_t[1]) 
                                         )[0]
                
                if inds_in_box.size == 0:
                    continue
                
                data_in_this_interval = data[inds_in_box, :]

                if verbose == True:
                    print '(i,j) = (%02d, %02d) %06d points in box z = [%.2f, %.2f ], t = [%4.1f,%4.1f ] ' % (i, j, 
                                                                                                                  inds_in_box.size, 
                                                                                                                  interval_z[0], interval_z[1], 
                                                                                                                  interval_t[0], interval_t[1]),
                t0 = time.time()
                v = self.all_fInterps[i,j](data_in_this_interval)
                dt = time.time()-t0
            
                #setting the values in the array to be returned
                all_values[inds_orig[inds_in_box]] = v

                if verbose == True:
                    print 'interpolated in %5.2f sec, %05d points at %e points/sec' % (dt, inds_in_box.size, inds_in_box.size/dt)
        #
        
        return all_values

    def save(self, fpath):
        '''saves the object data needed to set it up again to a file'''
        
        numpy.savez_compressed(fpath, 
                               method='sectioned',
                               data=[self.intervals_z, 
                                     self.ghost_z, 
                                     self.intervals_t,
                                     self.ghost_t,
                                     self.all_fInterps]
                               )
        print 'saved interpolation function to : %s' % fpath
        
    
    def load(self, fpath):
        '''loads saved object data needed to set it up again to a file'''
        
        # loading the data
        obj = numpy.load(fpath)['data']
        
        # setting the data to the attributes
        self.intervals_z = obj[0]
        self.ghost_z = obj[1]
        self.intervals_t = obj[2]
        self.ghost_t = obj[3]
        self.all_fInterps = obj[4]

def fill_nans(z):
    
    res = z.shape
    
    x, y = numpy.meshgrid(numpy.arange(0, res[0]), numpy.arange(0, res[1]))
    
    points = numpy.vstack((x.flatten(), y.flatten())).T
    values = z.flatten()
    
    inds = numpy.where( numpy.isnan(values) == False)[0]
    
    z_filled = gdata(points[inds,:], values[inds], (x, y), method='cubic')
    
    return z_filled
#