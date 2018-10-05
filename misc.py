import numpy
import time
import logging
import sys
import pylab
import collections
import pickle
    

def scale(array, bottomTmp, topTmp, minm, maxm, *args1, **args2):
    """Copied from the original IDL implementation scale.pro

       :retuens: a scaled array,this function modifies the entered array (call by
        reference),it retuenes the same array scaled between bottom and top
        such that the maximum value is top and the minimum value is bottom
       :param ARRAY: can have an arbitrary dimension (mxn)
       :param BOTTOM: minimum value at which the data will be scaled (ie the minimum
         returned value will be not smaller than 'bottom')
       :param TOP: maximum value at which the data will be scaled (ie the maximum
         returned value will be not larger than 'top')
       :param MAXIMUM: Is the maximum value the entries of "array" can take.
       :param MINIMUM: Is the mainimum value the entries of "array" can take.
       :param integer: ( keyword = True|False(default) ) if this keyword is set the then the entries of array are scaled
         between bot and top-1 and returned as long integers
       :param eps: ( keyword = double ) this number is decremented from the scaled numbers such that
         the maximum is not reached. i.e the numbers are scaled to TOP-1
                                                                            
       example:
    
       .. code-block:: python

          x = [20,30,40,50, 80,81,82,83,84,85]
          mn = numpy.min(x)
          mx = numpy.max(x)
         
          # scaling between 0 and 1, assuming the interval of the data [0,100]
          y = scale(x, 0.0, 1.0, 0.0, 100.0)
          # the output wud be
          #      y = 0.2 ,  0.3 ,  0.4 ,  0.5 ,  0.8 ,  0.81,  0.82,  0.83,  0.84,  0.85

          # scaling between 0 and 1, assuming the interval of the data [mn,mx]
          y = scale(x, 0.0, 1.0, mn, mx)
          # the output wud be
          #      y =  0,  0.153, 0.307, 0.461, 0.923, 0.938, 0.953, 0.969, 0.984, 1.  
        
        
       Changelog: 2009-03-06 : changed the keywords min and max to regular parameters.
    """

    bottom = numpy.double(bottomTmp); top = numpy.double(topTmp);
    mn = numpy.double(minm); mx = numpy.double(maxm)
    factor1 = mx  - mn
    factor2 = top - bottom
    xs = ( (array - mn)/factor1 )*factor2 + bottom
    
    if 'integer' in args2:
        if 'eps' not in args2:
            eps = 1e-15
        else:
            eps = args2['eps']
    
        xs = numpy.int64( (1.0-eps)*xs )

    return xs


def fetchNestedDtypeValue(var, namesLst):
    """Returns the value of a nested dtype given the names to look for as a list. 
       For example, to return the value of var['xx']['yy']['zz'], u can use :
        
            value = fetchNestedDtypeValue[ var, ['xx', 'yy', 'zz']
        
       an exception is raised if any of the names is not found
        
       .. note:: it might be cool to implement a feature where give zz, 
           without xx and yy it can return the value of zz (assuming it
           is uniqe, i.e no other name is zz). Hint: a tree lokup maybe??
        
       see testDtypeUtils.py
    """

    currVar = var
    for name in namesLst:
        #print 'looking for current name %s' % name
        currVarNames = currVar.dtype.names
        #print 'names in the current var are : ', currVarNames
        #print '-----------------------------------------------'
    
        if name in currVarNames:
            #print 'found the name *%s*' % name
            currVar = currVar[name]
        else:
            raise NameError('name *%s* not found in the variable "currVar"' % name)
        
    return currVar


def default_logger():
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

    return logger
