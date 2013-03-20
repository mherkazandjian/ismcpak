# looks up and restors a 2D grided of data from a file name
#------------------------------------------------------------------------------------
import numpy
import pickle
import os

def fetchRadexGrid( dirname = None, specStr = None, gmechSec = None, 
                    transition = None, Av_max = None, verbose = None):

    dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr
    parmsFile     = dirname + 'parms.out'
    fileInfoFile  = dirname + 'filesInfo.out'

    #parameters dictionary used by ismcpak to produce each grid 
    parms = pickle.load(open(parmsFile))
    
    #loading the picke file holding the information of the file (a list of dicts)
    filesInfo = pickle.load(open(fileInfoFile))
    
    if verbose != None and verbose == True:
        print '-------------------------------------------'
        print '   zSec       transition    filename '
        print '-------------------------------------------'
        for fileInfo in filesInfo:
            print '%1.3e       %s         %s' % (10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename'])
    
    #looking for the reference grid (zero gmech) for this transition
    for fileInfo in filesInfo:
        
        cond1 = numpy.fabs(10.0**fileInfo['zSec'] - gmechSec) < 1e-14
        cond2 = fileInfo['transition'] == transition
        cond3 = numpy.fabs(fileInfo['Av_max'] - Av_max) < 1e-10
         
        if cond1 and cond2 and cond3:
            
            print 'found the following grid file:',
            print 10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename'], fileInfo['Av_max'] 
            fname = fileInfo['filename']
            #if os.path.isfile(fname+'-cleaned'):
            #    fname = fname+'-cleaned'
            grd = numpy.loadtxt(fname, dtype = numpy.float64)
            return grd, fname
               
    #no grid found with matching specs
    print 'could not find : '
    print 'specStr = ', specStr, 'gmechSec = ', gmechSec, 'transition = ', transition
    raise ValueError('no grid found matching specs')
        
    return None