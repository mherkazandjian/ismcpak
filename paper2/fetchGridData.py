# looks up and restors a grid from a file name
#------------------------------------------------------------------------------------
import numpy as np
import pickle, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
import matplotlib.cm as cm


def fetchRadexGrid( dirname = None, specStr = None, gmechSec = None, transition = None, verbose = None):

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
        if (np.fabs(10.0**fileInfo['zSec'] - gmechSec) < 1e-14) and fileInfo['transition'] == transition:
            print 'found the following grid file:',
            print 10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename']
            fname = fileInfo['filename']
            grd = np.loadtxt(fname, dtype = np.float64)
            return grd
               
    #no grid found with matching specs
    print 'could not find : '
    print 'specStr = ', specStr, 'gmechSec = ', gmechSec, 'transition = ', transition

    asdasd
    return None