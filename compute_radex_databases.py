import numpy
import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')
#matplotlib.use('PS')
import multiprocessing

import pylab
import meshUtils
#########################################parameters##########################################################
home = '/home/mher'

specStr_PDR   = 'CO'
specStr_Radex = 'CO'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMesh-z-0.5/',
         #'dirPath'      : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',
         
         'relativeGmech' : False,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         #'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plot'          : False, 
         'showGrids'     : False,
         'nThreads'      : 1,
          
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : False,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  #range which will be used in extracting data needed by radex from the PDR models
                                                                    #(only relevent to constructing databases)
                             'path'                 : home + '/ism/code/radex/Radex/bin-gcc/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : specStr_Radex,
                             'freqRange'            : [0, 50000],
                             #'xH2_Min'              : 2*0.0000000001
                             'xH2_Min'              : -1.0,
                             #'collisionPartners'    : ['H2','H+','H','e-','He'],
                             #'collisionPartners'    : ['H2','H','H+','e-'],
                             'collisionPartners'    : ['H2'],
                             'use_pdr_gas_den_H2'   : True,   #<----------
                             'tBack'                : 2.73,
                             'lineWidth'            : 1.0,
                             'verbose'              : False,
                             'maxDisplayTranistion' : 20,
                             ###----------extra convergence params-----------------------
                             'checkOutputIntegrity' : True,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 3e-2,
                             'changeFracTrial'      : 0.01,
                             'strict'               : True,
                             'nMaxTrial'            : 100,
                            },
        }
#############################################################################################################

#
def worker(Av):

    # reading the archive
    print 'setting up the archive'
    t0 = time.time()
    arxv = meshUtils.meshArxv(readDb = True, **parms)
    print 'time reading data %f' % (time.time() - t0)

    parms['radex']['Av_range'][1] = Av

    """construct radex databases for a bunch of species for a single Av"""
    #species = ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS']
    #species = ['CN']  #the pop dense do not add to 1...so this is done saperatly (need to set 'checkOutputIntegrity' to False)
    #species = ['HNC', 'HCO+']
    species = ['CN']#, 'HCN']
    for specStr in species:
        parms['radex']['specStr'] = specStr
        arxv.constructRadexDatabase(writeDb = True)
    
#


pool = multiprocessing.Pool(2)
pool.map(worker, [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0])
#pool.map(worker, [3.0, 4.0, 5.0, 6.0])

        
#worker([3.0])

'''    
if False:
    """construct radex databases for a bunch of species for a bunch of Avs"""
    for Av in numpy.arange(3.0, 30.0, 1.0):
        parms['radex']['Av_range'][1] = Av
        species = ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS', 'SiO']
        for specStr in species:
            parms['radex']['specStr'] = specStr
            arxv.constructRadexDatabase(writeDb = True)
            
        species = ['CN']  #the pop dense do not add to 1...so this is done saperatly (need to set 'checkOutputIntegrity' to False)
        parms['radex']['checkOutputIntegrity'] = False
        for specStr in species:
            parms['radex']['specStr'] = specStr
            arxv.constructRadexDatabase(writeDb = True)
'''
    
    
print 'done'