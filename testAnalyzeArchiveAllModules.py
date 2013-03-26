import numpy
import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')

import pylab
import meshUtils
#########################################parameters##########################################################
home = '/home/mher'

specStr_PDR   = 'O'
specStr_Radex = 'CO'

parms = {
         #path to the database files
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-2.0/',
         #'dirPath'     : home + '/ism/runs/oneSided/uniformSweep2-z-2-no-mech/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweepNew-1and2/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-2/',         
         #'dirPath'      : home + '/ism/runs/oneSided/singleModels-z-2.0/',
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-1.0-high-res-no-gmech/',
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-1.0/',
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plotRanges'    : [[-1,7],[-1,7  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[-3,7],[-3,7],[-51, -15]],  # uniform gmech
         
         'plot'          : True, 
         'showGrids'     : True,
         'gridsInfo'     : { '00' : {#some quantity
                                    'show'     : True,
                                    'quantity' : ['state', 'gasT'],
                                    #'quantity' : ['therm', 'heating'],
                                    'slabIdx'  : 0,
                                    },
                             '01' : {# abundance 
                                    'show'     : True, 
                                    'quantity' : ['state', 'abun'],
                                    'slabIdx'  : 0,
                                    'specStr'  : specStr_PDR,
                                    },
                             '10' : { # column density
                                    'show'     : True,
                                    'maxAv'    : 10,
                                    'specStr'  : specStr_PDR,
                                    },
                             '11' : { # line intensitities
                                     'show'           : True,
                                     #------------------comment those if radex parms is 'pdr' is selected below this------------                                    
                                     'type'           : 'pdr', #if type = pdr, quantity should point to a valid destination in the dtype in arxv.meshes[i]
                                     'quantity'      : ['fineStructureCoolingComponents','C','rate','1-0'], # for use with 'pdr'
                                     'specStr'        : 'C',     # database to be restored/computed
                                     #-----------comment those radex parms if 'pdr' is selected above this--------------
                                     #'type'           : 'radex',
                                     #'specStr'        : specStr_Radex,     # database to be restored/computed
                                     #'transitionIndx' : 0,
                                     #'quantity'       : 'fluxcgs',
                                     #----------------end radex parms---------------------------------------------------
                                     'showContours'   : True,
                                     'Av_max'         : 10.0,  #the maximum Av to be used  
                                    },
                           },
         'gridsRes'      : 100,
         
         'meshPltAvRng'  : [0, 30.0], #plotting range as a function of Av
          
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : True,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  #range which will be used in extracting data needed by radex from the PDR models
                                                                    #(only relevent to constructing databases)
                             'path'                 : home + '/ism/code/radex/Radex/bin/radex',  
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
                             'checkOutputIntegrity' : False,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 1e-2,
                             #'popDensSumTol'        : 10,
                             'changeFracTrial'      : 0.001,
                             'nMaxTrial'            : 10,
                            },
        }
#############################################################################################################

# reading the archive
print 'setting up the archive'
t0 = time.time()
arxv = meshUtils.meshArxv(readDb = True, **parms)
print 'time reading data %f' % (time.time() - t0)

if parms['plot']:
    # plotting stuff
    pylab.ioff()
    arxv.plotGrids()
    pylab.show()

if False:
    """construct radex databases for a bunch of species"""
    #species = ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS', 'CN']
    species = ['CN']  #the pop dense do not add to 1...so this is done saperatly (need to set 'checkOutputIntegrity' to False)
    for specStr in species:
        parms['radex']['specStr'] = specStr
        arxv.constructRadexDatabase(writeDb = parms['radex']['writeDb'])
    
if False:
    arxv.save_radex_grids(
                          relativeDirPath = 'analysis/%s/' % parms['radex']['specStr'],
                          #relativeDirPath = 'analysis/%s/colFmt/' % parms['radex']['specStr'],
                          basename        = 'radexGrid',
                          transitionInds  = [0,1,2,3,4,5], #CO transitions
                          #transitionInds  = range(40), 
                          #transitionInds  = [0,1],
                          quantity        = 'fluxcgs',
                          fileFormat      = 'numpytxt')  # 'numpytxt' or '3columns'

if False:
    arxv.save_PDR_emission_grids(relativeDirPath = 'analysis/%s/' % parms['gridsInfo']['11']['quantity'][1],
                                 basename        = 'pdrGrid',
                                 transitions     = ['1-0'],
                                 quantity        = 'rate',)

if False:
    arxv.save_PDR_quantity_grids(relativeDirPath = 'analysis/surfaceHeating/',
                                 basename        = 'grid',
                                 quantity        = ['therm', 'heating'],
                                 slabIdx         = 0,
                                )
if False:
    arxv.do_something_for_all_meshes()
    
    
print 'done'