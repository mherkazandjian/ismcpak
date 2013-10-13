'''
test making a copy of the meshArxv object using the methods .as_dict and .set_data_from_dict.
'''
import numpy
import time
import sys
import os

import meshUtils
#########################################parameters##########################################################
home = '/home/mher'

specStr_PDR   = 'HCN'
specStr_Radex = 'HCN'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMesh-z-1.0/',
                   
         #------------------------------------------
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
         #'relativeGmech' : False,  # False => 3rd dim is gMech 
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plotRanges'    : [[-1,7],[-1,7  ],[-12, 1.2]],     # adaptive gMech 
         #'plotRanges'     : [[-4,6],[-4,6],[-51, -15]],  # uniform gmech
         #'plotRanges'    : [[-1,7],[-1,7  ],[-18, -12]],     # variable CR rate 

         #useful for variable CR runs
         #'grid_qx'       : ['hdr','nGas'],
         #'grid_qy'       : ['hdr','G0'],
         #'grid_qz'       : ['from_meshes_info', 'parms', 3, 'CR_rate'],  # 3 indicates the 4th column in self.infoAll['parms']
         
         #-----------------------------------         
         'plot'          : True, 
         'showGrids'     : True,
         'gridsInfo'     : { '00' : {#some quantity
                                    'show'     : True,
                                    #'quantity' : ['state', 'gasT'],
                                    'quantity' : ['therm', 'heating'],
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
                                    'maxAv'    : 30,
                                    'specStr'  : specStr_PDR,
                                    },
                             '11' : { # line intensitities
                                     'show'           : True,
                                     ##------------------comment those if radex parms is 'pdr' is selected below this------------                                    
                                     #'type'           : 'pdr', #if type = pdr, quantity should point to a valid destination in the dtype in arxv.meshes[i]
                                     #'quantity'      : ['fineStructureCoolingComponents','C+','rate','1-0'], # for use with 'pdr'
                                     #'specStr'        : 'C',     # database to be restored/computed
                                     ##-----------comment those radex parms if 'pdr' is selected above this--------------
                                     'type'           : 'radex',
                                     'specStr'        : specStr_Radex,     # database to be restored/computed
                                     'transitionIndx' : 0,
                                     'quantity'       : 'fluxcgs',
                                     #----------------end radex parms---------------------------------------------------
                                     'showContours'   : True,
                                     'Av_max'         : 5.0,  #the maximum Av to be used  
                                    },
                           },
         'gridsRes'      : 100,
         'nThreads'      : 1,
         'meshPltAvRng'  : [0, 30.0], #plotting range as a function of Av
          
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : False,
                             'plot_model_from_Db'   : False,   #plot the emission ladder from the database [do not run radex]
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 30.0],  #range which will be used in extracting data needed by radex from the PDR models
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
                             'checkOutputIntegrity'       : True,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'         : 1.0, 
                             'popDensSumTol'              : 1e-2,
                             'changeFracTrial'            : 0.05,
                             'strict'                     : True,
                             'nMaxTrial'                  : 100,
                            },
        }
#############################################################################################################

# reading the archive
print 'setting up the archive'
t0 = time.time()
arxv = meshUtils.meshArxv(readDb = True, **parms)
print 'time reading data %f' % (time.time() - t0)


#getting the data which are 'pickalable' as a dict 
data_as_dict = arxv.as_dict()

#making a copy of arxv from the data
arxv_new = arxv.make_copy_from_data(data_as_dict)
print 'made a copy of the data'
