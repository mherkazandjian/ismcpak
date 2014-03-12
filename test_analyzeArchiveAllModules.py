import numpy
import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')
#matplotlib.use('PS')

import pylab
import meshUtils
#########################################parameters##########################################################
#home = '/home/mher'
#home = '/home/mher/data2/mher'
home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

specStr_PDR   = 'HCN'
specStr_Radex = 'HCN'

parms = {
         #path to the database files
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-2.0/',
         #'dirPath'     : home + '/ism/runs/oneSided/uniformSweep2-z-2-no-mech/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweepNew-1and2/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-2/',         
         #'dirPath'      : home + '/ism/runs/oneSided/singleModels-z-2.0/',
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-1.0-high-res-no-gmech/',
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1-copy/',
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1-copy2/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-1.0/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-1.0/',
         'dirPath'      : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',
         #'dirPath'      : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',
         #'dirPath'      : home + '/ism/runs/oneSided/sph-db-z-0.2/',
         #'dirPath'      : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',
         #'dirPath'      : home + '/ism/runs/oneSided/sph-db-test/',
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMesh-z-1.0/',
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMesh-z-1.0-750-mw-CR/',
                   
         #------------------------------------------
         #'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
         'relativeGmech' : False,  # False => 3rd dim is gMech 
         #'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         #'plotRanges'    : [[-1,7],[-1,7  ],[-12, 1.2]],     # adaptive gMech 
         'plotRanges'     : [[-4,7],[-4,7],[-51, -15]],  # uniform gmech
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
                                     'quantity'       : 'fluxcgs', #,'fluxcgs', 'fluxKkms'
                                     #----------------end radex parms---------------------------------------------------
                                     'showContours'   : True,
                                     'Av_max'         : 10.0,  #the maximum Av to be used  
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
                             'Av_range'             : [0.0, 0.1],  #range which will be used in extracting data needed by radex from the PDR models
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
                             'quantity'             : 'fluxKkms',  # 'fluxKkms', 'fluxcgs' 
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

if parms['plot']:
    # plotting stuff
    pylab.ioff()
    arxv.plotGrids()
    pylab.show()

if False:
    """construct radex databases for a bunch of species for a single Av"""
    #species = ['13CO', 'HCN', 'HNC', 'HCO+', 'CS']
    #species = ['CN']  #the pop dense do not add to 1...so this is done saperatly (need to set 'checkOutputIntegrity' to False)
    #species = ['HNC', 'HCO+']
    species = ['CN']
    for specStr in species:
        parms['radex']['specStr'] = specStr
        arxv.constructRadexDatabase(writeDb = True)

if False:
    """construct radex databases for a bunch of species for a bunch of Avs"""

    parms['radex']['compute'] = True
    parms['radex']['writeDb'] = True
    
    #for Av in numpy.arange(10.0, 30.0+0.0001, 2.0):
    for Av in numpy.arange(1.0, 10.0+0.0001, 1.0):
    #for Av in numpy.array([0.01, 0.1, 1.0, 2.0]):
    #for Av in numpy.array([0.01]):
        parms['radex']['Av_range'][1] = Av
        #species = ['HCN'] #, 'HNC', 'HCO+',]
        #species = ['SiO'] #, 'HNC', 'HCO+',]
        species = ['CN'] #, 'HNC', 'HCO+',]
        #species = ['CO', '13CO']
        #species = ['CO', '13CO'] 
        for specStr in species:
            parms['radex']['specStr'] = specStr
            arxv.constructRadexDatabase(writeDb = True)
        
        '''    
        species = ['CN']  #the pop dense do not add to 1...so this is done saperatly (need to set 'checkOutputIntegrity' to False)
        parms['radex']['checkOutputIntegrity'] = False
        for specStr in species:
            parms['radex']['specStr'] = specStr
            arxv.constructRadexDatabase(writeDb = True)
        '''
            
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
    
if False:
    
    #makes a copy by selecting every other slice in the grid in each dimension
    #(i.e half the resolution)
    arxv.make_copy(dirName='/home/mher/ism/runs/oneSided/sph-db-z-1.0-low-res/', 
                   x = arxv.grid_x_unique[::2],
                   y = arxv.grid_y_unique[::2],
                   z = arxv.grid_z_unique
                  )

    #makes a copy by selecting every other slice in the grid in each dimension
    """
    arxv.make_copy(dirName='/home/mher/ism/runs/oneSided/sph-db-test/', 
                   x = [0, 2, 4], 
                   y = [0, 2, 4], 
                   z = [-50, -25, -20]
                   )
    """
    
if False:
    '''
    construct and save interpolation functions for various quantities and save these 
    interpolation function onto the disk.
    '''
    
    ## interpolation functions for emission in Kkms
    #lines  = ('13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3',)
    #lines += ('CO1-0', 'CO2-1', 'CO3-2', 'CO4-3',)
    
    for q in ['fluxKkms', 'fluxcgs']:
        
        lines = ['CN%d-%d' % (i+1, i) for i in range(8)]
    
        #lines = ['HCN1-0', 'HCN2-1', 'HCN3-2', 'HCN4-3', 'HCN5-4', 'HCN6-5', 'HCN7-6', ]
        #lines  = ('SiO1-0', 'CS1-0', 'HCO+1-0', 'HCN1-0', 'HNC1-0','13CO1-0', 'CO1-0',)
             
        for line in lines:
            print 'making interpolation function for %s' % line
            F = arxv.get_4D_interp_quantity(
                                            info={'source':'radex'},
                                            save=True,
                                            sectioned=True,
                                            line = line,
                                            Avs  = 'all',
                                            quantity = q, #'fluxcgs', # 'fluxcgs', #'fluxKkms'
                                           )
    
    ## interpolation functions for intergrated quantities (N(x), gamma, lamda..)
    pass

print 'done'
