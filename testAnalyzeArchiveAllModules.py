from numpy import *
from time import *
import sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
from meshUtils import *

#########################################parameters##########################################################
home = '/home/mher'

parms = {
         #path to the database files
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-2.0/',
         #'dirPath'     : home + '/ism/runs/oneSided/uniformSweep2-z-2-no-mech/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweepNew-1and2/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-2/',         
         'dirPath'      : home + '/ism/runs/oneSided/singleModels-z-0.1/',
         #'dirPath'      : home + '/ism/runs/oneSided/surfaceGrid-z-0.1/',
         
         # reference database
         'runDirPath2'   : home + '/ism/runs/oneSided/surfaceGrid-z-0.1-high-res-no-gmech/',
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'plotRanges'    : [[0,6],[0,6  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[0,6],[0,6],[-51, -15]],  # uniform gmech
         'metallicity'    : 0.1,

         'plotGrids'     : False,
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
                                    'specStr'  : 'H+',
                                    },
                             '10' : { # column density
                                    'show'     : True,
                                    'maxAv'    : 30,
                                    'specStr'  : 'CO',
                                    },
                             '11' : { # line intensitities
                                     'show'           : True,
                                     #'type'           : 'pdr', #if type = pdr, quantity should point to a valid destination in the dtype in arxv.meshes[i]
                                     #'quantity'      : ['fineStructureCoolingComponents','O','rate','1-0'], # for use with 'pdr'
                                     'type'           : 'radex',
                                     'specStr'        : 'CO',     # database to be restored/computed
                                     'transitionIndx' : 0,
                                     'quantity'       : 'fluxcgs',
                                     'showContours'   : True,
                                    },
                           },
         'gridsRes'      : 100,
         
         'radex'         : { 'use'                  : False,
                             ###-----------radex database parms-----------------
                             'compute'              : True, #if true, runns radex on all meshes
                             'writeDb'              : True, #if true, writes the computed stuff to a db
                             'path'                 : home + '/ism/code/radex/Radex/bin/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : 'CO',
                             'freqRange'            : [0, 50000],
                             #'xH2_Min'              : 2*0.0000000001
                             'xH2_Min'              : -1.0,
                             'collisionPartners'    : ['H2','H+','H','e-','He'],
                             #'collisionPartners'    : ['H2','H','H+','e-'],
                             #'collisionPartners'    : ['H2'],
                             'tBack'                : 2.73,
                             'lineWidth'            : 1.0,
                             'verbose'              : True, 
                             'maxDisplayTranistion' : 20,
                             ###----------extra convergence params-----------------------
                             'checkOutputIntegrity' : False,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 1e-2,
                             #'popDensSumTol'        : 10,
                             'changeFracTrial'      : 0.001,
                             'nMaxTrial'            : 100,
                            },
         #-----------------chemical network parameters------------------------
         'chemistry'     : {
                            'rxnFile'       : home + '/ism/code/ismcpak/data/rate99Fixed.inp',
                            'specNumFile'   : home + '/ism/code/ismcpak/data/species.inp',
                            'underAbunFile' : home + '/ism/code/ismcpak/data/underabundant.inp',
                            'removeManual'  : ['13CH3'],
                            'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
                            'umistVer'      : 'umist99',
                           }
        }
#############################################################################################################

# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv( **parms )
arxv.readDb( check = True)
print 'time reading data %f' % (time() - t0)

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech         = parms['relativeGmech'], 
                                   referenceDatabasePath = parms['runDirPath2'] )

if parms['radex']['use'] and parms['gridsInfo']['11']['show']:
    if parms['radex']['compute']:
        arxv.constructRadexDatabase(writeDb = parms['radex']['writeDb'])
    else:
        arxv.readDbRadex(parms['radex']['specStr'], check = True)


if parms['plotGrids']:
    # plotting stuff
    pyl.ioff()
    arxv.plotGrids()
    pyl.show()

if False:
    arxv.save_radex_grids(
                          relativeDirPath = 'analysis/%s/' % parms['radex']['specStr'],
                          #relativeDirPath = 'analysis/%s/colFmt/' % parms['radex']['specStr'],
                          basename        = 'radexGrid',
                          transitionInds  = [0,1,2,3,4,5,6],
                          #transitionInds  = range(40), 
                          #transitionInds  = [0,1],
                          quantity        = 'fluxcgs',
                          fileFormat      = 'numpytxt')  # 'numpytxt' or '3columns'

if False:
    arxv.save_PDR_emission_grids(relativeDirPath = 'analysis/%s/' % parms['gridsInfo']['11']['quantity'][1],
                                 basename        = 'pdrGrid',
                                 transitions     = ['1-0'],
                                 quantity        = 'rate')

if False:
    arxv.save_PDR_quantity_grids(relativeDirPath = 'analysis/surfaceHeating/',
                                 basename        = 'grid',
                                 quantity        = ['therm', 'heating'],
                                 slabIdx         = 0,
                                )
if False:
    arxv.do_something_for_all_meshes()
    
    
print 'done'