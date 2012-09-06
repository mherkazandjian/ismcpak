from numpy import *
from time import *
import sys
import pylab as pyl

from meshUtils import *

#########################################parameters##########################################################
home = '/home/mher'

parms = {
         #path to the database files
         #'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         #'dirPath'     : home + '/ism/runs/oneSided/uniformSweep2-z-2-no-mech/',
         'dirPath'      : home + '/ism/runs/oneSided/uniformSweepNew-1and2/',
         #'dirPath'      : home + '/ism/runs/oneSided/uniformSweep2-z-2/',

         # reference database
         'runDirPath2'   : home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/',

         'relativeGmech' : False,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'zSec'          : -30,  # section in the 3D dimension to be used for generating 
                                  # grids. usuall this is the log10 of the mechanical heating
                                  # it can be used as the ratio of mechanical heating to the
                                  # surface heating(gMech = 0)
         'plotLog10zAxs' : False, # set this to true to take the log of arxv.grid_z when plotting
         
         #plotRanges    : [[0,6],[0,6],[-30, -15]],
         'plotRanges'     : [[0,6],[0,6],[-30, -15]],
         'metallicity'    : 1.0,

         'gridsInfo'     : { '00' : {#some quantity
                                    'quantity' : ['state', 'gasT'],
                                    'slabIdx'  : 0,
                                    },
                             '01' : {# abundance 
                                    'quantity' : ['state', 'abun'],
                                    'slabIdx'  : 0,
                                    'specStr'  : 'H+',
                                    },
                             '10' : { # column density
                                    'maxAv'    : 10,
                                    'specStr'  : 'H2O',
                                    },
                             '11' : { # line intensitities
                                     'quantity' : ['state', 'abun'],
                                     'slabIdx'  : 0,
                                     'specStr'  : 'H+',
                                    },
                           },
         'gridsRes'       : 100,   
         'radexParms'     : { 'use'                  : False,
                              'radexPath'            : home + '/ism/code/radex/Radex/bin/radex',  
                              'molDataDirPath'       :  home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                              'specStr'              : 'CO',
                              #'xH2_Min'              : 2*0.0000000001
                              'xH2_Min'              : -1.0,
                              'collisionPartners'    : ['H2','H+','H'],
                              #'collisionPartners'    : ['H2','H','H+','e-']
                              #'collisionPartners'    : ['H2']
                              #'collisionPartners'    : ['H2','H+','e-','H']
                              'plotTransitionInGrid' : 0,
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

# read and setting up the chemical network used in the 
t0 = time()
arxv.setupChemistry()
print 'time setting up the chemistry %f' % (time() - t0)

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech         = parms['relativeGmech'], 
                                   referenceDatabasePath = parms['runDirPath2'] )

# plotting stuff
pyl.ioff()
arxv.plotGrid(parms['gridsRes'], 
              parms['zSec'], 
              radex = parms['radexParms'], 
              ranges = parms['plotRanges'], 
              gridsInfo = parms['gridsInfo'], 
              plotLog10zAxs = parms['plotLog10zAxs'],
              parms = parms)

pyl.show()
"""
#arxv.saveGridsToFiles(gridsRes, lgammaMechSec, radexParms)
pyl.show()
"""
print 'done'
