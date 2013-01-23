# purpose : produce bar plots for the line emissions for different lines, column
#           column densities...etc..
# keywords: plot, line ratio, intensity 
#--------------------------------------------------------------------------------

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
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         
         # reference database
         'runDirPath2'   : home + '/ism/runs/oneSided/surfaceGridHighRes-z-1.0/',
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'plotRanges'    : [[0,6],[0,6  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[0,6],[0,6],[-51, -15]],  # uniform gmech
         'metallicity'    : 1.0,

         'plotGrids'     : False,
         
         'radex'         : { 'use'                  : True,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
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
##############################################setting up and reading the data###############################################################

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

#reading all the available precomputed radex databases
arxv.readDbsRadex(species=['13CO', 'CN', 'CO', 'CS',
                           'HCN', 'HCO+', 'HNC', 'SiO'])


############################extracting what we want and plotting####################################
arxv.use_radexDb('CO')
fig = pyl.figure()
ax = fig.add_subplot(111)

log10n, log10G0 = (3.0, 30,)
gm_v = np.log10(np.array([0.1, 1.0, 5.0, 10.0, 50.0, 100.0])/100.0)
#-------------------------------------------------------------------------
#get the fluxes for this transition for different mechanical heating rates 
transitionIdx = 0
values1 = []
inds1 = np.arange(gm_v.size)

#find the index of the particular model of interest
gm = -10 
ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm) 
v_zero = arxv.meshesRadex[ind][transitionIdx]['fluxcgs']

for gm in gm_v:
    ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm)
    values1.append(arxv.meshesRadex[ind][transitionIdx]['fluxcgs'])

values1 = np.log10(values1)
rects1 = ax.bar(inds1, values1, width=0.2, bottom=0, color = 'g')
#-------------------------------------------------------------------------
#get the fluxes for this transition for different mechanical heating rates 
transitionIdx = 1
values2 = []
inds2 = np.arange(gm_v.size)

#find the index of the particular model of interest
gm = -10 
ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm) 
v_zero = arxv.meshesRadex[ind][transitionIdx]['fluxcgs']

for gm in gm_v:
    ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm)
    values2.append(arxv.meshesRadex[ind][transitionIdx]['fluxcgs'])

values2 = np.log10(values2)
rects2 = ax.bar(inds2+0.2, values2, width=0.2, bottom=0, color = 'r')
#-------------------------------------------------------------------------
#get the fluxes for this transition for different mechanical heating rates 
transitionIdx = 2
values3 = []
inds2 = np.arange(gm_v.size)

#find the index of the particular model of interest
gm = -10 
ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm) 
v_zero = arxv.meshesRadex[ind][transitionIdx]['fluxcgs']

for gm in gm_v:
    ind = arxv.get_mesh_index(x = log10n, y = log10G0, z = gm)
    values3.append(arxv.meshesRadex[ind][transitionIdx]['fluxcgs'])

values3 = np.log10(values3)
rects2 = ax.bar(inds2+0.2+0.2, values3, width=0.2, bottom=0, color = 'b')

pyl.show()


print 'done'
