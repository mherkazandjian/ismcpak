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

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plotRanges'    : [[0,6],[0,6],[-12, 6]],     # adaptive gMech 
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
                                     'Av_max'         : 10.0,  #the maximum Av to be used  
                                    },
                           },
         'gridsRes'      : 100,
         
         'meshPltAvRng'  : [0, 30.0], #plotting range as a function of Av
          
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : False,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  #range which will be used in extracting data needed by radex from the PDR models
                                                                    #(only relevent to constructing databases)
                             'path'                 : home + '/ism/code/radex/Radex/bin/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : 'CO',
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
    
arxv.save_PDR_quantity_grids(relativeDirPath = 'analysis/surfaceHeating/',
                             basename        = 'grid',
                             quantity        = ['therm', 'heating'],
                             slabIdx         = 0,
                            )

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# plot a grid from a pdr quantity (see mesh data format) into one panel with
# colorbar, contour lines (values are user specified)

import numpy
import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')

import pylab
import meshUtils

import pickle
import matplotlib.cm as cm

subDirName    = 'surfaceHeating'
plotTitle     = r'$\Gamma_{surf}$ [erg cm$^{-3}$ s$^{-1}]$' 
quantity      = ['therm','heating']
slabIdx       = 0
cLevels       = numpy.arange(20) - 30
cbarTicks     = numpy.arange(-30.0, -15.0, 2.0)
dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % subDirName
parmsFile     = dirname + 'parms.out'
fileInfoFile  = dirname + 'filesInfo.out'
colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/%s-%d-base.eps' % (subDirName, slabIdx)
#======================================================================================

width  = 3    #figure width (non normalized) 
height = 3.4  #figure heigh (non normalized)
as_rat = width/height #aspect ratio of the figure

ax_xs  = 0.17 #axses x start (normalized)
ax_ys  = 0.17 #axses y start (normalized)
ax_sz  = 0.75 #axses size (normalized)

cbar_xs = ax_xs        #colorbar x start
cbar_ys = ax_sz + 0.13 #colorbar x start
cbar_sc = 0.99         #scale of the width of the cbar (relative to the width of ax)
cbar_w  = 0.02         #width of the cbar (normalized)

fig    = pylab.figure(0, figsize = (width, height) )
fig.set_facecolor('white')
ax1    = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])
axCbar = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_sz, cbar_ys, cbar_sc*ax_sz - (0.5*(1.0-cbar_sc))*ax_sz, cbar_w])
 
#loading the pickle file holding the parameters used to generate the data (a dict)
parms = pickle.load(open(parmsFile))
#loading the picke file holding the information of the file (a list of dicts)
filesInfo = pickle.load(open(fileInfoFile))

print '-------------------------------------------'
print '   zSec       transition    filename '
print '-------------------------------------------'
for fileInfo in filesInfo:
    print '%1.3e       %s         %s' % (10.0**fileInfo['zSec'], fileInfo['slabIdx'], fileInfo['filename'])

ranges = parms['plotRanges']

#looking for the reference grid (zero gmech) for this transition
for fileInfo in filesInfo:
    if (numpy.fabs(10.0**fileInfo['zSec'] - 1e-10) < 1e-14):
        print 10.0**fileInfo['zSec'], fileInfo['slabIdx'], fileInfo['filename']
        fname = fileInfo['filename']        
grd = numpy.loadtxt(fname, dtype = numpy.float64)

rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

im = ax1.imshow(grd, extent = rangesLst, origin='lower', cmap = colormap)

CS = ax1.contour(grd, cLevels, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')

pylab.clabel(CS, fmt = '%.1f' )
cbar = pylab.colorbar(im, cax = axCbar, orientation = 'horizontal')

#-----------------
#---i used xslect() and picked the points manually from the plot...i had to
#   do it just once, so it is ok to use primitive techniques :D
#
g0Line = numpy.array([ 1.72071948,  1.40111874,  1.11400923,  2.09678229,  1.26451065,
        1.70563702,  1.50500901,  1.54691774,  1.79392117,  1.51833355,
        1.95240084,  1.88088171,  1.61911698,  1.70208335,  2.42088828,
        2.13019158, -1.54607488,  2.74487783,  2.6605203 ,  2.1827508 ,
        2.53768582,  1.87725719,  1.9073518 ,  2.63187265,  2.97012525,
        1.70722868,  1.93239909,  2.89917134,  1.93696925,  0.89165463,
        2.03208654,  1.59173189,  1.26532584,  1.22406598])
nLine = numpy.array([-0.95343206, -0.74296246, -0.65047267, -0.48159976, -0.4178872 ,
       -0.38536022, -0.05497469, -0.10211952,  0.75582297,  0.61765387,
        0.79158139,  1.83037719,  1.79642624,  2.1968747 ,  2.48060574,
        2.93450754, -5.20260207,  3.48034674,  3.30123644,  2.90356134,
        2.87118067,  2.81722183,  2.81334312,  2.82651464,  2.75454441,
        2.63855186,  2.59596492,  2.59132914,  2.46318178,  2.39410328,
        2.2353289 ,  2.15386878,  2.12006321,  2.07260759])
ax1.plot(nLine, g0Line, 'ko')
ax1.set_xlim([0,6])
ax1.set_ylim([0,6])
#-----------------
ax1.set_xlabel('$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='large')
ax1.set_ylabel('$\log_{10}$ [G$_0$]', size='large')

axCbar.set_title(plotTitle)
cbar.set_ticks(cbarTicks)

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pylab.show()

print 'min-max log10 intensities = ', numpy.nanmin(grd), numpy.nanmax(grd) 
print 'done'
