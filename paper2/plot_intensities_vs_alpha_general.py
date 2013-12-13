# plots the relative changes in the emissions grids. Each grid represnts a certain
#relative gmech. 6 of those are plotted
#------------------------------------------------------------------------------------
import os
import matplotlib

if 'particle3' in os.uname():
    matplotlib.use('Qt4Agg')
import pylab
import numpy
    
import meshUtils
import plot_utils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
imageSaveDir = '/home/mher/tmp'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
'''
parms = {
        'line' : {
                   'code' : 'C+158',
                   'type' : 'pdr'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-1, 1], # range of the values, also that of the cbar
        'res'          : [20.0, 20.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
#        'f_interp_dim' : '2D',     
#        'interp'       : 'linear',
#        'zoom'         : 20,
        }
'''
parms = {
        'line' : {
                   'code' : 'CN1_0.5-0_0.5',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [10.0, 10.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        }

plot_utils.plot_relative_change_intensity_grid(arxvPDR, imageSaveDir=imageSaveDir, **parms)

print 'done'