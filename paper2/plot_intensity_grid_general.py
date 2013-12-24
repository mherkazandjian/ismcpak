# plots the emissions grid
#------------------------------------------------------------------------------------
import os
import matplotlib
import numpy

if 'particle3' in os.uname():
    matplotlib.use('Qt4Agg')
    
import meshUtils
import plot_utils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
imageSaveDir = '/home/mher/tmp'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : numpy.log10(1e-10), #-10.0,
        'res'       : [200.0, 200.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cbar_range': [-12, -6], 
        'cLevels'   : [-10.0, -8.0 , -7.0, -6.5, -6.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],
        'clip'      : None,
        'removeNans': False,
        'get_section_data_using' : '3d-interp',
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
#        'zoom'      : 10,
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
print 'done'

