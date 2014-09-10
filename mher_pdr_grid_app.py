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

'''
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
#        'line' : {
#                   'code' : 'C+158',
#                   'type' : 'pdr'
#                 },

        'line' : {
                   'code' : 'CO1-0',
                   'type' : 'radex-lvg'
                 },

        'Av_use'    : 10.0,
        'z_sec'     : numpy.log10(1e-10),
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : numpy.log10([1e-7, 2e-7, 2.5e-7]),
        #'cLevels'   : [-8, -7, -6],
        'cbarTicks' : [-6.0, -5.0, -4.0, -3.0],
        'clip'      : None,
        #clip       : [-15.0, 0.0],
        #clip_max_n': 4.0,
        'removeNans': False,
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
#grd = numpy.log10(10.0**grd / (2.0*numpy.pi))
plot_utils.plot_intensity_grid(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
'''
print 'done'
