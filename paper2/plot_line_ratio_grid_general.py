# plots the emissions grid
#------------------------------------------------------------------------------------
import os
import matplotlib
import numpy
import sys
from PyQt4 import QtGui
matplotlib.use('Qt4Agg')
    
import meshUtils
import plot_utils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
image_save_dir = '/home/mher/tmp'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line1' : {
                   'code' : 'O63',
                   'type' : 'pdr'
                  },
        'line2' : {
                   'code' : 'C+158',
                   'type' : 'pdr'
                  },
        
        'Av_use'    : 10.0,
        'z_sec'     : numpy.log10(1e-1),
        'res'       : [100, 100],
        'ranges'    : [[3.0, 6.0], [3.0, 6.0]],
        'cmap'      : {'obj'     : matplotlib.cm.jet,
                       'v_range' : [-2, 2],
                       'log'     : True,
                      }, 
        'cbar'      : {
                       'ticks'  : numpy.linspace(-2, 2, 5),
                       'format' : '10^x',  # None, '10^x', '%e' (or any format)
                      },
        'c_levels'  : {
                       #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                       'values' : numpy.log10([0.01,   0.1,  1,    10,    30,   100]),
                       'format' : '%.2f',  # None, '10^x', '%e' (or any format such as '%.2f'), 'strs'
                       'pow10'  : True,
                       #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                      },
        'clip'      : None,
        #clip       : [-15.0, 0.0],
        #clip_max_n': 4.0,
        'removeNans': False,
        'f_interp_dim' : '3D',
        }

#--------------------------------------------------------------------------------------------------------------------------

if False:
    line_ratio = plot_utils.line_ratio(arxvPDR, image_save_dir=image_save_dir, **parms)

if True:
    def main():
        app = QtGui.QApplication(sys.argv)
        gui = plot_utils.line_ratio_app(arxvPDR=arxvPDR, initial_parms=parms)
        #gui = plot_utils.line_ratio_app(arxvPDR=arxvPDR)
        sys.exit(app.exec_())    
    
    if __name__ == '__main__':
        main()    

print 'done'
