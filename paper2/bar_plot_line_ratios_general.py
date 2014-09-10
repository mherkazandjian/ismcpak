# purpose : - produce bar plots for the line emissions for different lines, column
#             column densities...etc..
#           - here focusing on CO lines
# keywords: plot, line ratio, intensity, molecular, CO, fine structure
#--------------------------------------------------------------------------------

import numpy
import sys, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import meshUtils
import plot_utils


#########################################parameters##########################################################
home           = '/home/mher'
dirPath        = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
image_save_dir = '/home/mher/tmp'


Av_max      = 10.0
#line_ratios = ['O63/C369', 'O63/C609', 'O63/C+158', 'C+158/C369', 'C+158/C609', 'C369/C609'] 
line_ratios = ['CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO16-15/CO1-0', 'CO7-6/CO3-2', 'CO10-9/CO7-6', 'CO16-15/CO10-9'] 
ranges      = [ [-3.0, 2.0], [-1.0, 3.0], [-1.0, 3.0], [-1.0, 2.0], [-1.0, 4.0], [-1.0, 4.0] ]
names       = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4']
grid_coords = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ]


#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges=ranges, 
                                names=names, 
                                grid_coords=grid_coords,
                                line_ratios = line_ratios,
                                Av_max = Av_max,
                                image_save_dir = '/home/mher/tmp',
                                plot_title = 'CO-ratios'
                                )
