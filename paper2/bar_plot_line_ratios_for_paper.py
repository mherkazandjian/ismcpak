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

home           = '/home/mher'
image_save_dir = os.path.join(home, 'ism/docs/paper02/src/figs')
#image_save_dir = '/home/mher/tmp/figs'

#########################################parameters##########################################################
dirPath        = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)
#############################################################################################################

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-0.5, 2.0], [-0.5, 2.0], [-0.5, 2.5], [0.0, 3.0], [0.0, 5.0], [0.0, 5.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['O63/C369', 'O63/C609', 'O63/C+158', 'C+158/C369', 'C+158/C609', 'C369/C609'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-atomic-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 2.0], [-1.0, 3.0], [-1.0, 3.0], [-1.0, 2.0], [-1.0, 4.0], [-1.0, 4.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO16-15/CO1-0', 'CO7-6/CO3-2', 'CO10-9/CO7-6', 'CO16-15/CO10-9'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-CO-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 2.0], [-1.0, 3.0], [-1.0, 3.0], [-1.0, 2.0], [-1.0, 4.0], [-1.0, 4.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO16-15/13CO1-0', '13CO7-6/13CO3-2', '13CO10-9/13CO7-6', '13CO16-15/13CO10-9'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-13CO-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 0.0], [-2.0, 0.0], [-2.0, 0.0], [-2.0, 0.0], [-2.0, 0.0], [-2.0, 0.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['13CO1-0/CO1-0', '13CO2-1/CO2-1', '13CO3-2/CO3-2', '13CO4-3/CO4-3', '13CO7-6/CO7-6', '13CO16-15/CO16-15'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-13CO-CO-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-1, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-2.0, 2.0], [-2.0, 2.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HNC4-3/HNC1-0', 'HCN4-3/HCN1-0', 'HNC1-0/HCN1-0', 'HNC4-3/HCN4-3', 'HNC4-3/HCN1-0'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HNC-HCN-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 1.0], [-3.0, 1.0], [-3.0, 2.0], [-3.0, 2.0], [-3.0, 3.0], [-3.0, 3.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+4-3/HCO+1-0', 'CN2_1.5-1_1.5/CN1_0.5-0_0.5', 'CS4-3/CS1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HCO+-CN-CS-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2, 1.0], [-3.0, 2.0], [-2.0, 2.0], [-2.0, 2.0], ], 
                                names          = [ 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+1-0/CO1-0', 'HCO+1-0/13CO1-0', 'HCN1-0/CO1-0', 'HNC1-0/HCO+1-0', 'HCN4-3/HCO+4-3', 'CN1_0.5-0_0.5/HCN1-0', 'CN2_1.5-1_1.5/HCN1-0', 'CS4-3/13CO1-0'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-misc1-z-%.1f' % arxvPDR.metallicity,
                                gm_v           = numpy.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])/100.0,
                                colors         = [            'k', 'g', 'b', 'c' , 'y',   'm',  'r']
                                )


################################################################################################
# Z = 0.1 Z_sun bar plots
################################################################################################

dirPath = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-0.1/')
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2.0, 2.0], [-2.0, 2.0], [-0.5, 3.0], [-0.5, 3.0], [0.0, 4.0], [0.0, 5.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['O63/C369', 'O63/C609', 'O63/C+158', 'C+158/C369', 'C+158/C609', 'C369/C609'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-atomic-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-1, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 2.0], [-1.0, 2.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HNC4-3/HNC1-0', 'HCN4-3/HCN1-0', 'HNC1-0/HCN1-0', 'HNC4-3/HCN4-3', 'HNC4-3/HCN1-0'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HNC-HCN-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 1.0], [-3.0, 1.0], [-3.0, 2.0], [-3.0, 2.0], [-3.0, 3.0], [-3.0, 3.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+4-3/HCO+1-0', 'CN2_1.5-1_1.5/CN1_0.5-0_0.5', 'CS4-3/CS1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HCO+-CN-CS-z-%.1f' % arxvPDR.metallicity,
                                )


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2, 1.0], [-3.0, 2.0], [-2.0, 2.0], [-2.0, 2.0], ], 
                                names          = [ 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+1-0/CO1-0', 'HCO+1-0/13CO1-0', 'HCN1-0/CO1-0', 'HNC1-0/HCO+1-0', 'HCN4-3/HCO+4-3', 'CN1_0.5-0_0.5/HCN1-0', 'CN2_1.5-1_1.5/HCN1-0', 'CS4-3/13CO1-0'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-misc1-z-%.1f' % arxvPDR.metallicity,
                                gm_v           = numpy.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])/100.0,
                                colors         = [            'k', 'g', 'b', 'c' , 'y',   'm',  'r']
                                )

################################################################################################
# Z = 0.5 Z_sun bar plots
################################################################################################

dirPath = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-0.5/')
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2.0, 2.0], [-2.0, 2.0], [-0.5, 3.0], [-0.5, 3.0], [0.0, 4.0], [0.0, 5.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['O63/C369', 'O63/C609', 'O63/C+158', 'C+158/C369', 'C+158/C609', 'C369/C609'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-atomic-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-1, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 2.0], [-1.0, 2.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HNC4-3/HNC1-0', 'HCN4-3/HCN1-0', 'HNC1-0/HCN1-0', 'HNC4-3/HCN4-3', 'HNC4-3/HCN1-0'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HNC-HCN-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 1.0], [-3.0, 1.0], [-3.0, 2.0], [-3.0, 2.0], [-3.0, 3.0], [-3.0, 3.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+4-3/HCO+1-0', 'CN2_1.5-1_1.5/CN1_0.5-0_0.5', 'CS4-3/CS1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HCO+-CN-CS-z-%.1f' % arxvPDR.metallicity,
                                )


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2, 1.0], [-3.0, 2.0], [-2.0, 2.0], [-2.0, 2.0], ], 
                                names          = [ 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+1-0/CO1-0', 'HCO+1-0/13CO1-0', 'HCN1-0/CO1-0', 'HNC1-0/HCO+1-0', 'HCN4-3/HCO+4-3', 'CN1_0.5-0_0.5/HCN1-0', 'CN2_1.5-1_1.5/HCN1-0', 'CS4-3/13CO1-0'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-misc1-z-%.1f' % arxvPDR.metallicity,
                                gm_v           = numpy.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])/100.0,
                                colors         = [            'k', 'g', 'b', 'c' , 'y',   'm',  'r']
                                )

################################################################################################
# Z = 2.0 Z_sun bar plots
################################################################################################

dirPath = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-2.0/')
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2.0, 2.0], [-2.0, 2.0], [-0.5, 3.0], [-0.5, 3.0], [0.0, 4.0], [0.0, 5.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.8],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['O63/C369', 'O63/C609', 'O63/C+158', 'C+158/C369', 'C+158/C609', 'C369/C609'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-atomic-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-1, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0], [-1.0, 2.0], [-1.0, 2.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HNC4-3/HNC1-0', 'HCN4-3/HCN1-0', 'HNC1-0/HCN1-0', 'HNC4-3/HCN4-3', 'HNC4-3/HCN1-0'],
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HNC-HCN-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 1.0], [-3.0, 1.0], [-3.0, 2.0], [-3.0, 2.0], [-3.0, 3.0], [-3.0, 3.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+4-3/HCO+1-0', 'CN2_1.5-1_1.5/CN1_0.5-0_0.5', 'CS4-3/CS1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HCO+-CN-CS-z-%.1f' % arxvPDR.metallicity,
                                )

plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-3, 1.0], [-3.0, 1.0], [-3.0, 2.0], [-3.0, 2.0], [-3.0, 3.0], [-3.0, 3.0] ], 
                                names          = [ 'MA1', 'MA2', 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [1.0, 2.0], [2.0, 1.0], [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+4-3/HCO+1-0', 'CN2_1.5-1_1.5/CN1_0.5-0_0.5', 'CS4-3/CS1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-HCO+-CN-CS-z-%.1f' % arxvPDR.metallicity,
                                )


plot_utils.bar_plot_line_ratios(
                                arxvPDR, 
                                ranges         = [ [-2, 1.0], [-3.0, 2.0], [-2.0, 2.0], [-2.0, 2.0], ], 
                                names          = [ 'M1', 'M2', 'M3', 'M4'], 
                                names_pos      = [0.6, 0.1],
                                grid_coords    = [ [3.0, 3.0], [3.0, 5.0], [5.5, 3.0], [5.5, 5.0] ],
                                line_ratios    = ['HCO+1-0/CO1-0', 'HCO+1-0/13CO1-0', 'HCN1-0/CO1-0', 'HNC1-0/HCO+1-0', 'HCN4-3/HCO+4-3', 'CN1_0.5-0_0.5/HCN1-0', 'CN2_1.5-1_1.5/HCN1-0', 'CS4-3/13CO1-0'], 
                                Av_max         = 10.0,
                                image_save_dir = image_save_dir,
                                plot_title     = 'bar-plots-lineRatios-misc1-z-%.1f' % arxvPDR.metallicity,
                                gm_v           = numpy.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])/100.0,
                                colors         = [            'k', 'g', 'b', 'c' , 'y',   'm',  'r']
                                )
