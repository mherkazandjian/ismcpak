# plots 4D plots (grid of grids) of line ratios as a functio of gmech and transitions
#------------------------------------------------------------------------------------
import os
import matplotlib
import numpy
import sys
from PyQt4 import QtGui
matplotlib.use('Qt4Agg')
import pylab

import meshUtils
import plot_utils
#------------------------------------------------------------------------------------

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
#image_save_dir = '/home/mher/ism/docs/paper02/src/figs'
image_save_dir = '/home/mher/tmp/figs/'

#---------------------------------------------------------------------------------------
#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
if True:
    parms = {
            'fname'  : 'CO-CO_13CO-13CO_CO-13CO_AV10.eps',            
            'ratios' : (#CO/CO line ratios
                        'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO16-15/CO1-0',
                        'CO7-6/CO3-2', 'CO10-9/CO7-6', 'CO16-15/CO10-9',
                        #13CO/13CO line ratios
                        '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO16-15/13CO1-0',
                        '13CO7-6/13CO3-2', '13CO10-9/13CO7-6', '13CO16-15/13CO10-9',
                        #CO/13CO line ratios
                        'CO1-0/13CO1-0', 'CO2-1/13CO2-1', 'CO3-2/13CO3-2', 'CO4-3/13CO4-3',
                        'CO7-6/13CO6-5', 'CO10-9/13CO10-9', 'CO16-15/13CO16-15',                    
                        ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            #-----------------------------------------------------------------------------------------------------------
            
            'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01,   0.1,  1,    10,    100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            'fig'       :  {'kwargs':{
                                      'figsize' : (5, 9),
                                      },
                            },
            'axs'       :   {
                             'left':0.09, 'bottom':0.05, 'w':0.5, 'h':0.8
                            },
            'cbar'      :  {
                            'range' : [-2, 2],
                            'scale' : 0.7,
                            'sz'    : 0.02,
                            'space' : 0.04,
                            'ticks' : [-2, -1, 0, 1, 2],
                           }, 
            'xticks'    : [1, 3, 5],
            'yticks'    : [1, 3, 5],            
            }
    
    p = plot_utils.line_ratio_grid_of_grid2(arxvPDR, image_save_dir=image_save_dir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
if True:
    parms = {
            'fname'  : 'HCO+-CO_HCN-CO_CS-HCO+-CN-HCN_AV10.eps',                         
            'ratios' : (#HCO+/CO line ratios
                        'HCO+1-0/CO1-0', 'HCN1-0/CO1-0', 'CS1-0/HCO+1-0', 'CS4-3/HCO+4-3',
                        'CN1_0.5-0_0.5/HCN1-0', 'CN2_1.5-1_1.5/HCN1-0',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            #-----------------------------------------------------------------------------------------------------------
            
            'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01,   0.1,  1,    10,    100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            'fig'       :  {'kwargs':{
                                      'figsize' : (5, 5),
                                      },
                            },
            'axs'       :   {
                             'left':0.145, 'bottom':0.135, 'w':0.6, 'h':0.67
                            },
            'cbar'      :  {
                            'range' : [-2, 2],
                            'scale' : 0.7,
                            'sz'    : 0.02,
                            'space' : 0.08,
                            'ticks' : [-2, -1, 0, 1, 2],
                           }, 
            'xticks'    : [1, 3, 5],
            'yticks'    : [1, 3, 5],
            }
    
    p = plot_utils.line_ratio_grid_of_grid2(arxvPDR, image_save_dir=image_save_dir, **parms)

#--------------------------------------------------------------------------------------------------------------------------
if True:
    parms = {
            'fname'  : 'HCN-HCO+_AV5.eps',                         
            'ratios' : (#HCN/HCO+ line ratios
                        'HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3',   
                        'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  5.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HCO+'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HCN-HCO+_AV10.eps',                         
            'ratios' : (#HCN/HCO+ line ratios
                        'HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3',   
                        'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HCO+'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HCN-HCO+_AV30.eps',                         
            'ratios' : (#HCN/HCO+ line ratios
                        'HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3',   
                        'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  30.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HCO+'), **parms)
#--------------------------------------------------------------------------------------------------------------------------
if True:
    parms = {
            'fname'  : 'HCN-HNC_AV5.eps',                         
            'ratios' : (#HCN/HNC line ratios
                        'HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3',   
                        'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  5.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HNC'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HCN-HNC_AV10.eps',                         
            'ratios' : (#HCN/HNC line ratios
                        'HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3',   
                        'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'    : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HNC'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HCN-HNC_AV30.eps',                         
            'ratios' : (#HCN/HNC line ratios
                        'HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3',   
                        'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),            
            'ranges'    : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'    :  30.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HCN-HNC'), **parms)

#--------------------------------------------------------------------------------------------------------------------------
if True:
    parms = {
            'fname'  : 'HNC-HCO+_AV5.eps',                         
            'ratios' : (#HNC/HCO+ line ratios
                        'HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3',   
                        'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',
                       ),
            'zsecs'  : (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges' : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'   : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'    : [100, 100],
            'cbar'   : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  5.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HNC-HCO+'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HNC-HCO+_AV10.eps',                         
            'ratios' : (#HNC/HCO+ line ratios
                        'HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3',   
                        'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'    :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HNC-HCO+'), **parms)
    #--------------------
    parms = {
            'fname'  : 'HNC-HCO+_AV30.eps',                         
            'ratios' : (#HNC/HCO+ line ratios
                        'HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3',   
                        'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',
                       ),
            'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),            
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  30.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','gmech','HNC-HCO+'), **parms)

#-------------------------------------------------------------------------------------------------------------------------
if True:
    home         = '/home/mher'
    dirPath      = os.path.join(home, 'ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-sweep/')
    
    #---------------------------------------------------------------------------------------
    #reading and setting up the pdr database
    arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True,
                                 grid_qx = ['hdr','nGas'],
                                 grid_qy = ['hdr','G0'],
                                 grid_qz = ['from_meshes_info', 'parms', 3, 'CR_rate'],  # 3 indicates the 4th column in self.infoAll['parms']                             
                                 )
    
    #--------------------------------------------------------------------------------------------------------------------------
    parms = {
            'fname'  : 'HCN-HNC_AV10.eps',
            'ratios' : ('HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3', 
                        'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',),
            'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
            'zsecType' : 'zeta',            
            #-----------------------------------------------------------------------------------------------------------
            
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','CR') , **parms)
    #--------------------------------------------------------------------------------------------------------------------------
    parms = {
            'fname'  : 'HNC-HCO+_AV10.eps',
            'ratios' : ('HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3', 
                        'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',),
            'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
            'zsecType' : 'zeta',            
            #-----------------------------------------------------------------------------------------------------------
            
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','CR') , **parms)
    #--------------------------------------------------------------------------------------------------------------------------
    parms = {
            'fname'  : 'HCN-HCO+_AV10.eps',             
            'ratios' : ('HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3',   
                        'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',),
            'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
            'zsecType' : 'zeta',                        
            #-----------------------------------------------------------------------------------------------------------
            
            'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
            'cmap'      : {'obj'     : matplotlib.cm.jet,
                           'v_range' : [-2, 2],
                           'log'     : True,
                          }, 
            'res'       : [100, 100],
            'cbar'      : {
                           'ticks'  : numpy.linspace(-2, 2, 5),
                           'format' : '10^x',  # None, '10^x', '%e' (or any format)
                          },
            'Av_use'       :  10.0,
            'c_levels'  : {
                           #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                           'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          },
            'clip'      : None,
            }
    
    p = plot_utils.line_ratio_grid_of_grid(arxvPDR, image_save_dir=os.path.join(image_save_dir,'line-ratio-grid-grids','CR') , **parms)

#------------------------------------------------------------------------------------------------------
pylab.show()
print 'done' 

