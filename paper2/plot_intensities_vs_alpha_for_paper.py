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
#imageSaveDir = '/home/mher/ism/docs/paper02/src/figs'
imageSaveDir = '/home/mher/tmp'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'C+158',
                   'type' : 'pdr'
                 },         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-1, 1], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-1, 1],  
        'cbarTicks' : [-1, -0.5, 0.0, 0.5, 1.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, False, True], [True, True, True] ],
#        'f_interp_dim' : '2D',
#        'interp'    : 'cubic',
        }
plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'C609',
                   'type' : 'pdr'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-1, 1], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range'   : [-1, 1],  
        'cbarTicks'    : [-1, -0.5, 0.0, 0.5, 1.0],
        'xticks'       : [0, 1, 2, 3, 4, 5],
        'yticks'       : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',   
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'C369',
                   'type' : 'pdr'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-1, 1], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range'   : [-1, 1],  
        'cbarTicks'    : [-1, -0.5, 0.0, 0.5, 1.0],
        'xticks'       : [0, 1, 2, 3, 4, 5],
        'yticks'       : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',   
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'O63',
                   'type' : 'pdr'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range'   : [-1, 1],  
        'cbarTicks'    : [-1, -0.5, 0.0, 0.5, 1.0],
        'xticks'       : [0, 1, 2, 3, 4, 5],
        'yticks'       : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, False, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',           
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
#        'get_section_data_using' : '3d-interp',        
#        'get_section_data_using' : 'select',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO2-1',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO3-2',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO4-3',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO6-5',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO7-6',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO10-9',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1.0, 0.0, 1.0, 2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CO16-15',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [0, 6], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [0, 6],  
        'cbarTicks' : [0, 2, 4, 6],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : '13CO1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-3, 3], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-3, 3],  
        'cbarTicks' : [-3, -2, -1, 0, 1, 2, 3],
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : '13CO2-1',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-3, 3], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-3, 3],  
        'cbarTicks' : [-3, -2, -1, 0, 1, 2, 3],       
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : '13CO3-2',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-3, 3], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-3, 3],  
        'cbarTicks' : [-3, -2, -1, 0, 1, 2, 3], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : '13CO6-5',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-3, 3], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-3, 3],  
        'cbarTicks' : [-3, -2, -1, 0, 1, 2, 3], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HCN1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',       
        }
plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HCN4-3',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',       
        }
plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HNC1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HNC4-3',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HCO+1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                
        }
plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HCO+4-3',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                        
        }
plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CN1_0.5-0_0.5',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                                
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CN2_1.5-1_1.5',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[False, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                        
        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CS1-0',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'CS4-3',
                   'type' : 'radex-lvg'
                 },
         
        'Av_use'       :  10.0,
        'relGmech_ref' : -10.0,
        'relGmech'     : [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ],
        'v_range'      : [-2, 2], # range of the values, also that of the cbar
        'res'          : [100.0, 100.0],
        'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'         : matplotlib.cm.jet,
        'cLevels'      : [0],
        'clip'         : None,
        'removeNans'   : False,
        'cbar_range': [-2, 2],  
        'cbarTicks' : [-2, -1, 0, 1, 2], 
        'xticks'    : [0, 1, 2, 3, 4, 5],
        'yticks'    : [0, 1, 2, 3, 4, 5],     
        'show_contours' : [[True, True, True], [True, True, True] ],
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',                        
        
        }

plot_utils.plot_relative_change_intensity_grid2(arxvPDR, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------

print 'done'