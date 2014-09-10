# plots the emissions grid
#------------------------------------------------------------------------------------
import os
import matplotlib

if 'particle3' in os.uname():
    matplotlib.use('Qt4Agg')
    
import meshUtils
import plot_utils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
#imageSaveDir = '/home/mher/ism/docs/paper02/src/figs'
imageSaveDir = '/home/mher/tmp/foo'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'C+158',
                   'type' : 'pdr'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cbar_range': [-6, -2], 
        'cLevels'   : [-5.0, -4.0, -3.0, -2.5],
        'cbarTicks' : [-6.0, -5.0, -4.0, -3.0, -2.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],
        'clip'      : None,
        'removeNans': False,
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'C609',
                   'type' : 'pdr'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cbar_range': [-6, -2], 
#        'cbar_range': [-5.5, -4.5], 
        'cLevels'   : [-5.0, -4.8, -4.7, -4.6],
        'cbarTicks' : [-6.0, -4.0, -2.0, 0.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],
        'clip'      : None,
        'removeNans': False,
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'C369',
                   'type' : 'pdr'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cbar_range': [-6, -2],  
        'cLevels'   : [-5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8],
        'cbarTicks' : [-6.0, -4.0, -2.0, 0.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],
        'clip'      : None,
        'removeNans': False,
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'O63',
                   'type' : 'pdr'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cmap'      : matplotlib.cm.jet,
        'cbar_range': [-6, 0.0],  
        'cLevels'   : [-6.0, -5.0, -4.0, -3.0, -2.0, -1.0],
        'cbarTicks' : [-6.0, -4.0, -2.0, 0.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],
        'clip'      : None,
        'removeNans': False,
        'f_interp_dim' : '2D',
        'interp'    : 'cubic',
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -7.0, -6.5, -6.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : '13CO1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -7.0, -6.8, -6.6],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO2-1',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -7.0, -6.0, -5.5, -5.2],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO3-2',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.5, -5.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO4-3',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.5, -5.0, -4.5],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO6-5',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.5, -5.0, -4.5, -4.0],        
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO7-6',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.5, -5.0, -4.5, -4.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO10-9',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -4.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,
        'clip_max_n': 4.0,
        'removeNans': True,                  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CO16-15',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -4.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,
        'clip_max_n': 3.5,
        'removeNans': True,                  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : '13CO2-1',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -7.0, -6.0, -5.5],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }

grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : '13CO3-2',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.5, -5.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : '13CO6-5',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-10, -3.0],
        'cbarTicks' : [-10.0, -8.0, -6.0, -4.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-10.0, -8.0, -6.0, -5.0, -4.0],
        'clip'      : [-10.0, 0.0],
        'removeNans': False,  
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'HCN1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.8, -6.6],
        'clip'      : [-12.0, 0.0],
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'HCN4-3',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.0, -5.0],        
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'HNC1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.5, -8.0, -7.0, -7.5],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'line' : {
                   'code' : 'HNC4-3',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.0, -5.0],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'HCO+1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.7],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'HCO+4-3',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.0, -5.0],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CN1_0.5-0_0.5',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.5],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CN2_1.5-1_1.5',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -8.0, -7.0, -6.0],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CS1-0',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -9.0, -8.5, -8.0, -7.0, -6.0],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        'clip_max_n': 4.0,
        'removeNans': True,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
parms = {
         
        'line' : {
                   'code' : 'CS4-3',
                   'type' : 'radex-lvg'
                 },
        'Av_use'    : 10.0,
        'z_sec'     : -10.0,
        'res'       : [100.0, 100.0],
        'ranges'    : [[0.0, 6.0], [0.0, 6.0]],
        'cbar_range': [-12, -5.0],
        'cbarTicks' : [-12.0, -10.0, -8.0, -6.0],
        'xticks'    : [0, 1, 2, 3, 4, 5, 6],
        'yticks'    : [0, 1, 2, 3, 4, 5, 6],                          
        'cmap'      : matplotlib.cm.jet,
        'cLevels'   : [-12.0, -10.0, -9.0, -8.0, -7.0, -6.5],
        'clip'      : [-12.0, 0.0],
        'removeNans': False,        
        'clip_max_n': 4.0,
        'removeNans': True,        
        }
grd = arxvPDR.get_emission_grid_from_databases(**parms)
plot_utils.plot_intensity_grid2(grd, imageSaveDir=imageSaveDir, **parms)
#--------------------------------------------------------------------------------------------------------------------------
'''
print 'done'