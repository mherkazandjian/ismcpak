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
#dirPath      = os.path.join(home, 'ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-sweep/')
#imageSaveDir = '/home/mher/ism/docs/paper02/src/figs'
imageSaveDir = '/home/mher/tmp/foo'

#---------------------------------------------------------------------------------------
#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True,
#                             grid_qx = ['hdr','nGas'],
#                             grid_qy = ['hdr','G0'],
#                             grid_qz = ['from_meshes_info', 'parms', 3, 'CR_rate'],  # 3 indicates the 4th column in self.infoAll['parms']                             
                             )

#--------------------------------------------------------------------------------------------------------------------------
parms = {
        #'ratios' : ('CO1-0/13CO1-0', 'CO2-1/13CO2-1', 'CO3-2/13CO3-2', 'CO3-2/13CO3-2', 'CO4-3/13CO4-3', 
        #            'CO5-4/13CO5-4', 'CO6-5/13CO6-5', 'CO7-6/13CO7-6', 'CO8-7/13CO8-7', 'CO9-8/13CO9-8',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #----------------------------------------------------         
        #'ratios' : ('CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO6-5/CO1-0', 'CO7-6/CO1-0', 
        #            'CO8-7/CO1-0', 'CO9-8/CO1-0', 'CO10-9/CO1-0', 'CO11-10/CO1-0', 'CO12-11/CO1-0',
        #            'CO13-12/CO1-0',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #----------------------------------------------------         
        #'ratios' : ('HCN1-0/CO1-0', 'HCN2-1/CO2-1', 'HCN3-2/CO3-2', 'HCN4-3/CO4-3', 
        #            'HCN5-4/CO5-4', 'HCN6-5/CO6-5', 'HCN7-6/CO7-6', 'HCN8-7/CO8-7', 'HCN9-8/CO9-8',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #------------------------------------------------
        #'ratios' : ('HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3', 
        #            'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #------------------------------------------------        
        #'ratios' : ('HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3', 
        #            'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #------------------------------------------------                
        #'ratios' : ('HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3', 
        #            'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #-----------------------------------------------------------
        #'ratios' : ('CS1-0/13CO1-0', 'CS2-1/13CO1-0', 'CS3-2/13CO1-0', 'CS3-2/13CO1-0', 'CS4-3/13CO1-0', 
        #            'CS5-4/13CO1-0', 'CS6-5/13CO1-0', 'CS7-6/13CO1-0', 'CS8-7/13CO1-0', 'CS9-8/13CO1-0',),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #-----------------------------------------------
        #'ratios' : ('HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2',),     
        #'zsecs': (1e-10, 0.01,),
        #-----------------------------------------------
        #'ratios' : (
        #            'CO7-6/HCN1-0', 'CO8-7/HCN1-0', 'CO9-8/HCN1-0',
        #            ),
        #                    'CO7-6/HCN2-1', 'CO8-7/HCN2-1', 'CO9-8/HCN2-1',
        #                    'CO7-6/HCN3-2', 'CO8-7/HCN3-2', 'CO9-8/HCN3-2',
        #                    'CO7-6/HCN4-3', 'CO8-7/HCN4-3', 'CO9-8/HCN4-3',     
        #                   ),
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #-------------------------------------------------------------------------------------------------------
        #'ratios' : ('CO1-0/HCO+1-0', 'CO2-1/HCO+2-1', 'CO3-2/HCO+3-2', 'CO4-3/HCO+4-3',
        #            'CO1-0/HNC1-0' , 'CO2-1/HNC2-1' , 'CO3-2/HNC3-2' , 'CO4-3/HNC4-3',
        #            'CO1-0/HCN1-0' , 'CO2-1/HCN2-1' , 'CO3-2/HCN3-2' , 'CO4-3/HCN4-3',
        #            ),     
        #'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
        #-------------------------------------------------------------------------------------------------------     
        'ratios' : 
        (
          #'CO1-0/HCO+1-0', 'CO4-3/HCO+4-3',
          'HCO+1-0/CO1-0',
          'HCN1-0/CO1-0',
          'HNC1-0/CO1-0',
          'CO1-0/HNC1-0',
          'CO1-0/HCN1-0' , 'CO3-2/HCN3-2' , 'CO4-3/HCN4-3',
          'HCN1-0/HNC1-0',
          'HCN1-0/HCO+1-0',
          'HNC1-0/HCO+1-0',
          'HCN4-3/HCO+4-3',
          'HCO+4-3/HCO+1-0',
          #'HCO+7-6/HCO+1-0',
          #'CO2-1/CO1-0', 
          #'CO3-2/CO1-0', 
          #'CO4-3/CO3-2', 
          #'CO6-5/CO3-2', 
          #'CO7-6/CO3-2', 
          #'CO8-7/CO3-2', 
          #'CO9-8/CO3-2', 
          #'CO10-9/CO3-2', 
          #'CO11-10/CO3-2', 
          #'CO12-11/CO3-2',
          #'CO13-12/CO3-2',
          'HCO+4-3/CO1-0',
         ),
        'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
         
        ################################COSMIC RAY PARMS############################################################
        #'ratios' : ('HCN1-0/HCO+1-0', 'HCN2-1/HCO+2-1', 'HCN3-2/HCO+3-2', 'HCN4-3/HCO+4-3', 
        #            'HCN5-4/HCO+5-4', 'HCN6-5/HCO+6-5', 'HCN7-6/HCO+7-6', 'HCN8-7/HCO+8-7', 'HCN9-8/HCO+9-8',),
        #'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
        #-----------------------------------------------------------------------------------------------------------
        #'ratios' : ('HNC1-0/HCO+1-0', 'HNC2-1/HCO+2-1', 'HNC3-2/HCO+3-2', 'HNC4-3/HCO+4-3', 
        #            'HNC5-4/HCO+5-4', 'HNC6-5/HCO+6-5', 'HNC7-6/HCO+7-6', 'HNC8-7/HCO+8-7', 'HNC9-8/HCO+9-8',),
        #'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
        #-----------------------------------------------------------------------------------------------------------
        #'ratios' : ('HCN1-0/HNC1-0', 'HCN2-1/HNC2-1', 'HCN3-2/HNC3-2', 'HCN4-3/HNC4-3', 
        #            'HCN5-4/HNC5-4', 'HCN6-5/HNC6-5', 'HCN7-6/HNC7-6', 'HCN8-7/HNC8-7', 'HCN9-8/HNC9-8',),
        #'zsecs': (5e-17, 5e-16, 5e-15, 5e-14, 5e-13),
        #'zsecType' : 'zeta',
        #-----------------------------------------------------------------------------------------------------------
        
        'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
        #'ranges'       : [[0.0, 6.0], [0.0, 6.0]],
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

p1 = plot_utils.line_ratio_grid_of_grid(arxvPDR, imageSaveDir=imageSaveDir, **parms)

parms['Av_use'] = 10.0
p2 = plot_utils.line_ratio_grid_of_grid(arxvPDR, imageSaveDir=imageSaveDir, **parms)

parms['Av_use'] = 30.0
p3 = plot_utils.line_ratio_grid_of_grid(arxvPDR, imageSaveDir=imageSaveDir, **parms)

pylab.show()

########################################################################################

import line_ratio_utils


''' 
#----------------------reading line ration---------------------------------------------
#path to the excel sheet containing the data
data_sheet = '/home/mher/ism/marissa/NGC_253_data/NGC_253_Fluxes_Master_newsheet.xlsx' 
ratios = (
          'CO1-0/HCO+1-0', 
          #'CO4-3/HCO+4-3',
          
          #'CO1-0/HNC1-0',
          
          #'CO1-0/HCN1-0' , 
          #'CO3-2/HCN3-2' , 
          #'CO4-3/HCN4-3',
          
          'HCN1-0/HNC1-0',
          
          #'HCN1-0/HCO+1-0',
          'HCN4-3/HCO+4-3',
          )

#reading the lines from the excel sheet
obs = line_ratio_utils.read_observations(data_sheet)

#making the observed line ratio
obs_ratios = line_ratio_utils.ratios()
obs_ratios.make_ratios(obs, ratios)
#---------------------------------------------------------------------------------
'''


#setting line ratios manually
#---------------------------


if True:
    #Leonen 2008 narrow sample (most populated reigon) (best contraining case)    
    obs_ratios = line_ratio_utils.ratios()
    obs_ratios.add_ratio('HCN1-0/HNC1-0'  , mn=1.5 , mx=4.0)    
    obs_ratios.add_ratio('HCN1-0/HCO+1-0' , mn=0.6 , mx=3.2)   
    obs_ratios.add_ratio('HNC1-0/HCO+1-0' , mn=0.3 , mx=1.0)   
    #NGC235
    obs_ratios.add_ratio('HCO+4-3/CO1-0'  , mn=0.1, mx=0.5)       #good
    

if False:
    #Leonen 2008 narrow sample (most populated reigon)    
    obs_ratios = line_ratio_utils.ratios()
    obs_ratios.add_ratio('HCN1-0/HNC1-0'  , mn=1.5 , mx=4.0)    
    obs_ratios.add_ratio('HCN1-0/HCO+1-0' , mn=0.6 , mx=3.2)   
    obs_ratios.add_ratio('HNC1-0/HCO+1-0' , mn=0.3 , mx=1.0)   
    obs_ratios.add_ratio('HCO+1-0/CO1-0'  , mn=0.03, mx=0.1)
    obs_ratios.add_ratio('HNC1-0/CO1-0'   , mn=0.01, mx=0.1)
    obs_ratios.add_ratio('HCN1-0/CO1-0'   , mn=0.03, mx=0.1)

if False:
    #Leonen 2008 narrow sample (most populated reigon)    
    obs_ratios = line_ratio_utils.ratios()
    
    obs_ratios.add_ratio('HCN1-0/HNC1-0'   , mn=1.5 , mx=4.0)
    obs_ratios.add_ratio('HCN1-0/CO1-0'    , mn=0.03 , mx=0.3)
    obs_ratios.add_ratio('HCN1-0/CO1-0'    , mn=0.01 , mx=0.1)
    
       
    #obs_ratios.add_ratio('HNC1-0/HCO+1-0' , mn=0.3 , mx=1.0)   
    #obs_ratios.add_ratio('HCO+1-0/CO1-0'  , mn=0.03, mx=0.1)
    #obs_ratios.add_ratio('HNC1-0/CO1-0'   , mn=0.01, mx=0.1)
    #obs_ratios.add_ratio('HCN1-0/CO1-0'   , mn=0.03, mx=0.1)


#Leonen 2008 full sample range 
#obs_ratios.add_ratio('HCN1-0/HNC1-0'  , mn=0.6, mx=7.5)
#obs_ratios.add_ratio('HCN1-0/HCO+1-0' , mn=0.3 , mx=4.5)
#obs_ratios.add_ratio('HNC1-0/HCO+1-0' , mn=0.2 , mx=3.2)

#obs_ratios.add_ratio('CO4-3/CO3-2', mn=0.1 , mx=1.0)
#obs_ratios.add_ratio('CO7-6/CO3-2', mn=1.0 , mx=2.0)



#-------------
p1.set_obs_ratios(obs_ratios)
fig1, axs1 = p1.fit_observations_visual()

p2.set_obs_ratios(obs_ratios)
fig2, axs2 = p2.fit_observations_visual()

p3.set_obs_ratios(obs_ratios)
fig3, axs3 = p3.fit_observations_visual()

#merge the axes into one figure (hstack axes)
#--------------------------------------------
axs = numpy.vstack((axs1, axs2, axs3))

fig, axs_new = pylab.subplots(axs.shape[0], axs.shape[1], sharex = True, sharey = True, figsize = (9, 4),
                              subplot_kw = {'xlim':[parms['ranges'][0][0], parms['ranges'][0][1] ],  #keywords passed to figure.add_subplot()      
                                            'ylim':[parms['ranges'][1][0], parms['ranges'][1][1] ],
                                            #'aspect':'equal',
                                            #'adjustable':'datalim',
                                            'autoscale_on' : False, 
                                            }
                              )
pylab.subplots_adjust(left=0.2, bottom=0.25, right=0.8, top=0.7, wspace=0.0, hspace=0.0)       
#ax.set_autoscale_on(False)

for r, axs_r in enumerate(axs_new):
    
    for c, ax in enumerate(axs_r):
        
        ax_old = axs[r,c]
        
        #copying the ticks 
        ax.set_xticks( ax_old.get_xticks() )
        ax.set_yticks( ax_old.get_yticks() )
        
        #copying the xtick labels
        ax.set_xticklabels( [l.get_text() for l in ax_old.get_xticklabels()] )
        ax.set_yticklabels( [l.get_text() for l in ax_old.get_xticklabels()] )
        
        #copying the image data
        ax.imshow(ax_old.get_images()[0].get_array().T, origin='lower', extent = numpy.array(parms['ranges']).flatten())

        
#fig.savefig('/home/mher/foo.eps')

pylab.draw()
print 'done' 


