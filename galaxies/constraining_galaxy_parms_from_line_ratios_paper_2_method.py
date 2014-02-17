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
from mylib.utils.misc import scale
from galaxies import fi_utils
#------------------------------------------------------------------------------------

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
#imageSaveDir = '/home/mher/ism/docs/paper02/src/figs/constraining'
imageSaveDir = '/home/mher/tmp/'

#---------------------------------------------------------------------------------------
#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#--------------------------------------------------------------------------------------------------------------------------
parms = {
        'ratios' : 
        (
         'CO2-1/CO1-0' ,
         'CO3-2/CO1-0' ,
         'CO4-3/CO1-0' ,
         'CO5-4/CO1-0' ,
         'CO6-5/CO1-0' ,
         ),
        'zsecs': (1e-10, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
         
        #'ranges'       : [[3.0, 6.0], [3.0, 6.0]],
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
        'Av_use'     :  5.0,
        'c_levels'  : {
                       #'values' : numpy.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50]),
                       'values' : numpy.log10([0.01, 0.05, 0.1, 0.5,  1,  5,  10,   50,   100]),
                       'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                       #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                      },
        'em_unit'   : 'Kkms',
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

obs_ratios_all = []

obs_ratios = line_ratio_utils.ratios()
obs_ratios.add_ratio('CO2-1/CO1-0'  , mn=0.88 , mx=0.92)
obs_ratios_all.append(obs_ratios) 

obs_ratios = line_ratio_utils.ratios()
obs_ratios.add_ratio('CO3-2/CO1-0' , mn=0.85 , mx=0.87)
obs_ratios_all.append(obs_ratios) 

obs_ratios = line_ratio_utils.ratios()
obs_ratios.add_ratio('CO4-3/CO1-0' , mn=0.78 , mx=0.82)  
obs_ratios_all.append(obs_ratios) 

obs_ratios = line_ratio_utils.ratios()
obs_ratios.add_ratio('CO6-5/CO1-0' , mn=0.35, mx=0.5)
obs_ratios_all.append(obs_ratios) 

def make_figure_for_Av(p, fname):
        
    #-------------
    p.set_obs_ratios(obs_ratios_all[0])
    fig1, axs1 = p.fit_observations_visual()
    
    p.set_obs_ratios(obs_ratios_all[1])
    fig2, axs2 = p.fit_observations_visual()
    
    p.set_obs_ratios(obs_ratios_all[2])
    fig3, axs3 = p.fit_observations_visual()
    
    p.set_obs_ratios(obs_ratios_all[3])
    fig4, axs4 = p.fit_observations_visual()
    
    #merge the axes into one figure (hstack axes)
    #--------------------------------------------
    axs = numpy.vstack((axs1, axs2, axs3, axs4))
    
    fig, axs_new = pylab.subplots(axs.shape[0], axs.shape[1], sharex = True, sharey = True, figsize = (8, 4),
                                  subplot_kw = {'xlim':[parms['ranges'][0][0], parms['ranges'][0][1] ],  #keywords passed to figure.add_subplot()      
                                                'ylim':[parms['ranges'][1][0], parms['ranges'][1][1] ],
                                                #'aspect':'equal',
                                                #'adjustable':'datalim',
                                                'autoscale_on' : False, 
                                                }
                                  )
    
    xs, xe = 0.07, 0.7
    ys, ye = 0.15, 0.9
    pylab.subplots_adjust(left=xs, bottom=ys, right=xe, top=ye, wspace=0.05, hspace=0.1)       
    
    for r, axs_r in enumerate(axs_new):
        
        for c, ax in enumerate(axs_r):
            
            ax_old = axs[r,c]
            
            #copying the ticks 
            ax.set_xticks( ax_old.get_xticks() )
            ax.set_yticks( ax_old.get_yticks() )
            
            #copying the xtick labels
            ax.set_xticklabels( [l.get_text() for l in ax_old.get_xticklabels()] )
            ax.set_yticklabels( [l.get_text() for l in ax_old.get_xticklabels()] )
            
            #copying and setting the image data
            #data_this = ax_old.get_images()[0].get_array()        
            #ax.imshow(data_this, origin='lower', extent = numpy.array(parms['ranges']).flatten(), vmin=0, vmax=4)
    
            #computing and setting the image data
            if r == 0:
                data_show = ax_old.get_images()[0].get_array()
            else:
                data_this_col_upto_this_row = []
                
                for k, ax_this_col in enumerate(axs[0:r+1, c]):
                    data = ax_this_col.get_images()[0].get_array()
                    #data *= (k+1)  #if this is uncommented, then the maximum values of all the intesections would be factorial(n_ratios_ 
                    data_this_col_upto_this_row.append(data)
                #
                
                data_show = numpy.array(data_this_col_upto_this_row).sum(axis=0)
            
            print numpy.unique(data_show)
            #scaling the data to be displayed to the highes value corresponding to the color
            #which is visually distict from the one which will be assigned to the maximum
            data_show_scaled = scale(data_show, 0.0, 0.75, 0.0, 4.0)
     
            inds = numpy.where(data_show == 4)
            if inds[0].size > 0:
                data_show_scaled[inds] = 1.0
 
            ax.imshow(data_show_scaled, origin='lower', extent = numpy.array(parms['ranges']).flatten(), 
                      vmin=0, vmax=1.0, cmap=matplotlib.cm.jet)
           
            #displaying the title at the top of the subplots in the top row 
            if r == 0:
                
                title_str = r'  $\alpha$' + '\n' + '%.2f' % parms['zsecs'][c] #alpha string
                    
                ax.set_title(title_str, 
                             size='small', rotation=0, 
                             horizontalalignment='center', 
                             verticalalignment='bottom')
    
            #printing the labels (n) at the bottom row 
            if r==axs.shape[0]-1:
                ax.set_xlabel(r'$\log_{10}(n)$')
            if c==0:
                ax.set_ylabel(r'$\log_{10}(G_0)$')
            
            if c==axs.shape[1]-1:
                
                def make_string(r):
                    if r != axs.shape[0]-1:
                        return obs_ratios_all[r].keys()[0]
                    else:
                        strng = ''
                        for i, key in enumerate(obs_ratios_all[-1]):
                            strng += key 
                            if (i-1) % 2 == 0 and i > 0:
                                strng += '\n'
                            else:
                                strng += ','
                        #
                        return strng
                    #
                #
                 
                bb = ax.get_position().bounds
                
                pylab.figtext(xe + 0.01, bb[1] + 0.02, make_string(r), size='small')
    #        
    
    fig.savefig(os.path.join(imageSaveDir,fname))
    
    pylab.draw()
    print 'done' 
#


make_figure_for_Av(p1, 'Av5.eps')
make_figure_for_Av(p2, 'Av10.eps')
make_figure_for_Av(p3, 'Av30.eps')
