#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
from galaxies import fi_utils
import lineDict
import pickle

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'snap_index': numpy.arange(4, 4 + 1, 1),
          'ranges'    : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                         'sph':{
                               'min_log_n_use'  : -3.0,      
                               'min_log_G0_use' : -3.0,
                               'min_log_gm_use' : -50.0,
                               'Av_use'         :  [0.0, 20000000.0],
                               'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                               },
                      
                        #the size of the box to be displayed (particles outside the range are discarded)
                        'box_size' : [-8.0, 8.0] | units.kpc, 
                        },
          'line_ratios' : [
                            [
                              'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO3-2/CO2-1', #low-j ratios 
                              'CO4-3/CO1-0', 'CO6-5/CO1-0', 'CO6-5/CO4-3', #high-j ratios'
                            ], 
                            [
                              '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO3-2/13CO2-1', #low-j ratios 
                               '13CO4-3/13CO1-0', '13CO6-5/13CO1-0', '13CO6-5/13CO4-3', #high-j ratios'
                            ],
                            [
                              '13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO3-2/CO2-1', #low-j ratios 
                              '13CO4-3/CO1-0', '13CO6-5/CO1-0', '13CO6-5/CO4-3', #high-j ratios'
                            ],
                          ],
          'fig_save_path' : '/home/mher/ism/docs/paper02.5/src/figs/line_ratios_vs_beam_size.eps'
        }
################################## Plotting the line for different beam sizes #####################################

## loading the pickle file which holds the luminisity info
fname = os.path.join(params['rundir'],'analysis','line_luminosity_CO_13CO_snap_%d.pkl' % params['snap_index'][0])

fObj = open(fname, 'r')
    
luminosity = pickle.load(fObj )

fObj.close()

print 'read luminisities file :\n\t\t %s' % fname

################################# plotting the line ratios as a function of beam size ######################

fig, axs = pylab.subplots(1, 3, sharex=False, sharey=False, figsize=(12.0, 4.0))

## for ax in axs[:,0] : ax.set_ylabel('y(kpc)')
for ax in axs: ax.set_xlabel('beam size (kpc)', size=10)
axs[0].set_ylabel(r'$\log_{10}$ [line ratio]')

pylab.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.95, wspace=0.25, hspace=0.25)
titles = ['CO/CO', '13CO/13CO', '13CO/CO']
syms = ['r-', 'g-', 'b-', 'r--', 'g--', 'b--']
yranges = [[0.1, 3.0], [0.03, 1.0], [0.005, 0.3]]

for i, bunch in enumerate(params['line_ratios']):

    for j, line_ratio in enumerate(bunch):
        
        line1, line2 = line_ratio.split('/')

        print i, j, line_ratio, line1, line2
        
        beam_size = luminosity['beam_r']*2.0
        
        curve_data = luminosity['lum_r'][line1]/luminosity['lum_r'][line2]

        axs[i].plot(beam_size, curve_data, syms[j])
        
        npts = beam_size.size 
        
        axs[i].text(beam_size[npts/2], curve_data[npts/2], line_ratio, size=9)
        
    print '-----------------'
    
pylab.show()

if params['fig_save_path'] != None:
    fig.savefig(params['fig_save_path'])
    print 'saved image file to :\n\t\t\t %s' % params['fig_save_path']

print 'done'