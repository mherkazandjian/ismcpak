#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab
import scipy

from numpy import *
from pylab import *

from amuse.units import units
from mylib.utils.misc  import default_logger 
from galaxies import fi_utils
import meshUtils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'
#home = os.path.join('/net', os.environ['HOST'], 'data2', 'mher')

params = {
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100',  # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : [],#'CO'], #'13CO', 'HCN', 'HNC', 'HCO+'],
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'original-only', #by-number', #'matched',  #'original-only' ,#None ,#by-number          
          'snap_index' : numpy.arange(4, 4 + 1, 1),
          'ranges'      : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                          'sph':{
                                 'min_log_n_use'  : -3.0,      
                                 'min_log_G0_use' : -3.0,
                                 'min_log_gm_use' : -50.0,
                                 'Av_use'         :  [0.0, 20000000.0],
                                 'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                                },
                          
                          #the size of the box to be displayed (particles outside the range are discarded)
                          'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                         },
          'check'   : 'default',
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
        }

#fig_save_path1 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/luminosities.eps'
#fig_save_path2 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/pdfs.eps'
#fig_save_path3 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/line_ratio_grid_HCN1-0_HNC1-0.svg'
#fig_save_path4 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/line_ratio_grid_HCN1-0_HCO+1-0.svg'

#fig_save_path1 = None
#fig_save_path2 = None
#fig_save_path3 = None
#fig_save_path4 = None

fig_save_path_app_1 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/pdfs_app1.eps'
fig_save_path_app_2 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/pdfs_app2.eps'
fig_save_path_app_3 = '/home/mher/ism/docs/paper04/src/figs/results/synthetic_pdfs/luminosities_app3.eps'

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

################################################################################################
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

F = arxvPDR.load_interp_func(info={'source':'radex'}, line='HCN1-0', quantity='fluxKkms')

########################

nmin, nmax = -3, 6.0

import paper_plots

#paper_plots.synthetic_flux_contribution_from_pdf(arxvPDR, F, params, -3.0, 6.0, yrng=[1e-8, 1e2], fig_save_path=fig_save_path1)
#paper_plots.plot_PDFs_of_synthetic_luminosity_distribututions(-3.0, 8.0, fig_save_path=fig_save_path2)
#paper_plots.line_ratio_grid_from_pdf_sweep(arxvPDR, params, nmin, nmax, fig_save_path=fig_save_path3)
#paper_plots.line_ratio_grid_from_pdf_sweep_HCOP_HCN(arxvPDR, params, nmin, nmax, fig_save_path=fig_save_path4)


## appendix plots
#paper_plots.plot_PDFs_of_synthetic_luminosity_distribututions_app1(-3.0, 8.0, fig_save_path=fig_save_path_app_1)
#paper_plots.plot_PDFs_of_synthetic_luminosity_distribututions_app2(-3.0, 8.0, fig_save_path=fig_save_path_app_2)

paper_plots.synthetic_flux_contribution_from_pdf_app3(arxvPDR, F, params, -3.0, 6.0, yrng=[1e-8, 1e2], fig_save_path=fig_save_path_app_3)


print 'done'
