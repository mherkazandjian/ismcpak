#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from galaxies import fi_utils

home = '/home/mher'

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

snaps = numpy.arange(4, 4 + 1, 1)

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'latex'    : True,
          'latex_dir': os.path.join(home, 'ism','docs', 'paper02.5', 'src', 'figs'),
          'all_maps' : {
                  'map1' : {
                            'attr' : 'em_fluxKkms_CO1-0',
                            'title': r'$\log_{10}[ Flux ({^{12}{\rm CO}(1-0)} / K.km.s^{-1}) ]$',
                            'v_rng'   : [-10.0, 4.0],
                           },
                  'map2' : {
                            'attr' : 'em_fluxKkms_13CO1-0',
                            'title': r'$\log_{10}[ Flux ({^{13}{\rm CO}(1-0)} / K.km.s^{-1}) ]$',
                            'v_rng'   : [-10.0, 4.0],                             
                           },
                      },
          'save_image' : True,          
        }

for snapIndex in snaps:    
    fi_utils.plot_all_maps_for_snapshot_from_saved_data(snapIndex, params)

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

snaps = numpy.arange(20, 20 + 1, 1)


params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'latex'    : True,
          'latex_dir': os.path.join(home, 'ism','docs', 'paper02.5', 'src', 'figs'),
          'all_maps' : {
                  'map1' : {
                            'attr' : 'em_fluxKkms_CO1-0',
                            'title': r'$\log_{10}[ Flux({^{12}{\rm CO}(1-0)} / K.km.s^{-1}) ]$',
                            'v_rng'   : [-10.0, 4.0],
                           },
                  'map2' : {
                            'attr' : 'em_fluxKkms_13CO1-0',
                            'title': r'$\log_{10} Flux({^{13}{\rm CO}(1-0)} / K.km.s^{-1})$',
                            'v_rng'   : [-10.0, 4.0],                             
                           },
                      },
          'save_image' : True,
        }

for snapIndex in snaps:
    fi_utils.plot_all_maps_for_snapshot_from_saved_data(snapIndex, params)

