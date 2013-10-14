#########################################################################################################
import time, sys, os
from IPython.parallel import Client

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
from numpy import log10
import pylab

from amuse.units import units, constants, nbody_system
import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

snaps = numpy.arange(4, 4 + 1, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'maps'   : {
                      'attr' : 'em_fluxKkms_13CO1-0', #'mass', 'G0', 'gmech', 'Av'
                      'title': r'$L({^{13}{\rm CO}(1-0)} / K.km.s^{-1})$', 
                     },
        }
#############################################################################################################
#############################################################################################################
#############################################################################################################

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

def plot_maps(snapIndex, params):
    
    #path to processed fi snapshot  
    filename = os.path.join(params['rundir'],'analysis', 'fiout.%06d.%s.npz' % (snapIndex, params['maps']['attr']))  
    
    loaded_data = numpy.load(filename)
    
    map_data, map_params = loaded_data['map_data'], loaded_data['params'].tolist()

    #getting the time of the snapshot 
    runinfo = fi_utils.parse_old_runinfo_file(params['rundir'] + '/firun/runinfo')
    snap_time = (float(runinfo['dtime'])*timeUnit.number) * (float(runinfo['noutbod']) * snap)

    bs_min, bs_max = map_params['ranges']['box_size'].number
    
    #displaying all the maps in a single plot
    fig = pylab.figure(0, figsize=(8,8))
    ax = fig.add_axes([0.15, 0.085, 0.75, 0.75])
    
    ax.set_xlabel('x(kpc)', size='large')    
    ax.set_ylabel('y(kpc)', size='large')
    ax.set_xlim([bs_min, bs_max])
    ax.set_ylim([bs_min, bs_max])
    
    im = ax.imshow(map_data,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=map_params['maps']['v_rng'][0],  
                   vmax=map_params['maps']['v_rng'][1], 
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

    ax.set_title(params['maps']['title'], size='large')
        
    cbar_ax = fig.add_axes([0.2, 0.9, 0.6, 0.02]) 

    ax.tick_params(axis='both', which='major', labelsize=20)
    #ax.tick_params(axis='both', which='minor', labelsize=20)
        
    pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal')

    ax.text(-8.0 + 1.5, 8.0 - 1.5, '%.2f' % snap_time + 'Gyr', color='white', size='xx-large')
    
    filename_fig = filename + '.eps' 
    fig.savefig(filename_fig)
    print 'saved image file to :\n\t\t\t %s' % filename_fig
    
    pylab.show()
#

for snap in snaps:    
    plot_maps(snap, params)

pylab.show()