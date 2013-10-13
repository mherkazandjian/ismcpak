import time
import sys

import matplotlib
import mylib
import mylib.utils
matplotlib.use('Qt4Agg')

import numpy
from numpy import log10
import pylab
from mpl_toolkits.mplot3d import Axes3D

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants, nbody_system
from mylib.utils.misc  import xselect
from mylib.utils.histogram import hist_nd
import fi_utils
#===========================================================================================================
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',   # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',   # the path of the dir containing the simulation
          'imres' : 20,                                                  # resolution of the maps to be produced imres x imres
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                            },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-1, 1] | units.kpc,
                      },          
          }

snap_index   = 4
vMin, vMax   = [0, 5]
plot_map     = True
plot_n_gm_g0 = True
#===========================================================================================================
 
#extracting/guessing the metallicity from the name of the directory of the run
metallicity = None

if '-std' in params['rundir']:
    metallicity = 0.2
if '-sol' in params['rundir']:
    metallicity = 1.0

suffix = '%06d' % snap_index
#path to fi snapshots
snapName = 'fiout.%s' % suffix 
filename = params['rundir'] + '/firun/' + snapName 
    
#########################setting up the snapshot data######################
#########################setting up the snapshot data######################
#########################setting up the snapshot data######################
#########################setting up the snapshot data######################
#loading the sph simulation data
print 'loading snapshot %s : ' % filename
gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

#computing some quantities of the sph particles and converting some other
#quantites to units compatible with the SPH particles
gas = fi_utils.convert_units_to_pdr_units(gas, metallicity)

print 'done reading fi snapshot : %s' % filename
print 'number of sph particles in snapshot = %d' %  len(gas)

nBins  = params['imres']
bs_min = params['ranges']['box_size'][0].number
bs_max = params['ranges']['box_size'][1].number

'''
inds = numpy.where(
                   (numpy.log10(n_gas_cgs) >= params['ranges']['sph']['min_log_n_use'])*
                   (numpy.log10(G0) >= params['ranges']['sph']['min_log_G0_use'])*
                   (numpy.log10(g_mech) >= params['ranges']['sph']['min_log_gm_use'])                   
                  )[0]

x1, y1, z1 = xkpc[inds], ykpc[inds], zkpc[inds]
n1, g01, g_mech1 = n_gas_cgs[inds], G0[inds], g_mech[inds]
Av1 = Av_mean[inds]
m1 = m_cgs[inds]
'''

hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn=bs_min, mx=bs_max, nbins=nBins, reverse_indicies=True) 
gas  = gas[hist.inds_in]

#q = gas.temperat.value_in(units.K)
#q = gas.T.number
#q = (gas.rho / (1.3|units.amu)).value_in(units.cm **-3)
#q = 6.54*gas.fuvheat.value_in( units.none )
q_str = 'gas.G0'
#q = gas.G0
q = eval(q_str) #gas.G0


map_f = numpy.zeros(hist.f.shape, 'f8')

for i in numpy.arange(nBins):
    for j in numpy.arange(nBins):
        inds_in_bin = hist.get_indicies([i,j])

        n_in_bin = inds_in_bin.size
        
        if n_in_bin > 0:
            map_f[i,j] = numpy.mean(q[inds_in_bin])

print 'map min,max = ', map_f.min(), map_f.max()

if plot_map:
    
    #clipping the map values outside the specified colorbar ranges
    map_f = log10(map_f)
    print 'map in log scale before clipping = ', map_f.min(), map_f.max()
    
    map_f = numpy.clip(map_f, vMin, vMax)
    
    im = pylab.imshow(
                      map_f,
                      extent=[bs_min, bs_max, bs_min, bs_max],
                      vmin=vMin, 
                      vmax=vMax, 
                      interpolation='bessel', #intepolation used for imshow
                      origin='lower'
                     )   
    pylab.colorbar(im, orientation='vertical')
        
    pylab.xlim([bs_min, bs_max])
    pylab.ylim([bs_min, bs_max])
    pylab.title(q_str)
    pylab.show()

if plot_n_gm_g0:
    
    #plotting the particles states
    x, y, z = log10(gas.n), log10(gas.G0), log10(gas.gmech)
    
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, '.', alpha=0.1)
    ax.set_xlim(-3, 6)
    ax.set_ylim(-3, 6)
    ax.set_zlim(-35, -18)
    
    ax.set_xlabel('n')
    ax.set_ylabel('G0')
    ax.set_zlabel('g_mech')
    
    z_rng    = ax.get_zlim()
    n_z_secs = 15
    z_secs   = numpy.linspace(z_rng[0], z_rng[1], n_z_secs)

    #plotting the n,g0 distributions for different sections in gmech
    n_rows = numpy.int32(numpy.sqrt(n_z_secs)) + 1
    
    fig, axes = pylab.subplots(n_rows, n_rows, sharex=False, sharey=False, figsize=(10, 10))
    pylab.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.01, hspace=0.01)
    
    #setting the plotting ranges and labels
    for ax in axes.flatten():
        ax.set_xlim([-3, 6])
        ax.set_ylim([-3, 6])
    for ax in axes[:,0]: ax.set_ylabel('G0')
    for ax in axes[-1,:]: ax.set_xlabel('n')
    
    #plotting the distribution functions in terms of n,g0 for sections in gmech
    for s in numpy.arange(n_z_secs-1):
        
        j = int(s) / int(n_rows)
        i = int(s) % n_rows 
        ax = axes[j,i]
         
        print s, z_secs[s], z_secs[s+1]

        z_sec_min, z_sec_max = z_secs[s], z_secs[s+1]

        #keeping the elements withing this section of z
        inds = numpy.where(  (z >= z_sec_min) * (z <= z_sec_max) )[0]

        if inds.size > 0:
                    
            x_sec, y_sec, z_sec = x[inds], y[inds], z[inds]
            #pylab.plot( x_sec[::10], y_sec[::10], '.', alpha=0.1, markersize=3)

            setattr(ax,'inds' , inds)
            setattr(ax,'gmech', (10.0**z_sec).mean())
            
            data = numpy.vstack((x_sec, y_sec))
            hist = mylib.utils.histogram.hist_nd(data, bs = 0.25, mn = -3.0, mx = 6.0)
            
            ax.contour(hist.f.T, origin='lower', extent=[-3, 6, -3, 6])
        
            ax.text(ax.get_xlim()[0] + 1, ax.get_ylim()[1] - 2, 
                    r'<$\Gamma_{mech}$> = %.2e' % (10.0**z_sec).mean(),
                    size='small')
            ax.text(ax.get_xlim()[0] + 1, ax.get_ylim()[1] - 3, 
                      'N = %.2e' % x_sec.size, size='small')
        else:
            
            ax.text(0, 0, 'no particles', size='small')

    mylib.utils.removeAxesLabels.removeSubplotLabels(axes, keep_first_last_labels=False)
    
    pylab.show()



f = pylab.figure()

j, i = 2, 2
pylab.plot(gas.x[axes[j,i].inds], gas.y[axes[j,i].inds], 'r.')

j, i = 2, 3
pylab.plot(gas.x[axes[j,i].inds], gas.y[axes[j,i].inds], 'b.')

pylab.xlim([bs_min, bs_max])
pylab.ylim([bs_min, bs_max])
pylab.show()
