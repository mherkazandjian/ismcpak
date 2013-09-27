import time
import sys

import matplotlib
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

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',   # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',   # the path of the dir containing the simulation
          'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                            },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-4, 4] | units.kpc,
                      },          
          }

snap_index   = 20
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

#plotting the temperature vs the gas density
meanmwt=1.3|units.amu
n_gas_cgs = (gas.rho/meanmwt).value_in( units.cm **-3 )
gasT = gas.temperat.value_in(units.K)
G0 = 6.54*gas.fuvheat.value_in( units.none )
g_mech = gas.dethdt.as_quantity_in( units.erg / (units.g * units.s ) )
# converting the mechanical heating rate from per unit mass (of H gas)     
# to per unit volume (see notesISM.odt)
g_mech = g_mech.value_in( g_mech.unit ) # gMech in erg / (g s )
g_mech = g_mech * 1.6474e-24 * n_gas_cgs # gMech in erg / (cm^3 s) 
m_cgs = gas.mass.value_in(units.g)

sigma_v = gas.vdisp.value_in(units.km / units.s)  #sph particle velocity despertion in km/s    
Pe = (1.085 * gasT + 54.0 * sigma_v**2.0)*n_gas_cgs
Av_mean = 0.22 * metallicity * ( 1520.0 / 100.0) * numpy.sqrt(Pe/1e4)

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
#gas  = gas[hist.inds_in]
gas  = gas[hist.inds_in]

#q = gas.temperat.value_in(units.K)
#q = gas.T.number
#q = (gas.rho / (1.3|units.amu)).value_in(units.cm **-3)
#q = 6.54*gas.fuvheat.value_in( units.none )
q = gas.G0

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
    
    im = pylab.imshow(map_f, 
                      extent=[bs_min, bs_max, bs_min, bs_max],
                      vmin=vMin, 
                      vmax=vMax, 
                      interpolation='bessel', #intepolation used for imshow
                      origin='lower')   
    pylab.colorbar(im, orientation='vertical')
        
    pylab.xlim([bs_min, bs_max])
    pylab.ylim([bs_min, bs_max])
        
    pylab.show()


if plot_n_gm_g0:
    
    x, y, z = log10(n_gas_cgs), log10(G0), log10(g_mech)
    
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, '.', alpha=0.1)
    ax.set_xlim(-3, 6)
    ax.set_ylim(-3, 6)
    ax.set_zlim(-35, -15)
    
    ax.set_xlabel('n')
    ax.set_ylabel('G0')
    ax.set_zlabel('g_mech')
    
    z_rng    = ax.get_zlim()
    n_z_secs = 20
    z_secs   = numpy.linspace(z_rng[0], z_rng[1], n_z_secs)
        
    pylab.show()