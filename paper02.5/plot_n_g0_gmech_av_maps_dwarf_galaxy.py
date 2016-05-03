import time
import sys

import matplotlib
from IPython.core.debugger import Tracer

matplotlib.use('Qt4Agg')

import numpy
from numpy import log10
import pylab
from mpl_toolkits.mplot3d import Axes3D

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants, nbody_system
from mylib.utils.misc import xselect
from mylib.utils.histogram import hist_nd
sys.path.append('../galaxies')
import fi_utils

#===========================================================================================================
home = '/home/mher'

params = {
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',   # the path of the dir containing the simulation
          'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                            },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8, 8] | units.kpc,
                      },          
          }

snap_index   = 20
# vMin, vMax   = [0, 5]
# vMin, vMax   = [-30, -16]
# vMin, vMax   = [0, 30]
plot_map     = True
plot_n_gm_g0 = False
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
if True:
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
    gas = gas[hist.inds_in]

################################################################################
# setup the figure
fig, axs = pylab.subplots(1, 4, sharex=True, sharey=True, figsize=(7, 3.1),
                          subplot_kw = {'xlim':[-2, 2],  #keywords passed to figure.add_subplot()
                                        'ylim':[-2, 2],
                                        'aspect':'auto',
                                        'adjustable':'datalim',
                                        'autoscale_on' : False,   # enforces the set limits (ranges) of the subplot axes
                                        'xticks': pylab.linspace(-1, 1, 3),
                                        'yticks': pylab.linspace(-1, 1, 3),
                                        }
                          )

pylab.subplots_adjust(left=0.1, bottom=0.18,
                      right=0.98, top=0.65,
                      wspace=0.05, hspace=0.0)

axs[0].set_ylabel('y(kpc)', size=8)
for ax in axs:
    ax.set_xlabel('x(kpc)', size=8)
    ax.tick_params(axis='both', labelsize=8)

fontsz = 9
cbar1 = fig.add_axes([0.1, 0.75, 0.21, 0.03])
cbar1.tick_params(axis='both', which='major', labelsize=fontsz)
cbar1.set_title(r'$\log_{10}$<n/${\rm cm}^{-3}$>', size=fontsz)
cbar1.set_yticks([])
cbar1.set_xlim(-1, 4)
cbar1_ticks = numpy.linspace(0, 3, 4)
cbar1.set_xticks(cbar1_ticks)
cbar1.set_xticklabels(cbar1_ticks, size=10)

cbar2 = fig.add_axes([0.323, 0.75, 0.21, 0.03])
cbar2.tick_params(axis='both', which='major', labelsize=fontsz)
cbar2.set_title(r'$\log_{10}$<$G/G_0$>', size=fontsz)
cbar2.set_yticks([])
cbar2.set_xlim(-1, 4)
cbar2_ticks = numpy.linspace(0, 3, 4)
cbar2.set_xticks(cbar2_ticks)
cbar2.set_xticklabels(cbar2_ticks, size=10)

cbar3 = fig.add_axes([0.545, 0.75, 0.21, 0.03])
cbar3.tick_params(axis='both', which='major', labelsize=fontsz)
cbar3.set_title(r'$\log_{10}$<$\Gamma_{\rm mech}$/${\rm erg}$ ${\rm cm}^{-3}$ ${\rm s}^{-1}$>', size=fontsz)
cbar3.set_yticks([])
cbar3.set_xlim(-26, -20)
cbar3_ticks = numpy.linspace(-26, -20, 4)
cbar3.set_xticks(cbar3_ticks)
cbar3.set_xticklabels(cbar3_ticks, size=10)

cbar4 = fig.add_axes([0.77, 0.75, 0.21, 0.03])
cbar4.tick_params(axis='both', which='major', labelsize=fontsz)
cbar4.set_title(r'<A$_{\rm V}$/mag >', size=fontsz)
cbar4.set_yticks([])
cbar4.set_xlim(0, 5)
cbar4_ticks = numpy.linspace(1, 5, 5)
cbar4.set_xticks(cbar4_ticks)
cbar4.set_xticklabels(cbar4_ticks, size=10)

################################################################################
q = gas.n
vMin, vMax = -1, 4
map_f = numpy.zeros(hist.f.shape, 'f8')

for i in numpy.arange(nBins):
    for j in numpy.arange(nBins):
        inds_in_bin = hist.get_indicies([i,j])
        n_in_bin = inds_in_bin.size
        if n_in_bin > 0:
            map_f[i,j] = numpy.mean(q[inds_in_bin])

#clipping the map values outside the specified colorbar ranges
map_f = log10(map_f)
print 'map in log scale before clipping = ', map_f.min(), map_f.max()

map_f = numpy.clip(map_f, vMin, vMax)

im = axs[0].imshow(map_f,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=vMin,
                   vmax=vMax,
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')
pylab.colorbar(im, cax=cbar1, orientation='horizontal', ticks=cbar1_ticks)

################################################################################
q = gas.G0
vMin, vMax = -1, 4
map_f = numpy.zeros(hist.f.shape, 'f8')

for i in numpy.arange(nBins):
    for j in numpy.arange(nBins):
        inds_in_bin = hist.get_indicies([i,j])
        n_in_bin = inds_in_bin.size
        if n_in_bin > 0:
            map_f[i,j] = numpy.mean(q[inds_in_bin])

map_f = log10(map_f)
print 'map in log scale before clipping = ', map_f.min(), map_f.max()

map_f = numpy.clip(map_f, vMin, vMax)

im = axs[1].imshow(map_f,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=vMin,
                   vmax=vMax,
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

pylab.colorbar(im, cax=cbar2, orientation='horizontal', ticks=cbar2_ticks)

################################################################################
q = gas.gmech
vMin, vMax = -26, -20
map_f = numpy.zeros(hist.f.shape, 'f8')

for i in numpy.arange(nBins):
    for j in numpy.arange(nBins):
        inds_in_bin = hist.get_indicies([i,j])
        n_in_bin = inds_in_bin.size
        if n_in_bin > 0:
            map_f[i,j] = numpy.mean(q[inds_in_bin])

#clipping the map values outside the specified colorbar ranges
map_f = log10(map_f)
print 'map in log scale before clipping = ', map_f.min(), map_f.max()

map_f = numpy.clip(map_f, vMin, vMax)

im = axs[2].imshow(map_f,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=vMin,
                   vmax=vMax,
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

pylab.colorbar(im, cax=cbar3, orientation='horizontal', ticks=cbar3_ticks)

################################################################################
q = gas.Av
vMin, vMax = 0, 5
map_f = numpy.zeros(hist.f.shape, 'f8')

for i in numpy.arange(nBins):
    for j in numpy.arange(nBins):
        inds_in_bin = hist.get_indicies([i,j])
        n_in_bin = inds_in_bin.size
        if n_in_bin > 0:
            map_f[i,j] = numpy.mean(q[inds_in_bin])

# map_f = log10(map_f)
print 'map in log scale before clipping = ', map_f.min(), map_f.max()

map_f = numpy.clip(map_f, vMin, vMax)

im = axs[3].imshow(map_f,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=vMin,
                   vmax=vMax,
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

pylab.colorbar(im, cax=cbar4, orientation='horizontal', ticks=cbar4_ticks)

pylab.show()

# fig.savefig('/home/mher/ism/docs/paper02.5/src/figs/n_g_gmech_av_maps/dwarf.eps')

