'''
plot a grid of line ratios from a snapshot
'''
#########################################################################################################
import time, sys, os

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units
from mylib.utils.misc  import default_logger
from mylib.utils.histogram import hist_nd 
import fi_utils
import lineDict
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'r_res'   : 25, # resolution in the radia scale of the plots as a function of r
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snap_index': numpy.arange(4, 4 + 1, 1),
          'ranges'    : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                         'sph':{
                               'min_log_n_use'  : -3.0,      
                               'min_log_G0_use' : -3.0,
                               'min_log_gm_use' : -50.0,
                               'Av_use'         :  [0.0, 20000000.0],
                               'Av_clip'        :  [3.0, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                               },
                      
                        #the size of the box to be displayed (particles outside the range are discarded)
                        'box_size' : [-8.0, 8.0] | units.kpc, 
                        },
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

#setting up the line ratio map grid
curves = [[],[],[]]
attrs = {}
for i, bunch in enumerate(params['line_ratio_curves']):
    
    for j, ratio in enumerate(bunch):
        
        line1, line2 = ratio.split('/')

        attr1 = params['em_unit'] + '_' + line1
        attr2 = params['em_unit'] + '_' + line2
        
        curves[i].append([ratio, attr1, attr2])

        print i, j, ratio, line1, line2, attr1, attr2
        
        attrs[attr1] = True
        attrs[attr2] = True

#getting the unique transitions whose maps will be constructed
curve_attrs = attrs.keys()
print 'emissions to be extracted from the processed snapshopt'
for attr in curve_attrs: print '\t%s' % attr

#getting the species involved in those emissions
species = {}
for key in curve_attrs:
    species[lineDict.lines[key.replace(params['em_unit'] + '_', '')]['specStr']] = True
print 'Species invloved = ', species.keys()

#path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

#loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, species, load_pdr=params['pdr_sph'])    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

#keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

#keeping the gas particles within a radius R = bs_max
bs_min, bs_max = params['ranges']['box_size'].number

#making the 1D histogram as a function of R
x, y = gas.x, gas.y
r = numpy.sqrt(x**2 + y**2)

print 'getting the spatial distrubutions....'
hist = hist_nd(r.reshape(1, r.size), mn = 0.0, mx=bs_max, nbins=params['r_res'], reverse_indicies=True, loc=True)
hist.info()

#keeping the gas particles which are withing the ranges of the histogram
gas = gas[hist.inds_in]

print '\t\tdone getting the spatial distributuions'

fig = pylab.figure()
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

x, y = gas.x, gas.y 


#plotting the points
for i in numpy.arange(hist.nBins):
    inds_in_bin = hist.get_indicies(i)
    ax.plot(x[inds_in_bin], y[inds_in_bin], '.')

pylab.show()
 
################################################computing the emission###########################################   
#computing the emission maps of the lines involved in the ratios
print 'computing the emission and the mean excitation temperatures in all the sectors...'

#curve which will hold the computed data the sectors 
emission = numpy.zeros(hist.f.size, 'f8')
Tex_mean = numpy.zeros(hist.f.size, 'f8')
gmech_per_n_mean = numpy.zeros(hist.f.size, 'f8')

#emission data for this line for all the gas particles
data_em  = getattr(gas, 'em_fluxKkms_CO1-0')
data_Tex = getattr(gas, 'em_Tex_CO1-0')
data_gm  = getattr(gas, 'gmech')
data_n   = getattr(gas, 'n')

#looping over the sectors and getting the emissions for this sector    
for i in numpy.arange(hist.nBins):

    inds_in_bin = hist.get_indicies(i)
    
    emission[i] = numpy.sum(data_em[inds_in_bin])
    Tex_mean[i] = numpy.average(data_Tex[inds_in_bin], weights=data_em[inds_in_bin]) 
    gmech_per_n_mean[i] = numpy.average(data_gm[inds_in_bin]/data_n[inds_in_bin], weights=data_em[inds_in_bin]) 
    
    print '\t\t... done'
        
print '\tfinished computing the emissions'

#plotting the emission as a function of radius
fig1 = pylab.figure()
ax1 = fig1.add_axes([0.15, 0.15, 0.8, 0.8])

fig2 = pylab.figure()
ax2 = fig2.add_axes([0.15, 0.15, 0.8, 0.8])

ax1.plot(hist.f.cntrd, Tex_mean)
ax1.set_xlabel(r'$r$ (kpc)', size=20)
ax1.set_ylabel(r'<T$_{ex}$>', size=20)

ax2.plot(hist.f.cntrd, numpy.log10(gmech_per_n_mean))
#ax1.set_xlabel(r'$r$ (kpc)', size=20)
#ax1.set_ylabel(r'<T$_{ex}$>', size=20)

pylab.show()

