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
          
          'r_res'   : 50, # resolution in the radia scale of the plots as a function of r
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
          'em_unit'           : 'em_fluxKkms',
          'line_ratio_curves' : [
                                ['CO2-1/CO1-0', 'CO3-2/CO1-0' , 'CO7-6/CO1-0' , 'CO10-9/CO1-0',
                                 'CO3-2/CO2-1', 'CO7-6/CO2-1' , 'CO10-9/CO2-1',
                                 'CO7-6/CO3-2', 'CO10-9/CO3-2',
                                 'CO10-9/CO7-6',
                                ],
                                ['13CO2-1/13CO1-0', '13CO3-2/13CO1-0' , '13CO7-6/13CO1-0' , '13CO10-9/13CO1-0',
                                 '13CO3-2/13CO2-1', '13CO7-6/13CO2-1' , '13CO10-9/13CO2-1',
                                 '13CO7-6/13CO3-2', '13CO10-9/13CO3-2',
                                 '13CO10-9/13CO7-6',
                                ],
                                ['13CO2-1/CO1-0' , '13CO3-2/CO1-0' , '13CO7-6/CO1-0' , '13CO10-9/CO1-0',
                                 '13CO3-2/CO2-1' , '13CO7-6/CO2-1' , '13CO10-9/CO2-1',
                                 '13CO7-6/CO3-2' , '13CO10-9/CO3-2',
                                 '13CO10-9/CO7-6',
                                ],
                               ],
          'plot' : {
                    'v_rng' : [-3.0, 1.0],
                    }
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
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

x, y = gas.x, gas.y

#plotting the points
for i in numpy.arange(hist.nBins):
    inds_in_bin = hist.get_indicies(i)
    ax.plot(x[inds_in_bin], y[inds_in_bin], '.')

pylab.show()
 
################################################computing the emission###########################################   
emissions_data = {}
#computing the emission maps of the lines involved in the ratios
print 'computing the emission from all the sectors...'

area_secotrs = numpy.pi*(hist.f.epos**2 - hist.f.spos**2)

for i, this_attr in enumerate(curve_attrs):
    
    print 'computing the emission as a function of radius for %s' % this_attr
    
    #curve which will hold the computed emission intensity in that sector
    emissions_data[this_attr] = numpy.zeros(hist.f.size, 'f8')

    #emission data for this line for all the gas particles
    data_attr = getattr(gas, this_attr)

    #looping over the sectors and getting the emissions for this sector    
    for i in numpy.arange(hist.nBins):

        inds_in_bin = hist.get_indicies(i)
        
        emissions_data[this_attr][i] = (data_attr[inds_in_bin].sum())
    
    print '\t\t... done'
        
print '\tfinished computing the emissions'

#plotting the emission as a function of radius
fig = pylab.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

plts = []
for i, this_attr in enumerate(curve_attrs):
    plt, = ax.plot(hist.f.cntrd, emissions_data[this_attr]/area_secotrs)
    plts.append(plt)
ax.legend(plts, curve_attrs, loc=0)

pylab.show()


############################################Plotting the line ratios versus r###########################################
   
v_min, v_max = params['plot']['v_rng']

fig, axs = pylab.subplots(1, 3, sharex=False, sharey=False, figsize=(12.0, 4.0))
                                
#for ax in axs[:,0] : ax.set_ylabel('y(kpc)')
for ax in axs: ax.set_xlabel('R(kpc)', size=10)
axs[0].set_ylabel(r'$\log_{10}$ [line ratio] ')

pylab.subplots_adjust(left=0.10, bottom=0.15, right=0.9, top=0.9, wspace=0.15, hspace=0.15)
titles = ['CO/CO', '13CO/13CO', '13CO/CO']
syms = ['r-', 'g-', 'r--', 'g--', 'k-', 'b--', 'k--', 'r-.', 'g-.', 'b-.', 'k-.']

for i, bunch in enumerate(curves):

    for j, curve in enumerate(bunch):
        print i, j, curve 
        curve_data = numpy.log10(emissions_data[curve[1]] / emissions_data[curve[2]])
        axs[i].plot(hist.f.cntrd, curve_data, syms[j])
    
    print '-----------------' 
    axs[i].set_ylim(-5, 0)
    axs[i].set_title(titles[i])
pylab.show()
