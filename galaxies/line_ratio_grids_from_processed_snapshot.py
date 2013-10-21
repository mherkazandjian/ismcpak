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
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
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
          'em_unit'          : 'em_fluxKkms',
          'line_ratio_grids' : [
                                ['CO2-1/CO1-0'     , 'CO3-2/CO1-0'    ,      'CO7-6/CO1-0', 'CO10-9/CO1-0'],
                                ['13CO3-2/13CO1-0' , 'CO3-2/CO2-1'    ,      'CO7-6/CO2-1', 'CO10-9/CO2-1'],
                                ['13CO7-6/13CO1-0' , '13CO7-6/13CO2-1',      'CO7-6/CO3-2', 'CO10-9/CO3-2'],
                                ['13CO10-9/13CO1-0', '13CO7-6/13CO2-1', '13CO10-9/13CO2-1', 'CO10-9/CO7-6'],
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
n_rows = len(params['line_ratio_grids'])
n_col  = len(params['line_ratio_grids'][0])

maps_grid = numpy.zeros((n_rows, n_col), '3S20')

for row_ind, row in enumerate(params['line_ratio_grids']):
    
    for nc, ratio_map_info in enumerate(row):

        col_ind = (n_rows - len(row)) + nc
            
        #print row_ind, col_ind, ratio_map_info

        line1, line2 = ratio_map_info.split('/') 
        
        attr1 = params['em_unit'] + '_' + line1
        attr2 = params['em_unit'] + '_' + line2

        print attr1, attr2
        print '-----------------------------------'
        
        maps_grid[row_ind, col_ind][0] = ratio_map_info   
        maps_grid[row_ind, col_ind][1] = attr1   
        maps_grid[row_ind, col_ind][2] = attr2

#getting the unique transitions whose maps will be constructed
map_attrs = numpy.unique(numpy.hstack(maps_grid[:,:,1:3]))
print 'emissions to be extracted from the processed snapshopt'
for attr in map_attrs: print '\t%s' % attr

#getting the species involved in those emissions
species = {}
for key in map_attrs:
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

#making the 2D histogram
bs_min, bs_max = params['ranges']['box_size'].number

print 'getting the spatial distrubutions....'
hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
hist.info()
print '\t\tdone getting the spatial distributuions'

maps = {}
#computing the emission maps of the lines involved in the ratios
print 'making the maps of all the lines...'
for i, this_attr in enumerate(map_attrs):
    maps[this_attr] = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.sum)
print '\t\tfinished making the maps'

################################################PLOTTING########################################################
#plotting the line ratio maps

fig, axs = pylab.subplots(n_rows, n_col, sharex=True, sharey=True, 
                          figsize=(12.0*numpy.float(n_col)/numpy.float(n_rows), 12), 
                          subplot_kw = {'xlim':[bs_min, bs_max],
                                        'ylim':[bs_min, bs_max],
                                        'aspect':'equal',
                                        'adjustable':'datalim',
                                        })
                                 
for ax in axs[:,0] : ax.set_ylabel('y(kpc)')
for ax in axs[-1,:]: ax.set_xlabel('x(kpc)')

v_min, v_max = params['plot']['v_rng']

for r in numpy.arange(n_rows):
    for c in numpy.arange(n_col):
        
        print 'making the line ratio map for'
        print '\t\t', maps_grid[r,c][0]
        
        ratio_map = maps[maps_grid[r,c][1]]/maps[maps_grid[r,c][2]]
        
        ratio_map = numpy.log10(ratio_map)
        
        #ratio_map = ratio_map.clip(-3,3)
        
        im = axs[r,c].imshow(ratio_map.T,   
                             extent=[bs_min, bs_max, bs_min, bs_max],
                             vmin=v_min,     
                             vmax=v_max, 
                             interpolation='bessel', #intepolation used for imshow
                             origin='lower')
    
        axs[r, c].set_title(maps_grid[r,c][0])
            
cbar_ax = fig.add_axes([0.2, 0.95, 0.6, 0.01]) 
cbar_ax.tick_params(axis='both', which='major', labelsize=30)
cbar_ax.set_title(r'$\log_{10}$ [line_ratio]', size=25)
    
pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal', ticks=numpy.linspace(v_min, v_max, 5))

pylab.show()
