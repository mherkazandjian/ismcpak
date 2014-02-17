'''
    - compute the histogram
    - from the histogram of the projected particles on the sky 
    - from the histogram compute the intensity weighed intensity per bin (pixel) for each line (this is in erg.cm2.s-1) 
    - from this compute the luminosity within each pixel
    - Superpose the beam on each map and sum the luminosity covered by each beam to compute the total 
      lumonisity in that beam.  
    - plot the line ratios for different line as a function of beam size (to get an idea how they depend on r)
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
from galaxies import fi_utils
import lineDict
import pickle

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

fig_save_path = '/home/mher/ism/docs/paper02.5/src/figs/total_luminosity_vs_J.eps'
#fig_save_path = None #'/home/mher/ism/docs/paper02.5/src/figs/total_luminosity_vs_J.eps'

#########################################################################################################
######################DISK GALAXY DISK GALAXY DISK GALAXY DISK GALAXY DISK GALAXY #######################
#########################################################################################################

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
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
          'em_unit'   : 'em_fluxKkms',
          'lines'     : [
                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3'  ,# 'CO5-4'  , 'CO6-5'  ,
                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3',# '13CO5-4', '13CO6-5',
                        ],
          'beam_radii': numpy.linspace(0.2, 8.0, 30),
          'fig_save'  : False,
        }

## setting up the logger object
logger = default_logger()

## setting up the line attribute array whose emission will be retried from the snapshot
attrs = {}
for i, line in enumerate(params['lines']):
    
        attr = params['em_unit'] + '_' + line
        
        attrs[attr] = True
#

## getting the unique transitions whose maps will be constructed
curve_attrs = attrs.keys()
print 'emissions to be extracted from the processed snapshopt'
for attr in curve_attrs: print '\t%s' % attr

## getting the species involved in those emissions
species = {}
for key in curve_attrs:
    species[lineDict.lines[key.replace(params['em_unit'] + '_', '')]['specStr']] = True
print 'Species invloved = ', species.keys()

## path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

## loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, species, load_pdr=params['pdr_sph'])    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

## keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

## makingt the 2d histogram
bs_min, bs_max = params['ranges']['box_size'].number

print 'getting the spatial distrubutions....'
hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
hist.info()
print '\t\tdone getting the spatial distributuions'

## keeping the gas particles which are within the ranges of the histogram
gas = gas[hist.inds_in]

## getting the luminosity maps for each line
luminosity = {
              'lines' : numpy.array(params['lines'],'S'),
              'beam_r': params['beam_radii'],
              'lum_r' : {},
              'maps'  : {},
              'units' : {
                         'luminosity' : 'K.km.s^-1.kpc^2',
                         },
             }

print 'computing the luminosity from all the pixles in for each line map...'

print 'making the maps of all the lines...'
for i, this_attr in enumerate(attrs):
    
    ## the intensity map (intensity weight averaged intensity map)
    #this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.average, weights=this_attr)
    this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #ge, weights=this_attr)
    
    ## computing the luminsoty by mutiplying by the area of each pixel
    this_map_luminosity = this_map_intensity * hist.f.dl.prod()
    
    ## getting the line code from the attribute name
    line = this_attr.replace(params['em_unit']+'_','')
    
    luminosity['maps'][line] = this_map_luminosity

    luminosity['lum_r'][line] = numpy.zeros(params['beam_radii'].size,'f8')
    
print '\t\tfinished making the luminosity maps'

## looping over the annuli and getting the luminosity with each annulus as well as in the total luminosity
## from the beam whose radius is the outer radius of by each anulus

## getting the distance each the pixels from the center of the map
r_pix = numpy.sqrt(hist.f.cntrd[0]**2 + hist.f.cntrd[1]**2)

for i, r_beam in enumerate(luminosity['beam_r']):

    inds_pix_in_beam = numpy.where(r_pix < r_beam)
    
    for j, line in enumerate(luminosity['lines']):

        print i, j, r_beam, line
        
        luminosity['lum_r'][line][i] = luminosity['maps'][line][inds_pix_in_beam].sum()
             
print '\tfinished computing the emissions'

beam_radii = luminosity['beam_r']

use_beam_size = 8.0
ind_8kpc = numpy.argmin(numpy.fabs(beam_radii - use_beam_size))
print 'available beam size used : %e' % beam_radii[ind_8kpc]

use_beam_size = 1.0
ind_1kpc = numpy.argmin(numpy.fabs(beam_radii - use_beam_size))
print 'available beam size used : %e' % beam_radii[ind_1kpc]

ladder_disk_CO_8kpc = []  
ladder_disk_13CO_8kpc = [] 
ladder_disk_CO_1kpc = []  
ladder_disk_13CO_1kpc = [] 

for i in range(4):
    CO_line_str = 'CO%d-%d' % (i+1,i)
    CO13_line_str = '13CO%d-%d' % (i+1,i)

    ladder_disk_CO_8kpc.append(luminosity['lum_r'][CO_line_str][ind_8kpc])
    ladder_disk_13CO_8kpc.append(luminosity['lum_r'][CO13_line_str][ind_8kpc])
    ladder_disk_CO_1kpc.append(luminosity['lum_r'][CO_line_str][ind_1kpc])
    ladder_disk_13CO_1kpc.append(luminosity['lum_r'][CO13_line_str][ind_1kpc])

#########################################################################################################
######################DWARF GALAXY DWARF GALAXY DWARF GALAXY DWARF GALAXY DWARF GALAXY #######################
#########################################################################################################

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : True, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snap_index': numpy.arange(20, 20 + 1, 1),
          'ranges'    : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                         'sph':{
                               'min_log_n_use'  : -3.0,      
                               'min_log_G0_use' : -3.0,
                               'min_log_gm_use' : -50.0,
                               'Av_use'         :  [0.0, 20000000.0],
                               'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                               },
                      
                        #the size of the box to be displayed (particles outside the range are discarded)
                        'box_size' : [-2.0, 2.0] | units.kpc, 
                        },
          'em_unit'   : 'em_fluxKkms',
          'lines'     : [
                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3'  , #'CO5-4'  , 'CO6-5'  ,
                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', #'13CO5-4', '13CO6-5',
                        ],
          'beam_radii': numpy.linspace(0.2, 2.0, 30),
          'fig_save'  : False,
        }


## setting up the logger object
logger = default_logger()

## setting up the line attribute array whose emission will be retried from the snapshot
attrs = {}
for i, line in enumerate(params['lines']):
    
        attr = params['em_unit'] + '_' + line
        
        attrs[attr] = True
#

## getting the unique transitions whose maps will be constructed
curve_attrs = attrs.keys()
print 'emissions to be extracted from the processed snapshopt'
for attr in curve_attrs: print '\t%s' % attr

## getting the species involved in those emissions
species = {}
for key in curve_attrs:
    species[lineDict.lines[key.replace(params['em_unit'] + '_', '')]['specStr']] = True
print 'Species invloved = ', species.keys()

## path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % params['snap_index'] + '.states.npz'  

## loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, species, load_pdr=params['pdr_sph'])    
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

## keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

## makingt the 2d histogram
bs_min, bs_max = params['ranges']['box_size'].number

print 'getting the spatial distrubutions....'
hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
hist.info()
print '\t\tdone getting the spatial distributuions'

## keeping the gas particles which are within the ranges of the histogram
gas = gas[hist.inds_in]

## getting the luminosity maps for each line
luminosity = {
              'lines' : numpy.array(params['lines'],'S'),
              'beam_r': params['beam_radii'],
              'lum_r' : {},
              'maps'  : {},
              'units' : {
                         'luminosity' : 'K.km.s^-1.kpc^2',
                         },
             }

print 'computing the luminosity from all the pixles in for each line map...'

print 'making the maps of all the lines...'
for i, this_attr in enumerate(attrs):
    
    ## the intensity map (intensity weight averaged intensity map)
    this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #average, weights=this_attr)
    
    ## computing the luminsoty by mutiplying by the area of each pixel
    this_map_luminosity = this_map_intensity * hist.f.dl.prod()
    
    ## getting the line code from the attribute name
    line = this_attr.replace(params['em_unit']+'_','')
    
    luminosity['maps'][line] = this_map_luminosity

    luminosity['lum_r'][line] = numpy.zeros(params['beam_radii'].size,'f8')
    
print '\t\tfinished making the luminosity maps'

## looping over the annuli and getting the luminosity with each annulus as well as in the total luminosity
## from the beam whose radius is the outer radius of by each anulus

## getting the distance each the pixels from the center of the map
r_pix = numpy.sqrt(hist.f.cntrd[0]**2 + hist.f.cntrd[1]**2)

for i, r_beam in enumerate(luminosity['beam_r']):

    inds_pix_in_beam = numpy.where(r_pix < r_beam)
    
    for j, line in enumerate(luminosity['lines']):

        print i, j, r_beam, line
        
        luminosity['lum_r'][line][i] = luminosity['maps'][line][inds_pix_in_beam].sum()
             
print '\tfinished computing the emissions'

beam_radii = luminosity['beam_r']

use_beam_size = 2.0
ind_2kpc = numpy.argmin(numpy.fabs(beam_radii - use_beam_size))
print 'available beam size used : %e' % beam_radii[ind_2kpc]

ladder_dwarf_CO_2kpc = []  
ladder_dwarf_13CO_2kpc = [] 

for i in range(4):
    CO_line_str = 'CO%d-%d' % (i+1,i)
    CO13_line_str = '13CO%d-%d' % (i+1,i)

    ladder_dwarf_CO_2kpc.append(luminosity['lum_r'][CO_line_str][ind_2kpc])
    ladder_dwarf_13CO_2kpc.append(luminosity['lum_r'][CO13_line_str][ind_2kpc])
######################################################################################################
############################plotting the image to be saved #################################

fig = pylab.figure(figsize=(4,4))
ax = fig.add_axes([0.18, 0.15, 0.8, 0.8])

Ju_all = [1,2,3,4]#,5,6]

ax.set_xlim([1, 4])

ladder_disk_CO_8kpc = numpy.array(ladder_disk_CO_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2 
pltCO, = ax.semilogy(Ju_all, ladder_disk_CO_8kpc, 'k')
ax.text(Ju_all[2], ladder_disk_CO_8kpc[2], 'disk')

ladder_disk_13CO_8kpc = numpy.array(ladder_disk_13CO_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
plt13CO, = ax.semilogy(Ju_all, ladder_disk_13CO_8kpc, 'k--')
ax.text(Ju_all[2], ladder_disk_13CO_8kpc[2], 'disk')

ladder_disk_CO_1kpc = numpy.array(ladder_disk_CO_1kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
ax.semilogy(Ju_all, ladder_disk_CO_1kpc, 'r-')
ax.text(Ju_all[1], ladder_disk_CO_1kpc[1], 'center')

ladder_disk_13CO_1kpc = numpy.array(ladder_disk_13CO_1kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
ax.semilogy(Ju_all, ladder_disk_13CO_1kpc, 'r--')
ax.text(Ju_all[1], ladder_disk_13CO_1kpc[1], 'center')

ladder_dwarf_CO_2kpc = numpy.array(ladder_dwarf_CO_2kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
ladder_dwarf_13CO_2kpc = numpy.array(ladder_dwarf_13CO_2kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
ax.semilogy(Ju_all, ladder_dwarf_CO_2kpc, 'b-')
ax.semilogy(Ju_all, ladder_dwarf_13CO_2kpc, 'b--')

ax.text(Ju_all[1], ladder_dwarf_CO_2kpc[1], 'dwarf')
ax.text(Ju_all[1], ladder_dwarf_13CO_2kpc[1], 'dwarf')

ax.legend([pltCO, plt13CO], ['CO', r'$^{13}$CO'], loc=0)

ax.set_xlabel(r'J$_{\rm up}$', size=10)
ax.set_ylabel(r'luminosity [K.km.s$^{-1}$.pc$^2$]', size=10)

ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xticks(Ju_all)
ax.set_xticklabels(Ju_all, '%d')

pylab.draw()
pylab.show()
###############################################################################################################

if fig_save_path != None:
    fig.savefig(fig_save_path)
    print 'saved image file to :\n\t\t\t %s' % fig_save_path
