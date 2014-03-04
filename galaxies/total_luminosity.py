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

#fig_save_path = '/home/mher/ism/docs/paper02.5/src/figs/total_luminosity_vs_J.eps'
fig_save_path = None #'/home/mher/ism/docs/paper02.5/src/figs/total_luminosity_vs_J.eps'

#########################################################################################################
######################DISK GALAXY DISK GALAXY DISK GALAXY DISK GALAXY DISK GALAXY #######################
#########################################################################################################

params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext',  # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
          'weights' : 'original-only', #'by-number', #'matched',  #'original-only' ,#None 
           
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
#                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3'  ,'CO5-4'  , 'CO6-5'  ,
                          'CO2-1'  , 'CO3-2'  , 'CO4-3'  ,'CO5-4'  , 'CO6-5', 'CO7-6', 'CO8-7', 'CO9-8', 'CO10-9', 'CO11-10','CO12-11', 'CO13-12', 'CO14-13', 'CO15-14',
                          '13CO2-1'  , '13CO3-2'  , '13CO4-3'  ,'13CO5-4'  , '13CO6-5', '13CO7-6', '13CO8-7', '13CO9-8', '13CO10-9', '13CO11-10','13CO12-11', '13CO13-12', '13CO14-13', '13CO15-14',
#                          'HCN2-1'  , 'HCN3-2'  , 'HCN4-3'  ,'HCN5-4'  , 'HCN6-5', 'HCN7-6',
#                          'HNC2-1'  , 'HNC3-2'  , 'HNC4-3'  ,'HNC5-4'  , 'HNC6-5', 'HNC7-6',
#                          'HCO+2-1'  , 'HCO+3-2'  , 'HCO+4-3'  ,'HCO+5-4'  , 'HCO+6-5', 'HCO+7-6',

#                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', '13CO5-4', '13CO6-5',
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

## setting the weights
weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % params['snap_index']
gas.use_weights(weighting=params['weights'], weights_filename = weights_filename)

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
    
    print '\t%s' % this_attr
    
    ## the intensity map (intensity weight averaged intensity map)
    this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.average, weights='weights')
    #this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #ge, weights=this_attr)
    
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

#use_beam_size = 1.0
#ind_1kpc = numpy.argmin(numpy.fabs(beam_radii - use_beam_size))
#print 'available beam size used : %e' % beam_radii[ind_1kpc]

ladder_disk_CO_8kpc = []
ladder_disk_13CO_8kpc = [] 
ladder_disk_CO_1kpc = []
ladder_disk_13CO_1kpc = [] 

#### CO and 13CO
Ju_all = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]

for i in numpy.array(Ju_all) - 1:
    CO_line_str = 'CO%d-%d' % (i+1,i)
    CO13_line_str = '13CO%d-%d' % (i+1,i)

    ladder_disk_CO_8kpc.append(luminosity['lum_r'][CO_line_str][ind_8kpc])
    ladder_disk_13CO_8kpc.append(luminosity['lum_r'][CO13_line_str][ind_8kpc])
#    ladder_disk_CO_1kpc.append(luminosity['lum_r'][CO_line_str][ind_1kpc])
#    ladder_disk_13CO_1kpc.append(luminosity['lum_r'][CO13_line_str][ind_1kpc])

#### HCN, HNC, HCO+
Ju_all_HD = [2,3,4,5,6,7]

ladder_disk_HCN_8kpc = [] 
ladder_disk_HNC_8kpc = [] 
#ladder_disk_HCOP_8kpc = []

for i in numpy.array(Ju_all_HD) - 1:
    
    HCN_line_str  = 'HCN%d-%d' % (i+1,i)
    HNC_line_str  = 'HNC%d-%d' % (i+1,i)
#    HCOP_line_str = 'HCO+%d-%d' % (i+1,i)

#    ladder_disk_HCN_8kpc.append(luminosity['lum_r'][HCN_line_str][ind_8kpc])
#    ladder_disk_HNC_8kpc.append(luminosity['lum_r'][HNC_line_str][ind_8kpc])
#    ladder_disk_HCOP_8kpc.append(luminosity['lum_r'][HCOP_line_str][ind_8kpc])


fig = pylab.figure(figsize=(4,4))
ax = fig.add_axes([0.18, 0.15, 0.8, 0.8])

#Ju_all = [1,2,3,4,5,6]
Ju_all = [2,3,4,5,6,7,8,9,10,11,12,13,14,15]

ladder_disk_CO_8kpc = numpy.array(ladder_disk_CO_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2 
pltCO, = ax.semilogy(Ju_all, ladder_disk_CO_8kpc, 'k')
ax.text(Ju_all[2], ladder_disk_CO_8kpc[2], 'disk-CO')

ladder_disk_13CO_8kpc = numpy.array(ladder_disk_13CO_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
plt13CO, = ax.semilogy(Ju_all, ladder_disk_13CO_8kpc, 'k--')
ax.text(Ju_all[2], ladder_disk_13CO_8kpc[2], 'disk-13CO')

#ladder_disk_HCN_8kpc = numpy.array(ladder_disk_HCN_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
#pltHCN, = ax.semilogy(Ju_all_HD, ladder_disk_HCN_8kpc, 'k--')
#ax.text(Ju_all_HD[2], ladder_disk_HCN_8kpc[2], 'disk-HCN')

#ladder_disk_HNC_8kpc = numpy.array(ladder_disk_HNC_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
#pltHNC, = ax.semilogy(Ju_all_HD, ladder_disk_HNC_8kpc, 'k--')
#ax.text(Ju_all_HD[2], ladder_disk_HNC_8kpc[2], 'disk-HNC')

#ladder_disk_HCOP_8kpc = numpy.array(ladder_disk_HCOP_8kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
#pltHCOP, = ax.semilogy(Ju_all_HD, ladder_disk_HCOP_8kpc, 'k--')
#ax.text(Ju_all_HD[2], ladder_disk_HCOP_8kpc[2], 'disk-HCO+')

#ladder_disk_CO_1kpc = numpy.array(ladder_disk_CO_1kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
#ax.semilogy(Ju_all, ladder_disk_CO_1kpc, 'r-')
#ax.text(Ju_all[1], ladder_disk_CO_1kpc[1], 'center')

#ladder_disk_13CO_1kpc = numpy.array(ladder_disk_13CO_1kpc)*1e6 ## converting from K.km.s-1.kpc^2 to K.km.s-1.pc^2
#ax.semilogy(Ju_all, ladder_disk_13CO_1kpc, 'r--')
#ax.text(Ju_all[1], ladder_disk_13CO_1kpc[1], 'center')

#ax.legend([pltCO, plt13CO], ['CO', r'$^{13}$CO'], loc=0)

ax.set_xlabel(r'J$_{\rm up}$', size=10)
ax.set_ylabel(r'luminosity [K.km.s$^{-1}$.pc$^2$]', size=10)

ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xticks(Ju_all)
ax.set_xticklabels(Ju_all, '%d')
ax.set_xlim([1, 16])
ax.set_ylim([10, 1e10])

pylab.draw()
pylab.show()
###############################################################################################################

if fig_save_path != None:
    fig.savefig(fig_save_path)
    print 'saved image file to :\n\t\t\t %s' % fig_save_path
