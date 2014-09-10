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
import fi_utils
import lineDict
import pickle

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext',  # the path of the dir containing the simulation
          
          'imres'   : 100,   # resolution of the image (over which the beams will be ovelayed)
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
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
                         'CO1-0'  , 'CO2-1'  , 'CO3-2'  , 'CO4-3'  , 'CO5-4'  , 'CO6-5'  ,
                         '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', '13CO5-4', '13CO6-5',
                        ],
          'beam_radii': numpy.linspace(0.2, 8.0, 60),
          'save' : False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

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
              'lines' : numpy.array(params['lines'], 'S'),
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
    this_map_intensity = fi_utils.make_map(gas, hist, attr=this_attr, func=numpy.mean) #func=numpy.average, weights=this_attr)
    
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

################################## plotting the emission as a function of beam radius ############################
fig1 = pylab.figure()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])

fig2 = pylab.figure()
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

names1, names2 = [], []
plts1, plts2 = [], []
area_secotrs = numpy.pi*(hist.f.epos**2 - hist.f.spos**2)

for line in luminosity['lines']:
    
    if lineDict.lines[line]['specStr'] == 'CO':
        plt1, = ax1.plot(luminosity['beam_r'], numpy.log10(luminosity['lum_r'][line]))
        names1.append(line)
        plts1.append(plt1)
    if lineDict.lines[line]['specStr'] == '13CO':
        plt2, = ax2.plot(luminosity['beam_r'], numpy.log10(luminosity['lum_r'][line]))
        names2.append(line)
        plts2.append(plt2)
ax1.legend(plts1, names1, loc=0)
ax2.legend(plts2, names2, loc=0)

ax1.set_xlabel('beam radius')
ax1.set_ylabel('log10 [line luminosity / (K.km/s kpc2)]')

ax2.set_xlabel('beam radius')
ax2.set_ylabel('log10 [line luminosity / (K.km/s kpc2)]')
pylab.show()

if params['save'] == True:
    
    fname = os.path.join(params['rundir'],'analysis','line_luminosity_CO_13CO_snap_%d.pkl' % params['snap_index'][0])
    
    fObj = open(fname, 'w')
    
    pickle.dump(luminosity,fObj)
    
    fObj.close()
    
    print 'wrote the luminosity dict to :\n\t\t %s' % fname
