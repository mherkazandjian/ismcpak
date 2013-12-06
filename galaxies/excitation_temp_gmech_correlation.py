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
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'species' : ['CO'],#, '13CO'],
          'pdr_sph' : False, #if set to true looks for the file fiout.xxxxxx.states.npz.pdr.npz and tries to load it
           
          'snaps'   : numpy.arange(4, 4 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 20000000.0],
                             'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                      },
          'all_maps' : {  
                  'map1'   : {
                              'attr'    : 'em_fluxKkms_CO1-0', #'mass', 'G0', 'gmech', 'Av'
                              'v_rng'   : [-10.0, 4.0],
                              'title'   : r'$f(L_{CO(1-0} K.km.s-1))$', 
                              'as_log10': True,
                              'func'    : numpy.mean,
                             },
                        },
        'save_maps' : False,
        'save_image': False,
        }

#############################################################################################################
#############################################################################################################
#############################################################################################################

#setting up the logger object
logger = default_logger()

bs_min, bs_max = params['ranges']['box_size'].number

def generate_maps(snap_index, params):
    
    #path to processed fi snapshot  
    snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snap_index + '.states.npz'  
    
    #loading the processed sph simulation data with the emissions 
    logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
    gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'], load_pdr=params['pdr_sph'])    
    logger.debug('done reading fi snapshot : %s' % snap_filename)
    logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

    #keeping gas particles within the specified ranges
    gas = fi_utils.select_particles(gas, params['ranges'])
    logger.debug('got the sph particles in the required ranges')
    logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

    #making the 2D histogram
    print 'getting the spatial distrubutions'
    hist = hist_nd(numpy.vstack((gas.x, gas.y)), mn = bs_min, mx=bs_max, nbins=params['imres'], reverse_indicies=True, loc=True)
    hist.info()
    print 'done getting the spatial distributuions'

    #getting all the map
    for i, this_map in enumerate(params['all_maps']):
        
        this_map_info = params['all_maps'][this_map]
        
        map_data = fi_utils.make_map(gas, hist, **this_map_info)
        
        fi_utils.plot_map(map_data, params, this_map_info, 
                          fi_utils.get_snapshot_time(snap_index, params), 
                          params,
                          snap_filename
                          )
        
    return gas, hist
#

for snap in params['snaps']:        
    
    gas, hist = generate_maps(snap, params)

    #compute the mean excitation temperature map
    map_Tex   = numpy.zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)
    map_gmech = numpy.zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)
    map_gmech_per_H = numpy.zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)
    
    Tex   = getattr(gas, 'em_Tex_CO1-0')
    em    = getattr(gas, 'em_fluxKkms_CO1-0')
    ncm3  = getattr(gas, 'n')
    gmech = getattr(gas, 'gmech')
    g0    = getattr(gas, 'G0')
    
    for i in numpy.arange(hist.nBins[0]):
            
        for j in numpy.arange(hist.nBins[1]):
                
            inds_in_bin = hist.get_indicies([i,j])
                
            if inds_in_bin.size > 0:
                
                map_Tex[i,j] = numpy.average(Tex[inds_in_bin], weights=em[inds_in_bin])
                map_gmech_per_H[i,j] = numpy.mean(gmech[inds_in_bin]/ncm3[inds_in_bin])
                #map_Tex[i,j] = numpy.mean(Tex[inds_in_bin])#, weights=em[inds_in_bin])
                #map_gmech_per_H[i,j] = numpy.mean(gmech[inds_in_bin]/ncm3[inds_in_bin])
            #
        #
    #
    
    #compute the mean gmech map
    x_pix, y_pix = hist.f.cntrd

    r_pix = numpy.sqrt(x_pix**2 + y_pix**2)    
    #
#

########################################################################################

fig, axs = pylab.subplots(1, 3, sharex=False, sharey=False, figsize=(12.0, 4.0))

pylab.subplots_adjust(left=0.10, bottom=0.15, right=0.95, top=0.9, wspace=0.25, hspace=0.15)

axs[0].set_xlabel(r'R(kpc)')
axs[0].set_ylabel(r'<$\Gamma_m / n $>', size=10)

axs[1].set_xlabel(r'<$\Gamma_m / n $>', size=10)
axs[1].set_ylabel(r'$T_{ex}$ (K)', size=10)

axs[2].set_xlabel(r'R(kpc)', size=10)
axs[2].set_ylabel(r'$T_{ex}$ (K)', size=10)

r = r_pix.flatten()
gm_per_h = map_gmech_per_H.flatten()
Tex = map_Tex.flatten()

#######
inds = numpy.where( (r > 0.0)*(r < 1.0) )
axs[0].semilogy(r[inds][::2]       , gm_per_h[inds][::2], 'r.')
axs[1].semilogx(gm_per_h[inds], Tex[inds], 'r.')
axs[2].plot(r[inds], Tex[inds], 'r.')

inds = numpy.where( (r > 1.0)*(r < 2.0) )
axs[0].semilogy(r[inds][::2]       , gm_per_h[inds][::2], 'g.')
axs[1].semilogx(gm_per_h[inds], Tex[inds], 'g.')
axs[2].plot(r[inds], Tex[inds], 'g.')

inds = numpy.where( (r > 2.0)*(r < 3.0) )
axs[0].semilogy(r[inds][::3]       , gm_per_h[inds][::3], 'b.')
axs[1].semilogx(gm_per_h[inds], Tex[inds], 'b.')
axs[2].plot(r[inds], Tex[inds], 'b.')

inds = numpy.where( (r > 3.0)*(r < 10.0) )
axs[0].semilogy(r[inds][::5]       , gm_per_h[inds][::5], 'c.')
axs[1].semilogx(gm_per_h[inds], Tex[inds], 'c.')
axs[2].plot(r[inds], Tex[inds], 'c.')

#######
axs[0].set_xlim([0, 8])
axs[0].set_ylim([1e-28, 1e-25])
axs[0].set_yticks(axs[0].get_yticks()[::2])

axs[1].set_xlim([1e-28, 1e-25])
axs[1].set_ylim([10.0, 100.0])
axs[1].set_xticks(axs[1].get_xticks()[::2])

axs[2].set_xlim([0, 8])
axs[2].set_ylim([10.0, 100.0])


###########plotting the mean Tex as a function of distance from the center###########


r = r.reshape(1, r.size)
Tex = Tex.flatten()

print 'getting the spatial distrubutions....'
hist = hist_nd(r, nbins=25, mn=0.0, reverse_indicies=True, loc=True)
hist.info()

r = r[0, hist.inds_in]
Tex = Tex[hist.inds_in]

Tex_mean = numpy.zeros(hist.nBins, 'f8')
#looping over the sectors and getting the emissions for this sector    
for i in numpy.arange(hist.nBins):

    inds_in_bin = hist.get_indicies(i)
    
    Tex_mean[i] = numpy.mean(Tex[inds_in_bin])
    
print '\tfinished computing the emissions'

axs[2].plot(hist.f.cntrd, Tex_mean, 'k--', linewidth=3)


pylab.show()

