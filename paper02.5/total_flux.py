#########################################################################################################

import matplotlib
matplotlib.use('Qt4Agg')

import numpy
import pylab

from amuse.units import units, nbody_system
from mylib.utils.misc  import default_logger

from galaxies import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

home = '/home/mher'

fig_save_path = '/home/mher/ism/docs/paper02.5/src/figs/total_flux_vs_J.eps'

#################################################the dwarf galaxy#####################################################
params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'snaps' : numpy.arange(20, 20 + 1, 1),
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          'species' : ['CO', '13CO'],
          
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [3.0, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-2.0, 2.0] | units.kpc, #kpc
                      },

          'J_max': 8,  
        }

#setting up the logger object
logger = default_logger()

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

bs_min, bs_max = params['ranges']['box_size'].number

snapIndex = 4    
#path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snapIndex + '.states.npz'  

#loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'])   
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

#computing the time of the snapshot
snap_time = fi_utils.get_snapshot_time(snapIndex, params)

#keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

#the upper transitions
Ju_all  = numpy.arange(params['J_max']) + 1
flux_CO_dwarf = numpy.zeros(Ju_all.size, 'f8')
flux_13CO_dwarf = numpy.zeros(Ju_all.size, 'f8')

#getting the CO ladder for the whole galaxy 
for Ju in Ju_all:
    
    #emissions of CO of all the gas particles
    em_all_this_transition = getattr(gas, 'em_fluxKkms_CO%d-%d' % (Ju, Ju-1))
    flux_CO_dwarf[Ju-1] = em_all_this_transition.sum()

    #emissions of 13CO of all the gas particles
    em_all_this_transition = getattr(gas, 'em_fluxKkms_13CO%d-%d' % (Ju, Ju-1))
    flux_13CO_dwarf[Ju-1] = em_all_this_transition.sum()

#################################################the dwarf galaxy#####################################################
params = {'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol', # the path of the dir containing the simulation
          'snaps' : numpy.arange(4, 4 + 1, 1),
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          'species' : ['CO', '13CO'],
          
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200.0],
                             'Av_clip'        :  [3.0, 29.9],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-8.0, 8.0] | units.kpc, #kpc
                      },

          'J_max': 8,  
        }

#setting up the logger object
logger = default_logger()

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

bs_min, bs_max = params['ranges']['box_size'].number

snapIndex = 4    
#path to processed fi snapshot  
snap_filename = params['rundir'] + '/firun/' + 'fiout.%06d' % snapIndex + '.states.npz'  

#loading the processed sph simulation data with the emissions 
logger.debug('loading proccessed snapshot %s : ' % snap_filename) 
gas = fi_utils.load_gas_particle_info_with_em(snap_filename, params['species'])   
logger.debug('done reading fi snapshot : %s' % snap_filename)
logger.debug('number of sph particles in proccessed snapshot = %d' %  len(gas))

#computing the time of the snapshot
snap_time = fi_utils.get_snapshot_time(snapIndex, params)

#keeping gas particles within the specified ranges
gas = fi_utils.select_particles(gas, params['ranges'])
logger.debug('got the sph particles in the required ranges')
logger.debug('number of gas particles in the specified ranages = %d' %  len(gas))

#the upper transitions
Ju_all  = numpy.arange(params['J_max']) + 1
flux_CO_disk = numpy.zeros(Ju_all.size, 'f8')
flux_13CO_disk = numpy.zeros(Ju_all.size, 'f8')
flux_CO_disk_center = numpy.zeros(Ju_all.size, 'f8')
flux_13CO_disk_center = numpy.zeros(Ju_all.size, 'f8')

x = gas.x
y = gas.y 
R = numpy.sqrt(x**2.0 + y**2.0)
inds_center = numpy.where(R < 0.5)


#getting the CO ladder for the whole galaxy 
for Ju in Ju_all:
    
    #emissions of CO of all the gas particles
    em_all_this_transition = getattr(gas, 'em_fluxKkms_CO%d-%d' % (Ju, Ju-1))
    flux_CO_disk[Ju-1] = em_all_this_transition.sum()
    flux_CO_disk_center[Ju-1] = em_all_this_transition[inds_center].sum()

    #emissions of 13CO of all the gas particles
    em_all_this_transition = getattr(gas, 'em_fluxKkms_13CO%d-%d' % (Ju, Ju-1))
    flux_13CO_disk[Ju-1] = em_all_this_transition.sum()
    flux_13CO_disk_center[Ju-1] = em_all_this_transition[inds_center].sum()
    
    
###########################################plotting###########################################################


fig = pylab.figure(figsize=(4,4))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

ax.set_xlim([Ju_all.min(), Ju_all.max()])

pltCO, = ax.semilogy(Ju_all, flux_CO_disk, 'r-')
ax.text(Ju_all[6], flux_CO_disk[6], 'disk')

plt13CO, = ax.semilogy(Ju_all, flux_13CO_disk, 'b--')
ax.text(Ju_all[6], flux_13CO_disk[6], 'disk')

ax.semilogy(Ju_all, flux_CO_disk_center, 'r-o')
ax.text(Ju_all[3], flux_CO_disk_center[3], 'center')

ax.semilogy(Ju_all, flux_13CO_disk_center, 'b--o')
ax.text(Ju_all[1], flux_13CO_disk_center[1]/5.0, 'center')

ax.semilogy(Ju_all, flux_CO_dwarf, 'r-')
ax.semilogy(Ju_all, flux_13CO_dwarf, 'b--')
ax.text(Ju_all[1], flux_13CO_dwarf[1], 'dwarf')

ax.legend([pltCO, plt13CO], ['CO', '13CO'], loc=0)

ax.set_xlabel(r'J$_{\rm upper}$', size=10)
ax.set_ylabel(r'flux [K.km.s$^{-1}$]', size=10)

ax.tick_params(axis='both', which='major', labelsize=10)


pylab.draw()
pylab.show()
###############################################################################################################

fig.savefig(fig_save_path)
print 'saved image file to :\n\t\t\t %s' % fig_save_path
