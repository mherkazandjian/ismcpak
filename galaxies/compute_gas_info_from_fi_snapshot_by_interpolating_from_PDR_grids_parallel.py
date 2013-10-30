'''
- load pdr meshes DB
- make interpolation functions 
- compute the emissions of particles 
- store the emissions into files
'''
#########################################################################################################
import sys
from IPython.parallel import Client

import matplotlib
matplotlib.use('Qt4Agg')

import scipy
import numpy
import pylab

from amuse.units import units

from mylib.utils.misc  import default_logger

import fi_utils
import meshUtils

#======================================================parameters=================================================

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',  # the path to the dir containing the PDR database          
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
               
          'use_em'  : False, 
          'use_pdr' : True,          
          
          'snaps'  : numpy.arange(4, 4 + 1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 2000000.0],
                             'Av_clip'        :  [3.0, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-20.0, 20.0] | units.kpc, #kpc

                       #modes in the PDR arxv which are within those ranges will be used in constructing interpolation functions 
                      'interp'   : {'log_n'    : [-10.0,  10.0],
                                    'log_G0'   : [-2.0 ,  3.0 ],
                                    'log_gmech': [-50.0, -20.0],
                                    'Av'       : [ 3.0 ,  28.0],
                                    'Av_res'   : 1.0,   #resolution of constructing the interp funcs from the PDR data  
                                    }, #Av_res is used only for quantities interpolated upon from PDR meshes 
                      },
          
          'emission': {
                      'f_mean_em_<q>fluxcgs</q>_CO<template>'  : 15,
                      'f_mean_em_<q>fluxKkms</q>_CO<template>' : 15,
                      'f_mean_em_<q>Tex</q>_CO<template>'      : 15,
                      'f_mean_em_<q>tau</q>_CO<template>'      : 15,
##                      
                      'f_mean_em_<q>fluxcgs</q>_13CO<template>' : 15,
                      'f_mean_em_<q>fluxKkms</q>_13CO<template>': 15,
                      'f_mean_em_<q>Tex</q>_13CO<template>'     : 15,
                      'f_mean_em_<q>tau</q>_13CO<template>'     : 15,
                      },
          
          'maps'   : {
                      'f_n_part' : {'pos': [0,0], 'title': r'$f(N)$'        , 'v_rng': [0,6]      , 'log10': True}, 
                      'f_mass'   : {'pos': [0,1], 'title': r'$f(m)$'        , 'v_rng': [30,40]    , 'log10': True},
                      'f_mean_n' : {'pos': [0,2], 'title': r'$f(\bar{n})$'  , 'v_rng': [-3.0, 4.0], 'log10': True},
                      'f_mean_g0': {'pos': [0,3], 'title': r'$f(\bar{g_0})$', 'v_rng': [-3.0, 3.0], 'log10': True},
                      #--------
                      'f_mean_gm'              : {'pos': [1,0], 'title': r'$f(\bar{\Gamma_m})$'                , 'v_rng': [-35.0, -22.0], 'log10': True},
                      'f_mean_Av'              : {'pos': [1,1], 'title': r'$f(\bar{Av})$'                      , 'v_rng': [0.0  , 2.0] ,  'log10': True},
                      'T_mean'                 : {'pos': [1,2], 'title': r'$f(\bar{T})$'                        , 'v_rng': [0.0  , 5.0] ,  'log10': True},
                      #'f_mean_em_C+-1-0'       : {'pos': [1,1], 'title': r'$f(L_{C^+ 158 \mu m})$'               , 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_no_gm_C+-1-0' : {'pos': [1,2], 'title': r'$f(L_{C^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_CO-7-6'       : {'pos': [1,3], 'title': r'$f(L_{O^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #--------
                      'f_mean_em_<q>fluxcgs</q>_CO<template>'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>fluxKkms</q>_CO<template>'  : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>Tex</q>_CO<template>'       : {'pos': [2,1], 'title': r'$f(Tex_{CO(1-0})$'   , 'v_rng': [0.0, 10000.0], 'log10': False},
                      'f_mean_em_<q>tau</q>_CO<template>'       : {'pos': [2,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 1000.0] , 'log10': False},
                      ########
                      'f_mean_em_<q>fluxcgs</q>_13CO<template>'   : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>fluxKkms</q>_13CO<template>'  : {'pos': [3,2], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>Tex</q>_13CO<template>'       : {'pos': [3,1], 'title': r'$f(Tex_{13CO(1-0})$' , 'v_rng': [0.0, 10000.0], 'log10': False},                      
                      'f_mean_em_<q>tau</q>_13CO<template>'       : {'pos': [3,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 1000.0] , 'log10': False},
                      ########
                      #'f_mean_em_<q>fluxcgs</q>_CO1-0'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0 [erg.cm^2.s-1]})$'       , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_<q>fluxKkms</q>_CO1-0'  : {'pos': [2,1], 'title': r'$f(L_{CO(1-0}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
                      #'f_mean_em_<q>Tex</q>_CO1-0'       : {'pos': [2,2], 'title': r'$f(Tex_{CO(1-0} [K])$'                , 'v_rng': [0.0, 10000.0], 'log10': False},
                      #'f_mean_em_<q>tau</q>_CO1-0'       : {'pos': [2,3], 'title': r'$f(tau_{CO(1-0} [K])$'                , 'v_rng': [0.0, 1000.0], 'log10': False},
                      #--------
                      #'f_mean_em_<q>fluxcgs</q>_CO3-2'   : {'pos': [3,0], 'title': r'$f(L_{CO(3-2 [erg.cm^2.s-1]})$'       , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_<q>fluxcgs</q>_13CO1-0'       : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
                      'f_mean_pdr_NCO'   : {'pos': [3,0], 'title': r'$log10 N(CO) [cm^{-2}]$'  , 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      'f_mean_pdr_NH2'   : {'pos': [3,1], 'title': r'$log10 N(H2) [cm^{-2}]$'  , 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      'f_mean_pdr_N13CO' : {'pos': [3,2], 'title': r'$log10 N(13CO) [cm^{-2}]$', 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      },

          #'interpolator' : scipy.interpolate.NearestNDInterpolator,     
          'interpolator' : scipy.interpolate.LinearNDInterpolator, 
          'save_info'   : True,
          'save_secies' : ['CO', '13CO']
          }

#fluxcgs, fluxKkms, tau, Tex, T_R

#replacing the template maps info with the ones of the actual lines
em = params['emission']
for em_key in em:
    
    if em_key in params['maps']:
    
        template = params['maps'].pop(em_key)
        print 'using this template to make the range in maps', template
        
        for trans_num in numpy.arange(em[em_key]):
            em_key_this = em_key.replace('<template>','%d-%d' % (trans_num+1, trans_num)) 
            
            print em_key_this
            
            params['maps'][em_key_this] = template
        #
    #
#
######################################################################################################## 
#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

#setting up the logger object
logger = default_logger()

#######################################setting up the parallel environment##############################
#all the engines
clients = Client()

#a subset of the engines can be used by slicing is returned as a 'direct view'
dview = clients[:]

#syncing the imports
with dview.sync_imports():
    import sys
dview.push(dict(syspath=sys.path))
dview.execute('sys.path=syspath', block=True)

with dview.sync_imports():
    import fi_utils
    import sys, time
    from amuse.units import units
    from mylib.utils.misc  import default_logger
    import meshUtils

#pushing the main parameter dict to all the engines
dview.push(dict(params=params), block=True)

dview.execute('logger = default_logger()', block=True)

########################################reading the data of the PDR modes###############################
#reading and setting up the pdr database
print 'reading the PDR databases in parallel'
dview.execute('arxvPDR = meshUtils.meshArxv(dirPath = params["pdrDb"], readDb=True)', block=True)

#=============================================
#loading the molecular specie radex database
#arxvPDR.readDbsRadex(allDbs=True)
#arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'], Av=params['ranges']['interp']['Av'])
#arxvPDR.readDbsRadex(species=['CO','HCN', 'HCO+'], in_Av_rng=params['ranges']['interp']['Av'])
dview.execute('arxvPDR.readDbsRadex(species=["CO", "13CO"], in_Av_rng=params["ranges"]["interp"]["Av"])', block=True)
########################################################################################################
print 'read the DBs on all the engines'

#the keys of the emission maps only
if params['use_em'] == True:
    em_keys = fi_utils.get_emission_only_keys(params['maps'])
else:
    print '\n\n\n\n\n\nflag to compute emissions is set to false, NOT computing emissions'
    time.sleep(0.5)
    em_keys = {}

#scattering the emission keys (the maps corresponding to those will use interpolation to be computed)
dview.scatter('em_keys', em_keys)

#the keys of the maps whose quantites will be extraceted from the PDR models
if params['use_pdr'] == True:
    pdr_keys = fi_utils.get_pdr_only_keys(params['maps'])
else:
    print '\n\n\n\n\n\nflag to compute info from PDR datase is set to false, NOT computing pdr info'
    time.sleep(0.5)
    pdr_keys = {}

#scattering the pdr info keys (the maps corresponding to those will use interpolation to be computed)
dview.scatter('pdr_keys', pdr_keys)

#making all the interpolation functions
print 'constructing the interpolation function on all engines'
t0 = time.time()
dview.execute('em_interp_funcs, pdr_interp_funcs = fi_utils.make_interp_funcs_from_arxv(arxvPDR, em_keys, pdr_keys, params, logger)', block=True)
print 'done contructing interpolation functions in %.2e seconds' % (time.time() - t0)


'''computing the emissions of the sph particles from the interpolation functions in parallel (Distributing 
the tasks over the engines).
'''
print 'interpolating info of the gas particles using the interpolation functions constructed from the PDR grids....'
t0 = time.time()

for snap in params['snaps']:
    
    dview.execute('snap=%d' % snap, block=True)
    dview.execute('output=fi_utils.snapshot_interpolated_data(snap, arxvPDR, em_interp_funcs, pdr_interp_funcs, em_keys, pdr_keys, params, logger)', block=True)
    dview.execute('em_sph=output[0]', block=True)
    dview.execute('pdr_sph=output[1]', block=True)
    dview.execute('gas=output[2]', block=True)
    print 'done interpolating all the emissions in %.2e seconds' % (time.time() - t0)

    #collecting all the computed emissions to this process    
    em_sph = {}
    for em_this_view in dview['em_sph']: #looping over the list of data from the views
        for em in em_this_view: #looping over the dict items in each em info computed in this view 
            em_sph[em] = em_this_view[em]

    #collecting all the computed pdr info to this process    
    pdr_sph = {}
    for pdr_sph_this_view in dview['pdr_sph']: #looping over the list of data from the views
        for pdr_info in pdr_sph_this_view: #looping over the dict items in each em info computed in this view 
            pdr_sph[pdr_info] = pdr_sph_this_view[pdr_info]
    
    #reading the sph gas info from the zeroth engine
    gas = fi_utils.get_useful_gas_attr_from_dview(dview, 'gas', 0)
    
    
    #setting the emission info and the pdr info as attributes to the 'gas' particle set
    all_info = dict(em_sph.items() + pdr_sph.items()) 
    for key in all_info:
        print key
        attr_name = key.replace('f_mean_','').replace('</q>','').replace('<q>','')
        
        setattr(gas, attr_name, all_info[key])
    
    #saving the emissions into files
    if params['save_info'] == True:        
        snap_filename = params['rundir'] + '/firun/fiout.%06d' % snap  
        fi_utils.save_gas_particle_info_saperate_files(snap_filename, gas, params['save_secies'])
    
#plotting the emissions
#fi_utils.plot_maps(snap, em_sph, gas, params, logger)
#pylab.show()

print 'done'
