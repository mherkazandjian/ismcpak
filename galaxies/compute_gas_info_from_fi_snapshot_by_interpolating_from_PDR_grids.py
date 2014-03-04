'''
- load pdr meshes DB
- make interpolation functions 
- compute the emissions of particles 
- store the emissions into files
'''
#########################################################################################################
import matplotlib
matplotlib.use('Qt4Agg')

import scipy
import numpy
import time

from amuse.units import units

from mylib.utils.misc  import default_logger

import fi_utils
import meshUtils

#======================================================parameters=================================================

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std-test', # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          #'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-test',  # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-ext-test',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',          # the path to the dir containing the PDR database
          
          'use_em'  : False, 
          'use_pdr' : False,

          'use_sampled_set' : True,
          
          'snaps'  : numpy.arange(4, 4+1, 1),
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av_use'         :  [0.0, 200000.0],
                             'Av_clip'        :  [0.01, 28.0],  #sph particles with Av higher than this are clipped to this value                             
                            },
                      
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-20.0, 20.0] | units.kpc, #kpc

                       #models in the PDR arxv which are within those ranges will be used in constructing interpolation functions 
                      'interp'   : {'log_n'        : [-10.0,  10.0],
                                    'log_G0'       : [-2.0 ,  3.0 ],
                                    'log_gmech'    : [-50.0, -20.0],
                                    'Av'           : [ 1.0,  28.0],
                                    'Av_res'       : 1.0,   #resolution of constructing the interp funcs from the PDR data
                                    'extra_Av_sec' : [0.01, 0.1],
                                    },   #Av_res is used only for quantities interpolated upon from PDR meshes 
                      },
          
          'emission': {
#                      'f_mean_em_<q>fluxcgs</q>_CO<template>'  : 0,
#                      'f_mean_em_<q>fluxKkms</q>_CO<template>' : 1,
#                      'f_mean_em_no_gm_<q>fluxKkms</q>_CO<template>' : 1,                      
#                      'f_mean_em_<q>fluxKkms</q>_CO<template>' : 15,
#                      'f_mean_em_<q>fluxKkms</q>_HCN<template>' : 1,
#                      'f_mean_em_<q>Tex</q>_CO<template>'      : 0,
#                      'f_mean_em_<q>tau</q>_CO<template>'      : 0,
##                      
#                      'f_mean_em_<q>fluxcgs</q>_13CO<template>' : 0,
#                      'f_mean_em_<q>fluxKkms</q>_13CO<template>': 0,
#                      'f_mean_em_<q>Tex</q>_13CO<template>'     : 0,
#                      'f_mean_em_<q>tau</q>_13CO<template>'     : 0,
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
#                      'f_mean_em_<q>fluxcgs</q>_CO<template>'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CO<template>'  : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_HCN<template>'  : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_no_gm_<q>fluxKkms</q>_CO<template>'  : {'pos': [2,0], 'title': r'$f(L_{CO(1-0}_{no_gm})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_fluxKkms_CO1-0'  : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_<q>Tex</q>_CO<template>'       : {'pos': [2,1], 'title': r'$f(Tex_{CO(1-0})$'   , 'v_rng': [0.0, 10000.0], 'log10': False},
#                      'f_mean_em_<q>tau</q>_CO<template>'       : {'pos': [2,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 1000.0], 'log10': False},
                      ########
#                      'f_mean_em_<q>fluxcgs</q>_13CO<template>'   : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_13CO<template>'  : {'pos': [3,2], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
#                      'f_mean_em_<q>Tex</q>_13CO<template>'       : {'pos': [3,1], 'title': r'$f(Tex_{13CO(1-0})$' , 'v_rng': [0.0, 10000.0], 'log10': False},                      
#                      'f_mean_em_<q>tau</q>_13CO<template>'       : {'pos': [3,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 1000.0] , 'log10': False},
                      ########       
                      #'f_mean_em_<q>fluxcgs</q>_CO1-0'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0 [erg.cm^2.s-1]})$'       , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_<q>fluxKkms</q>_CO1-0'  : {'pos': [2,1], 'title': r'$f(L_{CO(1-0}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN1_0.5-0_0.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(1-0}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN1_1.5-0_0.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(2-1}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN2_1.5-1_1.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(3-2}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN2_1.5-1_0.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(4-3}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN2_2.5-1_1.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(5-4}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN3_2.5-2_2.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(6-5}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN3_2.5-2_1.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(7-6}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
#                      'f_mean_em_<q>fluxKkms</q>_CN3_3.5-2_2.5'  : {'pos': [2,1], 'title': r'$f(L_{CN(8-7}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
                      #'f_mean_em_<q>Tex</q>_CO1-0'       : {'pos': [2,2], 'title': r'$f(Tex_{CO(1-0} [K])$'                , 'v_rng': [0.0, 10000.0], 'log10': False},
                      #'f_mean_em_<q>tau</q>_CO1-0'       : {'pos': [2,3], 'title': r'$f(tau_{CO(1-0} [K])$'                , 'v_rng': [0.0, 1000.0], 'log10': False},
                      #'f_mean_em_<q>fluxcgs</q>_13CO1-0'       : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_pdr_NH'   : {'pos': [3,1], 'title': r'$log10 N(H2) [cm^{-2}]$'  , 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      #'f_mean_pdr_NH2'  : {'pos': [3,1], 'title': r'$log10 N(H2) [cm^{-2}]$'  , 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},                      
                      #'f_mean_pdr_NCO'  : {'pos': [3,0], 'title': r'$log10 N(CO) [cm^{-2}]$'  , 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      #'f_mean_pdr_N13CO': {'pos': [3,2], 'title': r'$log10 N(13CO) [cm^{-2}]$', 'v_rng': [10.0, 27.0], 'log10': True, 'log10_interp_func':True},
                      
                      ## temperature at the deepest part of the cloud
                      'f_mean_pdr_T_molec': {'pos': [3,2], 'title': r'$log10 T_${\rm molecular gas}$ [K]', 'v_rng': [1.0, 4.0], 'log10': True, 'log10_interp_func':True},

                      #the empdr thing is not implemented yet!!!
                      },

          #'interpolator' : scipy.interpolate.NearestNDInterpolator, 
          'interpolator' : scipy.interpolate.LinearNDInterpolator, 
          'save_info'    : False,
#          'save_secies'  : ['CO', '13CO', 'HCN'],
          'save_secies'  : [],
          }

#fluxcgs, fluxKkms, tau, Tex, T_R

#replacing the template maps info with the ones of the actual lines
if 'emission' in params:

    em = params['emission']
    
    for em_key in em:

        if em_key in params['maps']:
        
            template = params['maps'].pop(em_key)
            print 'using this template to make the range in maps', template
            
            for trans_num in numpy.arange(em[em_key]):
                em_key_this = em_key.replace('<template>','%d-%d' % (trans_num+1, trans_num)) 
                
                params['maps'][em_key_this] = template
            #
        #
    #
#
######################################################################################################## 
#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

#setting up the logger object
logger = default_logger()

########################################reading the data of the PDR modes###############################
#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

#=============================================
#loading the molecular specie radex database
#arxvPDR.readDbsRadex(allDbs=True)
#arxvPDR.readDbsRadex(species=['CO','HCN', 'HCO+'], in_Av_rng=params['ranges']['interp']['Av'])
arxvPDR.readDbsRadex(species=['CO'], in_Av_rng=params['ranges']['interp']['Av'])
########################################################################################################

#the keys of the emission maps only
if params['use_em'] == True:
    em_keys = fi_utils.get_emission_only_keys(params['maps'])
else:
    print '\n\n\n\n\n\nflag to compute emissions is set to false, NOT computing emissions'
    time.sleep(0.5)
    em_keys = {}

#the keys of the maps whose quantites will be extraceted from the PDR models
if params['use_pdr'] == True:
    pdr_keys = fi_utils.get_pdr_only_keys(params['maps'])
else:
    print '\n\n\n\n\n\nflag to compute info from PDR datase is set to false, NOT computing pdr info'
    time.sleep(0.5)
    pdr_keys = {}

#making all the interpolation functions
em_interp_funcs, pdr_interp_funcs = fi_utils.make_interp_funcs_from_arxv(arxvPDR,
                                                                         em_keys,
                                                                         pdr_keys, 
                                                                         params,
                                                                         logger)

for snap in params['snaps']:

    em_sph, pdr_sph, gas = fi_utils.snapshot_interpolated_data(
                                                               snap, arxvPDR,
                                                               em_interp_funcs, pdr_interp_funcs, 
                                                               em_keys, pdr_keys,
                                                               params, logger
                                                               )
       
    ## setting the emission info and the pdr info as attributes to the 'gas' particle set
    all_info = dict(em_sph.items() + pdr_sph.items()) 
    for key in all_info:
        print key
        attr_name = key.replace('f_mean_','').replace('</q>','').replace('<q>','')
        
        setattr(gas, attr_name, all_info[key])

    ## saving the info interpolated from the PDR grids
    if params['save_info'] == True:
        snap_filename = params['rundir'] + '/firun/fiout.%06d' % snap
        fi_utils.save_gas_particle_info_saperate_files(snap_filename, gas, params['save_secies'])
#

####add a plotting routine here for the masp
