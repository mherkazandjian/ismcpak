'''
 in this test file, we load a database of PDR meshes and extract information
 at conditions determined by the SPH simulations

  - exclude sph particles which have densities below which the LVG model codes can handel 
    (for example, one with n_gas < 1e-3))
  - exlcude particle with alpha > 1 (unless for very high densities (get the densitiy of 
    the specie whose  emission is to be computed and see if its abundance is too low, if 
    if is too low, just set the emission of the sph particle corresponding to that specie
    to zero.
  - abundances of CO for PDR modesls with n < 0.1 cm^-3 are ~ 1e-10 (which is too low compared to the
    high density modesl

'''
#########################################################################################################
import time
import sys
import os
from IPython.parallel import Client

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate
import numpy
from numpy import log10
import pylab
import logging

from amuse.io import read_set_from_file, fi_io
from amuse.units import units, constants, nbody_system

from mylib.utils.misc  import xselect, default_logger
from mylib.utils.histogram import hist_nd
from mylib.utils.interpolation import sectioned_4D_interpolator 
from fi_utils import parse_old_runinfo_file
import fi_utils
import lineDict
import meshUtils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

pool_size = 1
snaps = numpy.arange(4, 4 + 1, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 50,                                                 # resolution of the maps to be produced imres x imres
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',      # the path to the dir containing the PDR database
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2/',          # the path to the dir containing the PDR database          
          'use_em' : True, 
          
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

                       #modes in the PDR arxv which are within those ranges will be used in constructing interpolation functions 
                      'interp'   : {'log_n'    : [-10.0,  10.0],
                                    'log_G0'   : [-2.0 ,  3.0 ],
                                    'log_gmech': [-50.0, -20.0],
                                    'Av'       : [ 3.0 ,  30.0],
                                    },
                      },
          
          'emission': {'f_mean_em_<q>fluxcgs</q>_CO<template>' : 3,
                      'f_mean_em_<q>fluxKkms</q>_CO<template>' : 3,
                      'f_mean_em_<q>Tex</q>_CO<template>'      : 3,
                      'f_mean_em_<q>tau</q>_CO<template>'      : 3,
                      
                      'f_mean_em_<q>fluxcgs</q>_13CO<template>' : 3,
                      'f_mean_em_<q>fluxKkms</q>_13CO<template>': 3,
                      'f_mean_em_<q>Tex</q>_13CO<template>'     : 3,
                      'f_mean_em_<q>tau</q>_13CO<template>'      : 3,

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
                      'f_mean_em_<q>fluxKkms</q>_CO<template>'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'     , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>Tex</q>_CO<template>'       : {'pos': [2,1], 'title': r'$f(Tex_{CO(1-0})$'   , 'v_rng': [0.0, 10000.0], 'log10': False},
                      'f_mean_em_<q>tau</q>_CO<template>'      : {'pos': [2,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 100.0], 'log10': False},
                      ########
                      'f_mean_em_<q>fluxcgs</q>_13CO<template>'   : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>Tex</q>_13CO<template>'       : {'pos': [3,1], 'title': r'$f(Tex_{13CO(1-0})$' , 'v_rng': [0.0, 10000.0], 'log10': False},
                      'f_mean_em_<q>fluxKkms</q>_13CO<template>'  : {'pos': [3,2], 'title': r'$f(L_{13CO(1-0})$'   , 'v_rng': [-10.0, -4.0] , 'log10': True},
                      'f_mean_em_<q>tau</q>_13CO<template>'      : {'pos': [3,3], 'title': r'$f(tau_{CO(1-0})$'   , 'v_rng': [0.0, 100.0], 'log10': False},
                      ########                                            
                      #'f_mean_em_<q>fluxcgs</q>_CO1-0'   : {'pos': [2,0], 'title': r'$f(L_{CO(1-0 [erg.cm^2.s-1]})$'       , 'v_rng': [-10.0, -4.0], 'log10': True},
                      #'f_mean_em_<q>fluxKkms</q>_CO1-0'  : {'pos': [2,1], 'title': r'$f(L_{CO(1-0}[K.km.s-1)$'             , 'v_rng': [-10.0, 4.0] , 'log10': True},
                      #'f_mean_em_<q>Tex</q>_CO1-0'       : {'pos': [2,1], 'title': r'$f(Tex_{CO(1-0} [K])$'                , 'v_rng': [0.0, 1000.0], 'log10': False},
                      #--------
                      #'f_mean_em_HCN-1-0'       : {'pos': [3,0], 'title': r'$f(L_{HCN(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_<q>fluxcgs</q>_13CO1-0'       : {'pos': [3,0], 'title': r'$f(L_{13CO(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
                      },

          'interpolator' : interpolate.NearestNDInterpolator, 
          #'interpolator' : interpolate.LinearNDInterpolator, 
          'image_save'  : False,
          'image_ext'   : 'eps',
          'save_info'   : False,
          'save_secies' : ['CO', '13CO']
          }

#fluxcgs, fluxKkms, tau, Tex, T_R
#those are the lines which will be computed
lines_info = lineDict.lines

#replacing the template maps info with the ones of the actual lines
em = params['emission']
for em_key in em:
    
    if em_key in params['maps']:
    
        template = params['maps'].pop(em_key)
        print 'using this template to make the range in maps', template
        
        for trans_num in numpy.arange(em[em_key]):
            em_key_this = em_key.replace('<template>','%d-%d' % (trans_num+1, trans_num)) 
            
            params['maps'][em_key_this] = template
            #para
        #
    #
#
#adasd 
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

#setting up the logger object
logger = default_logger()

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

lgMechMin = arxvPDR.grid_z.min()
bs_min, bs_max = params['ranges']['box_size'].number

if params['use_em'] == True:    
    #=============================================
    #loading the molecular specie radex database
    #arxvPDR.readDbsRadex(allDbs=True)
    #arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'], Av=params['ranges']['interp']['Av'])
    #arxvPDR.readDbsRadex(species=['CO','HCN', 'HCO+'], in_Av_rng=params['ranges']['interp']['Av'])
    arxvPDR.readDbsRadex(species=['CO'], in_Av_rng=params['ranges']['interp']['Av'])
    
def funcRadex(mesh_radex, **kwargs):
    
    if mesh_radex == None:
        return numpy.nan
    else:
        transition_info = mesh_radex[kwargs['transitionIdx']] 
        quantity = transition_info[kwargs['quantity']]  #a radex quantity
        
        if kwargs['quantity'] == 'fluxcgs' or kwargs['quantity'] == 'fluxKkms':
            return log10(quantity) 
        else:
            return quantity

def funcPDR(meshObj, **kwargs):
    
    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = meshObj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    
    return log10(value)

def get_interpolation_function_pdr(**kwargs):
    """makes an interpolation function from a integrated PDR quantity"""
    
    v = arxvPDR.apply_function_to_all_meshes(funcPDR, func_kw = kwargs)

    v = numpy.array(v)
    inds_valid = numpy.isfinite(v)
    v = v[inds_valid] 
    
    xGrd, yGrd, zGrd = arxvPDR.grid_x[inds_valid], arxvPDR.grid_y[inds_valid], arxvPDR.grid_z[inds_valid]
    data = numpy.array([xGrd, yGrd, zGrd], dtype = numpy.float64).T
    
    t0 = time.time()
    
    if 'sectioned' in kwargs and kwargs['sectioned'] == True:
        fInterp = sectioned_4D_interpolator(data, v, params['interpolator'])
    else:
        fInterp = params['interpolator'](data, v)
        
    logger.debug('time contruncting interpolator for quantity %s %s in %.4f seconds.' % (kwargs['quantity'][1], kwargs['quantity'][3], time.time() - t0))

    return fInterp

def get_interpolation_function_radex(arxvPDR, params, **kwargs):
    """makes an interpolation function from the radex meshes from the database.
    :param string specStr: The specie whose info will be used (emission, or anything else radex outputs)
    :param int transitionIdx: The index of the transition (for example for CO: transitionIdx = 0, uses info of the J=1-0 transition 
    :param bool sectioned: if true uses a sectioned 4D interpolator. 
    """
    
    specStr  = kwargs.pop('specStr')
    Av_range = params['ranges']['interp']['Av']
    
    #reading only the radex databases which are within the specified range in Av
    available_Avs_for_radex_dbs = numpy.sort(numpy.array(arxvPDR.radexDbs.keys(), dtype='f8'))
    inds_read = numpy.where(
                            (available_Avs_for_radex_dbs >= Av_range[0])*
                            (available_Avs_for_radex_dbs <= Av_range[1])
                           )[0]
    Avs_read = available_Avs_for_radex_dbs[inds_read]

    #collecting the data from the database corresponding to all the data
    #in Av avaiable
    for i, Av in enumerate(Avs_read):
        
        #using the database of the specific Av
        arxvPDR.use_radexDb(specStr=specStr, Av=Av)
        
        #getting the data corresponding to this Av
        v = arxvPDR.apply_function_to_all_radex_meshes(funcRadex, func_kw = kwargs)
                    
        #keeping the points which are useful (finite ones) 
        v = numpy.array(v)
        inds_valid = numpy.isfinite(v)
        v = v[inds_valid]

        xGrd, yGrd, zGrd = arxvPDR.grid_x[inds_valid], arxvPDR.grid_y[inds_valid], arxvPDR.grid_z[inds_valid]
        AvGrd = numpy.ones(xGrd.shape, 'f')*Av
        data = numpy.array([xGrd, yGrd, zGrd, AvGrd], dtype = numpy.float64).T
                
        #soting the x,y,z,t, and v into arrays to be used later to construct the interpoaltion function
        if i == 0:
            data_all_Av = data
            v_all_Av = v
        else:
            data_all_Av = numpy.vstack( (data_all_Av, data) )
            v_all_Av = numpy.hstack( (v_all_Av, v) )
    #

    if kwargs['quantity'] == 'Tex' or kwargs['quantity'] == 'T_R':
        v_all_Av = v_all_Av.clip(0.0, 1e4)
    
    #using the subset of points in the radex databse specified by params['ranges']['interp']
    inds_keep = numpy.where( 
                             (data_all_Av[:,0] >= params['ranges']['interp']['log_n'][0])* 
                             (data_all_Av[:,0] <= params['ranges']['interp']['log_n'][1])*
                             
                             (data_all_Av[:,1] >= params['ranges']['interp']['log_G0'][0])*
                             (data_all_Av[:,1] <= params['ranges']['interp']['log_G0'][1])*
                             
                             (data_all_Av[:,2] >= params['ranges']['interp']['log_gmech'][0])* 
                             (data_all_Av[:,2] <= params['ranges']['interp']['log_gmech'][1])*
                             
                             (data_all_Av[:,3] >= params['ranges']['interp']['Av'][0])*
                             (data_all_Av[:,3] <= params['ranges']['interp']['Av'][1])
                           )[0]
    data_all_Av = data_all_Av[inds_keep,:] 
    v_all_Av    = v_all_Av[inds_keep] 
    #
    
    t0 = time.time()
    if 'sectioned' in kwargs and kwargs['sectioned'] == True:
        fInterp = sectioned_4D_interpolator(data_all_Av, v_all_Av, params['interpolator'])
    else:
        fInterp = params['interpolator'](data_all_Av, v_all_Av)
    logger.debug('time contruncting interpolation function for %s transition %d in %.4f seconds' % (specStr, kwargs['transitionIdx'], time.time() - t0))
    return fInterp

def make_emission_interp_funcs(params):
    '''returns a dictionary of functions with keys XYZ-UPPER-LOWER. Items in 'maps' which have
    _em_ are parsed and the transition is looked up (f_mean_em_C+-1-0 is parsed into C+-1-0) 
    and an interpolation function for C+-1-0 is computed.
    '''
        
    em_interp_funcs = {}
    
    for map_key in params['maps']:
        
        if '_em_' in map_key:

            line = map_key.split('_')[-1] # BLABLA_SPECSTRUPPER-LOWER
            
            if lines_info[line]['type'] == 'radex-lvg':
                em_interp_funcs[map_key] = get_interpolation_function_radex(
                                                                            arxvPDR,    
                                                                            params,
                                                                            specStr  = lines_info[line]['specStr'],
                                                                            quantity = map_key.split('<q>')[-1].split('</q>')[0],
                                                                            transitionIdx = lines_info[line]['radexIdx'],
                                                                            sectioned = True,
                                                                           )
            '''
            if lines_info[line]['type'] == 'pdr':
                em_interp_funcs[line] = get_interpolation_function_pdr(
                                                             quantity = lines_info[line]['quantity'],
                                                             sectioned = True,
                                                             )
            '''
            
    return em_interp_funcs

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

#constructing all the interpolation functions
if params['use_em']:
    logger.debug('getting interpolation functions of the emissions')
    em_interp_funcs = make_emission_interp_funcs(params)
    logger.debug('done getting interpolation functions of the emissions')


#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################

def plot_map_standalone(key, title, params):
    '''plots one of the generated maps. The data of the map should be generated a priori and set in 
    the apropriate params['maps'][]
    '''
    
    """
    plot_map_standalone('f_mean_gm'             , r'log$_{10}$ ${\Gamma_m}$ / erg cm$^{-3}$ s$^{-1}$')
    plot_map_standalone('f_mean_n'              , r'log$_{10}$ n$_{gas}$ / cm$^{-3}$')
    plot_map_standalone('f_mean_em_CO-1-0'      , r'log$_{10}$ CO(1-0) / erg cm$^{-2}$ s$^{-1}$')
    plot_map_standalone('f_mean_em_no_gm_CO-1-0', r'log$_{10}$ CO(1-0) / erg cm$^{-2}$ s$^{-1}$  , ${\Gamma_m}$ = 0')
    plot_map_standalone('f_mean_em_CO-3-2'      , r'log$_{10}$ CO(3-2) / erg cm$^{-2}$ s$^{-1}$')
    plot_map_standalone('f_mean_em_no_gm-CO-3-2', r'log$_{10}$ CO(3-2) / erg cm$^{-2}$ s$^{-1}$  , ${\Gamma_m}$ = 0')
    """
    
    fig = pylab.figure(figsize=(6,6))
    
    ax = fig.add_axes([0.10, 0.08, 0.8, 0.8])
    #ax_cbar = fig.add_axes([0.15, 0.8, 0.7, 0.03]) 
    
    #setting labels and ticklabels
    ax.set_xlabel('x(kpc)', size='x-large')
    ax.set_ylabel('y(kpc)', size='x-large')

    bs_min, bs_max = params['ranges']['box_size'].number     
    ax.set_xlim([bs_min, bs_max])
    ax.set_ylim([bs_min, bs_max])
    
    ticksv = numpy.linspace(bs_min, bs_max, 5)
    ax.set_xticks(ticksv)
    ax.set_yticks(ticksv)

    ax.set_xticklabels(ticksv, size='x-large')
    ax.set_yticklabels(ticksv, size='x-large')

    im_map = params['maps'][key] 

    im_data = log10(im_map['data']).T
    
    #clipping the map values outside the specified colorbar ranges
    im_data = numpy.clip(im_data, im_map['v_rng'][0], im_map['v_rng'][1])

    im = ax.imshow(im_data, 
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=im_map['v_rng'][0],     
                   vmax=im_map['v_rng'][1],     
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')
    
    cbar = pylab.colorbar(im, orientation = 'horizontal', format = '%.1f', shrink=0.8, 
                          ax = ax, 
                          ticks=numpy.linspace(im_map['v_rng'][0], im_map['v_rng'][1], 5))

    if title == None:
        title=''
        
    ax.set_title(title, size='x-large')
    fig.savefig('/home/mher/' + key + '.eps')
    
    return ax, cbar


def generate_maps(snap, params):
    #snapshot index
    snapIndex = snap
    
    suffix = '%06d' % snapIndex
    #path to fi snapshots
    snapName = 'fiout.%s' % suffix 
    snap_filename = params['rundir'] + '/firun/' + snapName 
    
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #loading the sph simulation data
    logger.debug('loading snapshot %s : ' % snap_filename)
    gas, dark, stars = read_set_from_file(snap_filename, format = fi_io.FiFileFormatProcessor)

    logger.debug('done reading fi snapshot : %s' % snap_filename)
    logger.debug('number of sph particles in snapshot = %d' %  len(gas))
    
    path = params['rundir'] + '/firun/runinfo'
    runinfo = parse_old_runinfo_file(path)

    snap_time = (float(runinfo['dtime'])*timeUnit.number) * (float(runinfo['noutbod']) * snap)

    #getting a new particle set for the gas particles with the attributes we will be using later
    #and converting the units to the same as those of the PDR models    
    gas = fi_utils.convert_units_to_pdr_units(gas, metallicity)
    logger.debug('converted and computed stuff from sph data (in units compatible with the pdr models)')
    
    print 'quantity          min              max  '
    print 'n_gas         %e     %e   '  % (gas.n.min()     , gas.n.max())
    print 'G0_gas        %e     %e   '  % (gas.G0.min()    , gas.G0.max())
    print 'g_mech        %e     %e   '  % (gas.gmech.min() , gas.gmech.max())
    print 'Av_gas        %e     %e   '  % (gas.Av.min()    , gas.Av.max())
    
    #keeping gas particles within the specified ranges
    gas = fi_utils.select_particles(gas, params['ranges'])
    logger.debug('got the sph particles in the required ranges')

    #setting the lower bound of the mechanical heating of the sph particles of the 
    #minimum of the pdr database which is almost negligable (since some of the 
    #sph particles have a zero gmech, their log would be -inf, so we clip it to the
    #mininmum of the pdr database
    gas.gmech = numpy.clip(gas.gmech, 10.0**lgMechMin, numpy.Infinity)

    #h = hist_nd(numpy.vstack((gas.x, gas.y)), nbins=100, plot=True)
    #h = hist_nd(numpy.vstack((log10(gas.n), log10(gas.G0))), nbins=50, plot=True, mn=-3.0, mx=3.0)
    #h = hist_nd(numpy.vstack((log10(gas.n), log10(gas.G0))), nbins=50, plot=True, mn=-3.0, mx=3.0)
    #h = hist_nd(numpy.vstack((log10(gas.n), log10(gas.Av))), nbins=100, plot=True, plot_log=True, mn=[-3.0, -2.0], mx=[4.0, 2.0]) 
    #########################done setting up the snapshot data######################    
    
    #computing alpha for each sph particle (amount of gm wrt surface heating for the PDR model)
    #and selecting particles which have alpha < 1.0
    fi_utils.compute_alpha(gas, arxvPDR)
    logger.debug('computed the alpha for the sph particles')
    inds = numpy.where( gas.alpha < 1.0 )
    gas = gas[inds]
     
    #h = hist_nd(numpy.vstack((log10(gas.n), log10(gas.T))), nbins=100, plot=True, plot_log=True, mn=[-3.0, 0.0], mx=[4.0, 5.0])
    #return gas #inds = xselect(log10(gas.n), log10(gas.T))
    #return h
    #asdasdad
    
    #setting the maximum Av of the sph particles to a predifined value (usually the emissions do not increase after some Av bec of the optical depth of the lines)
    gas.Av = numpy.clip(gas.Av, params['ranges']['sph']['Av_clip'][0], params['ranges']['sph']['Av_clip'][1])
    
    maps = params['maps']
    
    #computing the emissions of the SPH particles from the PDR interpolation functions
    print 'getting the emissions of each sph particle'    
    if params['use_em']:
        em_sph = fi_utils.gas_particle_all_emissions(gas, maps, arxvPDR, em_interp_funcs)

        #setting the emission info as attributes to the 'gas' particle set
        for key in em_sph:
            attr_name = key.replace('f_mean_','').replace('</q>','').replace('<q>','')
            setattr(gas, attr_name, em_sph[key])
    #
    
    #---------------------------------------------------------------------------------------
    # Done getting the distribution of mass, particle count as a function of number density
    #---------------------------------------------------------------------------------------
    print 'getting the spatial distrubutions'
    
    nBins = params['imres']    #mesh resolution use in binning
    
    for map_idx_key in maps:
        maps[map_idx_key]['data'] = numpy.zeros((nBins, nBins), dtype=numpy.float64) 
        maps[map_idx_key]['f_logn'] = None 
    
    x_sph, y_sph = gas.x, gas.y
    m_sph, G0_sph, n_sph, Av_sph, gmech_sph, T_sph = gas.mass, gas.G0, gas.n, gas.Av, gas.gmech, gas.T
     
    hist = hist_nd(numpy.vstack((x_sph, y_sph)), mn = bs_min, mx=bs_max, nbins=nBins, reverse_indicies=True, loc=True)
    hist.info()

    #looping over the bins of the 2D histogram of the x,y coordinates and computing the averages of the maps
    for i in numpy.arange(nBins):
        
        for j in numpy.arange(nBins):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            maps['f_n_part']['data'][i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:
                maps['f_mass']['data'][i,j]    = m_sph[inds_in_bin].sum()
                maps['f_mean_n']['data'][i,j]  = numpy.mean(n_sph[inds_in_bin])
                maps['f_mean_g0']['data'][i,j] = numpy.mean(G0_sph[inds_in_bin])
                maps['f_mean_gm']['data'][i,j] = numpy.mean(gmech_sph[inds_in_bin])
                maps['f_mean_Av']['data'][i,j] = numpy.mean(Av_sph[inds_in_bin])
                maps['T_mean']['data'][i,j]    = numpy.mean(T_sph[inds_in_bin])

                for key in em_sph:
                    if '_em_' in key:
                        maps[key]['data'][i,j] = numpy.mean(em_sph[key][inds_in_bin])    
    #
    #
    def fetch_grid_data(map_key, do_log10 = None):
        '''returns the data of a map given its key in parms['maps']. Basically this returns params['maps']['MAP_KEY]['data']'''
        if do_log10 != None and do_log10 == True:
            return log10(maps[map_key]['data'])
        else:
            return maps[map_key]['data']
    #
    #
    
    #displaying all the maps in a single plot    
    print 'done getting the spatial distributuions'
    
    fig, axs = pylab.subplots(4, 4, sharex=True, sharey=True, figsize=(12, 12), 
                              subplot_kw = {'xlim':[bs_min, bs_max],
                                            'ylim':[bs_min, bs_max],
                                            'aspect':'equal',
                                            'adjustable':'datalim',
                                           })
                              
    for ax in axs[:,0]: ax.set_ylabel('y(kpc)')
    for ax in axs[3,:]: ax.set_xlabel('x(kpc)')
    
    pylab.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.9, wspace=0.15, hspace=0.15)
    
    for map_key in maps:
        
        i, j = maps[map_key]['pos']
        
        im_map = maps[map_key]
        
        map_data = fetch_grid_data(map_key, do_log10 = im_map['log10']).T
        
        print 'idx = %d,%d quantity = %-30s max,min = %s,%s' % (im_map['pos'][0], im_map['pos'][1], map_key, map_data.min(), map_data.max()) 

        #clipping the map values outside the specified colorbar ranges
        map_data = numpy.clip(map_data, im_map['v_rng'][0], im_map['v_rng'][1])
        
        im = axs[i,j].imshow(map_data, 
                             extent=[bs_min, bs_max, bs_min, bs_max],
                             vmin=im_map['v_rng'][0], 
                             vmax=im_map['v_rng'][1], 
                             interpolation='bessel', #intepolation used for imshow
                             origin='lower')
        axs[i,j].set_title(im_map['title'], size='large')
        
        pylab.colorbar(im, ax=axs[i,j], orientation='vertical')
        
        axs[i,j].set_xlim([bs_min, bs_max])
        axs[i,j].set_ylim([bs_min, bs_max])
    #
    
    #plotting the time string
    pylab.figtext(0.05, 0.97, '%.3f %s' % (snap_time, timeUnit.unit.symbol), size='xx-large', color='k')
    image_fname =  params['rundir'] + '/analysis/maps.%s.%s' % (suffix, params['image_ext'])
    
    #saving the image
    if params['image_save'] == True:
        fig.savefig(image_fname)
    print 'wrote image file : %s' % image_fname

    #plotting the mesh of the histogram
    pylab.Figure()
    pylab.plot(hist.f.cntrd[0].flatten(), hist.f.cntrd[1].flatten(), '+')

    #saving the computed info
    if params['save_info'] == True: 
        fi_utils.save_gas_particle_info_saperate_files(snap_filename, gas, params['save_secies'])
        #gas_loaded = fi_utils.load_gas_particle_info(filename) 

    return gas, hist
#
#


t0 = time.time()
if pool_size == 1:
    for snap in snaps:
        gas, hist = generate_maps(snap, params)
else:
    pool = multiprocessing.Pool(pool_size)
    pool.map(generate_maps, snaps, params)
print 'total time = %s s' % (time.time() - t0)




pylab.show()
    