# in this test file, we load a database of PDR meshes and extract information
# at conditions determined by the SPH simulations

#  - exclude sph particles which have densities below which the LVG model codes can handel 
#    (for example, one with n_gas < 1e-3))
#  - exlcude particle with alpha > 1 (unless for very high densities (get the densitiy of 
#    the specie whose  emission is to be computed and see if its abundance is too low, if 
#    if is too low, just set the emission of the sph particle corresponding to that specie
#    to zero.
#  - abundances of CO for PDR modesls with n < 0.1 cm^-3 are ~ 1e-10 (which is too low compared to the
#    high density modesl

import time
import sys
import os
import subprocess
import multiprocessing

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate
import numpy
from numpy import log10
import pylab
import logging

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants, nbody_system
import meshUtils
from mylib.utils.misc  import xselect, default_logger
from mylib.utils.histogram import hist_nd
from fi_utils import make_maps, parse_old_runinfo_file
import fi_utils
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

pool_size = 1
snaps = numpy.arange(4, 4+1, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
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
                      'box_size' : [-2.0, 2.0],

                       #modes in the PDR arxv which are within those ranges will be used in constructing interpolation functions 
                      'interp'   : {'log_n'    : [-10.0,  10.0 ],
                                    'log_G0'   : [-2.0 ,  3.0 ],
                                    'log_gmech': [-50.0, -20.0],
                                    'Av'       : [ 3.0 ,  30.0],
                                    },                          
                      },
          'maps'   : {
                      'f_n_part' : {'pos': [0,0], 'title': r'$f(N)$'        , 'v_rng': [0,6]      , 'log10': True}, 
                      'f_mass'   : {'pos': [0,1], 'title': r'$f(m)$'        , 'v_rng': [30,40]    , 'log10': True},
                      'f_mean_n' : {'pos': [0,2], 'title': r'$f(\bar{n})$'  , 'v_rng': [-3.0, 4.0], 'log10': True},
                      'f_mean_g0': {'pos': [0,3], 'title': r'$f(\bar{g_0})$', 'v_rng': [-3.0, 3.0], 'log10': True},
                      #--------
                      'f_mean_gm'              : {'pos': [1,0], 'title': r'$f(\bar{\Gamma_m})$'                  , 'v_rng': [-35.0, -22.0], 'log10': True},
                      #'f_mean_Av'              : {'pos': [1,1], 'title': r'$f(\bar{Av})$'                      , 'v_rng': [0.0  , 2.0] ,  'log10': True},
                      'T_mean'                 : {'pos': [1,1], 'title': r'$f(\bar{T})$'                        , 'v_rng': [0.0  , 5.0] ,  'log10': True},
                      #'f_mean_em_C+-1-0'       : {'pos': [1,1], 'title': r'$f(L_{C^+ 158 \mu m})$'               , 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_no_gm_C+-1-0' : {'pos': [1,2], 'title': r'$f(L_{C^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      'f_mean_em_CO-7-6'       : {'pos': [1,3], 'title': r'$f(L_{O^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #--------
                      'f_mean_em_CO-1-0'       : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'                  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      'f_mean_em_no_gm_CO-1-0' : {'pos': [2,1], 'title': r'$f(L_{CO(1-0)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      'f_mean_em_CO-3-2'       : {'pos': [2,2], 'title': r'$f(L_{CO(3-2})$'                  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      'f_mean_em_no_gm_CO-3-2' : {'pos': [2,3], 'title': r'$f(L_{CO(3-2)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #--------
                      'f_mean_em_HCN-1-0'       : {'pos': [3,0], 'title': r'$f(L_{HCN(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
                      },

          #'interpolator' : interpolate.NearestNDInterpolator, 
          'interpolator' : interpolate.LinearNDInterpolator, 
          'image_save': False,
          'image_ext' : 'eps',
          }

#those are the lines which will be computed
lines_info = { #line str  transIdx    trans latex str               pdr|lvg
              'C+-1-0'  :  {          'latex': r'C$^+$ 158 $\mum$', 'type':'pdr', 'quantity':['fineStructureCoolingComponents','C+','rate','1-0']},
              'CO-1-0'  :  {'idx':0,  'latex': r'CO(1-0)'         , 'type':'lvg', },
              'CO-3-2'  :  {'idx':2,  'latex': r'CO(3-2)'         , 'type':'lvg', },
              'CO-7-6'  :  {'idx':6,  'latex': r'CO(7-6)'         , 'type':'lvg', },
              'CO-16-15':  {'idx':15, 'latex': r'CO(7-6)'         , 'type':'lvg', },
                          
              'HCN-1-0':   {'idx':0,  'latex': r'HCN(1-0)'        , 'type':'lvg', },
              'HCN-3-2':   {'idx':2,  'latex': r'HCN(3-2)'        , 'type':'lvg', },
 
              'HNC-1-0':   {'idx':0,  'latex': r'HNC(1-0)'        , 'type':'lvg', },
              'HNC-3-2':   {'idx':2,  'latex': r'HNC(3-2)'        , 'type':'lvg', },

              'HCO+-1-0':  {'idx':0,  'latex': r'HCO$^{+}$(1-0)'  , 'type':'lvg', },
              'HCO+-3-2':  {'idx':2,  'latex': r'HCO$^{+}$(3-2)'  , 'type':'lvg', },
             }

#extracting/guessing the metallicity from the name of the directory of the run
metallicity = fi_utils.guess_metallicity(params['rundir'])

logger = default_logger()

#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################


#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

lgMechMin = arxvPDR.grid_z.min()     
bs_min, bs_max = params['ranges']['box_size'] 

if params['use_em'] == True:    
    #=============================================
    #loading the molecular specie radex database
    #arxvPDR.readDbsRadex(allDbs=True)
    #arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'], Av=params['ranges']['interp']['Av'])
    arxvPDR.readDbsRadex(species=['CO','HCN', 'HCO+'], in_Av_rng=params['ranges']['interp']['Av'])

def funcRadex(mesh_radex, **kwargs):
    
    if mesh_radex == None:
        return numpy.nan
    else:
        return log10(mesh_radex[kwargs['transitionIdx']]['fluxcgs'])  #a radex quantity

def funcPDR(meshObj, **kwargs):
    
    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = meshObj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    
    return log10(value)

def interpolator_sectioned_by_gm_Av(data_db, v):
    """Returns an interpolator function which makes use of sub-interpolation function for ranges in 
       mechanical heating. data_db is an array whose columns are : n, G0, gm, Av.
    """
    
    #sections in gm within which saperate interpolation functions will be built.
    # .. note:: make sure there are enough points in the database within each section
    intervals_gm = [ 
                     [-50.0, -35.0],  
                     [-35.0, -34.0], [-34.0, -33.0], [-33.0, -32.0], [-32.0, -31.0], [-31.0, -30.0],    
                     [-30.0, -29.0], [-29.0, -28.0], [-28.0, -27.0], [-27.0, -26.0], [-26.0, -25.0],    
                     [-25.0, -24.0], [-24.0, -23.0], [-23.0, -22.0], [-22.0, -21.0], [-21.0, -20.0],    
                     #[-20.0, -13.0],
                   ]
    ghost_gm = 0.0 #1.1
    
    #sections in Av within which saperate interpolation functions will be built.
    # .. note:: make sure there are enough points in the database within each section
    intervals_Av = [ 
                     [3.0, 4.0],   [4.0, 5.0],   [5.0, 6.0],   [6.0, 7.0],   [7.0, 8.0],   [8.0, 9.0], [9.0, 10.0],
                     [10.0, 12.0], [12.0, 14.0], [14.0, 16.0], [16.0, 18.0], [18.0, 20.0], 
                     [20.0, 22.0], [22.0, 24.0], [24.0, 26.0], [26.0, 28.0], [28.0, 30.0],
                   ] 

    ghost_Av = 0.0 #1.1
    
    #splitting data into sections specified by intervals_gm and constructing an interpolation
    #function for each interval (the intevals include an extra dex in each bound of the specified
    #bound to make sure interpolation is done correctly...i.e the scipy interpolator doesnt return
    # nan's)
    
    data_db_gm, data_db_Av = data_db[:, 2], data_db[:, 3]

    #numpy 2D array whose dimesnions are the number of intervals in gm and Av (in each dimensions)
    #for instance all_fInterps[0,3] corresponds to an interpolation function within the range 
    #intervals_gm[0] and intervals_Av[3]    
    all_fInterps = numpy.ndarray((len(intervals_gm), len(intervals_Av)), dtype='object')
    
    #lookup points which lie in each section of gm and Av and get the points in there
    #and construct an interpolation function using those points
    for i, interval_gm in enumerate(intervals_gm):
        
        if interval_gm[0] < params['ranges']['interp']['log_gmech'][0] or interval_gm[1] > params['ranges']['interp']['log_gmech'][1]:
            raise ValueError("specified interval [%.2f, %.2f] for gmech is outside the range specified in params['ranges']['interp']" % (interval_gm[0], interval_gm[1],)) 
        
        interval_gm_low  = interval_gm[0] - ghost_gm
        interval_gm_high = interval_gm[1] + ghost_gm        
        
        for j, interval_Av in enumerate(intervals_Av):

            if interval_Av[0] < params['ranges']['interp']['Av'][0] or interval_Av[1] > params['ranges']['interp']['Av'][1]:
                raise ValueError("specified interval [%.2f, %.2f] for Av is outside the range specified in params['ranges']['interp']" % (interval_Av[0],interval_Av[1],)) 
        
            interval_Av_low  = interval_Av[0] - ghost_Av
            interval_Av_high = interval_Av[1] + ghost_Av

            #finding the points in the DB in this section
            inds_in_this_interval = numpy.where(
                                               (data_db_gm >= interval_gm_low)*(data_db_gm <= interval_gm_high)*
                                               (data_db_Av >= interval_Av_low)*(data_db_Av <= interval_Av_high)
                                              )[0]
                               
            if inds_in_this_interval.size != 0:
                
                #if enough points are found, make the interpolation funcion
                data_db_this_box = data_db[inds_in_this_interval, :]
                v_in_this_box = v[inds_in_this_interval]
                
                print '(i,j) = (%02d, %02d) %06d DB points in box gm = [%.2f, %.2f ], Av = [%04.1f,%04.1f ] ' % (i, j, inds_in_this_interval.size, interval_gm_low, interval_gm_high, interval_Av_low, interval_Av_high),
                t0 = time.time()
                fInterp_this_interval = params['interpolator'](data_db_this_box, v_in_this_box)
                print 'interpolation function constructed in %.2f sec' % (time.time()-t0)
                all_fInterps[i,j] = fInterp_this_interval
            
            else:                
                print 'no points found in box gm = [%f, %f ], Av = [%f,%f ]' % (interval_gm_low, interval_gm_high, interval_Av_low, interval_Av_high)

            
    def fInterp(data):
        """devides data (values at which we would like to get interpolated values) into sections of 
        gm and Av (3rd and 4th colomn) and returns the interpolated values by passing them to the 
        sub-interpolation function in all_fInterps.
        """
        
        nPts = data.shape[0]
        
        #array which will hold the interpolated values from all the sections
        all_values = numpy.zeros(nPts, dtype=numpy.float64)
        
        data_gm, data_Av = data[:, 2], data[:, 3]
        inds_orig = numpy.arange(nPts) 
        
        #interpolating from the sectioned interpolation function. Lookup the points
        #which are in the bounds of each interpolation function and interpolate the
        #values and set them in the array to be returned.
        for i, interval_gm in enumerate(intervals_gm):
            
            for j, interval_Av in enumerate(intervals_Av):

            
                inds_in_box = numpy.where( 
                                           (data_gm >= interval_gm[0])*
                                           (data_gm <  interval_gm[1])* 
                                           (data_Av >= interval_Av[0])*
                                           (data_Av <  interval_Av[1]) 
                                         )[0]
                                                                                
                if inds_in_box.size == 0:
                    continue
                
                data_in_this_interval = data[inds_in_box, :]

                print '(i,j) = (%02d, %02d) %06d SPH points in box gm = [%.2f, %.2f ], Av = [%4.1f,%4.1f ] ' % (i, j, inds_in_box.size, interval_gm[0], interval_gm[1], interval_Av[0], interval_Av[1]),                
                t0 = time.time()
                v = all_fInterps[i,j](data_in_this_interval)
                #v = 0.0
                dt = time.time()-t0
            
                #setting the values in the array to be returned
                all_values[inds_orig[inds_in_box]] = v
                
                print 'interpolated in %5.2f sec, %05d points at %e points/sec' % (dt, inds_in_box.size, inds_in_box.size/dt)

        return all_values
    #
    
    return fInterp

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
        fInterp = interpolator_sectioned_by_gm_Av(data, v)
    else:
        fInterp = params['interpolator'](data, v)
        
    logger.debug('time contruncting interpolator for quantity %s %s in %.4f seconds.' % (kwargs['quantity'][1], kwargs['quantity'][3], time.time() - t0))

    return fInterp

def get_interpolation_function_radex(**kwargs):
    """makes an interpolation function from the radex meshes from the database."""
    
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
        fInterp = interpolator_sectioned_by_gm_Av(data_all_Av, v_all_Av)
    else:
        fInterp = params['interpolator'](data_all_Av, v_all_Av)
    logger.debug('time contruncting interpolation function for %s transition %d in %.4f seconds' % (specStr, kwargs['transitionIdx'], time.time() - t0))
    return fInterp

def make_emission_interp_funcs():
    '''returns a dictionary of functions with keys XYZ-UPPER-LOWER. Items in 'maps' which have
    _em_ are parsed and the transition is looked up (f_mean_em_C+-1-0 is parsed into C+-1-0) 
    and an interpolation function for C+-1-0 is computed.
    '''
        
    funcs = {}
    
    for map_key in params['maps']:
        
        if '_em_' in map_key:
            line = map_key.split('_')[-1] # the line XYZ-UPPER-LOWER
            
            specStr= line.split('-')[0]   # XYZ
            
            if lines_info[line]['type'] == 'lvg':
                funcs[line] = get_interpolation_function_radex(
                                                               specStr  = specStr,
                                                               transitionIdx = lines_info[line]['idx'],
                                                               sectioned = True,
                                                              )
            '''
            if lines_info[line]['type'] == 'pdr':
                funcs[line] = get_interpolation_function_pdr(
                                                             quantity = lines_info[line]['quantity'],
                                                             sectioned = True,
                                                             )
            '''
            
    return funcs

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

#constructing all the interpolation functions
if params['use_em']:
    logger.debug('getting interpolation functions of the emissions')
    funcs = make_emission_interp_funcs()
    logger.debug('done getting interpolation functions of the emissions')


def plot_map_standalone(key, title):
    
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

#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################
#############################done setting   up the PDR databaseses and the interpolation function##################

def generate_map(snap):
    #snapshot index
    snapIndex = snap
    
    suffix = '%06d' % snapIndex
    #path to fi snapshots
    snapName = 'fiout.%s' % suffix 
    filename = params['rundir'] + '/firun/' + snapName 
        
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #########################setting up the snapshot data######################
    #loading the sph simulation data
    logger.debug('loading snapshot %s : ' % filename)
    gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)

    logger.debug('done reading fi snapshot : %s' % filename)
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
    
    print 'getting the emissions of each sph particle'
    #getting emissions of the sph particles with mechanical heating    
    em_sph = {}
    data_with_gm  =  numpy.array([log10(gas.n), log10(gas.G0), log10(gas.gmech), gas.Av]).T
    data_no_gm    =  numpy.array([log10(gas.n), log10(gas.G0), numpy.ones(len(gas))*lgMechMin, gas.Av]).T

    time_interp = 0.0
    
    #getting the emission of the sph particles from the PDR models
    if params['use_em']:
        
        for map_key in maps:
            
            if '_em_' in map_key:
                
                #what points to use? with or witout gm
                if '_no_gm_' in map_key:
                    data_use = data_no_gm
                else:
                    data_use = data_with_gm
    
                line = map_key.split('_')[-1] # the line XYZ-UPPER-LOWER

                logger.debug('interpolating emissions of %s of the sph particles for %d particles' % (map_key, data_use.shape[0]) )
                
                t0 = time.time()
                em_sph[map_key] = 10.0**funcs[line](data_use)
                dt = time.time() - t0
                
                logger.debug('time interpolating emission for map %s = %.3f seconds for %d points at %e points/sec' % (map_key, dt, data_use.shape[0], data_use.shape[0]/dt) )

                #setting emissions of particles with nan interpolated values to zero             
                indsNan = numpy.where(numpy.isnan(em_sph[map_key]))[0]
                if indsNan.size != 0:
                    em_sph[map_key][indsNan] = 1e-30
                    
                gas.CO10 = em_sph[map_key]
                
                time_interp += dt
            
    logger.debug('done getting the emissions of each sph particle in %.2f seconds' % time_interp)
    
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
     
    hist = hist_nd(numpy.vstack((x_sph, y_sph)), mn = bs_min, mx=bs_max, nbins=nBins, reverse_indicies=True)
    hist.info()
     
    for i in numpy.arange(nBins):
        
        for j in numpy.arange(nBins):
            
            inds_in_bin = hist.get_indicies([i,j])
    
            maps['f_n_part']['data'][i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:
                maps['f_mass']['data'][i,j]    = m_sph[inds_in_bin].sum()
                maps['f_mean_n']['data'][i,j]  = numpy.mean(n_sph[inds_in_bin])
                maps['f_mean_g0']['data'][i,j] = numpy.mean(G0_sph[inds_in_bin])
                maps['f_mean_gm']['data'][i,j] = numpy.mean(gmech_sph[inds_in_bin])
                #maps['f_mean_Av']['data'][i,j] = numpy.mean(Av_sph[inds_in_bin])
                maps['T_mean']['data'][i,j]     = numpy.mean(T_sph[inds_in_bin])

                for key in em_sph:
                    if '_em_' in key:
                        maps[key]['data'][i,j] = numpy.mean(em_sph[key][inds_in_bin])
    #
    #
    def get_histogram(x, y, nb, mn, mx):
        '''
        given a map key, it generates and the distribution of y as a function of the quantity x
        '''

        #computing the histogram
        hist = hist_nd(x.reshape((1, x.shape[0])),  
                       loc=True,
                       nbins = nb, 
                       mn = mn, 
                       mx = mx, 
                       reverse_indicies = True)

        dist = numpy.zeros((nb), dtype=numpy.float64) 
        
        for i in numpy.arange(nb):
            inds_in_bin = hist.get_indicies([i])
        
            if inds_in_bin.size > 0:
                dist[i] = y[inds_in_bin].sum()
        
        hist.f[:] = dist[:]
        return hist
    #
    #    

    """
    #getting the distributions as a functio of log_n (the gas density)
    for map_idx_key in maps:
        if map_idx_key == 'f_n_part':
            maps['f_n_part']['f_logn'] = get_histogram(log10(n1), numpy.ones(n1.shape), 100.0, -3.0, 3.0)
        if map_idx_key == 'f_mass':
            maps['f_mass']['f_logn'] = get_histogram(log10(n1), m1, 100.0, -3.0, 3.0) 
        if map_idx_key == 'f_mean_n':
            maps['f_mean_n']['f_logn'] = get_histogram(log10(n1), n1, 100.0, -3.0, 3.0)
        if map_idx_key == 'f_mean_g0':
            maps['f_mean_g0']['f_logn'] = get_histogram(log10(n1), g01, 100.0, -3.0, 3.0) 
        if map_idx_key == 'f_mean_gm':
            maps['f_mean_gm']['f_logn'] = get_histogram(log10(n1), g_mech1, 100.0, -3.0, 3.0)
        for key in em_sph:
            if '_em_' in key:
                maps[key]['f_logn'] = get_histogram(log10(n1), em_sph[key], 100.0, -3.0, 3.0)
    #
    """
    
    #
    #
    def fetch_grid_data(map_key, do_log10 = None):
        
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
        
    pylab.figtext(0.05, 0.97, '%.3f %s' % (snap_time, timeUnit.unit.symbol), size='xx-large', color='k')
    image_fname =  params['rundir'] + '/analysis/maps.%s.%s' % (suffix, params['image_ext'])
    
    if params['image_save'] == True:
        fig.savefig(image_fname)
    print 'wrote image file : %s' % image_fname


    #displaying the distributions as a function of n
    print 'done getting the spatial distributuions'
    
    """
    fig, axs = pylab.subplots(4, 4, sharex=True, sharey=False, figsize=(16,16), 
                              subplot_kw = {'xlim':[-3.0, 3.0],
                                            'ylim':[-5.0, 0.5],
                                            'aspect':'equal',
                                            'adjustable':'datalim',
                                           }
                               )
    for ax in axs[:,0]: ax.set_ylabel(r'$\log_{10} f / f_{max}$')
    for ax in axs[3,:]: ax.set_xlabel(r'$\log_{10} n$')
    
    pylab.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.95, wspace=0.3, hspace=0.3)
    
    for map_key in maps:
        
        i, j = maps[map_key]['pos']
    
        im_map = maps[map_key] 
        
        xPlt = maps[map_key]['f_logn'].f.cntrd
        yPlt = log10(maps[map_key]['f_logn'].f/maps[map_key]['f_logn'].f.max())
        axs[i,j].plot(xPlt, yPlt)
        
        axs[i,j].set_title(im_map['title'], size='large')
    
    pylab.figtext(0.05, 0.97, '%.3f %s' % (snap_time, timeUnit.unit.symbol), size='xx-large')
    image_fname =  params['rundir'] + '/analysis/maps-logn.%s.%s' % (suffix, params['image_ext'])

    if params['image_save'] == True:
        fig.savefig(image_fname)
    print 'wrote image file : %s' % image_fname
    """
        
    #subprocess.Popen(['xv', image_fname])
    pylab.show()
    
    return gas
#
#

t0 = time.time()
if pool_size == 1:
    for snap in snaps:
        gas = generate_map(snap)
else:
    pool = multiprocessing.Pool(pool_size)
    pool.map(generate_map, snaps)
print 'total time = %s s' % (time.time() - t0)


def test_plot():
    
    '''
    nPts = 10000.0
    n_arr  = log10(numpy.ones(nPts)*1e2)
    G0_arr = log10(numpy.ones(nPts)*1e2)
    #gm_arr = log10(numpy.ones(nPts)*1e-24)
    gm_arr = numpy.linspace(-30.0, -16, nPts)
    #Av_arr = numpy.linspace(5.0, 10.0, nPts)
    Av_arr = numpy.ones(nPts)*5.0
    '''
    
    inds = numpy.where( 
                        (log10(gas.gmech) >= -26.00)*(log10(gas.gmech) <= -24.00)*
                        (gas.Av >= 4.00)*(gas.Av <= 5.0)
                      )  

    #inds = [ 63, 64,  65, 184]
    
    n_arr  = log10(gas[inds].n)
    G0_arr = log10(gas[inds].G0)
    gm_arr = log10(gas[inds].gmech)
    Av_arr = gas[inds].Av
    #Av_arr = gas[inds].Av.clip(0.0, 9.9)

    data  =  numpy.array([n_arr, G0_arr, gm_arr, Av_arr,]).T
    
    print data.shape
    vi = funcs['CO-1-0'](data)
    
    print numpy.where(numpy.isnan(vi))[0].size
    """

    pylab.plot(gm_arr, vi, '.')
    pylab.show()
    """


    