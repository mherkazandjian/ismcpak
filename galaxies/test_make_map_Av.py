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
from numpy import log10

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants, nbody_system
import meshUtils
from mylib.utils.misc  import xselect
from mylib.utils.histogram import hist_nd
from fi_utils import make_maps, parse_old_runinfo_file
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

pool_size = 1
snaps = numpy.arange(4, 4+1, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',   # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',   # the path of the dir containing the simulation
          'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',        # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2/',        # the path to the dir containing the PDR database          
          'use_em' : False, 
           
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                             'Av'             :  [0.0, 19],
                            },
                      
                      #ranges of the pdr models to be used in constructing the interpolation functions
                      'pdr': {
                              'log_n' : [0,0],
                             },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-4,4] | units.kpc,
                      'Av'       : [3.0, 20.0],  #range in Av over which the interpolation function will be constructed (for radex, Dbs with an increment of 1 will be queried)
                      },          
          'maps'   : {
                      'f_n_part' : {'pos': [0,0], 'title': r'$f(N)$'        , 'v_rng': [0,6]      , 'log10': True}, 
                      'f_mass'   : {'pos': [0,1], 'title': r'$f(m)$'        , 'v_rng': [30,40]    , 'log10': True},
                      'f_mean_n' : {'pos': [0,2], 'title': r'$f(\bar{n})$'  , 'v_rng': [-3.0, 4.0], 'log10': True},
                      'f_mean_g0': {'pos': [0,3], 'title': r'$f(\bar{g_0})$', 'v_rng': [-3.0, 3.0], 'log10': True},
                      #--------
                      'f_mean_gm'              : {'pos': [1,0], 'title': r'$f(\bar{\Gamma_m})$'                  , 'v_rng': [-35.0, -22.0], 'log10': True},
                      'f_mean_Av'              : {'pos': [1,1], 'title': r'$f(\bar{Av_m})$'                      , 'v_rng': [0.0  , 2.0] , 'log10': False},
                      #'f_mean_em_C+-1-0'       : {'pos': [1,1], 'title': r'$f(L_{C^+ 158 \mu m})$'               , 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_no_gm_C+-1-0' : {'pos': [1,2], 'title': r'$f(L_{C^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #'f_mean_em_CO-7-6'       : {'pos': [1,3], 'title': r'$f(L_{O^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      #--------
                      'f_mean_em_CO-1-0'       : {'pos': [2,0], 'title': r'$f(L_{CO(1-0})$'                  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #'f_mean_em_no_gm_CO-1-0' : {'pos': [2,1], 'title': r'$f(L_{CO(1-0)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #'f_mean_em_CO-3-2'       : {'pos': [2,2], 'title': r'$f(L_{CO(3-2})$'                  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #'f_mean_em_no_gm_CO-3-2' : {'pos': [2,3], 'title': r'$f(L_{CO(3-2)})$ $\Gamma_m = 0$'  , 'v_rng': [-10.0, -2.0], 'log10': True},
                      #--------
                      #'f_mean_em_HCN-1-0'       : {'pos': [3,0], 'title': r'$f(L_{HCN(1-0})$'                 , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCN-1-0' : {'pos': [3,1], 'title': r'$f(L_{HCN(1-0)})$ $\Gamma_m = 0$' , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_HCO+-1-0'      : {'pos': [3,2], 'title': r'$f(L_{HCO+(1-0})$'                , 'v_rng': [-10.0, -7.0], 'log10': True},
                      #'f_mean_em_no_gm_HCO+-1-0': {'pos': [3,3], 'title': r'$f(L_{HCO+(1-0)})$ $\Gamma_m = 0$', 'v_rng': [-10.0, -7.0], 'log10': True},
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

lnMin, lG0Min, lgMechMin = params['ranges']['sph']['min_log_n_use'], params['ranges']['sph']['min_log_G0_use'], params['ranges']['sph']['min_log_gm_use']      
bsMin, bsMax = params['ranges']['box_size'].number 


#extracting/guessing the metallicity from the name of the directory of the run
metallicity = None

if '-std' in params['rundir']:
    metallicity = 0.2
if '-sol' in params['rundir']:
    metallicity = 1.0

#------------------------------------------

def setupLogger():
    """sets up the logger which will prepend info about the printed stuff. Assignes a value to self.logger."""
    
    # setting up the logger                                                                                                                                                                                                              
    # create logger                                                                                                                                                                                                                      
    logger = logging.getLogger('simple_example')
    if not len(logger.handlers):
        logger.setLevel(logging.DEBUG)

        # create console handler and set level to debug                                                                                                                                                                                      
        ch = logging.StreamHandler( sys.stdout )  # setting the stream to stdout                                                                                                                                                             
        ch.setLevel(logging.DEBUG)

        # create formatter                                                                                                                                                                                                                   
        formatter = logging.Formatter('[%(asctime)s %(funcName)s() %(filename)s:%(lineno)s] %(message)s') # this was the original in the example                                                                              

        # add formatter to ch                                                                                                                                                                                                                
        ch.setFormatter(formatter)

        # add ch to logger                                                                                                                                                                                                                   
        logger.addHandler(ch)

    return logger

logger = setupLogger()

#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################
#############################setting up the PDR databaseses and the interpolation function##################


#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)
    
if params['use_em'] == True:    
    #=============================================
    #loading the molecular specie radex database
    arxvPDR.readDbsRadex(allDbs=True)
    #arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'], Av=params['ranges']['Av'])

def funcRadex(mesh_radex, **kwargs):
    
    if mesh_radex == None:
        return numpy.nan
    else:
        return numpy.log10(mesh_radex[kwargs['transitionIdx']]['fluxcgs'])  #a radex quantity

def funcPDR(meshObj, **kwargs):
    
    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = meshObj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    value /= (2.0*numpy.pi)
    
    return numpy.log10(value)

def interpolator_sectioned_by_gm_Av(data_db, v):
    """Returns an interpolator function which makes use of sub-interpolation function for ranges in 
       mechanical heating. data_db is an array whose columns are : n, G0, gm, Av.
    """
    
    intervals_gm = [ [-50.0, -40.0],  #sections in gm within which saperate interpolation functions will be built.
                     [-45.0, -35.0],  # .. note:: make sure there are enough points in the database within each section
                     [-35.0, -32.0],
                     [-32.0, -30.0],
                     [-30.0, -28.0],
                     [-28.0, -26.0],
                     [-26.0, -24.0],
                     [-24.0, -22.0],
                     [-22.0, -20.0],
                     [-20.0, -18.0],
                     [-18.0, -14.0],
                   ]
    ghost_gm = 1.1
    
    intervals_Av = [ [3.0, 4.0],  #sections in Av within which saperate interpolation functions will be built.
                     [4.0, 5.0],  # .. note:: make sure there are enough points in the database within each section
                     [5.0, 6.0],
                     [6.0, 7.0],
                     [7.0, 8.0],
                     [8.0, 9.0],
                     [9.0, 10.0],
                   ]
    ghost_Av = 1.0
    
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
        
        interval_gm_low  = interval_gm[0] - ghost_gm
        interval_gm_high = interval_gm[1] + ghost_gm

        for j, interval_Av in enumerate(intervals_Av):
        
            interval_Av_low  = interval_Av[0] - ghost_Av
            interval_Av_high = interval_Av[1] + ghost_Av

            #finding the points in the DB in this section
            inds_in_this_box = numpy.where( 
                                            (data_db_gm >= interval_gm_low)*(data_db_gm <= interval_gm_high)*
                                            (data_db_Av >= interval_Av_low)*(data_db_Av <= interval_Av_high)
                                          )[0]
                               
            if inds_in_this_box.size != 0:
                
                #if enough points are found, make the interpolation funcion
                data_db_this_box = data_db[inds_in_this_box, :]
                v_in_this_box = v[inds_in_this_box]
                  
                print '(i,j) = (%02d, %02d) %06d DB points in box gm = [%.2f, %.2f ], Av = [%04.1f,%04.1f ] ' % (i, j, inds_in_this_box.size, interval_gm_low, interval_gm_high, interval_Av_low, interval_Av_high),
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

                print '(i,j) = (%02d, %02d) %06d SPH points in box gm = [%.2f, %.2f ], Av = [%4.1f,%4.1f ] ' % (i, j, inds_in_box.size, interval_gm_low, interval_gm_high, interval_Av_low, interval_Av_high),                
                t0 = time.time()
                #v = all_fInterps[i,j](data_in_this_interval)
                v = 0.0
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
    Av_range = kwargs.pop('Av_range')
    
    #collecting the data from the database corresponding to all the data
    #in Av avaiable
    for i, Av in enumerate(numpy.arange(Av_range[0], Av_range[1])):
        
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
                                                               Av_range = params['ranges']['Av'],
                                                               transitionIdx = lines_info[line]['idx'],
                                                               sectioned = True,
                                                              )
            '''
            if lines_info[line]['type'] == 'pdr':
                funcs[line] = get_interpolation_function_pdr(
                                                             quantity = lines_info[line]['quantity'],
                                                             Av_range = params['ranges']['Av'],
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

    bsMin, bsMax = params['ranges']['box_size'].number     
    ax.set_xlim([bsMin, bsMax])
    ax.set_ylim([bsMin, bsMax])
    
    ticksv = numpy.linspace(bsMin, bsMax, 5)
    ax.set_xticks(ticksv)
    ax.set_yticks(ticksv)

    ax.set_xticklabels(ticksv, size='x-large')
    ax.set_yticklabels(ticksv, size='x-large')

    im_map = params['maps'][key] 

    im_data = numpy.log10(im_map['data']).T
    
    #clipping the map values outside the specified colorbar ranges
    im_data = numpy.clip(im_data, im_map['v_rng'][0], im_map['v_rng'][1])

    im = ax.imshow(im_data, 
              extent=[bsMin, bsMax, bsMin, bsMax],
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

xx = {} #.. note:: just a temporary variable

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
    
    #plotting the temperature vs the gas density
    meanmwt=1.3|units.amu
    n_gas_cgs = (gas.rho/meanmwt).value_in( units.cm **-3 )
    gasT = gas.temperat.value_in(units.K)
    G0 = 6.54*gas.fuvheat.value_in( units.none )
    g_mech = gas.dethdt.as_quantity_in( units.erg / (units.g * units.s ) )
    # converting the mechanical heating rate from per unit mass (of H gas)     
    # to per unit volume (see notesISM.odt)
    g_mech = g_mech.value_in( g_mech.unit ) # gMech in erg / (g s )
    g_mech = g_mech * 1.6474e-24 * n_gas_cgs # gMech in erg / (cm^3 s) 
    m_cgs = gas.mass.value_in(units.g)
    
    #computing the mean Av of the sph particles using Inti's recipe
    def Av_sph_particles():
        Tk = gasT #sph gas particle temperature in K
        sigma_v = gas.vdisp.value_in(units.km / units.s)  #sph particle velocity despertion in km/s
        
        Pe = (1.085 * gasT + 54.0 * sigma_v**2.0)*n_gas_cgs
                
        Av_mean = 0.22 * metallicity * ( 1520.0 / 100.0) * numpy.sqrt(Pe/1e4)
        
        return Av_mean
    #
    Av_mean = Av_sph_particles()

    xkpc, ykpc, zkpc = gas.position.value_in(params['ranges']['box_size'].unit).T
    
    #h = hist_nd(numpy.vstack((log10(n_gas_cgs), log10(G0)), nbins=10, plot=True, mn=-3.0, mx=3.0)
    #h = hist_nd(numpy.vstack((log10(xx['n']), xx['av'])), nbins=10, plot=True)
    #h = hist_nd(numpy.vstack((log10(xx['n']), xx['av'])), nbins=100, plot=True, mn=[0,0], mx=[4,30], w = numpy.ones(xx['n'].shape,'f'), wfunc = lambda x: numpy.log10(numpy.sum(x)))

    xx['av'] = Av_mean
    xx['n']  = n_gas_cgs
    xx['g0'] = G0
    xx['x']  = xkpc
    asdasdasd
    
    #using  a minimum Av according to the speciecied range
    Av_mean = Av_mean.clip(params['ranges']['sph']['Av'][0], params['ranges']['sph']['Av'][1])
    xx['Av_parts'] = Av_mean # ;;;
    logger.debug('compute the mean Av of the sph particles')

    logger.debug('converted and computed stuff from sph data (in units compatible with the pdr models)')
    
    """
    pylab.plot(n_gas_cgs, Av_mean,'.')
    pylab.show()
    asdasd
    """
    #########################done setting up the snapshot data######################
    #########################done setting up the snapshot data######################
    #########################done setting up the snapshot data######################
    #########################done setting up the snapshot data######################
    
    
    #setting the lower bound of the mechanical heating of the sph particles of the 
    #minimum of the pdr database which is almost negligable (since some of the 
    #sph particles have a zero gmech, their log would be -inf, so we clip it to the
    #mininmum of the pdr database
    g_mech = g_mech.clip(10.0**lgMechMin)

    ##
    def select_sph_particles(ranges):
        '''
        returns the sph particle indicies which are of interest and can be used to get stuff from the 
        PDR database.
        
        :param ranges: the ranges which will be used to select the particles
        '''

        #selecting particles within a spatial range and within the 
        #density and g0 range of the database
        inds = numpy.where(
                           (numpy.log10(n_gas_cgs) >= ranges['sph']['min_log_n_use'])*
                           (numpy.log10(G0) >= ranges['sph']['min_log_G0_use'])*
                           (numpy.log10(g_mech) >= ranges['sph']['min_log_gm_use'])*                   
                           (xkpc > ranges['box_size'][0].number)*(xkpc < ranges['box_size'][1].number)*
                           (ykpc > ranges['box_size'][0].number)*(ykpc < ranges['box_size'][1].number)
                          )[0]
                          
        return inds
    #
    #
        
    inds = select_sph_particles(params['ranges'])
    x1, y1, z1 = xkpc[inds], ykpc[inds], zkpc[inds]
    n1, g01, g_mech1 = n_gas_cgs[inds], G0[inds], g_mech[inds]
    Av1 = Av_mean[inds]
    m1 = m_cgs[inds]
    logger.debug('got the sph particles in the required ranges')

    #computing alpha for each sph particle (amount of gm wrt surface heating for the PDR model)
    f_log_gamma_surf=arxvPDR.construct3DInterpolationFunction(quantity=['therm','heating'], slabIdx=0, log10=True, interpolator='nearest') 
    dataNew = numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.ones(n1.shape)*lgMechMin]).T
    gammaSurf_sph_from_pdr = 10.0**f_log_gamma_surf(dataNew)
    Tkin_sph_from_pdf = 10.0**f_log_gamma_surf(dataNew)
    alpha_sph = g_mech1/gammaSurf_sph_from_pdr
    
    indsNan = numpy.where(numpy.isnan(Tkin_sph_from_pdf))[0]
    if indsNan.size != 0:
        raise ValueError('interpolated variable Tkin_sph_from_pdf has nan values in it!!')
    
    #selecting particles which have alpha < 1.0
    logger.debug('compute the alpha for the sph particles')
    
    maps = params['maps']

    print 'getting the emissions of each sph particle'
    #getting emissions of the sph particles with mechanical heating    
    em_sph = {}
    data_with_gm  =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.log10(g_mech1), Av1]).T
    data_no_gm    =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.ones(n1.size)*lgMechMin, Av1]).T
    
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
                asdasdad
                #setting emissions of particles with nan interpolated values to zero             
                indsNan = numpy.where(numpy.isnan(em_sph[map_key]))[0]
                if indsNan.size != 0:
                    em_sph[map_key][indsNan] = 0.0
                time_interp += dt
            
    logger.debug('done getting the emissions of each sph particle in %.2f seconds' % time_interp)
    
    #global xx
    #xx = [data_with_gm, em_sph['f_mean_em_CO-1-0'], alpha_sph]

    """    
    m1      : masses of the particles we selected for the visualization  
    n1      : masses of the particles we selected for the visualization  
    g01     : masses of the particles we selected for the visualization  
    g_mech1 : masses of the particles we selected for the visualization  
    em_sph  : a dict of arrays each holding the emissions of a specific 
              line of the sph particles we selected for visualization. The keys
              are the same as those of params['maps'] i.e f_mean_C+-1-0 for example.
    """
     
    #---------------------------------------------------------------------------------------
    # Done getting the distribution of mass, particle count as a function of number density
    #---------------------------------------------------------------------------------------
    print 'getting the spatial distrubutions'
    
    nBins = params['imres']    #mesh resolution use in binning
    
    for map_idx_key in maps:
        maps[map_idx_key]['data'] = numpy.zeros((nBins, nBins), dtype=numpy.float64) 
        maps[map_idx_key]['f_logn'] = None 
    
    hist = hist_nd(numpy.vstack((x1,y1)), mn = bsMin, mx=bsMax, nbins=nBins, reverse_indicies=True) 
    for i in numpy.arange(nBins):
        for j in numpy.arange(nBins):
            inds_in_bin = hist.get_indicies([i,j])
    
            maps['f_n_part']['data'][i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:
                maps['f_mass']['data'][i,j]    = m1[inds_in_bin].sum()
                maps['f_mean_n']['data'][i,j]  = numpy.mean(n1[inds_in_bin])
                maps['f_mean_g0']['data'][i,j] = numpy.mean(g01[inds_in_bin])
                maps['f_mean_gm']['data'][i,j] = numpy.mean(g_mech1[inds_in_bin])
                maps['f_mean_Av']['data'][i,j] = numpy.mean(Av_mean[inds_in_bin])

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
            maps['f_n_part']['f_logn'] = get_histogram(numpy.log10(n1), numpy.ones(n1.shape), 100.0, -3.0, 3.0)
        if map_idx_key == 'f_mass':
            maps['f_mass']['f_logn'] = get_histogram(numpy.log10(n1), m1, 100.0, -3.0, 3.0) 
        if map_idx_key == 'f_mean_n':
            maps['f_mean_n']['f_logn'] = get_histogram(numpy.log10(n1), n1, 100.0, -3.0, 3.0)
        if map_idx_key == 'f_mean_g0':
            maps['f_mean_g0']['f_logn'] = get_histogram(numpy.log10(n1), g01, 100.0, -3.0, 3.0) 
        if map_idx_key == 'f_mean_gm':
            maps['f_mean_gm']['f_logn'] = get_histogram(numpy.log10(n1), g_mech1, 100.0, -3.0, 3.0)
        for key in em_sph:
            if '_em_' in key:
                maps[key]['f_logn'] = get_histogram(numpy.log10(n1), em_sph[key], 100.0, -3.0, 3.0)
    #
    """
    
    #
    #
    def fetch_grid_data(map_key, log10 = None):
        
        if log10 != None and log10 == True:
            return numpy.log10(maps[map_key]['data'])
        else:
            return maps[map_key]['data']
    #
    #

    #displaying all the maps in a single plot    
    print 'done getting the spatial distributuions'
    
    fig, axs = pylab.subplots(4, 4, sharex=True, sharey=True, figsize=(12, 12), 
                              subplot_kw = {'xlim':[bsMin, bsMax],
                                            'ylim':[bsMin, bsMax],
                                            'aspect':'equal',
                                            'adjustable':'datalim',
                                           })
                              
                                                                                                
    for ax in axs[:,0]: ax.set_ylabel('y(kpc)')
    for ax in axs[3,:]: ax.set_xlabel('x(kpc)')
    
    pylab.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.9, wspace=0.15, hspace=0.15)
    
    for map_key in maps:
        
        i, j = maps[map_key]['pos']
    
        im_map = maps[map_key] 
        
        map_data = fetch_grid_data(map_key, log10 = im_map['log10']).T
        
        print 'idx = %d,%d quantity = %-30s max,min = %s,%s' % (im_map['pos'][0], im_map['pos'][1], map_key, map_data.min(), map_data.max()) 

        #clipping the map values outside the specified colorbar ranges
        map_data = numpy.clip(map_data, im_map['v_rng'][0], im_map['v_rng'][1])
        
        im = axs[i,j].imshow(map_data, 
                             extent=[bsMin, bsMax, bsMin, bsMax],
                             vmin=im_map['v_rng'][0], 
                             vmax=im_map['v_rng'][1], 
                             interpolation='bessel', #intepolation used for imshow
                             origin='lower')
        axs[i,j].set_title(im_map['title'], size='large')
        
        pylab.colorbar(im, ax=axs[i,j], orientation='vertical')
        
        axs[i,j].set_xlim([bsMin, bsMax])
        axs[i,j].set_ylim([bsMin, bsMax])
        
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
        yPlt = numpy.log10(maps[map_key]['f_logn'].f/maps[map_key]['f_logn'].f.max())
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
#
#

t0 = time.time()
if pool_size == 1:
    for snap in snaps:
        generate_map(snap)
else:
    pool = multiprocessing.Pool(pool_size)
    pool.map(generate_map, snaps)
print 'total time = %s s' % (time.time() - t0)


def test_plot():
    nPts = 10000.0
        
    n_arr  = numpy.log10(numpy.ones(nPts)*1e2)
    
    G0_arr = numpy.log10(numpy.ones(nPts)*1e2)
    
    #gm_arr = numpy.log10(numpy.ones(nPts)*1e-24)
    gm_arr = numpy.linspace(-30.0, -16, nPts)
    
    #Av_arr = numpy.linspace(5.0, 10.0, nPts)
    Av_arr = numpy.ones(nPts)*5.0
    
    data  =  numpy.array([
                          n_arr, 
                          G0_arr, 
                          gm_arr,
                          Av_arr, 
                         ]
                        ).T
    
    vi = funcs['CO-1-0'](data)
    pylab.plot(gm_arr, vi, '.')
    pylab.show()


    