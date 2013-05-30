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
from mylib.utils.misc  import xselect
from mylib.utils.histogram import hist_nd
from fi_utils import make_maps, parse_old_runinfo_file
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

pool_size = 1
snaps = numpy.arange(4, 5, 1)

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std',   # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol-linked',   # the path of the dir containing the simulation
          'imres' : 100,                                                  # resolution of the maps to be produced imres x imres
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/',        # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2/',        # the path to the dir containing the PDR database          
           
          'ranges' : {#ranges in n,g0 and gm of the sph particles to be included in producing the maps
                      'sph':{
                             'min_log_n_use'  : -3.0,      
                             'min_log_G0_use' : -3.0,
                             'min_log_gm_use' : -50.0,
                            },
                      
                      #ranges of the pdr models to be used in constructing the interpolation functions
                      'pdr': {
                              'log_n' : [0,0]
                             },
                      #the size of the box to be displayed (particles outside the range are discarded)
                      'box_size' : [-4,4] | units.kpc,
                      'Av'       : 10.0, 
                      },          
          'maps'   : {
                      'f_n_part' : {'pos': [0,0], 'title': r'$f(N)$'        , 'v_rng': [0,6]      , 'log10': True}, 
                      'f_mass'   : {'pos': [0,1], 'title': r'$f(m)$'        , 'v_rng': [30,40]    , 'log10': True},
                      'f_mean_n' : {'pos': [0,2], 'title': r'$f(\bar{n})$'  , 'v_rng': [-3.0, 4.0], 'log10': True},
                      'f_mean_g0': {'pos': [0,3], 'title': r'$f(\bar{g_0})$', 'v_rng': [-3.0, 3.0], 'log10': True},
                      #--------
                      'f_mean_gm'              : {'pos': [1,0], 'title': r'$f(\bar{\Gamma_m})$'                  , 'v_rng': [-35.0, -22.0], 'log10': True},
                      'f_mean_em_C+-1-0'       : {'pos': [1,1], 'title': r'$f(L_{C^+ 158 \mu m})$'               , 'v_rng': [-6.0, -2.0]  , 'log10': True},
                      'f_mean_em_no_gm_C+-1-0' : {'pos': [1,2], 'title': r'$f(L_{C^+ 158 \mu m})$ $\Gamma_m = 0$', 'v_rng': [-6.0, -2.0]  , 'log10': True},
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
#=============================================
#loading the molecular specie radex database
arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'], Av=params['ranges']['Av'])

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
        fInterp = interpolator_sectioned_by_gm(data, v)
    else:
        fInterp = params['interpolator'](data, v)
        
    logger.debug('time contruncting interpolator for quantity %s %s in %.4f seconds.' % (kwargs['quantity'][1], kwargs['quantity'][3], time.time() - t0))

    return fInterp


def interpolator_sectioned_by_gm(data, v):
    """Returns an interpolator function which makes use of sub-interpolation function for ranges in 
       mechanical heating"""

    """    
    intervals = [[-50.0, -40.0],
                 [-45.0, -35.0],
                 [-35.0, -30.0],
                 [-30.0, -25.0],
                 [-25.0, -20.0],
                 [-20.0, -15.0],
                ] 
    ghost = 1.1
    """

    intervals = [[-50.0, -40.0],
                 [-45.0, -35.0],
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
    ghost = 1.1
                                        
    #splitting data into sections specified by gm_intervals and constructing an interpolation
    #function for each interval (the intevals include an extra dex in each bound of the specified
    #bound to make sure interpolation is done correctly...i.e the scipy interpolator doesnt return
    # nan's)
    
    gm_data = data[:,2]

    all_fInterps = []     
    for interval in intervals:
        
        gm_interval_low  = interval[0]-ghost
        gm_interval_high = interval[1]+ghost
        
        inds_in_this_interval = numpy.where( 
                                            (gm_data >= gm_interval_low)*
                                            (gm_data <= gm_interval_high) 
                                           )[0]
        
        logger.debug('       making interpolation from for gm range [%f %f]' % (gm_interval_low, gm_interval_high))                                     

        data_in_this_interval = data[inds_in_this_interval, :]
        v_in_this_interval = v[inds_in_this_interval]
        fInterp_this_interval = params['interpolator'](data_in_this_interval, v_in_this_interval)
        all_fInterps.append(fInterp_this_interval)
        
    def fInterp(data):
        """devides data into sections of gm (3rd colomn) and returns the interpolated values by passing
        them to the sub-interpolation function in all_fInterps"""
        
        nPts = data.shape[0]
        
        #array which will hold the interpolated values from all the sections
        all_values = numpy.zeros(nPts, dtype=numpy.float64)
        
        gm_data = data[:,2]
        inds_orig = numpy.arange(nPts) 
        
        #getting the values from each section
        for i, interval in enumerate(intervals):
            
            inds_in_this_interval = numpy.where( 
                                                (gm_data >= interval[0])*
                                                (gm_data <= interval[1]) 
                                               )[0]
                                               
            if inds_in_this_interval.size == 0:
                continue
                                                 
            data_in_this_interval = data[inds_in_this_interval, :]

            t0 = time.time()            
            v = all_fInterps[i](data_in_this_interval)
            dt = time.time()-t0
            
            #appending values from this section to the values to be returned
            all_values[inds_orig[inds_in_this_interval]] = v
            
            logger.debug('     time interpolating for sub-interval %d in %.3f seconds , %d points at %e points/sec' % (i, dt, inds_in_this_interval.size, inds_in_this_interval.size/dt) )

        return all_values
    
    return fInterp

def get_interpolation_function_radex(**kwargs):
    """makes an interpolation function from the radex meshes."""
    
    specStr  = kwargs.pop('specStr')
    up_to_Av = kwargs.pop('up_to_Av')
    
    arxvPDR.use_radexDb(specStr=specStr, Av=up_to_Av) #.. todo:: do this for different Avs too    

    v = arxvPDR.apply_function_to_all_radex_meshes(funcRadex, func_kw = kwargs)

    v = numpy.array(v)
    inds_valid = numpy.isfinite(v)
    v = v[inds_valid]
     
    xGrd, yGrd, zGrd = arxvPDR.grid_x[inds_valid], arxvPDR.grid_y[inds_valid], arxvPDR.grid_z[inds_valid]
    data = numpy.array([xGrd, yGrd, zGrd], dtype = numpy.float64).T
    
    t0 = time.time()    
    if 'sectioned' in kwargs and kwargs['sectioned'] == True:
        fInterp = interpolator_sectioned_by_gm(data, v)
    else:
        fInterp = params['interpolator'](data, v)
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
                                                               up_to_Av = params['ranges']['Av'],
                                                               transitionIdx = lines_info[line]['idx'],
                                                               sectioned = True,
                                                              )
            if lines_info[line]['type'] == 'pdr':
                funcs[line] = get_interpolation_function_pdr(
                                                             quantity = lines_info[line]['quantity'],
                                                             up_to_Av = params['ranges']['Av'],
                                                             sectioned = True,
                                                             )
            
    return funcs

#getting the time unit
conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

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

xdbg = None

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

    def Av_sph_particles():
        Tk = gasT #sph gas particle temperature in K
        sigma_v = gas.vdisp.value_in(units.km / units.s)  #sph particle velocity despertion in km/s
        
        Pe = (1.085 * gasT + 54.0 * sigma_v**2.0)*n_gas_cgs
                
        Av_mean = 0.22 * metallicity * ( 1520.0 / 100.0) * numpy.sqrt(Pe/1e4)
        
        return Av_mean
    
    Av_mean = Av_sph_particles()
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
    
    xkpc, ykpc, zkpc = gas.position.value_in(params['ranges']['box_size'].unit).T
    
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
    x1,y1, z1 = xkpc[inds], ykpc[inds], zkpc[inds]
    n1, g01, g_mech1 = n_gas_cgs[inds], G0[inds], g_mech[inds]
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
    data_with_gm  =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.log10(g_mech1)]).T
    data_no_gm    =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.ones(n1.size)*lgMechMin]).T

    time_interp = 0.0
    
    for map_key in maps:
        if '_em_' in map_key:
            if '_no_gm_' in map_key:
                data_use = data_no_gm
            else:
                data_use = data_with_gm

            line = map_key.split('_')[-1] # the line XYZ-UPPER-LOWER
            
            t0 = time.time()
            em_sph[map_key] = 10.0**funcs[line](data_use)
            dt = time.time() - t0
            
            logger.debug('time interpolating emission for map %s = %.3f seconds for %d points at %e points/sec' % (map_key, dt, data_use.shape[0], data_use.shape[0]/dt) )
            
            #setting emissions of particles with nan interpolated values to zero             
            indsNan = numpy.where(numpy.isnan(em_sph[map_key]))[0]
            if indsNan.size != 0:
                em_sph[map_key][indsNan] = 0.0
            time_interp += dt
            
    logger.debug('done getting the emissions of each sph particle in %.2f seconds' % time_interp)
    
    #global xdbg
    #xdbg = [data_with_gm, em_sph['f_mean_em_CO-1-0'], alpha_sph]

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
    
    fig, axs = pylab.subplots(4, 4, sharex=True, sharey=True, figsize=(12,12), 
                              subplot_kw = {'xlim':[bsMin, bsMax],
                                            'ylim':[bsMin, bsMax],
#                                            'aspect':'equal',
#                                            'adjustable':'datalim',
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
