import pdb
import os, sys, time
import subprocess, shlex

import numpy
from numpy import log10, where, argmin, fabs, linspace, log, exp, pi, sqrt
from numpy import isfinite, isnan, zeros

from scipy.ndimage import zoom
from scipy import stats
from scipy import interpolate

import pylab
from pylab import histogram
from matplotlib import ticker

import time

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
from amuse.io import read_set_from_file, fi_io

import meshUtils
import lineDict
import pdrDict
import line_ratio_utils
import constraining

from mylib.constants import M_SUN_SI
from mylib.utils.histogram import hist_nd
from mylib.utils.interpolation import sectioned_4D_interpolator 
import mylib.units
from mylib.constants import M_SUN_SI
from mylib.utils.ndmesh import ndmesh
from mylib.utils.misc import scale
from mylib.utils import templates

import collections
import multiprocessing

def parse_old_runinfo_file(path):
    """parses the old runinfo file and returns its contetnts as a dict"""
    
    dict_ret = dict()
    
    for i, line in enumerate(open(path)):
        
        line = line.strip()
        
        if len(line) == 0 or line[0] == 'C':
            continue
        
        if i == 1:
            continue
        
        if i == 2:
            key = 'datadir'
            value = line
            
        if i >= 3:
            
            lineSplt = line.split()
            
            if ' - ' in line:
                if  '..' not in line:
                    value = True
                    key = ''.join(lineSplt[3::])
                else:
                    continue
            else:
                value = lineSplt[0]
                key = '-'.join(lineSplt[1::])
                
            dict_ret[key] = value
    
    
    return dict_ret    

def convert_units_to_pdr_units(gas, metallicity):
    """convert some of the attributes of the gas particles to units that are compatible with the 
    units used in the PDR code. This addes the attributes :
    
        Av, G0, gmech, x, y, z, n, T, vdisp, Pe 

    the gas is returned as a gas_set object
    """
    
    gas_new = gas_set(len(gas))
    
    pos = gas.position.value_in(units.kpc).T

    gas_new.x = pos[0]
    gas_new.y = pos[1]
    gas_new.z = pos[2]

    gas_new.vx = gas.vx.value_in(units.kms) 
    gas_new.vy = gas.vy.value_in(units.kms)
    gas_new.vz = gas.vz.value_in(units.kms) 
    
    mean_spec_weight = 1.3|units.amu
    gas_new.n = ((gas.rho / mean_spec_weight).value_in( units.cm **-3 ))

    gas_new.G0 = 6.54*gas.fuvheat
    
    gas_new.T = gas.temperat.number
    
    # converting the mechanical heating rate from per unit mass (of H gas)     
    # to per unit volume (see notesISM.odt)
    dethdt  = gas.dethdt.value_in( units.erg / (units.g * units.s ) )
    gas_new.gmech = dethdt * 1.6474e-24 * gas_new.n     
            
    gas_new.mass = gas.mass.value_in(units.g)

    gas_new.vdisp = gas.vdisp.value_in(units.km / units.s)
    
    #computing the mean Av        
    gas_new.Pe = (1.085 * gas_new.T + 54.0 * gas_new.vdisp**2)*gas_new.n #how to set the unit to the same untis as  (units.K*gas.n.unit*constants.kB.unit)
    gas_new.Av = 0.22 * metallicity * ( 1520.0 / 100.0) * numpy.sqrt(gas_new.Pe/1e4)

    return gas_new
    
def select_particles(gas, ranges):
    '''
    returns the SPH particle set  which have attributes withing the specied radges'''

    #selecting particles within a spatial range and within the 
    #density and g0 range of the database
    inds = where(
                 (numpy.log10(gas.n)  >= ranges['sph']['min_log_n_use'])*
                 (numpy.log10(gas.G0) >= ranges['sph']['min_log_G0_use'])*
                 (numpy.log10(gas.gmech) >= ranges['sph']['min_log_gm_use'])*                   
                 (gas.x >= ranges['box_size'][0].number)*(gas.x <= ranges['box_size'][1].number)*
                 (gas.y >= ranges['box_size'][0].number)*(gas.y <= ranges['box_size'][1].number)*
                 (gas.Av >= ranges['sph']['Av_use'][0])*(gas.Av <= ranges['sph']['Av_use'][1])                       
                 )[0]

    #using a minimum Av according to the speciecied range
    #Av_mean = Av_mean.clip(params['ranges']['sph']['Av'][0], params['ranges']['sph']['Av'][1])
    
                    
    return gas[inds]
    
def compute_alpha(gas, arxvPDR):
    """Computes the ratio of the mechanical heating of an SPH particle to the surface heating of
    a pure PDR with the same n,G0.
    """
    
    lgMechMin = arxvPDR.grid_z.min() #the minim mechanical heating  (in log) which will be assumed to be that of a pure PDR

    f_log_gamma_surf=arxvPDR.construct3DInterpolationFunction(quantity=['therm','heating'], 
                                                              slabIdx=0, 
                                                              log10=True, 
                                                              interpolator='linear')

    dataNew = numpy.array([
                           numpy.log10(gas.n), 
                           numpy.log10(gas.G0), 
                           numpy.ones(gas.n.shape)*lgMechMin
                          ]).T
                          
    print 'interpolating the surface heating at the surface of the SPH particles'
    gammaSurf_sph_from_pdr = 10.0**f_log_gamma_surf(dataNew)
    gas.alpha = gas.gmech/gammaSurf_sph_from_pdr
    
    indsNan = where(numpy.isnan(gammaSurf_sph_from_pdr))[0]
    if indsNan.size != 0:
        err_strng = 'interpolated variable Tkin_sph_from_pdf has %d nans values in it!!' % indsNan.size
        err_strng += '\n\t\tgMechMin = %e' % lgMechMin                 
        raise ValueError(err_strng)
    
    
def guess_metallicity(runDir):
    if '-std' in runDir:
        metallicity = 0.2
    if '-sol' in runDir:
        metallicity = 1.0
    return metallicity

def gas_particle_info_from_interp_func(gas, map_key, arxvPDR, interp_funcs, take_exp=None):
    '''given a bunch of gas particles, returns the emission according to what is specied in the paramter
     map_key
    '''
    
    #what points to use? with or witout gm
    if '_no_gm_' in map_key:
        data_use = numpy.array([log10(gas.n), log10(gas.G0), numpy.ones(len(gas))*arxvPDR.grid_z.min(), gas.Av]).T
    else:
        data_use = numpy.array([log10(gas.n), log10(gas.G0), log10(gas.gmech), gas.Av]).T
    print 'interpolating emissions of %s of the SPH particles for %d particles' % (map_key, data_use.shape[0])

    t0 = time.time()
    v_ret = interp_funcs[map_key](data_use)
    dt = time.time() - t0

    if take_exp != None and take_exp == True:
        v_ret = 10.0**v_ret

    print 'time interpolating for %s = %.3f seconds for %d points at %e points/sec' % (map_key, dt, data_use.shape[0], data_use.shape[0]/dt)

    #setting emissions of particles with nan interpolated values to zero             
    indsNan = where(numpy.isnan(v_ret))[0]
    if indsNan.size != 0:
        v_ret[indsNan] = 1e-30
                
    return v_ret


def save_gas_particle_info_saperate_files(path_of_snapshot, gas, species):
    '''save the data of the SPH particles into a files. Only the attributes in the list are saved. The
     file is written into a npz file.
    '''

    #-----------------------------------------
    #saving the particle states
    filename = path_of_snapshot + '.states.npz'
    
    # the basic attributes to be saved
    attr_save = ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'n', 'G0', 'gmech', 'Av', 'T', 'alpha']
    
    # additional attributes related to the sampled particles
    for extra_attr in ['weights', 'parent', 'children']:
        if extra_attr in dir(gas):
            attr_save.append(extra_attr)
             
    save_gas_particle_info(filename, gas, attr_save)
    print '\t\t\twrote file : %s ' % filename
    #-----------------------------------------
    
    #the list of all the attributes    
    attr_list = dir(gas)
    
    #saving the emission info for each species into a saperate file
    for specStr in species:

        attr_save_this_spec = []
        
        for attr in attr_list:
            
            if ('em_' in attr or 'empdr_' in attr) and '__' not in attr:
                
                lineStr = attr.split('_')[-1]
         
                if specStr == lineDict.lines[lineStr]['specStr']:
                    attr_save_this_spec.append(attr)
                #
            #
        #
        
        
        if len(attr_save_this_spec) > 0:
            print attr_save_this_spec
            filename_spec = filename + '.em_%s.npz' % specStr
            save_gas_particle_info(filename_spec, gas, attr_save_this_spec)
            print '\t\t\twrote file : %s ' % filename_spec
        else:
            print 'no emission info found for this species : %s' % specStr
    #

    #saving the info related to the PDR (non emission stuff), i.e attributes which do not have 'em' in them
    #but have the '_pdr_' as a sub string
    attr_save_pdr = []
        
    for attr in attr_list:
            
        if ('pdr_' in attr and 'empdr_' not in attr) and '__' not in attr:
                
            attr_save_pdr.append(attr)
        #
    #
        
    if len(attr_save_pdr) > 0:
        print attr_save_pdr
        filename_pdr = filename + '.pdr.npz'
        save_gas_particle_info(filename_pdr, gas, attr_save_pdr)
        print '\t\t\twrote file : %s ' % filename_pdr
    else:
        print 'no pdr info found'
    #


def save_gas_particle_info(filename, gas, attr_name_list):
    '''save the data of the SPH particles into a file. Only the attributes in the list are saved. The
     file is written into a npz file.
    '''
    
    #putting the attribute names strings into a numpy array
    attr_arr = numpy.array(attr_name_list)
    
    #putting the attribute data into a tuple
    attr_data_list = ()
    for attr_name in attr_arr:
        attr_data_list += (getattr(gas, attr_name),)
        
    numpy.savez_compressed(filename, *attr_data_list, names=attr_arr)

def load_gas_particle_info(filename, load_pdr=None):
    '''load the particle states from the file'''
    
    data = numpy.load(filename)
    print '\tloaded file:\n\t\t%s' % filename    
    n = data['arr_0'].size

    gas = gas_set(n)

    for i, name in enumerate(data['names']):
        setattr(gas, name, data['arr_%d' % i])
        print '\t\tattr_name = %-30s [%d elements]' % (name, data['arr_%d' %i].size)
        
    all_names = data['names'] 
    
    if load_pdr == True:
        filename_pdr = filename + '.pdr.npz'
        data_pdr = numpy.load(filename_pdr)
        print '\tloaded file:\n\t\t%s' % filename_pdr
        for i, name in enumerate(data_pdr['names']):
            setattr(gas, name, data_pdr['arr_%d' % i])
            print '\t\tattr_name = %-30s [%d elements]' % (name, data['arr_%d' %i].size)
            
        all_names = numpy.hstack((all_names, data_pdr['names']))

    return gas, all_names

def load_gas_particle_info_with_em(filename, species, load_pdr=None):
    '''loads the specified file foo.states.npz, also looks for files
     foo.states.npz.em_XYZ.npz, where XYZ is an entry in the list of strings in species
    '''
    
    names_all = []
    
    #loading the states of the SPH particles
    gas, names = load_gas_particle_info(filename, load_pdr=load_pdr)
    
    names_all.append(names.tolist())
    
    for specStr in species:
        filename_this = filename + '.em_%s.npz' % specStr
        gas_this, names_this = load_gas_particle_info(filename_this)
        for name_this in names_this:
            attr_data_this = getattr(gas_this, name_this)
            setattr(gas, name_this, attr_data_this)
    
        names_all.append(names_this.tolist())
        
        del gas_this
    
    print 'att attributes loaded are : ', names_all
    
    print 'warning: make sure the computed optical depths make sense!!! it doesnt seem so..'
    return gas

def make_map(gas, hist, attr=None, as_log10=None, func=None, show=False, **kwargs):
    '''looping over the bins of the 2D histogram of the x,y coordinates and 
     computing the averages of the maps
    '''

    map_data = zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)

    if attr not in dir(gas[0]):  
        raise AttributeError('%s is not an attribute of the gas particles,' % attr)
        
    attr_data = getattr(gas, attr)
    if 'weights' in kwargs:
        w_data = getattr(gas, kwargs['weights'])
    
    for i in numpy.arange(hist.nBins[0]):
        
        for j in numpy.arange(hist.nBins[1]):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            map_data[i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:

                x = attr_data[inds_in_bin]
                
                if 'weights' in kwargs:
                    w = w_data[inds_in_bin]
                    
                if 'weights' in kwargs:
                    map_data[i,j] = func(x, weights=w)
                else:
                    map_data[i,j] = func(x)
            #
        #
    #
    
    if as_log10 != None and as_log10 == True:
        map_data = numpy.log10(map_data)

    if show!=False:
        
        if show == 'log10':
            data_show = numpy.log10(map_data).T
        else:
            data_show = map_data
             
        pylab.figure()
        pylab.imshow(data_show, origin='lower')
        pylab.colorbar()
        
    return map_data

def funcPDR(meshObj, **kwargs):
    
    quantity = kwargs['quantity']
    up_to_Av = kwargs['up_to_Av']
    
    value = meshObj.compute_integrated_quantity(quantity, Av_range = [0.0, up_to_Av])
    
    return log10(value)

def get_interpolation_function_pdr(arxvPDR, params, logger, **kwargs):
    """makes an interpolation function from a integrated PDR quantity"""

    Av_range = params['ranges']['interp']['Av']
    Av_res   = params['ranges']['interp']['Av_res']

    #computing the quantities at the desired location in Av, which will be used
    #to construct the interpolation functions
    Avs = linspace(Av_range[0], Av_range[1], (Av_range[1] - Av_range[0])/Av_res + 1) 

    if 'extra_Av_sec' in params['ranges']['interp']:
        Avs = numpy.hstack((params['ranges']['interp']['extra_Av_sec'], Avs))
    
    kwargs['up_to_Av'] = {}
    func = kwargs.pop('func')
    
    print '\t\tcomputing info from the PDR mesh at Av = : '
    for i, Av in enumerate(Avs):
        
        print '                                              %.2f' % Av
        kwargs['up_to_Av'] = Av
        
        #getting the data corresponding to this Av
        v = arxvPDR.apply_function_to_all_meshes(func, func_kw = kwargs)

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


    t0 = time.time()
    if 'sectioned' in kwargs and kwargs['sectioned'] == True:
        fInterp = sectioned_4D_interpolator(data_all_Av, v_all_Av, params['interpolator'])
    else:
        fInterp = params['interpolator'](data_all_Av, v_all_Av)
    logger.debug('time contruncting interpolation function for %s in %.4f seconds' % (kwargs['descr'], time.time() - t0))

    return fInterp


def get_interpolation_function_radex(arxvPDR, params, logger, **kwargs):
    """makes an interpolation function from the radex meshes from the database.
    :param string specStr: The specie whose info will be used (emission, or anything else radex outputs)
    :param int transitionIdx: The index of the transition (for example for CO: transitionIdx = 0, uses info of the J=1-0 transition 
    :param bool sectioned: if true uses a sectioned 4D interpolator. 
    """
    
    specStr  = kwargs.pop('specStr')
    Av_range = params['ranges']['interp']['Av']
    
    #reading only the radex databases which are within the specified range in Av
    available_Avs_for_radex_dbs = numpy.sort(numpy.array(arxvPDR.radexDbs.keys(), dtype='f8'))
    inds_read = where(
                      (available_Avs_for_radex_dbs >= Av_range[0])*
                      (available_Avs_for_radex_dbs <= Av_range[1])
                     )[0]
    Avs_read = available_Avs_for_radex_dbs[inds_read]
    
    if 'extra_Av_sec' in params['ranges']['interp']:
        Avs_read = numpy.hstack((params['ranges']['interp']['extra_Av_sec'], Avs_read))

    Av_min, Av_max = numpy.min(Avs_read), numpy.max(Avs_read)
        
    #collecting the data from the database corresponding to all the data
    #in Av avaiable
    for i, Av in enumerate(Avs_read):
        
        #using the database of the specific Av
        arxvPDR.use_radexDb(specStr=specStr, Av=Av)
        
        #getting the data corresponding to this Av
        v = arxvPDR.apply_function_to_all_radex_meshes(meshUtils.radex_mesh_quantity, func_kw = kwargs)
                    
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
        v_all_Av = v_all_Av.clip(0.0, 10000.0)
    
    #using the subset of points in the radex databse specified by params['ranges']['interp']
    inds_keep = where( 
                      (data_all_Av[:,0] >= params['ranges']['interp']['log_n'][0])* 
                      (data_all_Av[:,0] <= params['ranges']['interp']['log_n'][1])*
                             
                      (data_all_Av[:,1] >= params['ranges']['interp']['log_G0'][0])*
                      (data_all_Av[:,1] <= params['ranges']['interp']['log_G0'][1])*
                             
                      (data_all_Av[:,2] >= params['ranges']['interp']['log_gmech'][0])* 
                      (data_all_Av[:,2] <= params['ranges']['interp']['log_gmech'][1])*
                             
                      (data_all_Av[:,3] >= Av_min)*
                      (data_all_Av[:,3] <= Av_max)
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

def make_interp_funcs_from_arxv(arxvPDR, em_keys, pdr_keys, params, logger):
    '''returns a dictionary of functions with keys XYZ-UPPER-LOWER. Items in 'maps' which have
    _em_ are parsed and the transition is looked up (f_mean_em_C+-1-0 is parsed into C+-1-0) 
    and an interpolation function for C+-1-0 is computed.
    '''
    print 'contructing interpolation functions'
    
    arxvPDR.logger = logger

    em_interp_funcs  = {}
    pdr_interp_funcs = {}
    
    for map_key in em_keys:
        
        line = map_key.split('_')[-1] # an entry in the dictionary lineDict.lines
                    
        if lineDict.lines[line]['type'] == 'radex-lvg':
            em_interp_funcs[map_key] = get_interpolation_function_radex(
                                                                        arxvPDR,    
                                                                        params,
                                                                        logger,
                                                                        specStr  = lineDict.lines[line]['specStr'],
                                                                        quantity = map_key.split('<q>')[-1].split('</q>')[0],
                                                                        transitionIdx = lineDict.lines[line]['radexIdx'],
                                                                        sectioned = True,
                                                                       )
        ''' (not tested)
        if lines_info[line]['type'] == 'pdr':
            em_interp_funcs[map_key] = get_interpolation_function_pdr(
                                                                      quantity = lines_info[line]['quantity'],
                                                                      sectioned = True,
                                                                      )
        '''

    for map_key in pdr_keys:
        
        #getting interpolation function for column densities
        print map_key
        quantity = map_key.split('_')[-1] # an entry in the dictionary pdrDict.pdr
        
        pdr_interp_funcs[map_key] = get_interpolation_function_pdr(
                                                                   arxvPDR,
                                                                   params,
                                                                   logger,
                                                                   specStr = pdrDict.pdr[quantity]['specStr'],
                                                                   quantity = pdrDict.pdr[quantity]['ismcpak'],
                                                                   descr = pdrDict.pdr[quantity]['descr'],
                                                                   func = pdrDict.pdr[quantity]['func'],
                                                                   as_log10 = params['maps'][map_key]['log10_interp_func'], 
                                                                   sectioned = True,
                                                                  )

        
    print 'done making the interpolation functions'
    
    return em_interp_funcs, pdr_interp_funcs

def snapshot_interpolated_data(snap, arxvPDR, em_interp_funcs, pdr_interp_funcs, 
                               em_keys, pdr_keys, params, logger):
    '''get the emission info using supplied interpolation functions for a certain snapshot
     Returns the requested emissions as a dict and also returnes the gas particles with those
     set as attributes.
    '''
    
    ## snapshot index
    snapIndex = snap
    
    suffix = '%06d' % snapIndex
    
    ## path to fi snapshots
    snapName = 'fiout.%s' % suffix

    ## full path to the snapshot
    snap_filename = params['rundir'] + '/firun/' + snapName
    
    ## appending the extention if the parameters instruct using a sampled snapshot
    if 'use_sampled_set' in params and params['use_sampled_set'] == True:
        snap_filename += '.ext.npz'

    ## extracting/guessing the metallicity from the name of the directory of the run
    metallicity = guess_metallicity(params['rundir'])
    
    ## loading the SPH simulation data
    logger.debug('loading snapshot %s : ' % snap_filename)
    t0 = time.time()
    if 'use_sampled_set' in params and params['use_sampled_set'] == True:
        
        gas, attr_read = load_gas_particle_info(snap_filename)
        
    else:
        gas, dark, stars = read_set_from_file(snap_filename, format = fi_io.FiFileFormatProcessor)

        ## getting a new particle set for the gas particles with the attributes we will be using later
        ## and converting the units to the same as those of the PDR models    
        gas = convert_units_to_pdr_units(gas, metallicity)
        logger.debug('converted and computed stuff from SPH data (in units compatible with the pdr models)')
    
    file_size_bytes = os.stat(snap_filename).st_size
    dt = time.time() - t0
    logger.debug('done reading fi snapshot : %s \n\t\t\t time reading = %.2e seconds @ %.2e MB/s' % (snap_filename, dt, (file_size_bytes/dt)/1024.0**2) )
    logger.debug('number of SPH particles in snapshot = %d' %  len(gas))
        
    print 'quantity          min              max  '
    print 'n_gas         %e     %e   '  % (gas.n.min()     , gas.n.max())
    print 'G0_gas        %e     %e   '  % (gas.G0.min()    , gas.G0.max())
    print 'g_mech        %e     %e   '  % (gas.gmech.min() , gas.gmech.max())
    print 'Av_gas        %e     %e   '  % (gas.Av.min()    , gas.Av.max())
    
    ## keeping gas particles within the specified ranges
    gas = select_particles(gas, params['ranges'])
    logger.debug('got the SPH particles in the required ranges')

    lgMechMin = arxvPDR.grid_z.min()

    #setting the lower bound of the mechanical heating of the SPH particles of the 
    #minimum of the pdr database which is almost negligable (since some of the 
    #SPH particles have a zero gmech, their log would be -inf, so we clip it to the
    #mininmum of the pdr database
    gas.gmech = numpy.clip(gas.gmech, 10.0**lgMechMin, numpy.Infinity)
    #########################done setting up the snapshot data######################    
    
    #computing alpha for each SPH particle (amount of gm wrt surface heating for the PDR model)
    #and selecting particles which have alpha < 1.0
    compute_alpha(gas, arxvPDR)
    logger.debug('computed the alpha for the SPH particles')
    inds = where( gas.alpha < 1.0 )
    gas = gas[inds]
     
    #setting the maximum Av of the SPH particles to a predifined value (usually 
    #the emissions do not increase after some Av bec of the optical depth of the lines)
    gas.Av = numpy.clip(gas.Av, params['ranges']['sph']['Av_clip'][0], params['ranges']['sph']['Av_clip'][1])
    
    #computing the emissions of the SPH particles from the PDR interpolation functions
    print 'getting the emissions of each SPH particle'    

    em_sph = {}
    for em_key in em_keys: em_sph[em_key] = gas_particle_info_from_interp_func(gas, 
                                                                               em_key,  
                                                                               arxvPDR, 
                                                                               em_interp_funcs, 
                                                                               take_exp=params['maps'][em_key]['log10'])

    pdr_sph = {}        
    for pdr_key in pdr_keys: pdr_sph[pdr_key] = gas_particle_info_from_interp_func(gas, 
                                                                                   pdr_key, 
                                                                                   arxvPDR, 
                                                                                   pdr_interp_funcs, 
                                                                                   take_exp=params['maps'][pdr_key]['log10'])
    
    return em_sph, pdr_sph, gas

def plot_map(map_data, map_params, map_info, snap_time, params, processed_snap_filename):
    '''
    '''
    bs_min, bs_max = map_params['ranges']['box_size'].number
    
    fig = pylab.figure(figsize=(8,8))

    ax = fig.add_axes([0.15, 0.085, 0.75, 0.75])

    ax.set_xlabel('x(kpc)', size='large')    
    ax.set_ylabel('y(kpc)', size='large')
    ax.set_xlim([bs_min, bs_max])
    ax.set_ylim([bs_min, bs_max])
    
    inds_nan = where(numpy.isfinite(map_data) == False)
    map_data[inds_nan] = map_info['v_rng'][0]
    
    import matplotlib.cm as cm
    
    im = ax.imshow(map_data,  #i think the transopse of map_data should be taken here
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=map_info['v_rng'][0],
                   vmax=map_info['v_rng'][1],
                   interpolation='bessel', #intepolation used for imshow
                   cmap = cm.OrRd,
                   origin='lower')

    ax.set_title(map_info['title'], size=25, verticalalignment='bottom')
    ax.tick_params(axis='both', which='major')#, labelsize=20)
    
    cbar_ax = fig.add_axes([0.2, 0.97, 0.6, 0.02]) 
    xticks=linspace(map_info['v_rng'][0], map_info['v_rng'][1], 8)
    pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal',
                   ticks=xticks)
    
    xticks_strs = [r'$10^{%d}$' % v for v in xticks]
    cbar_ax.set_xticklabels(xticks_strs, size='large')

    
    pylab.figtext(0.01, 0.87, '%.2f' % snap_time + 'Gyr', 
            color='black', size='large', weight='bold')

    if 'save_image' in params and params['save_image'] == True:
        filename_fig = processed_snap_filename + '.eps' 
        fig.savefig(filename_fig)
        print 'saved image file to :\n\t\t\t %s' % filename_fig
    
    #copy the saved image to the latex direcotry
    if 'latex' in params and params['latex'] == True:

        copy_to = os.path.join(
                               params['latex_dir'], '%s-%s-%s' % (
                                                                  os.path.split(os.path.split(params['rundir'])[-2])[-1],    
                                                                  os.path.split(params['rundir'])[-1],
                                                                  os.path.split(filename_fig)[-1],
                                                                 )
                              )
        command = 'cp %s %s' % (filename_fig, copy_to)
        args = shlex.split(command)
        subprocess.call(args)
        
        print 'copied %s to \n\t %s' % (filename_fig, copy_to) 

    return fig

def plot_cube_sections(data_cube, params):
    '''given a cube, it plots the sections (in velocity channels) of the emission'''
    
    #parameters of a typical spectrum for the cube
    v_min, v_max = params['cube']['spec_rng']
    v_res = params['cube']['spec_res']
    
    #the velocity bins (.. todo:: xxx shift this to the centroid of the velcoty bins) 
    v = linspace(v_min, v_max, v_res)
    
    v_sec_wdith = (v_max - v_min)/params['cube']['plot']['n_vsec_plt']
    
    nx = params['cube']['plot']['n_per_row']
    ny = params['cube']['plot']['n_vsec_plt'] / nx
    
    fig_height = 2.5*ny
    fig_width = 7 #fig_height*numpy.float(nx)/numpy.float(ny)
    fig, axs = pylab.subplots(ny, nx, sharex=True, sharey=True, figsize=(fig_width, fig_height), 
                              subplot_kw = {'xlim':params['cube']['xy_rng'][0:2],
                                            'ylim':params['cube']['xy_rng'][2:4],
                                            'xticks' : params['cube']['plot']['xticks'],
                                            'yticks' : params['cube']['plot']['yticks'],
                                            'aspect':'auto',
                                            'adjustable': 'datalim',
                                            })
    axs = axs.reshape(ny, nx)
    for ax in axs[:,0] : ax.set_ylabel('y(kpc)')
    for ax in axs[-1,:]: ax.set_xlabel('x(kpc)')
        
    pylab.subplots_adjust(left=0.10, bottom=0.20, right=0.99, top=0.7, wspace=0.05, hspace=0.15)
    
    #plotting the maps for all the velocity channels
    for n in numpy.arange(params['cube']['plot']['n_vsec_plt']):
     
        #range in velocity for this velocity channel
        v_from = n*v_sec_wdith + v_min
        v_to   = v_from + v_sec_wdith
        #inds in the z direction of the cube (along the spectrum for this channel)
        v_inds = where( (v > v_from) * (v <= v_to ) )  
    
        #computing the map for this channel by adding all the spectra for this channel
        this_map = data_cube[:,:, v_inds[0]]
        this_map = numpy.log10( this_map.mean(axis=2) ) 
        
        #plotting the map
        j_plt = n / nx
        i_plt = n % nx
        
        
        print '(i,j) = %d, %d ' % (i_plt, j_plt), v_from, v_to
        print this_map.shape
        print '-------------------------'
        
        this_map = this_map.clip(*params['cube']['plot']['rng'])
        
        im = axs[j_plt, i_plt].imshow(this_map.T, 
                                      extent=params['cube']['xy_rng'],
                                      vmin=params['cube']['plot']['rng'][0], 
                                      vmax=params['cube']['plot']['rng'][1], 
                                      interpolation='bessel', #intepolation used for imshow
                                      origin='lower')
        
        axs[j_plt, i_plt].text(-7, 7, '%.0f' % ((n + 0.5)*v_sec_wdith + v_min) + 'km/s', color='w')
    #

    cbar_ax = fig.add_axes([0.2, 0.85, 0.6, 0.03]) 
    cbar_ax.tick_params(axis='both', which='major', labelsize=12)
    cbar_ax.set_title(r'$\log_{10}$ [ <T$_{\rm B}$> / K ]', size=12)
    

    cbar_ticks = numpy.linspace(params['cube']['plot']['rng'][0], 
                                params['cube']['plot']['rng'][1], 6)
    xticks_strs = [r'$10^{%d}$' % v for v in cbar_ticks]
    
    pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal', 
                   ticks=cbar_ticks)
    cbar_ax.set_xticklabels(xticks_strs)
    
    pylab.show()
    
    return fig 

def plot_map_from_saved_data(snapIndex, map_info, params):
    '''
    '''
    #path to processed fi snapshot  
    snap_data_filename = os.path.join(params['rundir'],'analysis', 'fiout.%06d.%s.npz' % (snapIndex, map_info['attr']))  

    #loading the data from disk
    loaded_data = numpy.load(snap_data_filename)

    #the map data and the parameters used to generate the data
    map_data, map_params = loaded_data['map_data'], loaded_data['params'].tolist()

    #the time of the snapshot
    snap_time = get_snapshot_time(snapIndex, params)

    fig = plot_map(map_data, map_params, map_info, snap_time, params, snap_data_filename)
    
    return fig
#

def plot_all_maps_for_snapshot_from_saved_data(snapIndex, params):
    '''
    '''
    
    figs = []
        
    #getting all the maps
    for this_map in params['all_maps']:
    
        this_map_info = params['all_maps'][this_map]
        
        fig = plot_map_from_saved_data(snapIndex, this_map_info, params)
        
        figs.append(fig)
        
    print '---------------------------------------------------------------'
    pylab.show()
    
    return figs
#

def get_useful_gas_attr_from_dview(dview, gas_var_name, view_num):
    '''
        #gets the attributes defined in this function from the 0th entry in the dview object
        
        data_dict = get_use_gas_attr_from_dview(dview, 'gas', 0)
        
    '''
    
    useful_gas_attr = ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'n', 'G0', 'gmech', 'Av', 'T', 'alpha']
    
    data = {}
    for attr in useful_gas_attr:
        data[attr] = dview[gas_var_name + '.%s' % attr][view_num]

    gas_new = Particles(data['mass'].size)
    
    for attr in useful_gas_attr:
        setattr(gas_new, attr, data[attr])
        
    return gas_new
    

def get_emission_only_keys(map_keys):
    '''from the names of the map keys, all the keys with an _em_ are filtered out 
     and returned as emission keys
    ''' 

    em_keys = []
    for map_key in map_keys:      
        if '_em_' in map_key or '_empdr_' in map_key:
            em_keys.append(map_key) 

    return em_keys

def get_pdr_only_keys(map_keys):
    '''picks keys which have the tag _pdr_ in them and not _em_. i.e one which correspond to info from the pdr
    which are not related to the emission.
    '''
    
    pdr_keys = []
    for map_key in map_keys:      
        if '_pdr_' in map_key:
            pdr_keys.append(map_key) 

    return pdr_keys

def get_snapshot_time(snapIndex, params):
    
    #getting the time unit
    conv = nbody_system.nbody_to_si(1 | units.kpc, 1e9 | units.MSun)
    timeUnit = conv.to_si(1 | nbody_system.time).in_(units.Gyr)

    #getting the time of the snapshot
    path = params['rundir'] + '/firun/runinfo'
    runinfo = parse_old_runinfo_file(path)
    snap_time = (float(runinfo['dtime'])*timeUnit.number) * (float(runinfo['noutbod']) * snapIndex)

    return snap_time

class galaxy_gas_mass_estimator(object):
    
    def __init__(self, luminosity_maps=None, gas=None, params=None, hist=None, 
                 line_ratios=None, arxvPDR=None, **kwargs):
    
        self.luminosity_maps = luminosity_maps #: the dict containing the luminosity maps
        self.gas = gas #: the object of gas particles
        self.params = params #: the parameters used in controlling everything
        self.hist = hist
        self.line_ratios = line_ratios
        self.arxvPDR = arxvPDR
        
        #####################
        self.hist_obs = None
        self.H2_mass_mesh = None
        self.H_mass_mesh = None
        self.H2_mass_no_gmech_mesh = None
        self.H_mass_no_gmech_mesh = None    
        self.fig = None        
        self.axs = None
        self.contraining = None
        self.means = None
        self.sdev = None
        
    def get_model_emission_from_involved_line_ratios(self):
        '''loading the emission info from all the models for all Avs (also we check for 
        the consistenscy of the number of models...i.e same number of models 
        for all the lines)
        '''
        
        ## make line ratios from the mock luminosities 
        obs_mock_ratios_template = line_ratio_utils.ratios()
        
        for line_ratio in self.line_ratios:
        
            line1, line2 = line_ratio_utils.lines_involved(line_ratio)
            
            v1, v2 = 1, 1
            
            obs_mock_ratios_template.make_ratios(
                                                 {
                                                   line1:{'fluxKkms': 1, 'err': 1.0}, 
                                                   line2:{'fluxKkms': 1, 'err': 1.0},
                                                 },
                                                 ratios = [line_ratio],
                                                 em_unit = 'fluxKkms'
                                                )
        
        ## determin the species and the line codes in the line ratios
        obs_mock_ratios_template.species_and_codes()
        
        ## getting the emission of the lines involved in the ratios from the PDR archive
        print 'getting the emission from the PDR models...'
        model_em = {}
        for i, line in enumerate(obs_mock_ratios_template.codes):
            v, grid_coords = self.arxvPDR.get_emission_from_all_radex_dbs_for_Av_range(
                                                                                       line = line,
                                                                                       Avs = 'all', 
                                                                                       quantity = 'fluxKkms',
                                                                                       keep_nans = True,
                                                                                      )
            model_em[line] = 10.0**v
            print '\t%-12s %d' % (line, v.size), grid_coords.shape
        
            ## some checks of the sizes
            if v.size != grid_coords.shape[0]:
                raise ValueError('number of elements in the emission values is different from the number of modesl.')
            
            if i == 0:
                nModels = v.size
            else:
                if nModels != v.size:
                    raise ValueError('the number of elements for this line differes at least from that of one of the other lines')
        
        self.mode_em = model_em
        self.grid_coords = grid_coords
        
    def setup_observed_grid(self):
        '''sets up the observed mesh. A lower resolution mesh than the synthetica map'''
        
        bs_min, bs_max = self.params['ranges']['box_size'].number

        fig = pylab.figure(0, (8,8))
        
        axs = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        self.fig, self.axs = fig, axs
        
        pylab.imshow(log10(self.luminosity_maps['maps']['CO1-0']), 
                     extent=[bs_min, bs_max, bs_min, bs_max], 
                     origin='lower')

        #pylab.plot(self.gas.x[::100], self.gas.y[::100], '.', markersize=1)
        
        #pylab.plot(self.hist.f.spos[0], self.hist.f.spos[1], 'w+', markersize=10)

        hist_obs = hist_nd(numpy.vstack((self.gas.x, self.gas.y)), 
                           mn = bs_min, mx = bs_max, 
                           nbins = self.params['obs_res'], 
                           reverse_indicies = True, loc = True)
        hist_obs.info()
        
        obs_mesh = hist_obs.f 
        
        pylab.plot(obs_mesh.spos[0], obs_mesh.spos[1], 'k+', markersize=100, linewidth=20)

        fig.canvas.mpl_connect('button_press_event', self.inspect_estimation)
            
        pylab.show()

        self.hist_obs = hist_obs
        self.H2_mass_mesh = zeros(hist_obs.f.shape, 'f8')
        self.H_mass_mesh = zeros(hist_obs.f.shape, 'f8')
        self.H2_mass_no_gmech_mesh = zeros(hist_obs.f.shape, 'f8')
        self.H_mass_no_gmech_mesh = zeros(hist_obs.f.shape, 'f8')
 
    def inspect_estimation(self, event):

        bs_min, bs_max = self.params['ranges']['box_size'].number
        obs_res = self.params['obs_res']
        
        # get the x and y coords, flip y from top to bottom
        xd, yd = event.xdata, event.ydata
        
        if event.button == 1:
    
            if self.axs.contains(event)[0]:
                
                print 'cursor x,y = ', xd, yd
 
                ## plotting the centroid of the pixel clicked in the coarse grid 
                x_ind = numpy.floor(scale(xd, 0.0, obs_res, bs_min, bs_max))
                y_ind = numpy.floor(scale(yd, 0.0, obs_res, bs_min, bs_max))
                
                self.estimate_mass_in_pixel(x_ind, y_ind, plot_fits=True)
                
    def estimate_mass_in_all_pixels(self):
        
        res_x, res_y = self.hist_obs.f.shape
        
        for i in numpy.arange(res_x):
            for j in numpy.arange(res_y):
                
                masses_in_pixel = self.estimate_mass_in_pixel(x_ind = i, y_ind = j)
                #mH_1e9_m_sun, mH2_1e9_m_sun, mH_no_gmech_1e9_m_sun, mH2_no_gmech_1e9_m_sun
                
                self.H_mass_mesh[i,j] = masses_in_pixel[0]
                self.H2_mass_mesh[i,j] = masses_in_pixel[1]
                self.H_mass_no_gmech_mesh[i,j] = masses_in_pixel[2]
                self.H2_mass_no_gmech_mesh[i,j] = masses_in_pixel[3]
                
    
    def estimate_mass_in_pixels_at_radius(self, r=None):
        '''constrains using line ratios the parameters of clouds with and without gmech for 
        pixels within r - dr, r + dr, where dr is the diagonal of a pixel
        
        pixels_info = estimator.estimate_mass_in_pixels_at_radius(r=2.0)
        
        #info of the first pixel in the list
        constraining_0, means_0 = pixels_info[0]
        #the fit parms for a certain pixel
        parms, xi, parms_no_gm, xi_no_gm = constraining_0.get_model_parms_for_min_Xi2()
        
        ''' 
        
        hist_obs = self.hist_obs

        x, y = hist_obs.f.cntrd

        # the radial distance of the pixel from the center of the map
        r_pixels = numpy.sqrt(x**2 + y**2)
        
        # the diagonal width of a pixel
        dr = numpy.sqrt(hist_obs.szBins[0]**2 + hist_obs.szBins[1]**2)
        
        # indicies of pixels within a radial range
        inds = numpy.where( (r_pixels < (r + dr/2)) * (r_pixels > (r - dr/2))  )
        
        pixels_info = []
        
        # doing the fits for all the pixels within the radial range
        for x_ind, y_ind in numpy.vstack(inds).T:
            
            x_this = hist_obs.f.cntrd[0][x_ind, y_ind]
            y_this = hist_obs.f.cntrd[1][x_ind, y_ind]
              
            pylab.plot(x_this, y_this, 'o')
            
            masses_in_pixel = self.estimate_mass_in_pixel(x_ind = x_ind, 
                                                          y_ind = y_ind)
            
            #mH_1e9_m_sun, mH2_1e9_m_sun, mH_no_gmech_1e9_m_sun, mH2_no_gmech_1e9_m_sun
                
            #self.H_mass_mesh[i,j] = masses_in_pixel[0]
            #self.H2_mass_mesh[i,j] = masses_in_pixel[1]
            #self.H_mass_no_gmech_mesh[i,j] = masses_in_pixel[2]
            #self.H2_mass_no_gmech_mesh[i,j] = masses_in_pixel[3]
        
            pixels_info.append([self.contraining, self.means])

        #getting the mean value of the fit parameters
        parms_all, xi_all, parms_all_no_gm, xi_all_no_gm = [], [], [], []
        means_all = []
        for constr, means, in pixels_info:
            parms, xi, parms_no_gm, xi_no_gm = constr.get_model_parms_for_min_Xi2()
            parms_all.append(parms)
            xi_all.append(xi)
            parms_all_no_gm.append(parms_no_gm)
            xi_all_no_gm.append(xi_no_gm)
            means_all.append(means.values())
        
        #degrees of freedom        
        DOF = numpy.float64(len(constr.model_em.keys()) - 4)
 
        print '#########################################################'
        print '####means of pixels at R = %-.2f kpc                   ###' % r 
        print '#########################################################'
        fit_means = tuple(numpy.array(parms_all).mean(axis=0))
        fit_sdev  = tuple(numpy.array(parms_all).std(axis=0))
        print 'mean fit parms          = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % fit_means 
        print 'std  fit parms          = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % fit_sdev
        print 'mean Xi2 per DOF        = %-+8.3f' % (numpy.array(xi_all).mean() / DOF)
        print 'std  Xi2 per DOF        = %-+8.3f' % (numpy.array(xi_all).std()  / DOF)
        print '-----------------------------------------------------------------------------------------'
        print 'mean fit parms   (pure) = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % tuple(numpy.array(parms_all_no_gm).mean(axis=0))
        print 'std  fit parms   (puer) = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % tuple(numpy.array(parms_all_no_gm).std(axis=0))
        print 'mean Xi2 per DOF (pure) = %-+8.3f' % (numpy.array(xi_all_no_gm).mean() / (DOF + 1.0))
        print 'std  Xi2 per DOF (pure) = %-+8.3f' % (numpy.array(xi_all_no_gm).std()  / (DOF + 1.0))
        print '-----------------------------------------------------------------------------------------'
        print 'actual mean cloud parameter values at '
        print '-----------------------------------------------------------------------------------------'
        print 'mean parms SPH          = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % tuple(numpy.array(means_all).mean(axis=0))
        print 'sdev parms SPH          = %-+8.3f %-+8.3f %-+8.3f %-+8.3f' % tuple(numpy.array(means_all).std(axis=0))
        
        return pixels_info
    
    def estimate_mass_in_pixel(self, x_ind=None, y_ind=None, plot_fits=False):
       
        hist = self.hist
        hist_obs = self.hist_obs
        obs_mesh = self.hist_obs.f
        line_ratios = self.line_ratios
        luminosity = self.luminosity_maps
        params = self.params
        model_em = self.mode_em 
        grid_coords = self.grid_coords 
        arxvPDR = self.arxvPDR
        m_gas = self.gas.mass
        
        ######################################################
        
        area_obs_pixel = numpy.product(obs_mesh.dl)
        
        if plot_fits == True:
            pylab.plot(obs_mesh.cntrd[0][x_ind, y_ind], 
                       obs_mesh.cntrd[1][x_ind, y_ind], 
                       'r+', markersize=10)

        ## plotting the centroids of the pixels of the pixel of the map inside the clicked pixel
        inds = where(
                     (hist.f.cntrd[0] >= obs_mesh.spos[0][x_ind, y_ind])*\
                     (hist.f.cntrd[0] <= obs_mesh.epos[0][x_ind, y_ind])*\
                     (hist.f.cntrd[1] >= obs_mesh.spos[1][x_ind, y_ind])*\
                     (hist.f.cntrd[1] <= obs_mesh.epos[1][x_ind, y_ind])
                    )
         
        if plot_fits == True:
            pylab.plot(hist.f.cntrd[0][inds], hist.f.cntrd[1][inds], 'c+')
        
        inds_gas = hist_obs.get_indicies([x_ind, y_ind])

        if plot_fits == True:
            
            pylab.plot(self.gas.x[inds_gas], self.gas.y[inds_gas], 'o')
        
            pylab.draw()
        
        ## make line ratios from the mock luminosities 
        obs_mock_ratios = line_ratio_utils.ratios()                

        for line_ratio in line_ratios:
        
            line1, line2 = line_ratio_utils.lines_involved(line_ratio)
            
            print line1, line2

            lum1, lum2 = luminosity['maps'][line1][inds].sum(), luminosity['maps'][line2][inds].sum() 

            int1, int2 = lum1/area_obs_pixel, lum2/area_obs_pixel
            
            print '--------------'
            
            obs_mock_ratios.make_ratios(
                                        {
                                          line1:{'fluxKkms': int1, 'err': params['error_bars']*int1}, 
                                          line2:{'fluxKkms': int2, 'err': params['error_bars']*int2}
                                        },
                                        ratios = [line_ratio],
                                        em_unit = 'fluxKkms',
                                        lum = {line1: lum1, line2: lum2}
                                       )
            '''
            print line_ratio, obs_mock_ratios[line_ratio]
            
            obs_mock_ratios.species_and_codes()
            '''
        #
        
        ## fitting
        f = constraining.Xi2_line_ratios_single_component(obs_data = obs_mock_ratios, 
                                                          model_data = model_em, 
                                                          model_parms = grid_coords,
                                                          line_ratios = line_ratios
                                                          )

        f.compute_model_line_ratios()

        f.compute_Xi2()
        
        f.print_minima()

        if plot_fits == True:
            
            fig = pylab.figure(figsize=(16,8))
            ax1 = fig.add_axes([0.1, 0.5, 0.15, 0.4])
            ax2 = fig.add_axes([0.3, 0.5, 0.15, 0.4])
            
            ax3 = fig.add_axes([0.1, 0.1, 0.15, 0.3])
            ax4 = fig.add_axes([0.3, 0.1, 0.15, 0.3])
            ax5 = fig.add_axes([0.5, 0.1, 0.15, 0.3])
            ax6 = fig.add_axes([0.7, 0.1, 0.15, 0.3])

            f.plot_results(fig = fig, ax = ax1)
            print '-------------------------'
            f.plot_results(fig = fig, ax = ax2, no_gmech=True)
            
            #plotting the emission distribution in the pixel 
            #g.em_fluxKkms_CO1-0, pdr_NH2, pdr_NH, pdr_NCO, pdr_N13CO
            gas_in_pixel = self.gas[inds_gas]
            
            # getting the emission distriubtion as a function of gas density
            #----------------------------------------------------------------------------------------
            log_n_gas = log10(gas_in_pixel.n).reshape((1, len(gas_in_pixel)))
            n_dist = hist_nd(log_n_gas, nbins=50.0, mn=-3.0, mx=4.0, loc=True, reverse_indicies=True)
            
            # the gas particles within the ranges of the distribution histogram 
            gas_in_pixel_in_hist = gas_in_pixel[n_dist.inds_in]
            
            #getting the emission distribution (a selection of lines)
            dist_0 = zeros(n_dist.totalBins, 'f8') 
            
            for i in numpy.arange(n_dist.totalBins):
                
                inds_in_bin = n_dist.get_indicies(i)
                
                if inds_in_bin.size != 0:
                    
                    gas_in_bin = gas_in_pixel_in_hist[inds_in_bin]
 
                    dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_CO1-0').sum()
                    #dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_13CO1-0').sum()
                #
            #
            
            x, y = n_dist.f.cntrd, n_dist.f
            
            ## ploting the locations in density where the cumilitive particle 
            ## mass (number) distribution sum is 10%, 50%, 90%
            ## also getting the CDF for the emission
            ax3.plot(zoom(x, 20, order=2), zoom(y / numpy.max(y), 20, order=3), 'b')
            ax3.plot(zoom(n_dist.f.cntrd, 20, order=3), 
                     zoom(dist_0 / numpy.max(dist_0), 20, order=3), 'r')

            # CDF of the mass 
            f_CDF_n = (y / y.sum()).cumsum()
            n_i = linspace(-3.0, 4.0, 1000.0)[::-1]
            f_CDF_n_i = numpy.interp(n_i, x, f_CDF_n)

            # CDF of the emission
            f_CDF_em = ((dist_0 / dist_0.sum())[::-1]).cumsum()[::-1]
            f_CDF_em_i = numpy.interp(n_i, x, f_CDF_em)
            
            # plotting the 10, 50, 90% indicators for the SPH mass CDF
            ind90percent = argmin(fabs(f_CDF_n_i - 0.9))
            ax3.plot([n_i[ind90percent], n_i[ind90percent]], [0, 2.5], 'b--')
            ax3.text(n_i[ind90percent] - 1.1, 2.35, '90%', color = 'b')
            ax3.text(n_i[ind90percent] - 1.1, 2.35 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind90percent]*100) + '%', 
                     color = 'r')

            ind50percent = argmin(fabs(f_CDF_n_i - 0.5))            
            ax3.plot([n_i[ind50percent], n_i[ind50percent]], [0, 2.2], 'b--', linewidth=1)
            ax3.text(n_i[ind50percent] - 1.1, 2.05, '50%' , color = 'b')
            ax3.text(n_i[ind50percent] - 1.1, 2.05 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind50percent]*100)  + '%', 
                     color = 'r')
            
            ind10percent = argmin(fabs(f_CDF_n_i - 0.1))
            ax3.plot([n_i[ind10percent], n_i[ind10percent]], [0, 1.85], 'b--', linewidth=1)
            ax3.text(n_i[ind10percent] - 1.1, 1.7 , '10%', color = 'b')
            ax3.text(n_i[ind10percent] - 1.1, 1.7 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind10percent]*100) + '%', 
                     color = 'r')

            
            # plotting the 10, 50, 90% indicators for the emission CDF
            ind10percent = argmin(fabs(f_CDF_em_i - 0.1))
            ax3.plot([n_i[ind10percent], n_i[ind10percent]], [0, 1.55], 'r--', linewidth=1)
            ax3.text(n_i[ind10percent], 1.4, '10%', color = 'r')
            ax3.text(n_i[ind10percent], 1.25, '%d' % (100.0 - f_CDF_n_i[ind10percent]*100) + '%', 
                     color = 'b')


            ind50percent = argmin(fabs(f_CDF_em_i - 0.5))
            ax3.plot([n_i[ind50percent], n_i[ind50percent]], [0, 1.25], 'r--', linewidth=1)
            ax3.text(n_i[ind50percent], 1.1, '50%', color = 'r')
            ax3.text(n_i[ind50percent], 0.95, '%d' % (100.0 - f_CDF_n_i[ind50percent]*100) + '%', 
                     color = 'b')
            
            ind90percent = argmin(fabs(f_CDF_em_i - 0.90))
            ax3.plot([n_i[ind90percent], n_i[ind90percent]], [0, 0.95], 'r--', linewidth=1)
            ax3.text(n_i[ind90percent], 0.8, '90%', color = 'r')
            ax3.text(n_i[ind90percent], 0.65, '%d' % (100.0 - f_CDF_n_i[ind90percent]*100) + '%', 
                     color = 'b')

            ax3.set_ylim([0, 2.5])
            ax3.set_xlabel(r'$\log_{10}$[n$_{\rm gas}$]')
            ax3.set_ylabel(r'$f$')
            
            
            # getting the emission distriubtion as a function of G0
            #----------------------------------------------------------------------------------------
            log_G0_gas = log10(gas_in_pixel.G0).reshape((1, len(gas_in_pixel)))

            G0_dist = hist_nd(log_G0_gas, nbins=50.0, mn=-3.0, mx=4.0, loc=True, reverse_indicies=True)
            
            # the gas particles within the ranges of the distribution histogram 
            gas_in_pixel_in_hist = gas_in_pixel[G0_dist.inds_in]
            
            #getting the emission distribution (a selection of lines)
            dist_0 = zeros(G0_dist.totalBins, 'f8') 
            
            for i in numpy.arange(G0_dist.totalBins):
                
                inds_in_bin = G0_dist.get_indicies(i) 
                
                if inds_in_bin.size != 0:
                    
                    gas_in_bin = gas_in_pixel_in_hist[inds_in_bin]
 
                    dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_CO1-0').mean()
                #
            #
            
            ax4.plot(G0_dist.f.cntrd, G0_dist.f / numpy.max(G0_dist.f))

            ax4.plot(G0_dist.f.cntrd, dist_0 / dist_0.sum(), 'k--')
                         
            ax4.set_ylim([0, 1.4])
            ax4.set_xlabel('log G0 gas')
            ax4.set_ylabel('f')
            
            
            # getting the emission distriubtion as a function of gmech
            #----------------------------------------------------------------------------------------
            log_gmech_gas = log10(gas_in_pixel.gmech).reshape((1, len(gas_in_pixel)))

            gmech_dist = hist_nd(log_gmech_gas, nbins=50.0, mn=-30.0, mx=-20.0, loc=True, reverse_indicies=True)
            
            # the gas particles within the ranges of the distribution histogram 
            gas_in_pixel_in_hist = gas_in_pixel[gmech_dist.inds_in]
            
            #getting the emission distribution (a selection of lines)
            dist_0 = zeros(gmech_dist.totalBins, 'f8') 
            
            for i in numpy.arange(gmech_dist.totalBins):
                
                inds_in_bin = gmech_dist.get_indicies(i) 
                
                if inds_in_bin.size != 0:
                    
                    gas_in_bin = gas_in_pixel_in_hist[inds_in_bin]
 
                    dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_CO1-0').mean()
                #
            #
            
            ax5.plot(gmech_dist.f.cntrd, gmech_dist.f / numpy.max(gmech_dist.f))

            ax5.plot(gmech_dist.f.cntrd, dist_0 / dist_0.sum(), 'k--')
                         
            ax5.set_ylim([0, 1.4])
            ax5.set_xlabel('log gmech gas')
            ax5.set_ylabel('f')

            # getting the emission distriubtion as a function of gmech
            #----------------------------------------------------------------------------------------
            Av_gas = (gas_in_pixel.Av).reshape((1, len(gas_in_pixel)))

            Av_dist = hist_nd(Av_gas, nbins=50.0, mn=0.0, mx=30.0, loc=True, reverse_indicies=True)
            
            # the gas particles within the ranges of the distribution histogram 
            gas_in_pixel_in_hist = gas_in_pixel[Av_dist.inds_in]
            
            #getting the emission distribution (a selection of lines)
            dist_0 = zeros(Av_dist.totalBins, 'f8') 
            
            for i in numpy.arange(Av_dist.totalBins):
                
                inds_in_bin = Av_dist.get_indicies(i) 
                
                if inds_in_bin.size != 0:
                    
                    gas_in_bin = gas_in_pixel_in_hist[inds_in_bin]
 
                    dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_CO1-0').mean()
                #
            #
            
            ax6.plot(Av_dist.f.cntrd, Av_dist.f / numpy.max(Av_dist.f))

            ax6.plot(Av_dist.f.cntrd, dist_0 / dist_0.max(), 'k--')
                         
            ax6.set_ylim([0, 1.4])
            ax6.set_xlabel('Av gas')
            ax6.set_ylabel('f')            

            
        #computing the averages of the physical parameters of the particles inside the pixel 
        #-----------------------------------------------------------------------------------
        gas_in_pixel = self.gas[inds_gas]
        
        # getting the emission distriubtion as a function of gas density
        #---------------------------------------------------------------
        log_n_gas = log10(gas_in_pixel.n).reshape((1, len(gas_in_pixel)))
        n_dist = hist_nd(log_n_gas, nbins=50.0, mn=-3.0, mx=4.0, loc=True, reverse_indicies=True)
        
        # the gas particles within the ranges of the distribution histogram 
        gas_in_pixel_in_hist = gas_in_pixel[n_dist.inds_in]
        
        means = collections.OrderedDict()
        
        means['n']  = log10(gas_in_pixel.n.mean())
        means['G0'] = log10(gas_in_pixel.G0.mean())
        means['gm'] = log10(gas_in_pixel.gmech.mean())
        means['Av'] = gas_in_pixel.Av.mean()
        #-----------------------------------------------------------------------------------
                    
        ##  computing the masses in the pixel
        
        ## estimating the mass of H and H2 from the fitted PDR model (with and without gmech)
        PDR_parms, Xi2_min, PDR_parms_no_gmech, Xi2_min_no_gmech = f.get_model_parms_for_min_Xi2()
        
        area_kpc2 = area_obs_pixel
        area_cm2 = area_kpc2 * 1000.0**2 * mylib.units.PC2CM**2.0
        
        ## getting the column densities of H and H2
        log10_n, log10_G0, log10_gmech, Av = PDR_parms

        mesh_index = arxvPDR.get_mesh_index(x = log10_n, y = log10_G0, z = log10_gmech)
        
        arxvPDR.mshTmp.setData(arxvPDR.meshes[mesh_index])
        
        NH, NH2 = arxvPDR.mshTmp.getColumnDensity(specsStrs = ['H', 'H2'], maxAv = Av)
        
        mH_1e9_m_sun = (NH * area_cm2 * 1.67e-27) / (1e9 * M_SUN_SI)
        mH2_1e9_m_sun = (NH2 * area_cm2 * (2.0*1.67e-27)) / (1e9 * M_SUN_SI)  
                             
        ## getting the column densities of H and H2 (when no gmech is assumed)
        log10_n_no_gmech, log10_G0_no_gmech, log10_gmech_no_gmech, Av_no_gmech = PDR_parms_no_gmech
        
        mesh_index_no_gmech = arxvPDR.get_mesh_index(x = log10_n_no_gmech,
                                                     y = log10_G0_no_gmech,
                                                     z = log10_gmech_no_gmech)
        
        arxvPDR.mshTmp.setData(arxvPDR.meshes[mesh_index_no_gmech])
        
        NH_no_gmech, NH2_no_gmech = arxvPDR.mshTmp.getColumnDensity(specsStrs = ['H', 'H2'], maxAv = Av_no_gmech)
        
        mH_no_gmech_1e9_m_sun = (NH_no_gmech * area_cm2 * 1.67e-27) / (1e9 * M_SUN_SI)
        mH2_no_gmech_1e9_m_sun = (NH2_no_gmech * area_cm2 * (2.0*1.67e-27)) / (1e9 * M_SUN_SI)  

        print 'm(H)             = %.5f  m(H2) = %.5f m(H + H2) = %.5f' % (mH_1e9_m_sun, mH2_1e9_m_sun, mH_1e9_m_sun + mH2_1e9_m_sun)
        print 'm(H)  (no gmech) = %.5f  m(H2) = %.5f m(H + H2) = %.5f' % (mH_no_gmech_1e9_m_sun, mH2_no_gmech_1e9_m_sun, mH_no_gmech_1e9_m_sun + mH2_no_gmech_1e9_m_sun)
        print '                                      m_gas_sim       = %.5f' % ((m_gas[inds_gas].sum()/100.0)  / (1e9 * M_SUN_SI))

        if plot_fits == True:
            strng = '     with gmech : '
            strng += r'$\log_{10}$[n$_{\rm gas}]$ = %-+.2f' % log10_n + ' '
            strng += r'$\log_{10}$[G$_0]$ = %-+.2f' % log10_G0 + ' '
            strng += r'$\log_{10}$[$\Gamma_{\rm mech}$] = %-+.2f' % log10_gmech + ' '
            strng += r'$A_V$ = %-+.2f' % Av + '\n'

            strng += 'without gmech : '
            strng += r'$\log_{10}$[n$_{\rm gas}]$ = %-+.2f' % log10_n_no_gmech + ' '
            strng += r'$\log_{10}$[G$_0]$ = %-+.2f' % log10_G0_no_gmech + ' '
            strng += r'$\log_{10}$[$\Gamma_{\rm mech}$] = %-+.2f' % log10_gmech_no_gmech + ' '
            strng += r'$A_V$ = %-+.2f' % Av_no_gmech + '\n'
             
            ax1.set_title('with gmech')
            ax2.set_title('no gmech')
            
            pylab.figtext(0.5, 0.8, strng)
            
        
        f.get_model_em_for_min_Xi2()
        self.contraining = f
        self.means = means
        
        return mH_1e9_m_sun, mH2_1e9_m_sun, mH_no_gmech_1e9_m_sun, mH2_no_gmech_1e9_m_sun


def set_weights_sampled_to_zero(gas):

    inds_not_sampled = where( isnan(gas.parent) )
    inds_sampled = where( isfinite(gas.parent) )

    new_weights = gas.weights.copy()
    
    new_weights[inds_not_sampled] = 1.0
    new_weights[inds_sampled] = 0.0
    
    gas.weights = new_weights
    
    return  gas

def set_weights_not_sampled_to_zero(gas):

    inds_not_sampled = where( isnan(gas.parent) )
    inds_sampled = where( isfinite(gas.parent) )

    new_weights = gas.weights.copy()
    
    new_weights[inds_not_sampled] = 0.0
    new_weights[inds_sampled] = 1.0
    
    gas.weights = new_weights
    
    return  gas

    return gas

class gas_set(Particles):

    def copy_attr_from_amuse_set(self, gas_amuse, attr_list):
        '''copies the attributes from the amuse set to the gas_set'''
        
        for attr in attr_list:
            
            ## getting the data of the attribute of the amuse set 
            attr_gas_amuse = getattr(gas_amuse, attr)
            
            ## setting the repeated attributes to the sampled gas particle set
            setattr(self, attr, attr_gas_amuse)
        
        
    def fit_interval(self, n_min, n_max):
        '''Fit a certain interval in densities. Returns a scipy.stats.norm function
        with the fit deviates in the specified range and the number of points in that
        range.
        
        returns the fit function and the number of points used in the fit in the range
        of the densities which were specified by the arguments.
        '''

        ln = log(self.n)
            
        ## getting the log the ranges
        lmin_fit, lmax_fit = log(n_min), log(n_max)
        
        ln_fit = ln[ numpy.where( (ln > lmin_fit)*( ln < lmax_fit) ) ]

        ## fit a curve to the variates
        mu_fit, sigma_fit = stats.norm.fit(ln_fit)

        ## generating the function with the fitted deviates
        rv = stats.norm(loc=mu_fit, scale=sigma_fit)
        
        return rv, ln_fit.size
    
    def plot_PDF_CDF(self, fit_func_rng=None):
        
        fig, axs = pylab.subplots(2, 3, figsize=(12,8))
    
        
        # getting the PDF and the CDF of all the particles
        l10min, l10max = -6, 6
        nBins = 100
        
        bin_sz = (l10max - l10min)/float(nBins) 
        
        l10n = log10(self.n) 
        hpdf, bins = pylab.histogram(l10n, range=[l10min, l10max], bins=nBins, normed=True)
        
        # plotting the PDF of log10(n)
        axs[0,0].plot(bins[1::], hpdf, 'b')
        axs[0,0].set_title('PDF log10(n)')
        
        axs[0,1].semilogy(bins[1::], (hpdf*bin_sz).cumsum())
        axs[0,1].set_title('CDF log10(n)')
        axs[0,1].set_ylim(ymax=10.0)

        # plotting PDF(log(x))
        ln = log(self.n) 
        hpdf, lbins = pylab.histogram(ln, range=[log(1e-6), log(1e6)], bins=nBins, normed=True)
        axs[0,2].semilogy(lbins[1::], hpdf, 'b', linewidth=2)
        axs[0,2].set_title('PDF log(n)')

        ########################################################
        #####################  just plotting ################### 
        ########################################################
        
        #ranges = [ [1e-1, 1e1], [1e-2, 1e2], [1e-1, 1e2], [1e-1, 1e3], [1e-1, 1e4] ]
        #colors = [    'y'     ,   'g',         'k',          'c',        'm']
        ranges = [ fit_func_rng,  ]
        colors = [    'y'     , ]
        x = numpy.linspace(log(1e-4), log(1e6), 100)
        
        plts = ()
        labels = ()
        for i, rng in enumerate(ranges):
            
            rv, N = self.fit_interval(*rng)
            plt, = axs[0,2].semilogy(x, rv.pdf(x), colors[i])
            
            plts += (plt,)
            labels += '%d,%d' % (log10(rng[0]), log10(rng[1])),
            
        axs[0,2].legend(plts, labels, loc=0)
        
        ## setting the appropriate labels on the x-axis
        axs[0,2].set_xticks( [log(1e-4), log(1e-2), log(1e-0), log(1e2), log(1e4), log(1e6) ])
        axs[0,2].set_xticklabels([-4, -2, 0, 2, 4, 6 ])
    
        axs[0,2].set_ylim(1e-10, 1)    
        #########################################################
        #########################################################        
        
        return fig, axs
    
    def sample_ln_densities(self, npp, rem, rv):
        '''returns a sample of npp size of sampled densities in ln() scale, rv is the scipy 
        generator. The sample size from which npp particles are picked are based on the remaining
        number of particles to be sampled 'rem'.
        '''
         
        if npp == 5:
            if rem > 10000:
                ln_trial = rv.rvs(1000)
            if rem <= 10000:
                ln_trial = rv.rvs(10000)
            if rem <= 200:
                ln_trial = rv.rvs(100000)
            if rem <= 10:
                ln_trial = rv.rvs(1000000)
                
        if npp == 20:
            if rem > 3000:
                ln_trial = rv.rvs(10000)
            if rem <= 3000:
                ln_trial = rv.rvs(1000000)
            if rem <= 200:
                ln_trial = rv.rvs(1000000)
            if rem <= 10:
                ln_trial = rv.rvs(10000000)
                
        if npp == 50:
            if rem > 500:
                ln_trial = rv.rvs(100000)
            if rem <= 500:
                ln_trial = rv.rvs(1000000)
            if rem <= 300:
                ln_trial = rv.rvs(1000000)
            if rem <= 100:
                ln_trial = rv.rvs(10000000)

        if npp == 100:
            if rem > 500:
                ln_trial = rv.rvs(1000000)
            if rem <= 500:
                ln_trial = rv.rvs(1000000)
            if rem <= 300:
                ln_trial = rv.rvs(1000000)
            if rem <= 100:
                ln_trial = rv.rvs(50000000)

        return ln_trial
    
    def plot_sampled(self, ln_new, w_new, rv, axs=None):
        ## plotting the PDF again, on top of which the PDF of the sampled particles
        ## will be plotted and compared to
        nBins = 100
        
        ln = log(self.n)
        
        lhpdf, lbins = pylab.histogram(ln, range=[log(1e-4), log(1e6)], bins=nBins, normed=True)
        axs[1,0].semilogy(lbins[1::], lhpdf, 'b', linewidth=2)
        
        ## setting the appropriate labels on the x-axis
        axs[1,0].set_xticks( [log(1e-4), log(1e-2), log(1e-0), log(1e2), log(1e4), log(1e6) ])
        axs[1,0].set_xticklabels([-4, -2, 0, 2, 4, 6 ])
    
        axs[1,0].set_ylim(1e-10, 1)
        
        hpdf_s, lbins_s = pylab.histogram(ln_new, range=[log(1e-4), log(1e6)], bins=nBins, normed=True,
                                          weights=w_new,
                                         )
        axs[1,0].semilogy(lbins_s[1::], hpdf_s, 'r', linewidth=1)
    
        x = numpy.linspace(log(1e-4), log(1e6), 100)
        plt, = axs[1,0].semilogy(x, rv.pdf(x), 'm--')
    
        ###############################
        # getting again CDF of all the particles in log10
        l10min, l10max = -6, 6
        nBins = 100
        
        l10bin_sz = (l10max - l10min)/float(nBins) 
        
        l10n = log10(self.n) 
        l10hpdf, l10bins = pylab.histogram(l10n, range=[l10min, l10max], bins=nBins, normed=True)
    
        l10ns_new = log10(exp(ln_new)) 
        l10hpdf_new, l10bins_new = pylab.histogram(l10ns_new, range=[l10min, l10max], bins=nBins, normed=True,  
                                                   weights=w_new)
        
        # plotting the CDF of log10(n)
        axs[1,1].semilogy(l10bins[1::], (l10hpdf*l10bin_sz).cumsum(), 'b--')
        axs[1,1].semilogy(l10bins_new[1::], (l10hpdf_new*l10bin_sz).cumsum(), 'r--')
        
        
        axs[1,1].set_title('CDF log10(n)')
        axs[1,1].set_ylim(ymax=10.0)

        ###############################
        # sampling a big set based on the fit function and plotting its PDF
        
        lnmin, lnmax = log(1e-6), log(1e6)
        nBins = 100

        x = numpy.linspace(lnmin, lnmax, 100)
        plt, = axs[1,2].semilogy(x, rv.pdf(x), 'm--')
        
        sample = rv.rvs(1000000)
        
        lPDFsample, lbins = pylab.histogram(sample, range=[log(1e-6), log(1e6)], bins=nBins, normed=True)
        axs[1,2].semilogy(lbins[1::], lPDFsample, 'b', linewidth=2)
        
        ## setting the appropriate labels on the x-axis
        axs[1,2].set_xticks( [log(1e-4), log(1e-2), log(1e-0), log(1e2), log(1e4), log(1e6) ])
        axs[1,2].set_xticklabels([-4, -2, 0, 2, 4, 6 ])        
    
    def sample_set_serial(self, ln_gt_X, keys, rv, npp):
        ## sampling one particle at a time

        ## array which will hold the sampled densities
        ln_s = numpy.array([],'f8')
        ## array which will hold the keys of the sampled densities
        keys_s = numpy.array([],'f8')
        ## array which will hold the weights sampled densities
        w_s = numpy.array([],'f8')
        ## weights of the parent particles 
        w_ln_gt_X = zeros(ln_gt_X.size, 'f8')
                
        n_to_be_sub_sampled = ln_gt_X.size

        for i, ln_this in enumerate(ln_gt_X):
            
            rem = n_to_be_sub_sampled - i
            
            t0 = time.time()
            
            ts = 0.0
            
            attempts = 0

            while True:
                
                ts1 = time.time()
                ln_trial = self.sample_ln_densities(npp, rem, rv)
                ts += (time.time() - ts1)
                
                ## keeping sampled particles whose density is higher than the current one             
                ln_trial_gt_ln_this = ln_trial[ ln_trial > ln_this]
                ## keeping all particles
                #ln_trial_gt_ln_this = ln_trial[:]
    
                attempts += 1
                
                if  ln_trial_gt_ln_this.size > npp:
                    
                    sample_keep = ln_trial_gt_ln_this[:npp]
    
                    ## setting weights by volume to the sampled particles
                    w_this_sample = (1.0/float(npp+1.0))*(exp(ln_this)/exp(sample_keep))
                    ## setting zero weights
                    #w_this_sample = zeros(npp, 'f8')
                    ## setting unit weights
                    #w_this_sample = numpy.ones(npp, 'f8')
                    ## setting weights by mass
                    #w_this_sample = numpy.ones(npp, 'f8')/float(npp+1)
                    #((1.0/float(npp))*(exp(ln_this)/exp(sample_keep)))**(2.0/3.0),
                    #(1.0/float(npp))*numpy.ones(npp, 'f8')/10.0,
                    
                    ln_s = numpy.hstack((ln_s, sample_keep))
                    keys_s = numpy.hstack((keys_s, numpy.ones(npp)*keys[i]))
                    w_s = numpy.hstack((w_s, w_this_sample))
                    w_ln_gt_X[i] = 1.0 - w_this_sample.sum()
                    #w_ln_gt_X[i] = 1.0   
                                    
                    #pdb.set_trace()
                    
                    print '%.3f%% %d | %e | %d | ' % (100.0*float(i)/float(n_to_be_sub_sampled), 
                                                     n_to_be_sub_sampled - i,   
                                                     log10(exp(ln_this)), 
                                                     attempts), log10(exp(sample_keep).max())
                    print '----------------------------'
                    break
    
            print '% time sampling = ', ts / (time.time() - t0)
    
        return ln_s, keys_s, w_s, w_ln_gt_X
    
    def sample_higher_densities(self, npp=None, n_min_sample=None, fit_func_rng=None, plot=True):

        if plot == True:
            fig, axs = self.plot_PDF_CDF(fit_func_rng=fit_func_rng)
        
        ## picking one gaussian to be used in the sampled
        rv, N = self.fit_interval(fit_func_rng[0], fit_func_rng[1])

        ## testing the PDF
        #x = 2.0
        #y = (1.0/(sqrt(2.0*pi)*rv.std()))*exp( -0.5*((x - rv.mean())/rv.std())**2.0 )
        #print rv.pdf(x), y

        ln = log(self.n)
        
        ## particles whose density is larger than X will be sampled 
        X = n_min_sample
    
        # selecting the particles and their keys
        inds = numpy.where(ln > log(X))
        keys = self.key[inds]
        ln_gt_X = ln[inds]
        gas_ln_gt_X = self[inds] 
        
        # selecting the partilces whose density is less than X
        inds = where(ln < log(X))
        ln_lt_X = ln[ inds ]
        gas_ln_lt_X = self[inds]
        
        # soring the densities and the corresponding keys in ascending order
        inds_sorted = numpy.argsort(ln_gt_X)
        keys = keys[inds_sorted]
        ln_gt_X = ln_gt_X[inds_sorted]
        gas_ln_gt_X = gas_ln_gt_X[inds_sorted]

        n_to_be_sub_sampled = ln_gt_X.size
        print 'number of SPH particles to be sampled = %d' % n_to_be_sub_sampled
        print 'these are %.2f percent of the original set.' % (float(n_to_be_sub_sampled)/float(ln.size))

        ## doing the actual sampling        
        ln_s, keys_s, w_s, w_ln_gt_X = self.sample_set_serial(ln_gt_X, keys, rv, npp)
        
        ## merging the data        
        ln_new = numpy.hstack((ln_lt_X                       , ln_gt_X  , ln_s))        
        w_new  = numpy.hstack((numpy.ones(ln_lt_X.size, 'f8'), w_ln_gt_X, w_s))

        
        ## doing some plots including the sampled densitieis
        if plot == True:        
            self.plot_sampled(ln_new, w_new, rv, axs)
        
        n_s = exp(ln_s)
        
        print 'log10 max n = ', log10(n_s.max())    
        
        
        return n_s, w_s, gas_ln_gt_X, w_ln_gt_X, gas_ln_lt_X

    def get_inds_children(self):
        '''returns the indicies of the particles which are sampled, i.e the children'''
        
        return numpy.where(numpy.isfinite(self.parent))[0]

    def get_inds_original_set(self):
        '''returns the indicies of the original set of the particles'''
        
        return numpy.where(numpy.isnan(self.parent))[0]

    def get_inds_has_children(self):
        '''returns the indicies of the original set of the particles which have been sampled,
        i.e the particles which have children
        '''
        
        return numpy.where(numpy.isfinite(self.children))[0]
    
    def plot_pdf(self, qx=None, qy=None, log10x=False, log10y=False, nbins=50, xrng=None):
        '''getting the PDF of one quantity vs the other'''

        N = len(self) 
        
        ## getting the data to be processed from the attributes        
        x, y = getattr(self, qx), getattr(self, qy)
        if log10x == True: x = log10(x)
        if log10y == True: y = log10(y)
        
        if xrng != None:
            xmin, xmax = xrng
        else:
            xmin, xmax = x.min(), x.max()
            
        ## constructing the PDF of the first variable 'x'
        x_dist = hist_nd(x.reshape((1, N)), 
                         nbins=nbins, mn=xmin, mx=xmax, loc=True, reverse_indicies=True)

        ## keeping the data within the bounds of the histogram
        x_in_hist = x[x_dist.inds_in]
        y_in_hist = y[x_dist.inds_in]

        #getting the distribution of the y variable as a function of the histogram in x
        y_dist = zeros(x_dist.totalBins, 'f8')

        for i in numpy.arange(x_dist.totalBins):
           
            inds_in_bin = x_dist.get_indicies(i)
           
            if inds_in_bin.size != 0:
                 
                y_dist[i] = y_in_hist[inds_in_bin].sum()  #.mean  .average
            #
        #
        
        pylab.plot(x_dist.f.cntrd, x_dist.f)
        pylab.show()
        '''        
        #----------------------------------------------------------------------------------------
        
        #getting the emission distribution (a selection of lines)
        dist_0 = zeros(n_dist.totalBins, 'f8') 
       
        for i in numpy.arange(n_dist.totalBins):
           
            inds_in_bin = n_dist.get_indicies(i)
           
            if inds_in_bin.size != 0:
                
                gas_in_bin = gas_in_pixel_in_hist[inds_in_bin]
 
                dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_CO1-0').sum()
                #dist_0[i] = getattr(gas_in_bin, 'em_fluxKkms_13CO1-0').sum()
            #
        #
       
        x, y = n_dist.f.cntrd, n_dist.f
       
        ## ploting the locations in density where the cumilitive particle 
        ## mass (number) distribution sum is 10%, 50%, 90%
        ## also getting the CDF for the emission
        ax3.plot(zoom(x, 20, order=2), zoom(y / numpy.max(y), 20, order=3), 'b')
        ax3.plot(zoom(n_dist.f.cntrd, 20, order=3), 
                 zoom(dist_0 / numpy.max(dist_0), 20, order=3), 'r')

        # CDF of the mass 
        f_CDF_n = (y / y.sum()).cumsum()
        n_i = linspace(-3.0, 4.0, 1000.0)[::-1]
        f_CDF_n_i = numpy.interp(n_i, x, f_CDF_n)

        # CDF of the emission
        f_CDF_em = ((dist_0 / dist_0.sum())[::-1]).cumsum()[::-1]
        f_CDF_em_i = numpy.interp(n_i, x, f_CDF_em)
       
        # plotting the 10, 50, 90% indicators for the SPH mass CDF
        ind90percent = argmin(fabs(f_CDF_n_i - 0.9))
        ax3.plot([n_i[ind90percent], n_i[ind90percent]], [0, 2.5], 'b--')
        ax3.text(n_i[ind90percent] - 1.1, 2.35, '90%', color = 'b')
        ax3.text(n_i[ind90percent] - 1.1, 2.35 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind90percent]*100) + '%', 
                 color = 'r')

        ind50percent = argmin(fabs(f_CDF_n_i - 0.5))            
        ax3.plot([n_i[ind50percent], n_i[ind50percent]], [0, 2.2], 'b--', linewidth=1)
        ax3.text(n_i[ind50percent] - 1.1, 2.05, '50%' , color = 'b')
        ax3.text(n_i[ind50percent] - 1.1, 2.05 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind50percent]*100)  + '%', 
                 color = 'r')
       
        ind10percent = argmin(fabs(f_CDF_n_i - 0.1))
        ax3.plot([n_i[ind10percent], n_i[ind10percent]], [0, 1.85], 'b--', linewidth=1)
        ax3.text(n_i[ind10percent] - 1.1, 1.7 , '10%', color = 'b')
        ax3.text(n_i[ind10percent] - 1.1, 1.7 - 0.15, '%d' % (100.0 - f_CDF_em_i[ind10percent]*100) + '%', 
                 color = 'r')

       
        # plotting the 10, 50, 90% indicators for the emission CDF
        ind10percent = argmin(fabs(f_CDF_em_i - 0.1))
        ax3.plot([n_i[ind10percent], n_i[ind10percent]], [0, 1.55], 'r--', linewidth=1)
        ax3.text(n_i[ind10percent], 1.4, '10%', color = 'r')
        ax3.text(n_i[ind10percent], 1.25, '%d' % (100.0 - f_CDF_n_i[ind10percent]*100) + '%', 
                 color = 'b')


        ind50percent = argmin(fabs(f_CDF_em_i - 0.5))
        ax3.plot([n_i[ind50percent], n_i[ind50percent]], [0, 1.25], 'r--', linewidth=1)
        ax3.text(n_i[ind50percent], 1.1, '50%', color = 'r')
        ax3.text(n_i[ind50percent], 0.95, '%d' % (100.0 - f_CDF_n_i[ind50percent]*100) + '%', 
                 color = 'b')
       
        ind90percent = argmin(fabs(f_CDF_em_i - 0.90))
        ax3.plot([n_i[ind90percent], n_i[ind90percent]], [0, 0.95], 'r--', linewidth=1)
        ax3.text(n_i[ind90percent], 0.8, '90%', color = 'r')
        ax3.text(n_i[ind90percent], 0.65, '%d' % (100.0 - f_CDF_n_i[ind90percent]*100) + '%', 
                 color = 'b')

        ax3.set_ylim([0, 2.5])
        ax3.set_xlabel(r'$\log_{10}$[n$_{\rm gas}$]')
        ax3.set_ylabel(r'$f$')
        '''
    def match_weights(self):
        
        pylab.figure()
        
        nbins = 100
        N = len(self)
        
        ln = log(self.n)
        
        w = zeros(N, 'f8') 
        
        
        inds_o = self.get_inds_original_set()
        inds_s = self.get_inds_children()

        w[inds_o] = 1.0
        w[inds_s] = 0.0

        ## selecting the original gas particles
        gas_orig_amuse = self[self.get_inds_original_set()]
        gas_orig = gas_set(len(gas_orig_amuse)) 
        gas_orig.copy_attr_from_amuse_set(gas_orig_amuse, ['n'])
         
        fit_func_rng = [1e-2, 1e4]
        rv, N = gas_orig.fit_interval(fit_func_rng[0], fit_func_rng[1])
        
        hpdf_s, lbins_s = pylab.histogram(ln, range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                          weights=w,
                                         )
        pylab.plot(lbins_s[1::], log10(hpdf_s))

        x = linspace(log(1e-3), log(1e+6), 100)
        pylab.plot(x, log10(rv.pdf(x)), 'r--', lw=1)
        
        
        ###############
        ln_o, w_o = ln[inds_o], w[inds_o]
        ln_s, w_s = ln[inds_s], w[inds_s]
        
        w_o[:] = 1.0
        w_o[ ln_o > 5.4 ] = 0.0
        
        w_s[:] = 1.0/20.0
        match = True
        
        if match == False:
            hpdf_s, lbins_s = pylab.histogram(numpy.hstack((ln_o, ln_s)), 
                                              range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                              weights=numpy.hstack((w_o, w_s)),
                                             )
        else:
            
            t0 = time.time()
            
            
            intervals =  numpy.vstack((numpy.linspace(5.4, 8.4, 16), 
                                       numpy.linspace(5.4, 8.4, 16) + 0.2)).T
    
    
            for interval in intervals:
                
                print '#'*100
                print interval
                print '---------------------------------------'
                                
                ln_rng = interval
                
                wTrialRng = numpy.array([0.0, 1.0])
                nTrial = 0
                ln_x = ln_rng.mean()
                
                while True:
                    
                    ln_min, ln_max =  ln_rng
                    
                    w_s[ ((ln_s > ln_min)*(ln_s < ln_max)) ] = wTrialRng.mean()            
                    
                    hpdf_s, lbins_s = pylab.histogram(numpy.hstack((ln_o, ln_s)), 
                                                      range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                                      weights=numpy.hstack((w_o, w_s)),
                                                     )
                    f = interpolate.interp1d(lbins_s[1::], hpdf_s)
                
                    if numpy.fabs(1.0 - f(ln_x)/rv.pdf(ln_x)) < 0.01:
                        break
                    else:
                        
                        if f(ln_x) < rv.pdf(ln_x):
                            wTrialRng[0] = wTrialRng.mean()
                        else:
                            wTrialRng[1] = wTrialRng.mean()                      
                    print nTrial, wTrialRng
                    #if True:
                    #    break
                #
            
            print time.time() - t0
            print '---------------------------------------'
        #
        #1.0 - f(log(1e4))/rv.pdf(log(1e4))
        pylab.step(lbins_s[1::], log10(hpdf_s), 'g-', lw=2)
        
        pylab.show()
        
        #asdasdW
'''
        intervals =  numpy.vstack((numpy.linspace(5.4, 6.0, 5), numpy.linspace(5.6, 6.2, 5))).T
        
        for interval in intervals:
            
            print interval
            print '------------------------------------------'
            
            wTrialRng = numpy.array([0.0, 1.0])
            nTrial = 0
            
            while True:
                
                wTry = wTrialRng.mean()
                
                w_s[ ((ln_s >= interval[0])*(ln_s < interval[1])) ] = wTry
                
                hpdf_s, lbins_s = pylab.histogram(numpy.hstack((ln_o, ln_s)), 
                                                  range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                                  weights=numpy.hstack((w_o, w_s)),
                                                 )
                f = interpolate.interp1d(lbins_s[1::], hpdf_s)
            
                print nTrial, f(wTry), rv.pdf(wTry), 
                
                if numpy.fabs(1.0 - f(wTry)/rv.pdf(wTry)) < 0.05:
                    break
                else:
                    
                    if f(wTry) < rv.pdf(wTry):
                        wTrialRng[1] = wTry
                    else:
                        wTrialRng[0] = wTry
                                              
                print wTry, wTrialRng
                #if True:
                #    break
            #
'''        
        
        
        
        
        