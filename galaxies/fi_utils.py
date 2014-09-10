from amuse.datamodel import Particles
from amuse.io import read_set_from_file, fi_io
from amuse.units import units, constants, nbody_system
import collections
import constraining
import lineDict
import line_ratio_utils
import meshUtils
import ismUtils
import multiprocessing
from mylib.constants import M_SUN_SI, M_PROTON_CGS, LSUN_ERG_S, M_SUN_CGS
from mylib.units import KPC2CM
import mylib.units
from mylib.utils import templates
from mylib.utils.histogram import hist_nd
from mylib.utils.interpolation import sectioned_4D_interpolator
from mylib.utils.misc import find_matching_indicies, scale
from mylib.utils.ndmesh import ndmesh
import os
import pdb
import pdrDict
import shlex
import subprocess
import sys
import time

from matplotlib import ticker
from numpy import *
import numpy
from pylab import *
import pylab
from scipy import interpolate, stats
import scipy
from scipy.ndimage import zoom
from scipy.stats import chisqprob
#######################
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
    
    # computing the mean Av        
    gas_new.Pe = (1.085 * gas_new.T + 54.0 * gas_new.vdisp**2)*gas_new.n #how to set the unit to the same untis as  (units.K*gas.n.unit*constants.kB.unit)
    gas_new.Av = 0.22 * metallicity * ( 1520.0 / 100.0) * numpy.sqrt(gas_new.Pe/1e4)

    # copying some other attributes
    gas_new.id = numpy.int32(gas.id.value_in(gas.id.unit))
    gas_new.radius = gas.radius.value_in(units.parsec)
    
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

def compute_alpha(gas, arxvPDR, if_nan_set_nearest=None):
    """Computes the ratio of the mechanical heating of an SPH particle to the surface heating of
    a pure PDR with the same n,G0.
    """
    
    lgMechMin = arxvPDR.grid_z.min() #the minim mechanical heating  (in log) which will be assumed to be that of a pure PDR

    dataNew = numpy.array([
                           numpy.log10(gas.n), 
                           numpy.log10(gas.G0), 
                           numpy.ones(gas.n.shape)*lgMechMin
                          ]).T
    
    def interpolate_surface_heating(interpolator, data):                          
        print 'constructing the surface heating interpolation function from the PDR arxv'
        f_log_gamma_surf=arxvPDR.construct3DInterpolationFunction(quantity=['therm','heating'], 
                                                                  slabIdx=0, 
                                                                  log10=True, 
                                                                  interpolator=interpolator)
    
        print 'interpolating the surface heating at the surface of the SPH particles'
        gammaSurf_sph_from_pdr = 10.0**f_log_gamma_surf(data)
    
        return gammaSurf_sph_from_pdr
    #
    
    # computing the surface heating of the sph particles using linear interpolation
    gammaSurf_sph_from_pdr = interpolate_surface_heating('linear', dataNew)
    
    indsNan = where(numpy.isnan(gammaSurf_sph_from_pdr))[0]
        
    if  if_nan_set_nearest != True:
        if indsNan.size != 0:
            err_strng = 'interpolated variable gammaSurf_sph_from_pdr has %d nans values in it!!' % indsNan.size
            err_strng += '\n\t\tgMechMin = %e' % lgMechMin                 
            raise ValueError(err_strng)
    else:
        print '''!!!!!!!!! fi_utils.compute_alpha: found an nan in the interpolated surface heating 
        setting the interpolated surface heating to the nearest using nearest neighbour interpolation'''

        print 'number of points with Nan surface heating is ', indsNan.size
        
        # computing the surface heating of the sph particles using nearest neighbour interpolation
        gammaSurf_sph_from_pdr_nearest = interpolate_surface_heating('nearest', dataNew[indsNan])
        
        gammaSurf_sph_from_pdr[indsNan] = gammaSurf_sph_from_pdr_nearest 
    #
    
    gas.alpha = gas.gmech/gammaSurf_sph_from_pdr
    
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

    #setting emissions of particles with nan interpolated values to a very low value        
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
    for extra_attr in ['weights', 'parent', 'children', 'id', 'radius', 'Pe', 'vdisp']:
        if extra_attr in dir(gas):
            attr_save.append(extra_attr)
             
    save_gas_particle_info(filename, gas, attr_save)
    print '\t\t\twrote file : %s ' % filename
    #-----------------------------------------
    
    #the list of all the attributes    
    attr_list = gas.all_attributes()
    
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

def load_gas_particle_info(filename, load_pdr=None, load_only=None):
    '''load the particle states from the file'''
    
    data = numpy.load(filename)
    print '\tloaded file:\n\t\t%s' % filename
    n = data['arr_0'].size

    gas = gas_set(n)

    if load_only != None:
        names_to_be_loaded = numpy.array(load_only)
    else:
        names_to_be_loaded = data['names']         
    
    for i, name in enumerate(data['names']):
        
        if name not in names_to_be_loaded:
            print '\t\tattr_name = %-30s [%d elements] ... SKIPPING' % (name, data['arr_%d' %i].size)
            continue
        
        setattr(gas, name, data['arr_%d' % i])
        print '\t\tattr_name = %-30s [%d elements]' % (name, data['arr_%d' %i].size)
    
    all_names = names_to_be_loaded
    
    if load_pdr == True:
        filename_pdr = filename + '.pdr.npz'
        data_pdr = numpy.load(filename_pdr)
        print '\tloaded file:\n\t\t%s' % filename_pdr
        for i, name in enumerate(data_pdr['names']):
            setattr(gas, name, data_pdr['arr_%d' % i])
            print '\t\tattr_name = %-30s [%d elements]' % (name, data['arr_%d' %i].size)
            
        all_names = numpy.hstack((all_names, data_pdr['names']))

    return gas, all_names

def load_gas_particle_info_with_em(filename, species, load_pdr=None, load_only_em=None):
    '''loads the specified file foo.states.npz, also looks for files
     foo.states.npz.em_XYZ.npz, where XYZ is an entry in the list of strings in species
     
     load_only_em = ['flux_cgs_CO1-0', 'flux_cgs_CO2-1']
        or
     load_only_em = ['CO1-0', 'CO2-1']        
    '''
    
    names_all = []
    
    #loading the states of the SPH particles
    gas, names = load_gas_particle_info(filename, load_pdr=load_pdr)
    
    names_all.append(names.tolist())

    for specStr in species:
        
        filename_this = filename + '.em_%s.npz' % specStr
        
        if load_only_em != None:
            '''keeping the emission attributes for this specific species'''
            load_only_em_this_species = []
            for em_attr_name in load_only_em:
                line = em_attr_name.replace('em_','').replace('fluxKkms_','').replace('fluxcgs_','')
                spec_this_em_attr = lineDict.lines[line]['specStr']
                if spec_this_em_attr == specStr:
                    load_only_em_this_species.append(em_attr_name)
        else:
            load_only_em_this_species = None
            
        gas_this, names_this = load_gas_particle_info(filename_this, load_only=load_only_em_this_species)
        for name_this in names_this:            
            attr_data_this = getattr(gas_this, name_this)
            setattr(gas, name_this, attr_data_this)
    
        names_all.append(names_this.tolist())
        
        del gas_this
    
    print 'att attributes loaded are : ', names_all
    
    if 'radius' in gas.all_attributes():
        print '!'*100; print '!'*100; print '!'*100;
        print '!!!! warning : converting loaded radii from pc to kpc'
        print '!'*100; print '!'*100; print '!'*100;
        gas.radius = gas.radius * 1e-3
    
    print 'warning: make sure the computed optical depths make sense!!! it doesnt seem so..'
    return gas


def make_map(gas, hist, attr=None, as_log10=None, func=None, show=False, in_ax=None, **kwargs):
    '''looping over the bins of the 2D histogram of the x,y coordinates and 
     computing the averages of the maps
    '''

    # the area of the binx (area of each pixel)
    bin_area = hist.szBins.prod()
    
    map_data = zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)

    if attr not in dir(gas[0]):  
        raise AttributeError('%s is not an attribute of the gas particles,' % attr)
        
    attr_data = getattr(gas, attr)
    
    mass = gas.mass
    radius = gas.radius
    
    for i in numpy.arange(hist.nBins[0]):
        
        for j in numpy.arange(hist.nBins[1]):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            map_data[i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:

                q = attr_data[inds_in_bin]
                r = radius[inds_in_bin]
                m = mass[inds_in_bin]
                 
                if func == total_luminosity:
                    bin_value = total_luminosity(fluxes=q, radii=r)
                if func == mean_flux:
                    bin_value = mean_flux(fluxes=q, radii=r, bin_area=bin_area)
                
                map_data[i,j] = bin_value
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
        
        if in_ax == None:
            pylab.figure()
            pylab.imshow(data_show, origin='lower', **kwargs)
            pylab.colorbar()
        else:
            in_ax.imshow(data_show, origin='lower', **kwargs)
        
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
    
    # the generated interpolation functions can be used in the following way:
        fInterp(array([[3.0, 3.0, -50.0, 10.0]]))

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
    compute_alpha(gas, arxvPDR, if_nan_set_nearest=True)
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
                 line_ratios=None, arxvPDR=None, em_unit='Kkms', **kwargs):
    
        self.luminosity_maps = luminosity_maps #: the dict containing the luminosity maps
        self.gas = gas #: the object of gas particles
        self.params = params #: the parameters used in controlling everything
        self.hist = hist
        self.line_ratios = line_ratios
        self.arxvPDR = arxvPDR
        self.em_unit = em_unit  # the unit of the emission of the gas to be used
        
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
        
        self.gas_in_pixel = None
        
    def get_model_emission_from_pdr_arxv_involved_line_ratios(self, em_unit='Kkms'):
        '''loading the emission info from all the models for all Avs (also we check for 
        the consistenscy of the number of models...i.e same number of models 
        for all the lines)
        '''
        
        
        model_em, grid_coords = self.arxvPDR.get_model_emission_from_involved_line_ratios(self.line_ratios,     
                                                                                          em_unit=em_unit)
        
        self.model_em = model_em
        self.grid_coords = grid_coords
        
    def setup_observed_grid(self, use_line_ratios=None, clean=False, inspect=''):
        '''sets up the observed mesh. A lower resolution mesh than the synthetica map'''
        
        self.clean = clean
        self.inspect_species = inspect
        
        if use_line_ratios != None:
            self.line_ratios = use_line_ratios
            
        bs_min, bs_max = self.params['ranges']['box_size'].number

        fig = pylab.figure(0, (8,8))
        
        axs = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        self.fig, self.axs = fig, axs
        
        pylab.imshow(log10(self.luminosity_maps['maps'][self.luminosity_maps['maps'].keys()[0]]),
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

        if self.axs.contains(event)[0]:
           
            print 'cursor x,y = ', xd, yd

            ## plotting the centroid of the pixel clicked in the coarse grid 
            x_ind = numpy.floor(scale(xd, 0.0, obs_res, bs_min, bs_max))
            y_ind = numpy.floor(scale(yd, 0.0, obs_res, bs_min, bs_max))
        
            if event.button == 1:
        
                print 'left click...'
        
                self.estimate_mass_in_pixel(x_ind, y_ind, plot_fits=True)
                
            if event.button == 3:

                print 'right click...'
                
                self.plot_gas_pdf_in_pixel(x_ind, y_ind)
                
    
    def estimate_mass_in_pixel(self, x_ind=None, y_ind=None, plot_fits=False):
        
        t0 = time.time()
        
        hist = self.hist
        hist_obs = self.hist_obs
        obs_mesh = self.hist_obs.f
        line_ratios = self.line_ratios
        luminosity = self.luminosity_maps
        params = self.params
        model_em = self.model_em 
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
            pylab.plot(hist.f.cntrd[0][inds], hist.f.cntrd[1][inds], 'w+', zorder=2)
        
        inds_gas = hist_obs.get_indicies([x_ind, y_ind])

        gas_in_pixel = self.gas[inds_gas]

        if plot_fits == True:
            
            pylab.plot(gas_in_pixel.x, gas_in_pixel.y, 'o', markersize=2, zorder=1)
         
            pylab.draw()
                
        ## make line ratios from the mock luminosities 
        obs_mock_ratios = line_ratio_utils.ratios()                

        for line_ratio in line_ratios:
        
            line1, line2 = line_ratio_utils.lines_involved(line_ratio)
            
            print line1, line2

            flux1 = gas_in_pixel.get_mean_flux(line1, em_unit=self.em_unit, bin_area=area_obs_pixel)
            flux2 = gas_in_pixel.get_mean_flux(line2, em_unit=self.em_unit, bin_area=area_obs_pixel)
            
            obs_mock_ratios.make_ratios(
                                        {
                                          line1:{'flux' + self.em_unit: flux1, 'err': params['error_bars']*flux1}, 
                                          line2:{'flux' + self.em_unit: flux2, 'err': params['error_bars']*flux2}
                                        },
                                        ratios = [line_ratio],
                                        em_unit = 'flux' + self.em_unit,
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
            
            fig = pylab.figure(figsize=(16,10))
            
            # plotting the line ratios in these axes
            ax1 = fig.add_axes([0.05, 0.55, 0.25, 0.4])
            ax2 = fig.add_axes([0.35, 0.55, 0.25, 0.4])

            # plotting the line ratios
            f.plot_results(fig = fig, ax = ax1)
            #print '-------------------------'
            f.plot_results(fig = fig, ax = ax2, no_gmech=True)
            
            # plotting the PDFs in these axes
            axes = numpy.ndarray((2,4), 'object')
            ys = 0.02
            xsz, ysz = 0.2, 0.2
            dx = 0.02
            
            axes[0,0] = fig.add_axes([0.05, 0.25, xsz, ysz]) 
            axes[1,0] = fig.add_axes([0.05, ys, xsz, ysz]) 

            axes[0,1] = fig.add_axes([0.05 + xsz + dx, 0.25, xsz, ysz]) 
            axes[1,1] = fig.add_axes([0.05 + xsz + dx, ys, xsz, ysz]) 

            axes[0,2] = fig.add_axes([0.05 + 2*(xsz + dx), 0.25, xsz, ysz]) 
            axes[1,2] = fig.add_axes([0.05 + 2*(xsz + dx), ys, xsz, ysz]) 

            axes[0,3] = fig.add_axes([0.05 + 3*(xsz + dx), 0.25, xsz, ysz]) 
            axes[1,3] = fig.add_axes([0.05 + 3*(xsz + dx), ys, xsz, ysz]) 
            pylab.show()
            
            #plotting the emission distribution in the pixel
            if self.inspect_species != '':
                gas_in_pixel.get_emission_pdfs(qxs=['n', 'G0', 'gmech', 'Av'],
                                               line=self.inspect_species,
                                               log10xs=[True, True, True, False],
                                               xrngs=[[-3.0, 6.0], [-3.0,6.0], [-30.0, -20.0], [0.0, 27.0]],
                                               in_axes=axes, 
                                               em_unit=self.em_unit.replace('flux',''),
                                               )
            else:
                print '#'*100
                print "self.inspect_species == '', nothing to inspect"
                print '#'*100
                                
        #computing the averages of the physical parameters of the particles inside the pixel 
        #-----------------------------------------------------------------------------------
        def get_pixel_mean_info(weights=None):
            
            means = collections.OrderedDict()
            sdevs = collections.OrderedDict()
            
            if weights == None:
                
                means['n']  = log10(gas_in_pixel.n.mean())
                means['G0'] = log10(gas_in_pixel.G0.mean())
                means['gm'] = log10(gas_in_pixel.gmech.mean())
                means['Av'] = gas_in_pixel.Av.mean()
                
            else:

                means['n']  = log10(average(gas_in_pixel.n, weights=gas_in_pixel.weights))
                means['G0'] = log10(average(gas_in_pixel.G0, weights=gas_in_pixel.weights))
                means['gm'] = log10(average(gas_in_pixel.gmech, weights=gas_in_pixel.weights))
                means['Av'] = average(gas_in_pixel.Av, weights=gas_in_pixel.weights)
                
                sdevs['n']  = log10(sqrt(average((gas_in_pixel.n     - 10.0**means['n'])**2 , weights=gas_in_pixel.weights)))
                sdevs['G0'] = log10(sqrt(average((gas_in_pixel.G0    - 10.0**means['G0'])**2, weights=gas_in_pixel.weights)))
                sdevs['gm'] = log10(sqrt(average((gas_in_pixel.gmech - 10.0**means['gm'])**2, weights=gas_in_pixel.weights)))
                sdevs['Av'] = log10(sqrt(average((gas_in_pixel.Av    - means['Av'])**2      , weights=gas_in_pixel.weights)))

            return means, sdevs
        
        means, sdevs = get_pixel_mean_info(weights=True)
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

        def gen_fit_string():
            
            strng = 'fit parameters\n' 
            strng += '     with gmech : '
            strng += r'$\log_{10}$[n$_{\rm gas}]$ = %-+.2f' % log10_n + ' '
            strng += r'$\log_{10}$[G$_0]$ = %-+.2f' % log10_G0 + ' '
            strng += r'$\log_{10}$[$\Gamma_{\rm mech}$] = %-+.2f' % log10_gmech + ' '
            strng += r'$A_V$ = %-+.2f' % Av + '\n'

            strng += 'without gmech : '
            strng += r'$\log_{10}$[n$_{\rm gas}]$ = %-+.2f' % log10_n_no_gmech + ' '
            strng += r'$\log_{10}$[G$_0]$ = %-+.2f' % log10_G0_no_gmech + ' '
            strng += r'$\log_{10}$[$\Gamma_{\rm mech}$] = %-+.2f' % log10_gmech_no_gmech + ' '
            strng += r'$A_V$ = %-+.2f' % Av_no_gmech + '\n'
            
            return strng
        
        def gen_pixel_means_string(means):
            
            strng  = '    pixel means : '

            strng += r'$\log_{10}$[n$_{\rm gas}]$ = %-+.2f' % means['n'] + ' '
            strng += r'$\log_{10}$[G$_0]$ = %-+.2f' % means['G0'] + ' '
            strng += r'$\log_{10}$[$\Gamma_{\rm mech}$] = %-+.2f' % means['gm'] + ' '
            strng += r'$A_V$ = %-+.2f' % means['Av'] + '\n'
            
            return strng
        
        if plot_fits == True:
            
            strng = gen_fit_string()            
            pylab.figtext(0.61, 0.8, strng, size=10)
            
            # the reduced Xi2
            n_lines = float(line_ratio_utils.all_lines_involved(self.line_ratios).size)
            
            DOF = n_lines - 4.0
            p_value = chisqprob(Xi2_min, DOF)
            ax1.set_title(r'with gmech ' + r'$\chi^2 = %.3f$ p-value = %.3f' % (Xi2_min, p_value))
            
            DOF = n_lines - 3.0
            p_value_no_gmech = chisqprob(Xi2_min_no_gmech, DOF)            
            ax2.set_title(r'no gmech ' + r'$\chi^2 = %.3f$ p-value = %.3f' % (Xi2_min_no_gmech, p_value_no_gmech))

            strng = 'pixel coordinates = [%d, %d]' % (x_ind - self.hist_obs.f.shape[0]/2, 
                                                      y_ind - self.hist_obs.f.shape[1]/2)
            strng += ' nLines = %d' % n_lines            
            pylab.figtext(0.61, 0.9, strng, size=10)                        
            
            strng = gen_pixel_means_string(means)
            pylab.figtext(0.61, 0.78, strng, size=10)
            
        print '-'*100
        print 'mean values in the pixel'    
        print means 
        print 'standard deviations in the pixel'    
        print sdevs 
        print '-'*100
        
        print 'line ratios of the observed and the best fit model (whichever with or without gmech)'
        f.get_model_em_for_min_Xi2()
        self.contraining = f
        self.means = means
        
        self.gas_in_pixel = gas_in_pixel
        
        if self.clean == True:
            del self.gas_in_pixel
        
        print 'time processing = %.2f seconds' % (time.time() - t0)
        return means, sdevs, PDR_parms, Xi2_min, PDR_parms_no_gmech, Xi2_min_no_gmech, n_lines

    def plot_gas_pdf_in_pixel(self, x_ind, y_ind):
        
        print 'disceting the gas'

        hist = self.hist
        hist_obs = self.hist_obs
        obs_mesh = self.hist_obs.f
        params = self.params
        line_ratios = self.line_ratios    # line ratio strings
        params = self.params
        gas = self.gas
        ######################################################
        
        area_obs_pixel = numpy.product(obs_mesh.dl)

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

        pylab.plot(hist.f.cntrd[0][inds], hist.f.cntrd[1][inds], 'c+')
        
        inds_gas = hist_obs.get_indicies([x_ind, y_ind])

        pylab.plot(gas.x[inds_gas], gas.y[inds_gas], 'o')
        
        pylab.draw()
        
        gas_in_pixel = gas[inds_gas].copy()

        specs = gas_in_pixel.get_species_with_em()
        
        #weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % params['snap_index']
        #gas_in_pixel.use_weights(weighting='original-only', weights_filename = weights_filename)
        #gas_in_pixel.use_weights(weighting='by-number', weights_filename = weights_filename)
        #gas_in_pixel.use_weights(weighting='mathced', weights_filename = weights_filename)
        
        if False:

            fig = pylab.figure()
            ax = pylab.subplot(111)
            
            for spec in specs:
                
                print 'getting the emission of %s' % spec
                
                idx, total_lum = gas_in_pixel.get_total_luminosity_ladder(spec, em_unit=self.em_unit)
            
                flux = total_lum / area_obs_pixel
                
                plt, = ax.semilogy(idx+1, flux, label=spec)
            
            ax.legend(loc=0)
            ax.set_xlabel(r'$J_{\rm up}$')
            ax.set_ylabel('Flux')
        
        self.gas_in_pixel = gas_in_pixel
        
        line = 'HCN1-0'
         
        #gas_in_pixel.get_emission_pdf(qx='n', line='CO1-0', log10x=True, nbins=100, xrng=[-3, 6])
        #ax1, ax2 = gas_in_pixel.get_emission_pdf(qx='n', line=line, log10x=True, nbins=50, xrng=[-3, 6], wtitle='sampled')
        #ax1, ax2 = gas_in_pixel.get_emission_pdf(qx='G0', line=line, log10x=True, nbins=50, xrng=[-3, 6], wtitle='sampled')

        gas_in_pixel.get_emission_pdfs(qxs=['n', 'G0', 'gmech', 'Av'],
                                       line=line,
                                       log10xs=[True, True, True, False],
                                       xrngs=[[-3.0, 6.0], [-3.0,6.0], [-30.0, -20.0], [0.0, 27.0]],
                                       )
        if False:
            self.gas_in_pixel_original = gas_in_pixel.copy()
            self.gas_in_pixel_original.set_radii(weighting='original-only', 
                                                 rundir=self.params['rundir'], 
                                                 snap_index=self.params['snap_index'])
            
            self.gas_in_pixel_original.get_emission_pdf(qx='n', line=line, log10x=True, nbins=50, 
                                                       xrng=[-3, 6], wtitle='original', in_ax1=ax1, in_ax2=ax2,
                                                       linestyle='--')
        #pylab.figure()
        #gs = self.gas_in_pixel
        #pylab.loglog(gs.n[gs.get_inds_has_children()], gs.radius[gs.get_inds_has_children()], '.')

    def estimate_mass_in_some_pixels(self):
        
        info_all = []
        
        #inds = [[0,0]]
        #inds = [[0,1], [0,-1], [1,0], [-1,0]]
        #inds = [[-1,1], [-1,-1], [1,-1], [1,1]]
        inds = [
                [-1, 2], [ 0, 2], [ 1, 2],
                [-1,-2], [ 0,-2], [ 1,-2],
                [ 2, 1], [ 2, 0], [ 2,-1],
                [-2, 1], [-2, 0], [-2,-1],
                ]

        xc, yc = array(self.hist_obs.f.shape)/2
        for i, j in inds:
            x_ind , y_ind = xc + i, yc + j
            info = self.estimate_mass_in_pixel(x_ind=x_ind, y_ind=y_ind, plot_fits=True)
            info_all.append(info)
            
        #means, sdevs, PDR_parms, Xi2_min, PDR_parms_no_gmech, Xi2_min_no_gmech, n_lines = info

        if False:    
            i_fit = 2  # select the mPDR fit info
            pfit = 4.0
            i_Xi2 = 3
        if True:
            i_fit = 4  # select the pure PDR fit info
            pfit  = 4.0
            i_Xi2 = 5
        
        print 'k    [i, j]     r      log10n  log10G0  log10Gmech  Av    Xi2      Xi2_red      p-value'
        for k, [i,j] in enumerate(inds):
            x_ind , y_ind = xc + i, yc + j
            r = sqrt(self.hist_obs.f.cntrd[0][x_ind,y_ind]**2 + 
                     self.hist_obs.f.cntrd[1][x_ind,y_ind]**2 )
            print '%d    %-+d,%-+d  %f' % (k, i, j, r),
            
            # 0       1       2         3           4                   5             6        
            #means, sdevs, PDR_parms, Xi2_min, PDR_parms_no_gmech, Xi2_min_no_gmech, n_lines = info_all[k]

            means = info_all[k][0]
            Xi2 = info_all[k][i_Xi2]
            Xi2_red = info_all[k][i_Xi2]//(info[6]-pfit)
            print '   %5.2f  %5.2f      %5.2f   %5.2f' % (means['n'],  means['G0'], means['gm'], means['Av']),
            print ' %5.2f    %5.2f' % (Xi2, Xi2_red)

        print '-'*100
        
        print 'averages pixel          %5.2f  %5.2f      %5.2f   %5.2f ' % (
                                                                            mean([info[0]['n']  for info in info_all]),     
                                                                            mean([info[0]['G0']  for info in info_all]), 
                                                                            mean([info[0]['gm']  for info in info_all]), 
                                                                            mean([info[0]['Av']  for info in info_all]), 
                                                                            )
    
            
        print 'averages fit            %5.2f  %5.2f      %5.2f   %5.2f  %5.2f    %5.2f' % (
                                                                                           mean([info[i_fit][0] for info in info_all]), 
                                                                                           mean([info[i_fit][1] for info in info_all]), 
                                                                                           mean([info[i_fit][2] for info in info_all]), 
                                                                                           mean([info[i_fit][3] for info in info_all]), 
                                                                                           mean([info[i_Xi2] for info in info_all]),
                                                                                           mean([info[i_Xi2]/(info[6]-pfit) for info in info_all]),
                                                                                           )
        return info_all
                
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
        
class gas_set(Particles):


    def __getitem__(self, items):
        '''overriding the slicing method of the amuse particle set class'''
        
        n = getattr(self, 'key')[items].size 
        
        gas_ret = gas_set(n)
        
        for attr_name in self.all_attributes():
            
            if attr_name in ['connected_components', 'mass_segregation_Gini_coefficient', 'total_mass',
                             'bound_subset', 'LagrangianRadii', 'acceleration', 'thermal_energy',
                             'total_angular_momentum', 'get_binaries', 'angular_momentum', 'potential',
                             'Qparameter', 'potential_energy', 'center_of_mass', 'potential_energy_in_field',
                             'moment_of_inertia', 'cluster_core', 'total_momentum', 
                             'new_particle_from_cluster_core', 'position', 'kinetic_energy',
                             'densitycentre_coreradius_coredens', 'oblateness', 'rotate', 'move_to_center',
                             'center_of_mass_velocity', 'virial_radius', 'total_radius', 'scale_to_standard',
                             'find_closest_particle_to', 'velocity', 'binaries', 'specific_kinetic_energy']:
                continue
            
            data = getattr(self, attr_name)
            
            if type(getattr(self, attr_name)) == type(numpy.ndarray([])):
                
                setattr(gas_ret, attr_name, data[items])
                        
        return gas_ret

    def copy(self):
        '''overriding the slicing method of the amuse particle set class'''
        
        n = len(self) 
        
        gas_ret = gas_set(n)
        
        for attr_name in self.all_attributes():
            
            if attr_name in ['connected_components', 'mass_segregation_Gini_coefficient', 'total_mass',
                             'bound_subset', 'LagrangianRadii', 'acceleration', 'thermal_energy',
                             'total_angular_momentum', 'get_binaries', 'angular_momentum', 'potential',
                             'Qparameter', 'potential_energy', 'center_of_mass', 'potential_energy_in_field',
                             'moment_of_inertia', 'cluster_core', 'total_momentum', 
                             'new_particle_from_cluster_core', 'position', 'kinetic_energy',
                             'densitycentre_coreradius_coredens', 'oblateness', 'rotate', 'move_to_center',
                             'center_of_mass_velocity', 'virial_radius', 'total_radius', 'scale_to_standard',
                             'find_closest_particle_to', 'velocity', 'binaries', 'specific_kinetic_energy']:
                continue
            
            data = getattr(self, attr_name).copy()
            
            if type(getattr(self, attr_name)) == type(numpy.ndarray([])):
                
                setattr(gas_ret, attr_name, data)
                        
        return gas_ret

    def delete(self):            
    
        for attr_name in self.all_attributes():
            
            if attr_name in ['connected_components', 'mass_segregation_Gini_coefficient', 'total_mass',
                             'bound_subset', 'LagrangianRadii', 'acceleration', 'thermal_energy',
                             'total_angular_momentum', 'get_binaries', 'angular_momentum', 'potential',
                             'Qparameter', 'potential_energy', 'center_of_mass', 'potential_energy_in_field',
                             'moment_of_inertia', 'cluster_core', 'total_momentum', 
                             'new_particle_from_cluster_core', 'position', 'kinetic_energy',
                             'densitycentre_coreradius_coredens', 'oblateness', 'rotate', 'move_to_center',
                             'center_of_mass_velocity', 'virial_radius', 'total_radius', 'scale_to_standard',
                             'find_closest_particle_to', 'velocity', 'binaries', 'specific_kinetic_energy']:
                continue
            
            data = getattr(self, attr_name)
            
            if type(getattr(self, attr_name)) == type(numpy.ndarray([])):
            
                del data
                print '\t deleting attribute %s' % attr_name
                        
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

        print 'log10(mu), log10(std) = ', log10(exp(mu_fit)),  log10(exp(sigma_fit))

        ## generating the function with the fitted deviates
        rv = stats.norm(loc=mu_fit, scale=sigma_fit)
        
        return rv, ln_fit.size
    
    def plot_PDF_CDF(self, fit_func_rng=None, in_ax=None, color='r'):
        
        if in_ax == None:
            fig, axs = pylab.subplots(2, 3, figsize=(12,8))
        else:
            axs = in_ax
            fig = pylab.gcf()
    
        
        # getting the PDF and the CDF of all the particles
        l10min, l10max = -6, 6
        nBins = 100
        
        bin_sz = (l10max - l10min)/float(nBins) 
        
        l10n = log10(self.n) 
        hpdf, bins = pylab.histogram(l10n, range=[l10min, l10max], bins=nBins, normed=True)
        
        # plotting the PDF of log10(n)
        axs[0,0].plot(bins[1::], hpdf, color)
        axs[0,0].set_title('PDF log10(n)')
        
        axs[0,1].semilogy(bins[1::], (hpdf*bin_sz).cumsum(), color)
        axs[0,1].set_title('CDF log10(n)')
        axs[0,1].set_ylim(ymax=10.0)

        # plotting PDF(log(x))
        ln = log(self.n) 
        hpdf, lbins = pylab.histogram(ln, range=[log(1e-6), log(1e6)], bins=nBins, normed=True)
        axs[0,2].semilogy(lbins[1::], hpdf, color, linewidth=2)
        axs[0,2].set_title('PDF log(n)')

        ########################################################
        #####################  just plotting ################### 
        ########################################################
        
        #ranges = [ [1e-1, 1e1], [1e-2, 1e2], [1e-1, 1e2], [1e-1, 1e3], [1e-1, 1e4] ]
        #colors = [    'y'     ,   'g',         'k',          'c',        'm']
        ranges = [ fit_func_rng,  ]
        colors = [    color+':'     , ]
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
    
    def sample_ln_densities(self, npp, ln_this, rem, rv):
        '''returns a sample of npp size of sampled densities in ln() scale, rv is the scipy 
        generator. The sample size from which npp particles are picked are based on the remaining
        number of particles to be sampled 'rem'.
        '''
        
        ## generating (1 - cdf(x))*npp points  
        n_generated = npp*(1.0/(1.0 - rv.cdf(ln_this)))

        ln_trial = rv.rvs(n_generated)
        
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
    
    def sample_set_serial(self, ln_gt_X, ln_Y, rv, npp):
        '''sample the particles whose densities are ln_gt_X one at a time. Keeping sampled only particles
        whose density is less than Y.
        '''

        #ln_gt_X = ln_gt_X[-600:-499:]  ## testing on a subset
        
        ## array which will hold the sampled densities
        ln_s = numpy.array([],'f8')
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
                ln_trial = self.sample_ln_densities(npp, ln_this, rem, rv)
                ts += (time.time() - ts1)

                if attempts == 0:
                    print '\t generated a pool of %e random numbers' % ln_trial.size

                ## keeping sampled particles whose density is higher than the current one             
                ln_trial_gt_ln_this = ln_trial[ (ln_trial > ln_this)*(ln_trial < ln_Y)]
                ## keeping all particles
                #ln_trial_gt_ln_this = ln_trial[:]
    
                attempts += 1

                if  ln_trial_gt_ln_this.size > npp:
                    
                    sample_keep = ln_trial_gt_ln_this[:npp]
    
                    ## setting weights by volume to the sampled particles
                    #w_this_sample = (1.0/float(npp+1.0))*(exp(ln_this)/exp(sample_keep))
                    ## setting zero weights
                    #w_this_sample = zeros(npp, 'f8')
                    ## setting unit weights
                    #w_this_sample = numpy.ones(npp, 'f8')
                    ## setting weights by mass
                    w_this_sample = numpy.ones(npp, 'f8')/float(npp+1)
                    #((1.0/float(npp))*(exp(ln_this)/exp(sample_keep)))**(2.0/3.0),
                    #(1.0/float(npp))*numpy.ones(npp, 'f8')/10.0,
                    
                    ln_s = numpy.hstack((ln_s, sample_keep))
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
    
        return ln_s, w_s, w_ln_gt_X
    
    def sample_higher_densities(self, npp=None, n_min_sample=None, n_max_sample=None, 
                                fit_func_rng=None, plot=True):

        if plot == True:
            fig, axs = self.plot_PDF_CDF(fit_func_rng=fit_func_rng)
        
        ## picking one gaussian to be used in the sampled
        rv, N = self.fit_interval(fit_func_rng[0], fit_func_rng[1])

        ## testing the PDF
        #x = 2.0
        #y = (1.0/(sqrt(2.0*pi)*rv.std()))*exp( -0.5*((x - rv.mean())/rv.std())**2.0 )
        #print rv.pdf(x), y

        ln = log(self.n)
        
        ## particles whose density is larger than X and less than Y will be sampled 
        X = n_min_sample
        Y = n_max_sample 
        
        # selecting the particles within the density range to be sampled
        inds = numpy.where((ln > log(X))*(ln < log(Y)))
        ln_gt_X = ln[inds]
        gas_ln_gt_X = self[inds] 
        
        # selecting the partilces whose density is less than X
        inds = where(ln < log(X))
        ln_lt_X = ln[ inds ]
        gas_ln_lt_X = self[inds]
        
        # soring the densities and the corresponding keys in ascending order
        inds_sorted = numpy.argsort(ln_gt_X)
        ln_gt_X = ln_gt_X[inds_sorted]
        gas_ln_gt_X = gas_ln_gt_X[inds_sorted]

        n_to_be_sub_sampled = ln_gt_X.size
        print 'number of SPH particles to be sampled = %d' % n_to_be_sub_sampled
        print 'these are %.2f percent of the original set.' % (float(n_to_be_sub_sampled)/float(ln.size))
        
        ## doing the actual sampling        
        ln_s, w_s, w_ln_gt_X = self.sample_set_serial(ln_gt_X, log(Y), rv, npp)
        
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

    def get_inds_of_has_children_with_corresponding_children(self):
        '''returns in the indicies of the parents, and the indicies of the children
        
        p, c = gas.get_inds_of_has_children_with_corresponding_children()
        
        ## chaning the values of the parents
        ids = gas.ids
        
        ids[p] = arange(p.size)
        
        ## setting the  
        '''
        
        asdadasd 
        ## the indicies of the particles which have children and the children in the curernt set
        inds_parents = self.get_inds_has_children()
        inds_children = self.get_inds_children()

        ## the ids of the parents and their children 
        ids_parents = self.id[inds_parents]
        ids_children = self.id[inds_children]
        
        ## getting the mask of the parents in the children id list 
        mask = numpy.in1d(ids_children, ids_parents)
        inds = inds_children[numpy.where(mask)]

        return inds_parents, inds
        
    def get_inds_original_set(self):
        '''returns the indicies of the original set of the particles'''
        
        return numpy.where(numpy.isnan(self.parent))[0]

    def get_inds_has_children(self):
        '''returns the indicies of the original set of the particles which have been sampled,
        i.e the particles which have children
        '''
        
        return numpy.where(numpy.isfinite(self.children))[0]

    def get_emission_pdf_populations(self, qx=None, line=None, log10x=False, nbins=50, xrng=None, 
                                     wtitle='', in_ax1=None, in_ax2=None, linestyle='-', em_unit='fluxKkms'):
        '''getting emission distributions of the sub-populations
        
        get_emission_pdf_populations(qx='n', line=line, log10x=True, nbins=50, xrng=[-3, 6])

        '''

        inds = self.get_inds_original_set()
        gas_pop = self.copy()[inds]
        ax1, ax2 = gas_pop.get_emission_pdf(qx=qx, line=line, log10x=log10x, nbins=nbins, xrng=xrng, 
                                            wtitle=wtitle, linestyle='-')

        inds = self.get_inds_has_children()
        gas_pop = self.copy()[inds]
        gas_pop.get_emission_pdf(qx=qx, line=line, log10x=log10x, nbins=nbins, xrng=xrng,       
                                 wtitle=wtitle, linestyle='-o', in_ax1=ax1, in_ax2=ax2)

        inds = self.get_inds_children()
        gas_pop = self.copy()[inds]
        ax1, ax2 = gas_pop.get_emission_pdf(qx=qx, line=line, log10x=log10x, nbins=nbins, xrng=xrng, 
                                            wtitle=wtitle, linestyle=':', in_ax1=ax1, in_ax2=ax2)

        return ax1, ax2

    def get_emission_pdfs(self, qxs=None, line=None, log10xs=False, nbins=50, xrngs=None, in_axes=None, em_unit='Kkms'):
        '''plots the PDFs as a function of multiple quantities
        
        gas.get_emission_pdfs(qxs=['n', 'G0', 'gmech', 'Av'],
                              line='CO1-0',
                              log10xs=[True, True, True, False],
                              xrngs=[[-3.0, 6.0], [-3.0,6.0], [-30.0, -20.0], [0.0, 30.0]],)  
        '''
        
        if in_axes == None:
            fig, axs = pylab.subplots(2, 4, figsize=(12, 6))
        else:
            axs = in_axes

        for i, qx in enumerate(qxs):
            
            self.get_emission_pdf(qx=qx, line=line, log10x=log10xs[i], nbins=nbins, xrng=xrngs[i], 
                                  in_ax1=axs[0,i], in_ax2=axs[1,i], em_unit=em_unit)
            
    def get_emission_pdf(self, qx=None, line=None, log10x=False, nbins=50, xrng=None, wtitle='', 
                         in_ax1=None, in_ax2=None, linestyle='-', em_unit='Kkms', legend=True,
                         color1='b', color2='r', reverse=True):
        
        '''getting emission PDF vs a quantity qx
        
        example.
        l10n, emPDF = gas.get_emission_pdf(qx='n', line='13CO1-0', log10x=True, nbins=50, xrng=[-3, 6])

        '''

        if in_ax1 == None and in_ax2 == None:
            pylab.figure()
        
        N = len(self) 
        
        print 'number of particles = ', N
        
        ## getting the data to be processed from the attributes        
        x = getattr(self, qx)
        y = self.get_luminosity(line, em_unit)
        if log10x == True: x = log10(x)
        
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
        w_in_hist = self.weights
        
        #getting the distribution of the y variable as a function of the histogram in x
        xw_dist = zeros(x_dist.totalBins, 'f8')
        y_dist  = zeros(x_dist.totalBins, 'f8')

        for i in numpy.arange(x_dist.totalBins):
           
            inds_in_bin = x_dist.get_indicies(i)
           
            if inds_in_bin.size != 0:
                
                #the weighted distribution in x
                xw_dist[i] = w_in_hist[inds_in_bin].sum() 
                
                # sum of the qy in that bin
                y =  y_in_hist[inds_in_bin].sum()
                
                y_dist[i] = y
            #
        #

        if in_ax1 == None and in_ax2 == None:
            ax1 = pylab.subplot(211)
            ax2 = pylab.subplot(212)
        else:
            ax1 = in_ax1
            ax2 = in_ax2

        ax1.plot(x_dist.f.cntrd, xw_dist/xw_dist.max(), color1, linestyle=linestyle, label=qx, 
                 drawstyle='steps-mid')
        ax1.plot(x_dist.f.cntrd, y_dist/y_dist.max(), color2, linestyle=linestyle, 
                 label=line + '\n max = %.5f' % y_dist.max(), 
                 drawstyle='steps-mid')
        ax1.set_ylim([0, 1.1])
        #ax1.set_title(' PDF')
        if legend == True: ax1.legend(loc=0, prop={'size':8})
        
        x = x_dist.f.cntrd
        y1 = (xw_dist.cumsum())/xw_dist.sum()
        y2 = (y_dist.cumsum())/y_dist.sum()

        if reverse == True:
            x = x[::-1]
            y1 = y1[::-1]
            y2 = y2[::-1]
              
        ax2.plot(x, y1, color1, linestyle=linestyle, label=qx, drawstyle='steps-mid')
        ax2.plot(x, y2, color2, linestyle=linestyle, label=line, drawstyle='steps-mid')
        ax2.set_ylim([0, 1.1])
        #ax2.set_title(' CDF')
        if legend == True: ax2.legend(loc=0, prop={'size':8})
        
        pylab.gcf().canvas.set_window_title(wtitle)
        pylab.show()
    
        return ax1, ax2
    
    def set_weights_children_zero_parents_one(self):
        
        w = numpy.zeros(len(self), 'f8')
        
        inds_o = self.get_inds_original_set()
        inds_s = self.get_inds_children()

        w[inds_o] = 1.0
        w[inds_s] = 0.0
        
        self.weights = w
    
    def set_weights_childrent_to_inverse_sampling_number(self, sampling_info_fname=None):

        w = numpy.zeros(len(self.weights), 'f8')

        ## loading the sampling info from the saved file
        info = self.load_weights_function(sampling_info_fname)
        
        ## indicies of the original particles    
        inds1 = self.get_inds_original_set()
        w[inds1] = 1.0
        
        ## indicies of the original particles which have children
        inds2 = self.get_inds_has_children()
        w[inds2] = 1.0 / (1.0 + info['npp'])
        
        ## indicies of the points which are sampled
        inds3 = self.get_inds_children()        
        w[inds3] = 1.0 / (1.0 + info['npp'])
        
        self.weights = w
        
        
    def match_weights(self, npp=None, fit_func_rng=None, n_min_sample=None, 
                      match_interval=None, nbins_in_match_interval=None,
                      matching_tol=None, nPasses=1, save_weights_info=None):
        '''set the weights of the parent particles to zero and set the weights
        of the sampled (child particles) such that a the fit lognormal PDF is 
        maintained'''
        
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
        
        rv, N = gas_orig.fit_interval(fit_func_rng[0], fit_func_rng[1])
        
        hpdf_s, lbins_s = pylab.histogram(ln, range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                          weights=w,
                                         )
        pylab.plot(lbins_s[1::], log10(hpdf_s))
        
        x = linspace(log(1e-3), log(1e+6), 100)
        pylab.plot(x, log10(rv.pdf(x)), 'r--', lw=1)

        pylab.xticks(log(10.0**linspace(-6, 6, 12+1)), linspace(-6, 6, 12+1))
        
        x_fit = exp(x)
        y_fit = rv.pdf(x)
        ###############
        ln_o, w_o = ln[inds_o], w[inds_o]
        ln_s, w_s = ln[inds_s], w[inds_s]
        
        w_o[:] = 1.0
        w_o[ ln_o > log(n_min_sample) ] = 0.0

        w_s[:] = 1.0/float(npp)
        match = True

        if match == False:
            hpdf_s, lbins_s = pylab.histogram(numpy.hstack((ln_o, ln_s)), 
                                              range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                              weights=numpy.hstack((w_o, w_s)),
                                             )
        else:
            
            t0 = time.time()
            
            for p in range(nPasses):
                
                intervals = numpy.linspace(log(match_interval[0]), log(match_interval[1]), nbins_in_match_interval)
                intervals = numpy.vstack((intervals, 
                                          intervals + (intervals[1] - intervals[0]))).T
                
                print 'the intervals where the densities will be matched'
                print intervals
                print 'interal size', intervals[0][1] - intervals[0][0]
                
                weights_in_interval = numpy.zeros(nbins_in_match_interval, 'f8')
                
                for i, interval in enumerate(intervals):
                    
                    print '#'*100
                    print interval
                    print log10(exp(interval))
                    print '----------------%d-----------------------' % i
                                    
                    ln_rng = interval
                    
                    wTrialRng = numpy.array([0.0, 1.0])
                    nTrial = 0
                    ln_x = ln_rng.mean()
    
                    ln_min, ln_max =  ln_rng
    
                    inds = where(((ln_s > ln_min)*(ln_s < ln_max)))[0]
    
                    if inds.size == 0:
                        continue
                    
                    def compute_combined_pdf(wTry):
                            
                        w_s[ inds ] = wTry             
                        
                        all_ln = numpy.hstack((ln_o, ln_s))
                        all_weights = numpy.hstack((w_o, w_s))
                        
                        hpdf_s, lbins_s = pylab.histogram(all_ln, 
                                                          range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                                          weights=all_weights,
                                                         )
                        bins_c = (lbins_s[:-1:] + lbins_s[1::])/2.0
                        f = interpolate.interp1d(bins_c, hpdf_s)
                          
                        return f
    
                    ####
                    # try a weights sweep (check)
                    def sweep_weights():
                        weightsTry = [0.0,
                                      0.001, 
                                      0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                                      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                                      1.0]     
                        print '  x       f(x)     f     |  f(x)-f   | 1 - f(x)/f'
                                                         
                        for wTry in weightsTry:
                            f = compute_combined_pdf(wTry)
                            print '%.3f  %.5f  %.5f  | %.5f  %.5f' % (wTry, f(ln_x), rv.pdf(ln_x), f(ln_x) - rv.pdf(ln_x), numpy.fabs(1.0 - f(ln_x)/rv.pdf(ln_x))  )
                    #sweep_weights()
                    ####
                    
                    while True:
                                            
                        wTry = (wTrialRng[0] + wTrialRng[1])/2.0
                        
                        
                        f = compute_combined_pdf(wTry)
                        
                        
                        print 'f_fit       f_hist\n%.2e     %.2e' % (rv.pdf(ln_x), f(ln_x)) 
                        
                        if numpy.fabs(1.0 - f(ln_x)/rv.pdf(ln_x)) < matching_tol:
                            break
                        else:
                            
                            if f(ln_x) < rv.pdf(ln_x):
                                wTrialRng[0] = wTry
                            else:
                                wTrialRng[1] = wTry
                                
                        print '\t', nTrial, wTrialRng, wTry
                        
                        nTrial += 1
                        
                        #if True:
                        #    break
                    #
                
                    weights_in_interval[i] = wTry
                #
            #    
            print time.time() - t0
            print '---------------------------------------'
        #
        #1.0 - f(log(1e4))/rv.pdf(log(1e4))
        all_ln = numpy.hstack((ln_o, ln_s))
        all_weights = numpy.hstack((w_o, w_s))
        
        hpdf_s, lbins_s = pylab.histogram(all_ln, 
                                          range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                          weights=all_weights,
                                         )
        
        pylab.step(lbins_s[1::], log10(hpdf_s), 'g--', lw=2)
        #pylab.step(lbins_s[1::][::2], log10(hpdf_s)[::2], 'go', markersize=5)
        
        pylab.show()
        
        if save_weights_info != None:
            self.save_weights_function(save_weights_info, exp(intervals), weights_in_interval,
                                       npp, fit_func_rng, n_min_sample, match_interval, 
                                       nbins_in_match_interval, matching_tol, nPasses, x_fit, y_fit)
            
        return numpy.hstack((w_o, w_s)), intervals, weights_in_interval
    
    def set_weights_from_saved_sampling_info(self, sampling_info_file_path, show_plot=True):
        '''sets the weights of the gas particles to the ones of the saved weight function.  If the
        gas set has some particles removed, e.g. the ones with gmech = gsurface, then the plotted
        PDF might not be lognormal.  However the weights and the matching should agree well when
        this is applied to the original snapshot.
        '''

        N = len(self)
        nbins = 100
        
        if show_plot:
            fig = pylab.figure()

        ln = log(self.n)
        
        w = zeros(N, 'f8') 
        
        inds_o = self.get_inds_original_set()
        inds_s = self.get_inds_children()

        w[inds_o] = 1.0
        w[inds_s] = 0.0
        
        ## loading the sampling info from the saved file
        info = self.load_weights_function(sampling_info_file_path)        
        
        ## setting the weights of the original particles above the sampling interval to zero
        gas_orig_amuse = self[self.get_inds_original_set()]
        gas_orig = gas_set(len(gas_orig_amuse)) 
        gas_orig.copy_attr_from_amuse_set(gas_orig_amuse, ['n'])

        if show_plot:        
    
            hpdf_s, lbins_s = pylab.histogram(ln, range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                              weights=w,
                                             )
    
            pylab.plot(lbins_s[1::], log10(hpdf_s))
    
            pylab.plot(log(info['x_fit']), log10(info['y_fit']), 'r--', lw=1)
    
            pylab.xticks(log(10.0**linspace(-6, 6, 12+1)), linspace(-6, 6, 12+1))

        ## setting the weights of the sampled gas particles from the saved function
        ln_o, w_o = ln[inds_o], w[inds_o]
        ln_s, w_s = ln[inds_s], w[inds_s]
        
        w_o[:] = 1.0
        w_o[ ln_o > log(info['n_min_sample']) ] = 0.0

        ###########################################
        for i, interval in enumerate(info['intervals']):
            
            print '#'*100
            print i, interval
            print '----------------%d-----------------------' % i

            ln_rng = log(interval)

            ln_min, ln_max =  ln_rng

            inds = where(((ln_s > ln_min)*(ln_s < ln_max)))[0]

            if inds.size == 0:
                continue

            w_s[ inds ] = info['weight_func'][i]             
        #
        ###########################################
        
        all_ln = numpy.hstack((ln_o, ln_s))
        all_weights = numpy.hstack((w_o, w_s))

        if show_plot:
        
            hpdf_s, lbins_s = pylab.histogram(all_ln, 
                                              range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                              weights=all_weights,
                                             )
            
            pylab.step(lbins_s[1::], log10(hpdf_s), 'g--', lw=2)

            pylab.show()
        
        self.weights = all_weights
        
    def set_weights_sampled_from_saved_weight_function(self, intervals, weight_func):
        '''sets the weights of the sampeld point from the weight function'''
        pass
    
    def save_weights_function(self, fpath, intervals, weights_in_interval, 
                              npp=None, fit_func_rng=None, n_min_sample=None, 
                              match_interval=None, nbins_in_match_interval=None,
                              matching_tol=None, nPasses=None, x_fit=None, y_fit=None):
        '''saves the weight function into a file and all the info that was used in producing it'''
        
        names = ['intervals', 'weight_func', 
                 'npp', 'fit_func_rng', 'n_min_sample', 'match_interval', 
                 'nbins_in_match_interval', 'matching_tol', 'nPasses', 'x_fit', 'y_fit']
        
        data = [intervals, weights_in_interval, 
                npp, fit_func_rng, n_min_sample, 
                match_interval, nbins_in_match_interval,
                matching_tol, nPasses, x_fit, y_fit]
        
        numpy.savez(fpath, data=data, names=names)
    
        print 'saved weight function and info to\n\t%s' % fpath
        
    def load_weights_function(self, fpath):
        '''loads the weight function and returns it as an interpolation function, also returns the
        info used in producing the weight function as a dict'''
        
        weights_file  = numpy.load(fpath)
        
        info_ret = {}
        
        for i, name in enumerate(weights_file['names']):
            
            info_ret[name] = weights_file['data'][i]
            
        return info_ret
    
    def use_weights(self, weighting=None, **kwargs):
        '''set the weights determined by the weighting keyword.
        
        weights_filename = fi_utils.gen_weights_filename(params['rundir'], params['snap_index'])
        gas.use_weights(weighting=params['weights'], weights_filename = weights_filename)

        '''
        
        if weighting == 'matched':
            ''' weights parents = 0, weight children match the PDF'''
            self.set_weights_from_saved_sampling_info(kwargs['weights_filename'], show_plot=False)
            print '\t--------> set mathched weights'

        elif weighting == 'original-only':
            ''' weights children = 0, weights original = 1 '''
            self.set_weights_children_zero_parents_one()
            print '\t---------> set weights such that sampled particles are ignored'
        
        elif weighting == 'by-number':
            ''' wegihts parent + weights children = 1'''
            self.set_weights_childrent_to_inverse_sampling_number(kwargs['weights_filename'])
            print '\t---------> set weights of childrent to ivnerse of number of sampling'
            
        else:
            raise ValueError('unknown weighting!!')
        
    def set_radii_children_to_original_parent_radii(self, rundir, snap_index, to_kpc=True):
        '''sets the radii to children to the the original radii of the parents'''
        
        cinds = self.get_inds_children()
        cid = int32(self.id[cinds])
        
        r = self.radius
        
        ## loading the arrays holding the original information of the parents
        n_o, r_o, ids_o = load_original_n_r_ids(rundir, snap_index)
        ids_o = int32(ids_o)
        
        ##         
        inds_list_of_c_in_o_arrys = find_matching_indicies(ids_o, cid, check=False) 

        ## computing the original radii assuing the radii were scaled using the 
        ## current weights 
        r[cinds] = r_o[inds_list_of_c_in_o_arrys] # these are in pc
        if to_kpc == True:
            r[cinds] *= 1e-3 # converting them to kpc
        
        self.radius = r
        
    def set_radii_children_parent_radii(self, rundir, snap_index, to_kpc=True):
        '''sets the radii to children to the the original radii of the parents'''
        
        cinds = self.get_inds_children()
        cid = int32(self.id[cinds])
        
        r = self.radius
        
        ## loading the arrays holding the original information of the parents
        n_o, r_o, ids_o = load_original_n_r_ids(rundir, snap_index)
        ids_o = int32(ids_o)

        ##         
        inds_list_of_c_in_o_arrys = find_matching_indicies(ids_o, cid, check=False) 

        ## computing the original radii assuing the radii were scaled using the 
        ## current weights 
        r[cinds] = r_o[inds_list_of_c_in_o_arrys] # these are in pc
        if to_kpc == True:
            r[cinds] *= 1e-3 # converting them to kpc
        
        self.radius = r
    
    def set_radii_using_original_weighting(self, rundir, snap_index):
        '''sets the radii and weights such that the sampled population is ignored. i.e the weights and radii
        of the children are set to zero, wheras the weights of the parents are set to one and their radii
        are restored to their original values.
        '''
        
        pinds = self.get_inds_has_children()
        pid = int32(self.id[pinds])
        
        r = self.radius
        
        ## loading the arrays holding the original information of the parents
        n_o, r_o, ids_o = load_original_n_r_ids(rundir, snap_index)
        ids_o = int32(ids_o)
        
        ##         
        inds_list_of_p_in_o_arrys = find_matching_indicies(ids_o, pid, check=False) 

        ## set the radii of the parent to their original value 
        r[pinds] = r_o[inds_list_of_p_in_o_arrys]  # these are in pc
        r[pinds] *= 1e-3 # converting to kpc
        
        ## setting the radii of the childrent to zero
        r[self.get_inds_children()] = 0.0
        
        self.radius = r

        ## setting the appropriate weights
        self.use_weights(weighting='original-only')

    def set_radii_from_Av_and_using_original_weighting(self):
        '''sets the radii and weights such that the sampled population is ignored. 
        and the radii are computed from the AV
        '''
        
        ## set the radii of all the particles to the scale length determined by the Av of the cloud
        r = self.Av / (0.016 * self.n) # computing the radii from the Av (in pc)
        r *= 1e-3 # converting to kpc
        
        ## setting the radii of the childrent to zero
        r[self.get_inds_children()] = 0.0
        
        self.radius = r
        print '\t---------> set the radii from the Av'
        ## setting the appropriate weights
        self.use_weights(weighting='original-only')   # radii are not modified

    def set_radii_from_Av(self):
        '''sets the radii from the Av of the particles'''
        
        ## set the radii of all the particles to the scale length determined by the Av of the cloud
        r = self.Av / (0.016 * self.n) # computing the radii from the Av (in pc)
        r *= 1e-3 # converting to kpc
        
        self.radius = r
        print '\t---------> set the radii from the Av'
        
    def set_radii_using_matched_weighting(self, rundir, snap_index):
        '''
        r_parent = 0
        r_children = matched weights
        '''
        
        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()
        
        ## setting the radii of the children to the original radii of their parents
        self.set_radii_children_to_original_parent_radii(rundir, snap_index) ## in kpc
        
        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)    
        self.use_weights(weighting='matched', weights_filename = weights_filename)

        w = self.weights
        r = self.radius
        n = self.n

        ## setting the radii of the parents to zero
        r[pinds] = 0.0
                
        ## the ids of the children
        cid = int32(self.id[cinds])
                
        ## loading the arrays holding the original information of the parents (the loaded r_o radii are in pc)
        n_o, r_o, ids_o = load_original_n_r_ids(rundir, snap_index)
        ids_o = int32(ids_o)
        
        ##         
        inds_list_of_c_in_o_arrys = find_matching_indicies(ids_o, cid, check=False) 

        ## computing the original radii assuing the radii were scaled using the current weights 
        r[cinds] = r[cinds] * (w[cinds] * n_o[inds_list_of_c_in_o_arrys]/n[cinds])**(1.0/3.0) # these are in pc
        # no need to convert the radii of the children to kpc since they were set and converted
        # in set_radii_children_to_original_parent_radii which are used as r[cinds] in the previous line
        
        self.radius = r
        self.weights = w
    
    def set_radii_using_by_number_weighting(self, rundir, snap_index):
        '''
        r_parent = r_parent_original * (1/(npp+1))**1/3
        w_parent = 1.0 / (npp + 1)
        r_child = r_parent_original * ( (1/(npp+1)) * n_parent / n_child )**1/3
        w_child =   1.0 / (npp + 1)
        '''

        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()
        
        ## first set the radii of the parent particles back to their oringal
        self.set_radii_using_original_weighting(rundir, snap_index)
        
        ## setting the radii of the children to the original radii of their parents       
        self.set_radii_children_to_original_parent_radii(rundir, snap_index)
        
        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)    
        self.use_weights(weighting='by-number', weights_filename = weights_filename)

        w = self.weights
        r = self.radius
        n = self.n        
        
        ## computing the new radii corresponding to the new weights of the parents
        r[pinds] = r[pinds] * (w[pinds])**(1.0/3.0) 

        ## the ids of the children
        cid = int32(self.id[cinds])
                
        ## loading the arrays holding the original information of the parents
        n_o, r_o, ids_o = load_original_n_r_ids(rundir, snap_index)
        ids_o = int32(ids_o)
        
        ##         
        inds_list_of_c_in_o_arrys = find_matching_indicies(ids_o, cid, check=False) 

        ## computing the original radii assuing the radii were scaled using the 
        ## current weights 
        r[cinds] = r[cinds] * (w[cinds] * n_o[inds_list_of_c_in_o_arrys]/n[cinds])**(1.0/3.0)
        
        self.radius = r

    def set_radii_using_by_conserving_total_area(self, rundir, snap_index):
        '''
        r_parent = r_child = r_parent_original * sqrt(N_sampled)
        w_parent = ??
        w_child =  ??
        '''

        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()

        ## first set the radii of the parent particles back to their oringal
        self.set_radii_using_original_weighting(rundir, snap_index)
        
        ## setting the radii of the children to the original radii of their parents       
        self.set_radii_children_to_original_parent_radii(rundir, snap_index)

        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)    
        self.use_weights(weighting='by-number', weights_filename = weights_filename)
        print '-'*100,'\n', '-'*100,'\n','-'*100,'\n'
        print 'weights set by number'
        print '-'*100,'\n', '-'*100,'\n','-'*100,'\n'
        
        w = self.weights
        r = self.radius
        
        ## computing the new radii corresponding to the new weights of the parents
        r = r * numpy.sqrt(w)
        
        self.radius = r

    def set_radii_from_Av_by_conserving_total_area(self, rundir, snap_index):
        '''
        r_parent = 
        r_parent = r_child = r_parent_original * sqrt(N_sampled)
        w_parent = ??
        w_child =  ??
        '''

        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()

        ## first set the radii of the particles to that computed from their Av
        ## the children have the same Av as the parent consequently the same radius
        self.set_radii_from_Av()
        
        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)
        self.use_weights(weighting='by-number', weights_filename = weights_filename)
        print '-'*100,'\n', '-'*100,'\n','-'*100,'\n'
        print 'weights set by number'
        print '-'*100,'\n', '-'*100,'\n','-'*100,'\n'
        
        w = self.weights
        r = self.radius
        
        ## computing the new radii corresponding to the new weights of the parents
        r = r * numpy.sqrt(w)
        
        self.radius = r

    def set_radii_custom1(self, rundir, snap_index):
        '''
        weights parents = 0, weight children match the PDF
        mass = weights_matched * mass
        Av = Av already set before calling this (not changed)
        Area = mass / (mu * mH * N(H) ) # N(H) computed from Av2NH 
        r = sqrt(area / pi) 
        
        !!! this assumes the metallicity is 1
        '''
        
        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()
        
        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)    
        self.use_weights(weighting='matched', weights_filename = weights_filename)

        m = self.mass
        w = self.weights
        r = self.radius
        n = self.n
        Av = self.Av
        
        ## setting the masses
        m *= w
        
        ## computing the column denisty of H
        NH = ismUtils.Av2NH(Av, 1.0)
        
        ## computing the effective area from the mass of the particle and NH
        area_cm2 = m / ( 1.6 * mylib.constants.M_PROTON_CGS * NH)
        
        ## converting the are to kpc2
        area_kpc2 = area_cm2 * (1.0 / mylib.units.KPC2CM)**2.0 

        ## computing the radii
        r = sqrt(area_kpc2 / pi) 

        ## setting the attributes                
        self.radius = r
        self.weights = w
        self.mass = m

    def set_radii_custom2(self, rundir, snap_index):
        '''
        weights parents = 1, weight children = 0
        mass = weights_matched * mass
        Av = Av already set before calling this (not changed)
        Area = mass / (mu * mH * N(H) ) # N(H) computed from Av2NH 
        r = sqrt(area / pi) 
        
        !!! this assumes the metallicity is 1
        '''
        
        metallicity = 1.0
        
        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()
        
        ## setting the new weights by number to the parents and the children
        weights_filename = gen_weights_filename(rundir, snap_index)    
        self.use_weights(weighting='original-only', weights_filename = weights_filename)

        m = self.mass
        m[pinds] *= 101.0
        m[cinds] = 0.0
        
        w = self.weights
        r = self.radius
        n = self.n
        Av = self.Av
        
        Pe = (1.085 * self.T + 54.0 * self.vdisp**2)*self.n #how to set the unit to the same untis as  (units.K*gas.n.unit*constants.kB.unit)
        Av = 0.22 * metallicity * ( 1520.0 / 100.0) * sqrt(Pe/1e4)
        
        ## computing the column denisty of H
        NH = ismUtils.Av2NH(Av, metallicity)
        
        ## computing the effective area from the mass of the particle and NH
        area_cm2 = m / ( 1.6 * mylib.constants.M_PROTON_CGS * NH)
        
        ## converting the area to kpc2
        area_kpc2 = area_cm2 * (1.0 / mylib.units.KPC2CM)**2.0 

        ## computing the radii
        r = sqrt(area_kpc2 / pi) 
        r[cinds] = 0.0
        
        ## setting the attributes
        self.radius = r
        self.weights = w
        self.mass = m
        self.Av = Av
        
    def set_radii_custom3(self, rundir, snap_index):
        '''
        weights parents = 1/101, weight children = 1/101
        mass = weights_matched * mass
        Av = Av already set before calling this (not changed)
        Area = mass / (mu * mH * N(H) ) # N(H) computed from Av2NH 
        r = sqrt(area / pi) 
        
        !!! this assumes the metallicity is 1
        '''

        metallicity = 1.0
                
        ## getting the indicies of the parents and the chidlren
        pinds = self.get_inds_has_children()
        cinds = self.get_inds_children()
        
        m = self.mass
        w = self.weights
        n = self.n

        Pe = (1.085 * self.T + 54.0 * self.vdisp**2)*self.n #how to set the unit to the same untis as  (units.K*gas.n.unit*constants.kB.unit)
        Av = 0.22 * metallicity * ( 1520.0 / 100.0) * sqrt(Pe/1e4)

        ## setting the new weights by number to the parents and the children
        w[pinds] = 1.0/101.0
        w[cinds] = 1.0/101.0
        
        ## computing the column denisty of H
        NH = ismUtils.Av2NH(Av, metallicity)
        
        ## computing the effective area from the mass of the particle and NH
        area_cm2 = m / ( 1.6 * mylib.constants.M_PROTON_CGS * NH)
        
        ## converting the are to kpc2
        area_kpc2 = area_cm2 * (1.0 / mylib.units.KPC2CM)**2.0 

        ## computing the radii
        r = sqrt(area_kpc2 / pi) 

        ## setting the attributes                
        self.radius = r
        self.weights = w
        self.mass = m
        self.Av = Av 
        
    def set_radii(self, weighting=None, rundir=None, snap_index=None):
        '''sets the radii and weights of the particles based on the suggested weighting'''
        
        if weighting == 'matched':
            ''' weights parents = 0, weight children match the PDF'''
            self.set_radii_using_matched_weighting(rundir, snap_index)
            print '\t--------> set radii and weights using mathched weights'

        elif weighting == 'original-only':
            ''' weights children = 0, weights original = 1 '''
            self.set_radii_using_original_weighting(rundir, snap_index)
            print '\t---------> set radii and weights such that sampled particles are ignored'
        
        elif weighting == 'by-number':
            ''' wegihts parent + weights children = 1'''
            self.set_radii_using_by_number_weighting(rundir, snap_index)
            print '\t---------> set radii and weights of childrent to ivnerse of number of sampling'
            
        elif weighting == 'conserve-area':
            ''' the radii are set such that area_parent + area_children = original parent area'''
            self.set_radii_using_by_conserving_total_area(rundir, snap_index)
            print '\t---------> set radii and weights of children such that total area is conserved'
            
            """
            elif weighting == 'by-Av-and-conserve-area':
                '''the radii are set based on the Av of the particle'''
                self.set_radii_from_Av(rundir, snap_index)
                print '\t---------> set radii and weights of children such that total area is conserved'
                
            elif weighting == 'by-Av-and-conserve-area':
                '''the radii are set based on the Av of the particle'''
                self.set_radii_from_Av(rundir, snap_index)
                print '\t---------> set radii and weights of children such that total area is conserved'
            """
        elif weighting == 'custom1':
            ''' customized weighting and setting radii'''
            self.set_radii_custom1(rundir, snap_index)
            print '\t--------> set radii and weights using customized weights (custom1)'
        elif weighting == 'custom2':
            ''' customized weighting and setting radii'''
            self.set_radii_custom2(rundir, snap_index)
            print '\t--------> set radii and weights using customized weights (custom2)'
        elif weighting == 'custom3':
            ''' customized weighting and setting radii'''
            self.set_radii_custom3(rundir, snap_index)
            print '\t--------> set radii and weights using customized weights (custom3)'
             
        else:
            raise ValueError('unknown weighting!!')
        
    def print_stats(self):
        pass

    def get_em_lines(self, specStr, em_unit='Kkms'):
        '''returns the list of attributes for a certain species given an emission unit. Also returns 
        the radex transition index.'''

        all_attr = self.all_attributes()
        
        line_codes = []
        
        for attr in all_attr:
            
            if 'em_flux' in attr and em_unit in attr:
                
                line_code_this = attr.replace('em_fluxKkms_','').replace('em_fluxcgs_','')
                
                spec_this_line = lineDict.lines[line_code_this]['specStr']
                  
                if spec_this_line == specStr:
                        
                    line_codes.append(line_code_this)
                     
        line_codes = numpy.array(line_codes)
        
        ## radex indicies
        idxs = numpy.zeros(line_codes.size, 'i4')
        
        for i, line_code in enumerate(line_codes):
            idxs[i] = lineDict.lines[line_code]['radexIdx']
        
        ## sorting by increasing radex index
        inds = numpy.argsort(idxs)
        idxs = idxs[inds]
        line_codes = line_codes[inds]
        
        return idxs, line_codes
    
    def get_species_with_em(self):
        '''returns a list of string of the species which have emission info'''
        
        all_attr = self.all_attributes()
        
        line_codes = []
        
        for attr in all_attr:
            
            if 'em_flux' in attr:

                line_codes.append(attr.replace('em_fluxKkms_','').replace('em_fluxcgs_','')) 
        
        spec_strs = []
        
        for line_code in numpy.unique(line_codes):
            
            spec_strs.append(lineDict.lines[line_code]['specStr'])
            
        return numpy.unique(spec_strs)
    
    def get_em_attr_name(self, line_str, em_unit='Kkms'):
        '''gets the name of the attribute given a line code and the emission unit'''
        
        if em_unit in ['Kkms', 'cgs']:
            
            attr_name = 'em_flux%s_%s' % (em_unit, line_str)
            
            if attr_name in self.all_attributes():
                
                return attr_name
    
    def get_flux(self, line_str, em_unit='Kkms'):
        '''returns the flux of the speciefied line in the specied unit. Deafult unit is in Kkms'''

        attr_name = self.get_em_attr_name(line_str, em_unit)
        
        return getattr(self, attr_name)         

    def get_luminosity_species(self, specs):
        pass
    
    def get_luminosity(self, line_str, em_unit='Kkms'):
        '''compute the luminosity of particles'''
    
        return luminosity(self.get_flux(line_str, em_unit), self.radius)
        
    def get_total_luminosity_ladder(self, specStr, em_unit='Kkms'):
        '''compute the luminosity of particles'''

        idx, line_codes = self.get_em_lines(specStr, em_unit=em_unit)
        
        total_lum = numpy.zeros(idx.size, 'f8')
        
        for i, line_code in enumerate(line_codes):

            total_lum[i] = self.get_total_luminosity(line_code, em_unit=em_unit)             
            
        return idx, total_lum
    
    def plot_total_luminosity_ladder(self, specStr, em_unit='Kkms', axs=None, *args, **kwargs):
        
        x, y = self.get_total_luminosity_ladder(specStr, em_unit=em_unit)
        
        if axs != None:
            ax_use = axs
        else:
            pylab.figure()
            ax_use = pylab.subplot(111)
            
        plt, = ax_use.plot(x, y, *args, **kwargs)
        
        return plt
    

    def get_total_luminosity(self, line_str, em_unit='Kkms'):
        '''compute the total luminosity of particles'''
         
        return self.get_luminosity(line_str, em_unit).sum() 
        
    def get_mean_flux(self, line_str, em_unit='Kkms', bin_area=1.0):
        '''compute the mean flux of particles given a box size'''
    
        return self.get_total_luminosity(line_str, em_unit)/ bin_area
    
    def convert_radii_to_kpc(self):
        '''multiplt the radii by 1e-3 to convert them from pc to kpc (assuming they are in pc)'''
        
        self.radius = self.radius * 1e-3  # converting the radii from pc to kpc
        print '\twarning converted the radii from pc to kpc'''
         
    def check_particles(self, what_to_check, logger):
        '''checks the attributes and takes care of irregularities'''
        
        logger.debug('='*100)
        logger.debug('checking particle attributes for irregularitie')
        
        if what_to_check == 'default':
            # checking for the know problems
            
            # checking for extremely high values of CO(1-0)
            if 'em_fluxKkms_CO1-0' in self.all_attributes():
                x = getattr(self, 'em_fluxKkms_CO1-0')
                inds = where(x > 500.0)[0]
                x[inds] = 0.0
                new_max = x.max()
                logger.debug('warning: found %d particles with Flux > 300.0 k km / s' % inds.size)
                logger.debug('         setting thier values of the max flux %f k km / s' % new_max)
                x[inds] = new_max
                setattr(self, 'em_fluxKkms_CO1-0', x)

            # checking for extremely high values of CO(1-0) in cgs
            if 'em_fluxcgs_CO1-0' in self.all_attributes():
                x = getattr(self, 'em_fluxcgs_CO1-0')
                inds = where(x > 1.0)[0]
                x[inds] = 0.0
                new_max = x.max()
                logger.debug('warning: found %d particles with Flux > 1.0 erg / cm2 / s' % inds.size)
                logger.debug('         setting thier values of the max flux %f erg / cm2 / s' % new_max)
                x[inds] = new_max
                setattr(self, 'em_fluxcgs_CO1-0', x)
    
def luminosity(fluxes=None, radii=None):
    '''compute the total total_luminosity of particles. F is the array of the fluxes, r 
    is the radius of thte particles'''

    areas = pi * radii**2 # the projected area of the SPH particles
    lum = fluxes*areas   # the luminosity of each particle

    return lum

def total_luminosity(fluxes=None, radii=None):
    '''compute the total total_luminosity of particles. F is the array of the fluxes, 
    r is the radius of thte particles'''
    
    return luminosity(fluxes, radii).sum()
    
def mean_flux(fluxes=None, radii=None, bin_area=None):
    '''compute the mean/total flux of particles. F is the array of the fluxes, r is the radius 
    of thte particles, area_pixel is the area of the pixel in which the particles are.'''
    
    
    lum_pixel = total_luminosity(fluxes, radii)
    
    return lum_pixel / bin_area

def gen_weights_filename(rundir, snap_index):
    '''returns the weights filename
    
        fname = gen_weights_filename(params['rundir'], 4)
    '''

    return rundir + '/firun/' + 'weights_func.%06d.npz' % snap_index

def load_original_n_r_ids(rundir, snap_index):
    
    fname = rundir + '/firun/' + 'fiout.%06d.original_n_r_ids.npz' % snap_index
    
    f = numpy.load(fname)
    
    n, r, ids = f['arr_0']
    
    return n, r, ids


def flux_contribution_vs_n(n_bins, cloud_parms, F, nPDF, r_func, N, label, in_ax=None, no_plot=False, color=None):
    '''
    '''
    n_bins = n_bins.copy()
    
    log10G0, log10Gmech, Av = cloud_parms
    
    dlog10n = log10(n_bins)[1] - log10(n_bins)[0]
    
    log10G0s = ones(n_bins.size, 'f8')*log10G0
    log10gmechs = ones(n_bins.size, 'f8')*log10Gmech
    Avs = ones(n_bins.size, 'f8')*Av
    
    data = numpy.array([log10(n_bins), log10G0s, log10gmechs, Avs]).T
    
    flux_bins = 10.0**F.get( data )

    prob_n_bins = nPDF(n_bins)
     
    flux_nbins = (prob_n_bins*dlog10n)*flux_bins
    
    inds = where( isnan(flux_nbins) == False )
    flux_nbins = flux_nbins[inds]
    n_bins = n_bins[inds]
    
    # scaling the fluxes (assuming there are N clouds)
    flux_nbins *= N     
        
    if no_plot == False:
        
        if in_ax == None:
            ax=gca()
        else:
            ax=in_ax
        
        x, y = n_bins, flux_nbins
        ax.loglog(x, y, color, label=label)
        
        y_norm = y / y.sum()
        y_norm_c = y_norm.cumsum()
        
        
        finterp = interpolate.interp1d(log10(x), y_norm_c)
        finterp2 = interpolate.interp1d(log10(x), y)
        x_interp = linspace(log10(x.min()), log10(x.max()), 1000.0)
        y_norm_c_interp = finterp(x_interp)
        
        inds = where(y_norm_c_interp > 0.9)[0]
        x, y = x_interp[inds[0]], y_norm_c_interp[inds[0]]
        ax.loglog(10.0**x, finterp2(x), color + 'o', markersize=4)
        
        #ind = where()

    return flux_nbins

def plot_luminosity_from_pdf(nPDF, r_func, arxvPDR, F, params, nmin, nmax,
                        yrng=[1e-6, 1e12]):
    '''
    '''
    
    ########################
    figure(figsize=(16,8))

    subplot(231)
    ## radius (in kpc) as a function of density (in cm-3)    
    xs = linspace(-3.0, 6.0, 100)
    loglog(10.0**xs, r_func(10.0**xs)*1e3, 'g--', linewidth=2)
    xlabel('n (cm-3)')
    ylabel('r (pc)')
    
    ## getting the fitted PDF of the density (reading it from the disk)
    subplot(232)
    loglog(10.0**xs, nPDF(10.0**xs), 'r--', lw=1)
    
    ################################################################################################


    # computing the total luminosity discretely over the whole PDF
    
    log10n, dlog10n  = linspace(nmin, nmax, 100, retstep=True, endpoint=False)
    
    ###################
    n_bins = 10.0**(log10n + dlog10n*0.5)
    #prob_n_bins = nPDF(n_bins)
    #r_nbins = r_func(n_bins)
    
    N_particles = 2e6

    subplot(2, 3, 4)
    for log10G0 in [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
        flux_contribution_vs_n(n_bins, [log10G0, -22, 10.0], F,  
                             nPDF, r_func, N_particles, log10G0)
    xlim(10.0**nmin, 10.0**nmax)
    ylim(*yrng)
    legend(loc=0, prop={'size':8})
    
    subplot(2, 3, 5)
    for log10Gmech in [-50.0, -30.0, -25.0, -23.0, -22.0, -21.0]:
        flux_contribution_vs_n(n_bins, [2.0, log10Gmech, 10.0], F, 
                             nPDF, r_func, N_particles, log10Gmech)
    xlim(10.0**nmin, 10.0**nmax)
    ylim(*yrng)    
    legend(loc=0, prop={'size':8})

    subplot(2, 3, 6)
    for Av in [0.1, 1.0, 5.0, 10.0, 20.0]:
        flux_contribution_vs_n(n_bins, [2.0, -22.0, Av], F,  
                             nPDF, r_func, N_particles, Av)
    xlim(10.0**nmin, 10.0**nmax)
    ylim(*yrng)
    legend(loc=0, prop={'size':8})
    
    pylab.show() 
    return info

def effective_area_tests(gas, gasp, params):

    def area_info(pid):
        
        r = gas.radius[ gas.id == pid ]
        
        area = pi * r**2
        total_area = area.sum()
        area_parent = pi*(gasp.radius[gasp.id == pid])**2
        
        return total_area, area_parent
    '''
    print 'number of particles which ave children', len(gasp)
    
    every = 1000
    n = gasp.n[::every]
    area_ratio = zeros(n.size, 'f8')
     
    for i, this_id in enumerate(gasp.id[::every]):
        
        r = gas.radius[ gas.id == this_id ]
        
        area = pi * r**2
        total_area = area.sum() 
        area_parent = pi*(gasp.radius[0])**2
        
        print r
        print i
        print 'total area = ', total_area
        print 'area parent = ', area_parent 
        print 'area ratios = ', total_area / area_parent
    
        area_ratio[i] = total_area / area_parent
    
    
    semilogx(n, area_ratio, 'o')
    '''
    
    '''
    ids of particles with n > 1e5
    1352118,  143025, 1850783,  636856, 1658545,  209078, 1956500,
    127946,  300592, 1604612, 1548951, 1522416,  234501, 1209118,
    1317109, 1212827, 1335993,  930197, 1321721, 1030499,  848041,
    802810, 1022996, 1767978, 1878603, 1491714,  786006, 1026492,
    1810486,  127299, 1318426, 1529623,  168216,  599553, 1601742,
    1384949,   88511,  189546,    5310,  592677,  222482, 1650119,
    1880076,  354930,  103668, 1256942,  212912,  858876, 1942415,
    612252, 1909999,  327141,  606273,  239376, 1192756, 1537287,
    589815, 1111410,  246552
    '''
    
    pid = 5310
    ind = where(gasp.id == pid)
    line = 'HCN1-0'
    unit = 'cgs'
    
    fluxp = gas.get_flux(line, em_unit=unit)
    flux = gas.get_flux(line, em_unit=unit)
    
    r = gas.radius[ gas.id == gasp.id[ind] ]
    n = gas.n[ gas.id == gasp.id[ind] ]
    f = flux[ gas.id == gasp.id[ind] ]
    
    loglog(n, f, 'o')
    loglog([1e2, 1e6], [fluxp[ind],fluxp[ind]],'--' )
    ylim(1e-10, 1e1)
    
    areas = area_info(pid)
     
    print areas[0]/areas[1]

#################################################
class density_distribution(object):
    
    def __init__(self, dist, *args, **kwargs):
        
        self.pdf = None
        
        if dist == 'lognormal':
            ''' '''
            self.pdf = self.make_lognormal_pdf(*args, **kwargs)
            
            
    def make_lognormal_pdf(self, *args, **kwargs):
        
        rv = stats.norm(
                        loc=log(kwargs['mu']), 
                        scale=log(kwargs['sigma'])
                        )
        
        return lambda n: rv.pdf(log(n))
                                

def mass_from_pdf(n_pdf, r_func, n_min, n_max, res=100):
    '''compute the total mass from a density PDF and the radii of the assumed spherical
    clouds within a density range.
    
    # radius density scaling function, r(kpc) = R(n[cm-3]) (kpc)     
    r_func = coset9_sol_info.r_sph_kpc

    ## the probability density at that density (not in log scale
    nPDF = coset9_sol_info.nPDF_sph.pdf
    #nPDF = coset9_sol_info.nPDF_sph_test.pdf

    fi_utils.mass_from_pdf(nPDF, r_func, 1.0, 1e6)
    '''
    
    # the locations on the denisty pdf in log10
    l10n, dlog10n = linspace(log10(n_min), log10(n_max), res, endpoint=False, retstep=True)
    dlogn = log(10.0**dlog10n)
    l10n += (dlog10n*0.5)

    # the density in cgs
    n_cgs = 10.0**l10n
    
    # the values of the PDF at these densities
    p = n_pdf(n_cgs)*dlogn
    subplot(221)
    loglog(n_cgs, p)
    ylabel('PDF')    
    print 'print total probability in this interval = ', p.sum()
    
    # volume of the clouds in the denisty bins (cm3)
    r = r_func(n_cgs)*KPC2CM  # radius in cm
    v = (4.0/3.0)*pi*r**3  # in cm^3
    subplot(222)
    loglog(n_cgs, r / KPC2CM)
    ylabel('r(kpc)')
    # mass of the clouds as a function of density
    rho = n_cgs*M_PROTON_CGS 
    m_cgs = rho*v # in grams
    m_msun = m_cgs / M_SUN_CGS
    print n_cgs[0], r[0] / KPC2CM, m_msun[0]

    # mass of the clouds in each density bin
    m_bin_msun = p * m_msun
    subplot(223)
    loglog(n_cgs, m_bin_msun, label='mass PDF')
    loglog(n_cgs, m_msun, label='mass cloud')
    legend(loc=0)
    print 'total mass of the distribution (M_sun) = %e' % m_bin_msun.sum()
    
    # mass of the clouds in each density bin
    #m_bin_msun = p * m_msun
    #subplot(413)
    #loglog(n_cgs, m_msun)

    # total mass
    m_total = m_bin_msun.sum()
    
    show()
    return m_total
        
        