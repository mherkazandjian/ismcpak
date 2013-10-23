import os, sys, time
import subprocess, shlex

import numpy
from numpy import log10

import pylab

import time

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
from amuse.io import read_set_from_file, fi_io

import meshUtils
import lineDict
import pdrDict

from mylib.utils.histogram import hist_nd
from mylib.utils.interpolation import sectioned_4D_interpolator 

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
        
    """
    
    gas_new = Particles(len(gas))
    
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
    returns the sph particle set  which have attributes withing the specied radges'''

    #selecting particles within a spatial range and within the 
    #density and g0 range of the database
    inds = numpy.where(
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
    """Computes the ratio of the mechanical heating of an sph particle to the surface heating of
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
                          
    print 'interpolating the corresponding surface heating at the surface of the SPH particles'
    gammaSurf_sph_from_pdr = 10.0**f_log_gamma_surf(dataNew)
    gas.alpha = gas.gmech/gammaSurf_sph_from_pdr
    
    indsNan = numpy.where(numpy.isnan(gammaSurf_sph_from_pdr))[0]
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
    map
    '''
    
    #what points to use? with or witout gm
    if '_no_gm_' in map_key:
        data_use = numpy.array([log10(gas.n), log10(gas.G0), numpy.ones(len(gas))*arxvPDR.grid_z.min(), gas.Av]).T        
    else:
        data_use = numpy.array([log10(gas.n), log10(gas.G0), log10(gas.gmech), gas.Av]).T

    print 'interpolating emissions of %s of the sph particles for %d particles' % (map_key, data_use.shape[0])

    t0 = time.time()
    v_ret = interp_funcs[map_key](data_use)
    dt = time.time() - t0

    if take_exp != None and take_exp == True:
        v_ret = 10.0**v_ret

    print 'time interpolating for %s = %.3f seconds for %d points at %e points/sec' % (map_key, dt, data_use.shape[0], data_use.shape[0]/dt)

    #setting emissions of particles with nan interpolated values to zero             
    indsNan = numpy.where(numpy.isnan(v_ret))[0]
    if indsNan.size != 0:
        v_ret[indsNan] = 1e-30
                
    return v_ret


def save_gas_particle_info_saperate_files(path_of_snapshot, gas, species):
    '''save the data of the sph particles into a files. Only the attributes in the list are saved. The
     file is written into a npz file.
    '''

    #-----------------------------------------
    #saving the particle states
    filename = path_of_snapshot + '.states.npz'
    save_gas_particle_info(filename, gas, ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'n', 'G0', 'gmech', 'Av', 'T', 'alpha'])
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
    '''save the data of the sph particles into a file. Only the attributes in the list are saved. The
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
    n = data['arr_0'].size

    gas = Particles(n)

    for i, name in enumerate(data['names']):
        setattr(gas, name, data['arr_%d' % i])
    
    all_names = data['names'] 
    
    if load_pdr == True:
        data_pdr = numpy.load(filename + '.pdr.npz')
        for i, name in enumerate(data_pdr['names']):
            setattr(gas, name, data_pdr['arr_%d' % i])
        
        all_names = numpy.hstack((all_names, data_pdr['names']))

    return gas, all_names

def load_gas_particle_info_with_em(filename, species, load_pdr=None):
    '''loads the specified file foo.states.npz, also looks for files
     foo.states.npz.em_XYZ.npz, where XYZ is an entry in the list of strings in species
    '''
    
    names_all = []
    
    #loading the states of the sph particles
    gas, names = load_gas_particle_info(filename, load_pdr=load_pdr)
    print 'loaded file: \n\t\t %s' % filename
    
    names_all.append(names.tolist())     
    
    for specStr in species:
        filename_this = filename + '.em_%s.npz' % specStr
        gas_this, names_this = load_gas_particle_info(filename_this)
        print '\tloaded file:\n\t\t\t %s' % filename_this
        for name_this in names_this:
            attr_data_this = getattr(gas_this, name_this)
            setattr(gas, name_this, attr_data_this)
    
        names_all.append(names_this.tolist())
        
        del gas_this
    
    print 'att attributes loaded are : ', names_all
    
    print 'warning: make sure the computed optical depths make sense!!! it doesnt seem so..'
    return gas

def make_map(gas, hist, attr=None, as_log10=None, func=None, **kwargs):
    '''looping over the bins of the 2D histogram of the x,y coordinates and 
     computing the averages of the maps
    '''

    map_data = numpy.zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)

    if attr not in dir(gas[0]):  
        raise AttributeError('%s is not an attribute of the gas particles,' % attr)
        
    attr_data = getattr(gas, attr)
    
    for i in numpy.arange(hist.nBins[0]):
        
        for j in numpy.arange(hist.nBins[1]):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            map_data[i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:
                
                map_data[i,j] = func(attr_data[inds_in_bin])
            #
        #
    #
    
    if as_log10 != None and as_log10 == True:
        map_data = numpy.log10(map_data)

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
    Avs = numpy.linspace(Av_range[0], Av_range[1], (Av_range[1] - Av_range[0])/Av_res + 1) 

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

def snapshot_interpolated_data(snap, arxvPDR, em_interp_funcs, pdr_interp_funcs, em_keys, pdr_keys, params, logger):
    '''get the emission info using supplied interpolation functions for a certain snapshot
     Returns the requested emissions as a dict and also returnes the gas particles with those
     set as attributes.
    '''
    
    #snapshot index
    snapIndex = snap
    
    suffix = '%06d' % snapIndex
    #path to fi snapshots
    snapName = 'fiout.%s' % suffix 
    snap_filename = params['rundir'] + '/firun/' + snapName 

    #extracting/guessing the metallicity from the name of the directory of the run
    metallicity = guess_metallicity(params['rundir'])
    
    #loading the sph simulation data
    logger.debug('loading snapshot %s : ' % snap_filename)
    t0 = time.time()
    gas, dark, stars = read_set_from_file(snap_filename, format = fi_io.FiFileFormatProcessor)
    
    file_size_bytes = os.stat(snap_filename).st_size
    dt = time.time() - t0
    logger.debug('done reading fi snapshot : %s \n\t\t\t time reading = %.2e seconds @ %.2e MB/s' % (snap_filename, dt, (file_size_bytes/dt)/1024.0**2) )
    logger.debug('number of sph particles in snapshot = %d' %  len(gas))
    
    #getting a new particle set for the gas particles with the attributes we will be using later
    #and converting the units to the same as those of the PDR models    
    gas = convert_units_to_pdr_units(gas, metallicity)
    logger.debug('converted and computed stuff from sph data (in units compatible with the pdr models)')
    
    print 'quantity          min              max  '
    print 'n_gas         %e     %e   '  % (gas.n.min()     , gas.n.max())
    print 'G0_gas        %e     %e   '  % (gas.G0.min()    , gas.G0.max())
    print 'g_mech        %e     %e   '  % (gas.gmech.min() , gas.gmech.max())
    print 'Av_gas        %e     %e   '  % (gas.Av.min()    , gas.Av.max())
    
    #keeping gas particles within the specified ranges
    gas = select_particles(gas, params['ranges'])
    logger.debug('got the sph particles in the required ranges')

    lgMechMin = arxvPDR.grid_z.min()

    #setting the lower bound of the mechanical heating of the sph particles of the 
    #minimum of the pdr database which is almost negligable (since some of the 
    #sph particles have a zero gmech, their log would be -inf, so we clip it to the
    #mininmum of the pdr database
    gas.gmech = numpy.clip(gas.gmech, 10.0**lgMechMin, numpy.Infinity)
    #########################done setting up the snapshot data######################    
    
    #computing alpha for each sph particle (amount of gm wrt surface heating for the PDR model)
    #and selecting particles which have alpha < 1.0
    compute_alpha(gas, arxvPDR)
    logger.debug('computed the alpha for the sph particles')
    inds = numpy.where( gas.alpha < 1.0 )
    gas = gas[inds]
     
    #setting the maximum Av of the sph particles to a predifined value (usually 
    #the emissions do not increase after some Av bec of the optical depth of the lines)
    gas.Av = numpy.clip(gas.Av, params['ranges']['sph']['Av_clip'][0], params['ranges']['sph']['Av_clip'][1])
    
    #computing the emissions of the SPH particles from the PDR interpolation functions
    print 'getting the emissions of each sph particle'    

    em_sph = {}        
    for em_key in em_keys: em_sph[em_key] = gas_particle_info_from_interp_func(gas, em_key, arxvPDR, em_interp_funcs, take_exp=params['maps'][em_key]['log10'])

    pdr_sph = {}        
    for pdr_key in pdr_keys: pdr_sph[pdr_key] = gas_particle_info_from_interp_func(gas, pdr_key, arxvPDR, pdr_interp_funcs, take_exp=params['maps'][pdr_key]['log10'])
    
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
    
    im = ax.imshow(map_data,
                   extent=[bs_min, bs_max, bs_min, bs_max],
                   vmin=map_info['v_rng'][0],  
                   vmax=map_info['v_rng'][1], 
                   interpolation='bessel', #intepolation used for imshow
                   origin='lower')

    ax.set_title(map_info['title'], size=30, verticalalignment='bottom')
    ax.tick_params(axis='both', which='major', labelsize=20)
    
    cbar_ax = fig.add_axes([0.2, 0.97, 0.6, 0.02]) 
    cbar_ax.tick_params(axis='both', which='major', labelsize=30)

    
    pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal', 
                   ticks=numpy.linspace(map_info['v_rng'][0], map_info['v_rng'][1], 5))
    
    pylab.figtext(0.01, 0.87, '%.2f' % snap_time + 'Gyr', 
            color='black', size='xx-large', weight='bold')

    if 'save_image' in params and params['save_image'] == True:
        filename_fig = processed_snap_filename + '.eps' 
        fig.savefig(filename_fig)
        print 'saved image file to :\n\t\t\t %s' % filename_fig
    
    #copy the saved image to the latex direcotry
    if 'latex' in params and params['latex'] == True:

        copy_to = os.path.join(params['latex'], '%s-%s-%s' % (os.path.split(os.path.split(params['rundir'])[-2])[-1],    
                                                              os.path.split(params['rundir'])[-1],
                                                              os.path.split(filename_fig)[-1]
                                                              )
                               )
        command = 'cp %s %s' % (filename_fig, copy_to)
        args = shlex.split(command)
        subprocess.call(args)
        
        print 'copied %s to \n\t %s' % (filename_fig, copy_to) 

def plot_cube_sections(data_cube, params):
    '''given a cube, it plots the sections (in velocity channels) of the emission'''
    
    #parameters of a typical spectrum for the cube
    v_min, v_max = params['cube']['spec_rng']
    v_res = params['cube']['spec_res']
    
    #the velocity bins (.. todo:: xxx shift this to the centroid of the velcoty bins) 
    v = numpy.linspace(v_min, v_max, v_res)
    
    v_sec_wdith = (v_max - v_min)/params['cube']['plot']['n_vsec_plt']
    
    nx = params['cube']['plot']['n_per_row']
    ny = params['cube']['plot']['n_vsec_plt'] / nx
    
    fig, axs = pylab.subplots(ny, nx, sharex=True, sharey=True, figsize=(12.0*numpy.float(nx)/numpy.float(ny), 12), 
                              subplot_kw = {'xlim':params['cube']['xy_rng'][0:2],
                                            'ylim':params['cube']['xy_rng'][2:4],
                                            'aspect':'equal',
                                            'adjustable':'datalim',
                                            })
                                  
    for ax in axs[:,0] : ax.set_ylabel('y(kpc)')
    for ax in axs[-1,:]: ax.set_xlabel('x(kpc)')
        
    pylab.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.9, wspace=0.15, hspace=0.15)
    
    
    #plotting the maps for all the velocity channels
    for n in numpy.arange(params['cube']['plot']['n_vsec_plt']):
     
        #range in velocity for this velocity channel
        v_from = n*v_sec_wdith + v_min
        v_to   = v_from + v_sec_wdith
        #inds in the z direction of the cube (along the spectrum for this channel)
        v_inds = numpy.where( (v > v_from) * (v <= v_to ) )  
    
        #computing the map for this channel by adding all the spectra for this channel
        this_map = data_cube[:,:, v_inds[0]]
        this_map = numpy.log10(this_map.sum(axis=2)) 
        
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

    cbar_ax = fig.add_axes([0.2, 0.95, 0.6, 0.01]) 
    cbar_ax.tick_params(axis='both', which='major', labelsize=30)
    cbar_ax.set_title(r'$\log_{10}$ [ T$_B$ / K.km.s$^{-1}$]', size=25)
    
    pylab.colorbar(im, ax=ax, cax=cbar_ax, orientation='horizontal', 
                   ticks=numpy.linspace(params['cube']['plot']['rng'][0], params['cube']['plot']['rng'][1], 5))

    pylab.show()
    
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

    plot_map(map_data, map_params, map_info, snap_time, params, snap_data_filename)
#

def plot_all_maps_for_snapshot_from_saved_data(snapIndex, params):
    '''
    '''    
    #getting all the maps
    for this_map in params['all_maps']:
    
        this_map_info = params['all_maps'][this_map]
        
        plot_map_from_saved_data(snapIndex, this_map_info, params)
        
    print '---------------------------------------------------------------'
    pylab.show()
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