from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
import numpy
import time
from numpy import log10
import lineDict

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
        raise ValueError('interpolated variable Tkin_sph_from_pdf has nan values in it!!')
    
    
def guess_metallicity(runDir):
    if '-std' in runDir:
        metallicity = 0.2
    if '-sol' in runDir:
        metallicity = 1.0
    return metallicity

def gas_particle_emission(gas, map_key, arxvPDR, interp_funcs):
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
    em_ret = 10.0**interp_funcs[map_key](data_use)
    dt = time.time() - t0

    print 'time interpolating emission for map %s = %.3f seconds for %d points at %e points/sec' % (map_key, dt, data_use.shape[0], data_use.shape[0]/dt)

    #setting emissions of particles with nan interpolated values to zero             
    indsNan = numpy.where(numpy.isnan(em_ret))[0]
    if indsNan.size != 0:
        em_ret[indsNan] = 1e-30
                
    return em_ret


def gas_particle_all_emissions(gas, map_keys, arxvPDR, em_interp_funcs):
    '''returns the interpolated emission of the gas particles whose states are gas_n_g0_gm_Av
     for the maps. 'maps' is a list of string specifying the emission lines for example: ('CO1-0', 'CO2-1'...etc..)

    '''
    
    em_sph = {}
    
    #filtering out the emissions which should be computed
    em_keys = []
    for map_key in map_keys: 
        if '_em_' in map_key:
            em_keys.append(map_key)
            
    if len(em_keys) > 0:
        for em_key in em_keys:
            em_sph[em_key] = gas_particle_emission(gas, em_key, arxvPDR, em_interp_funcs)

    return em_sph

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
            
            if 'em_' in attr and '__' not in attr:
                
                lineStr = attr.split('_')[-1]
         
                if specStr == lineDict.lines[lineStr]['specStr']:
                    attr_save_this_spec.append(attr)
                #
            #
        #
        print attr_save_this_spec
        filename_spec = filename + '.em_%s.npz' % specStr
        save_gas_particle_info(filename_spec, gas, attr_save_this_spec)
        print '\t\t\twrote file : %s ' % filename_spec
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

def load_gas_particle_info(filename):
    '''save the data of the sph particles into a file. Only the attributes in the list 
     are saved. The file is written into a npz file.
    '''
    
    data = numpy.load(filename)
    
    n = data['arr_0'].size

    gas = Particles(n)

    for i, name in enumerate(data['names']):
        setattr(gas, name, data['arr_%d' % i])
        
    return gas, data['names']

def load_gas_particle_info_with_em(filename, species):
    '''loads the specified file foo.states.npz, also looks for files
     foo.states.npz.em_XYZ.npz, where XYZ is an entry in the list of strings in species
    '''
    
    names_all = []
    
    gas, names = load_gas_particle_info(filename)
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

def make_map(map_key, gas, hist, as_log10=None):
    '''looping over the bins of the 2D histogram of the x,y coordinates and 
     computing the averages of the maps
    '''

    map_data = numpy.zeros((hist.nBins[0], hist.nBins[1]), dtype=numpy.float64)
    
    attr_data = getattr(gas, map_key)
    
    for i in numpy.arange(hist.nBins[0]):
        
        for j in numpy.arange(hist.nBins[1]):
            
            inds_in_bin = hist.get_indicies([i,j])
            
            map_data[i,j] = inds_in_bin.size
            
            if inds_in_bin.size > 0:
                
                #map_data[i,j] = attr_data[inds_in_bin].sum()
                map_data[i,j] = numpy.mean(attr_data[inds_in_bin])
            #
        #
    #
    
    if as_log10 != None and as_log10 == True:
        map_data = numpy.log10(map_data)

    return map_data
