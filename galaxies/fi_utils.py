from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles
import numpy

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
                       (gas.x >= ranges['box_size'][0])*(gas.x <= ranges['box_size'][1])*
                       (gas.y >= ranges['box_size'][0])*(gas.y <= ranges['box_size'][1])*
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
    #print '##########using nearest neigbour intepolator....use linear'
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
    
def make_maps(**kwargs):
    """
    """
    
    pass