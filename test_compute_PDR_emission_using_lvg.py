"""
<keywords>
PDR, LVG, radex, line width, line-width, velocity dispersion, emission
</keywords>
<description>
compute the emission of a PDR for varying line-widths
</description>
"""

import os
import numpy
import pylab
import time

import meshUtils

home = '/home/mher'

# path of the PDR database
pdr_dirPath  = os.path.join(home, 'ism/runs/oneSided/sph-db-z-1.0-low-res/')

# path of the radex excutable
radex_path = os.path.join(home, 'ism/code/radex/Radex/bin-gcc/radex')  
molData_path = os.path.join(home, 'ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles')

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = pdr_dirPath, readDb=True)

pdr_parms = [2.0, 4.0, -50.0] # log10n, log10G, log10Gamma_mech
line_widths = numpy.linspace(0.5, 10.0, 20.0)

radex_parms = {'path'                 : radex_path,   
               'molDataDirPath'       : molData_path,
               'specStr'              : 'CO',
               'freqRange'            : [0, 50000],
               #'xH2_Min'              : 2*0.0000000001
               'xH2_Min'              : -1.0,
               #'collisionPartners'    : ['H2','H+','H','e-','He'],
               #'collisionPartners'    : ['H2','H','H+','e-'],
               'collisionPartners'    : ['H2'],
               'use_pdr_gas_den_H2'   : True,   #<----------
               'tBack'                : 2.73,
               'lineWidth'            : 1.0,
               'verbose'              : False,
               'maxDisplayTranistion' : 20,
               'quantity'             : 'fluxcgs',  # 'fluxKkms', 'fluxcgs' 
               ###----------extra convergence params-----------------------
               'checkOutputIntegrity'       : True,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
               'popDensSumExpected'         : 1.0, 
               'popDensSumTol'              : 1e-2,
               'changeFracTrial'            : 0.05,
               'strict'                     : True,
               'nMaxTrial'                  : 100,
              }
Av_range = [0.0, 28.0]

## initializing the radex object
arxvPDR.setup_default_radex_instance(radex_parms=radex_parms)

def compute_pdr_emission_using_lvg_code(arxvPDR, pdr_parms=None, lvg_parms=None, Av_range=None):
        
    ## setting to data to the PDR mesh object to be used in computing the info
    ## needed by the LVG code in computing the emission
    mesh_data = arxvPDR.get_mesh_data(x=pdr_parms[0], y=pdr_parms[1], z=pdr_parms[2])
    
    
    ## running radex to compute the emission
    radex_flag, radex_pdr_parms = arxvPDR.run_radex_on_pdr_model(radex_parms = radex_parms, Av_range = Av_range)
    
    ## getting the transitions
    transitions = arxvPDR.radexObj.transitions

    print transitions
    print 'pdr mean weighted T', radex_pdr_parms[0] 
    print 'pdr mean weighted coll densities', radex_pdr_parms[1] 
    print 'pdr mean weighted column density', radex_pdr_parms[2] 

    return radex_flag, radex_pdr_parms, transitions
#

all_transitions = []

t0 = time.time()
for sigma in line_widths:
    
    arxvPDR.radexObj.setInFileParm('lineWidth', sigma**0.2)
    
    flag, pdr_lvg_parms, transitions = compute_pdr_emission_using_lvg_code(arxvPDR,     
                                                                           pdr_parms=pdr_parms,     
                                                                           lvg_parms=radex_parms,
                                                                           Av_range=Av_range)
    all_transitions.append(transitions)
print 'finished computing in %.2f seconds' % (time.time() - t0)

# plot the flux of a transitions
Jlow = 0
#pylab.semilogy(line_widths, [transition['fluxKkms'][Jlow] for transition in all_transitions], 'r-')
pylab.plot(line_widths, [transition['fluxKkms'][Jlow] for transition in all_transitions], 'r-')
#Jlow = 10
#pylab.semilogy(line_widths, [transition['fluxKkms'][Jlow] for transition in all_transitions], 'b-')

pylab.show()

print 'done'