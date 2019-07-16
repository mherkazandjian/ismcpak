"""
<keywords>
test, pdr, amuse
</keywords>
<description>
an example code to demonstrate how to execute the PDR code via the amuse 
interface.
</description>
"""
import os

from amuse.community.pdr import interface

nWorker = 1  # number of proccesses

dataDir     = '/ism/ismcpak/data/'
outputDir   = '../ism/data/oneSided/single_mesh/'

pdr = interface.pdrInterface(number_of_workers = nWorker, redirection='none')

log10_nGas = 1.0
log10_G0   = 3.0
lMech      = 1e-50

#---------------------amuse parameters------------------------------------------------
pdr.set_outputDir                  (outputDir + 'meshes/')
pdr.set_species_fName              (dataDir + 'pdr/species.inp')
pdr.set_underUbundant_fName        (dataDir + 'pdr/underabundant.inp')
pdr.set_rate99_fName               (dataDir + 'pdr/rate99.inp')
pdr.set_selfSheilding_CO_fName     (dataDir + 'pdr/self_shielding_CO.inp')
pdr.set_rotationalCooling_baseName (dataDir + 'pdr/rotationalcooling/rotcool')
pdr.set_vibrationalCooling_baseName(dataDir + 'pdr/vibrationalcooling/vibcool')
pdr.set_database_fName             (dataDir + 'pdr/z-1.0.dat')
pdr.set_zeta                       (5.0e-17)
pdr.set_S_depletion                (200.0)
pdr.set_TTol                       (1e-3)
pdr.set_CTol                       (1e-3)
pdr.set_metalicity                 (1.0)
pdr.set_AvMax                      (30.0)
pdr.set_slabSizeCrit               (0.5)
pdr.set_min_deltaAv                (0.01)
pdr.set_max_deltaAv                (0.5)
pdr.set_maxSlabs                   (200)
#------------------------------------------------------------------------------------

# initializing variables and running the models
pdr.initialize()
ids, err = pdr.add_mesh(log10_nGas, 1000.0, log10_G0, lMech)

pdr.calc_equilibrium()

mshErr, err = pdr.get_errorFlags(ids)

T, err = pdr.get_temperature(ids)

print 'done'
