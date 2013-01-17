from numpy import *
from time import *
import sys
import pylab as pyl

from amuse.community.pdr import interface
from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *
#---------------------------------------------------------------------------
home = '/home/mher'

outputDir   = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/'
metallicity = 0.1   # in terms of solar metallicity
alpha       = 1.0   #percent of the surface heating to be used as gmech
log10nGas   = 3.0
log10G0     = 5.0


#reference database path from which the surface heating in the case of zero gmech will be extracted
databasePath  = '/home/mher/ism/runs/oneSided/surfaceGrid-z-0.1/'

plotRangenG0  = [[0,6],[0,6]]
#-----------------chemical network parameters------------------------
chemParms  = {
              'rxnFile'       : home + '/ism/code/ismcpak/data/rate99Fixed.inp',
              'specNumFile'   : home + '/ism/code/ismcpak/data/species.inp',
              'underAbunFile' : home + '/ism/code/ismcpak/data/underabundant.inp',
              'removeManual'  : ['13CH3'],
              'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
              'umistVer'      : 'umist99',
              }

nWorker = 1  # number of proccesses
pdr     = interface.pdrInterface( number_of_workers = nWorker, redirection='none') 

 
#---------------------amuse parameters------------------------------------------------
pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              ("/home/mher/ism/speciesInfo/species.inp");
pdr.set_underUbundant_fName        ("/home/mher/ism/speciesInfo/underabundant.inp");
pdr.set_rate99_fName               ("/home/mher/ism/speciesInfo/rate99.inp");
pdr.set_selfSheilding_CO_fName     ("/home/mher/ism/speciesInfo/self_shielding_CO.inp");
pdr.set_rotationalCooling_baseName ("/home/mher/ism/speciesInfo/rotationalcooling/rotcool");
pdr.set_vibrationalCooling_baseName("/home/mher/ism/speciesInfo/vibrationalcooling/vibcool");
pdr.set_database_fName             ("/home/mher/ism/database/z-1.0.dat");
pdr.set_zeta                       (5.0e-17);
pdr.set_S_depletion                (200.0);
pdr.set_TTol                       (1e-3);
pdr.set_CTol                       (1e-3);
pdr.set_metalicity                 (metallicity);
pdr.set_AvMax                      (30.0);
pdr.set_slabSizeCrit               (0.5);
pdr.set_min_deltaAv                (0.01);
pdr.set_max_deltaAv                (0.5);
pdr.set_maxSlabs                   (200);
#------------------------------------------------------------------------------------

# getting the basic species defined in baseSpecies.py
import baseSpecies
baseSpecs = baseSpecies.baseSpecies()

# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv( dirPath = databasePath, metallicity = metallicity )
arxv.readDb( check = True )
print 'time reading %f' % (time() - t0)

#------------------------------------------------------------------
# read and setting up the chemical network used in the 
t0 = time()
arxv.setupChemistry( parms = chemParms)
#-------------------------------------------------------------------
# plotting stuff

quantity       = ['therm', 'heating']
plotRange_nG0  = [[0,6],[0,6]]
slabIdx        = 0
res            = [100,100]
lgammaMechSec  = -50.0
log            = True

f   = arxv.construct3DInterpolationFunction(quantity = quantity, slabIdx = slabIdx)
grd = arxv.computeInterpolated2DGrid(quantity = quantity, 
                                     slabIdx  = slabIdx, 
                                     ranges   = plotRange_nG0,
                                     res      = res, 
                                     zSec     = lgammaMechSec, 
                                     fInterp  = f)

#pyl.ion()

arxv.showGrid(quantity = quantity, 
              slabIdx  = slabIdx, 
              ranges   = plotRange_nG0,
              res      = res, 
              zSec     = lgammaMechSec,
              log      = log)

x = np.array([log10nGas])
y = np.array([log10G0])
z = np.array([-50])
data = np.array([x,y,z]).T

lgSurf = f(data) 
print lgSurf 

rho   = log10nGas
G0    = log10G0
Lmech = lgSurf*alpha

# initializing variables and running the models
pdr.initialize()
ids,err = pdr.add_mesh(rho, 1000.0, G0, Lmech)

pdr.calc_equilibrium()
mshErr,err = pdr.get_errorFlags(ids)
T,err = pdr.get_temperature(ids)

print 'done'