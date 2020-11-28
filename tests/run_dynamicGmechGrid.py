"""
<keywords>
example, pdr, amuse, dynamic, mechanical, heating
</keywords>
<description>
run a grid of PDR models using the amuse interface as a function of
n, G0 and gmech

run using the command:

    $ cd $ISMCPAK/test
    $ mpirun -np run_oneSidedPDRGrid.py

</description>
"""
#-------------------------------------------------------------------------------
import matplotlib
#matplotlib.use('Qt4Agg')
matplotlib.use('PS')

from numpy import *
from numpy.random import *
import pylab as pyl

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from enumSpecies import *
from ismUtils import *
from meshUtils import *
from time import *

HOME      =  os.environ['HOME']

dataDir   = HOME + '/ism/code/ismcpak/data/'
outputDir = HOME + '/ism/runs/tests/dynamicGrid/'

nWorker = 4  # number of proccesses
pdr     = interface.pdrInterface(channel_type = 'mpi', 
                                 number_of_workers = nWorker, 
                                 redirection='none') 

metallicity  =  1.0   # in terms of solar metallicity
plotRangenG0 = [[0,6],[0,6]]

# path of the database from which the surface mech heating 
# rates will be extracted
databasePath  = os.path.join(outputDir,'../oneSidedGrid_no_gmech/')

#---------chemical network parameters (not used in the modelling)-------------------
rxnFile       = HOME + '/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = HOME + '/ism/code/ismcpak/data/species.inp'
underAbunFile = HOME + '/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

#----------------amuse modeling parameters----------------------------------------
pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              (dataDir + 'pdr/species.inp');
pdr.set_underUbundant_fName        (dataDir + 'pdr/underabundant.inp');
pdr.set_rate99_fName               (dataDir + 'pdr/rate99.inp');
pdr.set_selfSheilding_CO_fName     (dataDir + 'pdr/self_shielding_CO.inp');
pdr.set_rotationalCooling_baseName (dataDir + 'pdr/rotationalcooling/rotcool');
pdr.set_vibrationalCooling_baseName(dataDir + 'pdr/vibrationalcooling/vibcool');
pdr.set_database_fName             (dataDir + 'pdr/z-1.0.dat');
pdr.set_zeta                       (5.0e-17);
pdr.set_S_depletion                (200.0);
pdr.set_TTol                       (1e-3);
pdr.set_CTol                       (1e-3);
pdr.set_metalicity                 (1.0);
pdr.set_AvMax                      (30.0);
pdr.set_slabSizeCrit               (0.5);
pdr.set_min_deltaAv                (0.01);
pdr.set_max_deltaAv                (0.5);
pdr.set_maxSlabs                   (200);

# getting the basic species defined in baseSpecies.py
import baseSpeciesDefault as baseSpecies
baseSpecs = baseSpecies.baseSpecies()

# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv(dirPath = databasePath, readDb = True)
print 'time reading %f' % (time() - t0)

#--------------------------grid point to be modelled-------------------------------
dx   = 3.0      # log10 density
xMin = 0.0
xMax = 6.01

dy   = 3.0     # log10 G0
yMin = 0.0  
yMax = 6.01

# factor of surface heating to be added as mechanical heating
#z = [1e-10, 0.001, 0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 1.0,]
z = [1e-10, 0.01, 0.05, 0.1,]
#z = [0.0001, 0.001]

#-----------------------getting the mech heating rates at grid points--------------
x, y = mgrid[ xMin:(xMax+1e-10):dx, yMin:(yMax+1e-10):dy]
grid_shape = x.shape

x = x.flatten()
y = y.flatten()

# getting the interpolation function for the surface heating
nPts = x.size
gMechZero = x.copy()

gMechZero[:] = -50.0 # lowest mechanical energy used (in log)
f         = arxv.construct3DInterpolationFunction(quantity = ['therm', 'heating'], slabIdx  = 0, log10 = True)
dataNew   = np.array( [x, y, gMechZero] ).T
gammaSurf = f(dataNew)

print np.array([x, y, gammaSurf]).T

inds = np.where( np.isnan(gammaSurf) == False  )
x = x[ inds ]
y = y[ inds ]
gammaSurf = gammaSurf[ inds ]

xTmp, yTmp, zTmp = [], [], []

for v in z:
    xTmp.append(x)
    yTmp.append(y)

    gMechGrid = np.log10(v * (10.0**gammaSurf) )
    zTmp.append( gMechGrid )
    
x = np.array(xTmp).flatten()
y = np.array(yTmp).flatten()
z = np.array(zTmp).flatten()

n = x.size

pyl.ion()
figModels = pyl.figure()
axsModels = figModels.add_axes([0.1, 0.1, 0.8, 0.8])

fig1 = pyl.Figure()
axsModels.plot( x, y, 'o' )
axsModels.set_xlim( [-1, 7] )
axsModels.set_ylim( [-1, 7] )
pyl.draw()

print 'number of models to run = ', n

rho   = x
G0    = y
Lmech = 10**z

# initializing variables and running the models
pdr.initialize()
ids,err = pdr.add_mesh(rho, 1000.0, G0, Lmech)

# writing the parameter into an ascii file
f = file(outputDir+'parameters.out', 'w')
for i in arange(n):
    f.write( '%06d  %.3f  %.3f  %.3f\n' % (ids[i], rho[i], G0[i], Lmech[i]))
f.close()

#running the models
pdr.calc_equilibrium()
mshErr,err = pdr.get_errorFlags(ids)
T,err = pdr.get_temperature(ids)

# writing the results to an ascii file
if 1:
    f = file(outputDir+'results.out', 'w')
    for i in arange(n):
        f.write( '%d %.3f %.3f %.3e %.d\n' % (ids[i], rho[i], G0[i], Lmech[i], mshErr[i]))
    f.close()

print 'done'