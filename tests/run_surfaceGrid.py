"""
<keywords>
example, pdr, amuse, static, surface, heating, grid
</keywords>
<description>
run a grid of PDR models to compute the surface heating. The purpose of
this grid is to be used as input for other grids where the mechanical
heating will be added as a fraction of the surface heating that is
computed in this grid.

The following is done in this script:

   - run the grid
   - pack the models into an models archive (database)

run using the command:

    $ cd $ISMCPAK/test
    $ mpirun -np 1 run_surfaceGrid.py

</description>
"""
#--------------------------------------------------------------------------------------------------
import os
import sys
import time
import pickle
import matplotlib
matplotlib.use('PS')
from numpy import *
from numpy.random import *
import pylab as pyl

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from enumSpecies import *
from ismUtils import *
from meshUtils import meshArxv

nWorker = 1  # number of proccesses
pdr     = interface.pdrInterface( channel_type = 'mpi', number_of_workers = nWorker, redirection='none')

metallicity = 1.0   # in terms of solar metallicity
dataDir   = '/ism/ismcpak/data/'
outputDir   = '../../data/oneSided/surfaceGrid-z-%.1f-no-gmech/' % metallicity

# create the directory structure of the output
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
meshes_dir = os.path.join(outputDir, 'meshes')
if not os.path.isdir(meshes_dir):
    os.makedirs(meshes_dir)

pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              (dataDir + "pdr/species.inp");
pdr.set_underUbundant_fName        (dataDir + "pdr/underabundant.inp");
pdr.set_rate99_fName               (dataDir + "pdr/rate99Fixed.inp");
pdr.set_selfSheilding_CO_fName     (dataDir + "pdr/self_shielding_CO.inp");
pdr.set_rotationalCooling_baseName (dataDir + "pdr/rotationalcooling/rotcool");
pdr.set_vibrationalCooling_baseName(dataDir + "pdr/vibrationalcooling/vibcool");
pdr.set_database_fName             (dataDir + "pdr/z-%.1f.dat" % metallicity);
pdr.set_zeta                       (5.0e-17);
pdr.set_S_depletion                (200.0);
pdr.set_TTol                       (1e-3);
pdr.set_CTol                       (1e-3);
pdr.set_metalicity                 (metallicity);
pdr.set_AvMax                      (0.01);
pdr.set_slabSizeCrit               (0.5);
pdr.set_min_deltaAv                (0.01);
pdr.set_max_deltaAv                (0.5);
pdr.set_maxSlabs                   (100);

# ======== dump parameter to be used in the analysis to a pickle file =========
# .. todo:: due to some formatting issues/inconvience the files used by ismcpak
# are slightly different from the ones used the C PDF code, hence the slightly
# different paths in the dict below
parms = {
    'dirPath'      : outputDir,
    'relativeGmech': False,
    'metallicity'  : metallicity,
    'chemistry'    : {
                      'rxnFile'       : dataDir + "rate99Fixed.inp",
                      'specNumFile'   : dataDir + "species.inp",
                      'underAbunFile' : dataDir + "underabundant.inp",  # .. todo:: fix this in chemicalNetwork.py / ignore first line
                      'removeManual'  : ['13CH3'],
                      'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
                      'umistVer'      : 'umist99',
                     }
}
pickle.dump(parms, open(os.path.join(outputDir, 'used_parms.pkl'), 'w'))
# ========================= done dumping parameter ============================

# define the grid parameters of the models
dx   = 0.5      # log10 density
xMin = 0.0
xMax = 6.0

dy   = 0.5     # log10 G0
yMin = 0.0
yMax = 6.0

zMin = -50.0    # log10 mechanical heating
zMax = -49.0
dz   =  1.0

# generating the parameter space
if zMin == zMax:
    x, y = mgrid[ xMin:(xMax+1e-10):dx, yMin:(yMax+1e-10):dy]
    z = ones(x.size)*zMin
else:
    x, y, z = mgrid[ xMin:(xMax+1e-10):dx, yMin:(yMax+1e-10):dy, zMin:(zMax+1e-10):dz]

x = x.flatten()
y = y.flatten()
z = z.flatten()

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
    f.write( format(int(ids[i]), '06') + ' ' + format(rho[i], '.3f') + ' ' + format(G0[i], '.3f') + ' ' + format(Lmech[i], '.3e') + '\n' )
f.close()


#running the models
pdr.calc_equilibrium()
mshErr,err = pdr.get_errorFlags(ids)
T,err = pdr.get_temperature(ids)

# writing the results to an ascii file
f = file(outputDir+'results.out', 'w')
for i in arange(n):
    f.write( '%d %.3f %.3f %.3e %.d\n' % (ids[i], rho[i], G0[i], Lmech[i], mshErr[i]))
f.close()

# constructing the archive
print 'pack the meshes into a database'
t0 = time.time()
arxvW = meshArxv(dirPath=outputDir)
arxvW.construct(meshNamePrefix='mesh', writeDb=True)
print 'time constructing %f' % (time.time() - t0)

print 'done'
