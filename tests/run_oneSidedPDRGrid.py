#---------------------------------------------------------------------------------------------------
import numpy
import pylab

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from ismUtils import *
import time, sys, os
import meshUtils

#---------------------------------------------------------------------------------------------------
# make sure to use the appropriate database when changing metallicities
#---------------------------------------------------------------------------------------------------
HOME      = os.environ['HOME']

dataDir   = HOME + '/ism/code/ismcpak/data/'
outputDir = HOME + '/ism/runs/tests/oneSidedGrid/'

nWorker   = 4  # number of proccesses
pdr       = interface.pdrInterface( channel_type = 'mpi', number_of_workers = nWorker, redirection='none') 

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

nx   =  3j   #5j   # log10 density
xMin =  0.0
xMax =  6.0

ny   =  3j  #5j   # log10 G0
yMin =  0.0  
yMax =  6.0

nz   =  2j    #51j
zMin = -25.0  #log10 mechanical heating
zMax = -23.0

###########################################################################################
###########################################################################################
###########################################################################################
# generating the parameter space 
x, y, z = numpy.mgrid[xMin:xMax:nx, yMin:yMax:ny, zMin:zMax:nz]
x, y, z = x.flatten(), y.flatten(), z.flatten()

print 'unique x'
print numpy.unique(x)
print 'unique y'
print numpy.unique(y)
print 'unique z'
print numpy.unique(z)

n = x.size

pylab.ioff()
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot( x, y, z, 'ro', markersize=4)
ax.set_xlim( [-4, 8] )
ax.set_ylim( [-4, 8] )
ax.set_zlim( [-51, -15] )

ax.set_xlabel('n')
ax.set_ylabel('g0')
ax.set_zlabel('gm')

pylab.draw()
pylab.show()
print 'number of models to run = ', n

rho   = x
G0    = y
Lmech = 10.0**z

# initializing variables and running the models
pdr.initialize()
ids,err = pdr.add_mesh(rho, 1000.0, G0, Lmech)

# writing the parameter into an ascii file
f = file(outputDir+'parameters.out', 'w')
for i in numpy.arange(n):
    f.write( format(int(ids[i]), '06') + ' ' + format(rho[i], '.3f') + ' ' + format(G0[i], '.3f') + ' ' + format(Lmech[i], '.3e') + '\n' )
f.close()

#running the models
pdr.calc_equilibrium()
mshErr,err = pdr.get_errorFlags(ids)
T,err = pdr.get_temperature(ids)

# writing the results to an ascii file
f = file(outputDir+'results.out', 'w')
for i in numpy.arange(n):
    f.write( '%d %.3f %.3f %.3e %.d\n' % (ids[i], rho[i], G0[i], Lmech[i], mshErr[i]))
f.close()

print 'done'
