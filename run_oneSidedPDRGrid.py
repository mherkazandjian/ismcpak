#---------------------------------------------------------------------------------------------------
import numpy
import pylab

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from enumSpecies import *
from ismUtils import *
import time, sys, os
import meshUtils

#---------------------------------------------------------------------------------------------------
# make sure to use the appropriate database when changing metallicities
#---------------------------------------------------------------------------------------------------
HOME      = '/home/mher'
outputDir = HOME + '/ism/runs/oneSided/foo-z-1.0/'
nWorker   = 7  # number of proccesses
pdr       = interface.pdrInterface( channel_type = 'mpi', number_of_workers = nWorker, redirection='none') 

exclude   = False  #False|True
exclude_meshes_from_db = HOME + '/ism/runs/oneSided/sph-db-z-1.0/'
#arxvExc = meshUtils.meshArxv(dirPath = exclude_meshes_from_db, readDb = True, set_defaults = True)
filter_gMech = None #None|[1e-3, 10.0]

metallicity = 1.0   # in terms of solar metallicity

pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              (HOME + "/ism/speciesInfo/species.inp");
pdr.set_underUbundant_fName        (HOME + "/ism/speciesInfo/underabundant.inp");
pdr.set_rate99_fName               (HOME + "/ism/speciesInfo/rate99.inp");
pdr.set_selfSheilding_CO_fName     (HOME + "/ism/speciesInfo/self_shielding_CO.inp");
pdr.set_rotationalCooling_baseName (HOME + "/ism/speciesInfo/rotationalcooling/rotcool");
pdr.set_vibrationalCooling_baseName(HOME + "/ism/speciesInfo/vibrationalcooling/vibcool");
pdr.set_database_fName             (HOME + "/ism/database/z-1.0.dat");
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

nx   =  10j   #5j   # log10 density
xMin = -3.0
xMax =  6.0

ny   =  10j  #5j   # log10 G0
yMin = -3.0  
yMax =  6.0

nz   =  26j   #51j
zMin = -35.0  #log10 mechanical heating
zMax = -10.0

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

if exclude:
    #excluding points which are up to 1e-5 in the l2 distance to meshes
    #in the previous database
    xex, yex, zex = arxvExc.grid_x, arxvExc.grid_y, arxvExc.grid_z
    
    cond = numpy.ones(x.shape)
    for i in numpy.arange(xex.size):
        diff = numpy.sqrt((x - xex[i])**2 + (y - yex[i])**2 + (z - zex[i])**2)
        
        inds = numpy.where(diff < 1e-5)[0]
        
        if inds.size > 0:    
            cond[inds] *= 0
    
    indsKeep = cond.nonzero()    
    x, y, z = x[indsKeep], y[indsKeep], z[indsKeep]

    if filter_gMech != None:
        #excluding points which are outside gmech bounds than the surfrace heating by a factors of 
        #filter_gMech. The surface heating with zero gmech is considered to be the one
        #of gmech = 1e-50
        indsKeep = []
        for i in numpy.arange(x.size):
            mesh_ref = arxvExc.get_mesh_data(x[i], y[i], -50.0)
            surfHeat = mesh_ref.data['therm']['heating'][0]
            gMechAdd = 10.0**z[i]
             
            if (gMechAdd/surfHeat >= filter_gMech[0]) and (gMechAdd/surfHeat <= filter_gMech[1]):
                indsKeep.append(i)
                print '%-+.2f %-+.2f %-+.2f %-+.2f' % (x[i], y[i], z[i], numpy.log10(surfHeat))
            
        x, y, z = x[indsKeep], y[indsKeep], z[indsKeep]
        
        
n = x.size

pylab.ioff()
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot( x, y, z, 'ro', markersize=4)
if exclude:
    ax.plot(xex, yex, zex, 'bo', markersize=7)
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
