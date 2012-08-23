# the models of paper1 (table-1)
#---------------------------------------------------------------------------------------------------
from numpy import *
from numpy.random import *
import pylab as pyl

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from enumSpecies import *
from ismUtils import *
import time

nWorker = 3  # number of proccesses
pdr     = interface.pdrInterface( number_of_workers = nWorker, redirection='none') 

outputDir   = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/'
metallicity = 1.0   # in terms of solar metallicity
 
pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              ("/home/mher/ism/speciesInfo/species.inp");
pdr.set_underUbundant_fName        ("/home/mher/ism/speciesInfo/underabundant.inp");
pdr.set_rate99_fName               ("/home/mher/ism/speciesInfo/rate99.inp");
pdr.set_selfSheilding_CO_fName     ("/home/mher/ism/speciesInfo/self_shielding_CO.inp");
pdr.set_rotationalCooling_baseName ("/home/mher/ism/speciesInfo/rotationalcooling/rotcool");
pdr.set_vibrationalCooling_baseName("/home/mher/ism/speciesInfo/vibrationalcooling/vibcool");
pdr.set_database_fName             ("/home/mher/ism/database/database2.dat");
pdr.set_zeta                       (5.0e-17);
pdr.set_S_depletion                (200.0);
pdr.set_TTol                       (1e-3);
pdr.set_CTol                       (1e-3);
pdr.set_metalicity                 (metallicity);
pdr.set_AvMax                      (20.0);
pdr.set_slabSizeCrit               (0.5);
pdr.set_min_deltaAv                (0.01);
pdr.set_max_deltaAv                (0.5);
pdr.set_maxSlabs                   (200);


dx   = 0.5      # log10 density
xMin = 0.0
xMax = 6.0

dy   = 0.5     # log10 G0
yMin = 0.0  
yMax = 6.0

# factor of surface heating to be added as mechanical heating
z = [0.01, 0.05, 0.1, 0.25, 0.50, 0.75, 1.0, 10.0, 1e2, 1e3, 1e4]

# generating the parameter space 
xg = arange(xMin, xMax, dx)
yg = arange(yMin, yMax, dy)
zg = arange(len(z))
x, y, z = mgrid[ xMin:(xMax+1e-10):dx, yMin:(yMax+1e-10):dy, 0:len(z):1]

###### replace the values of z with the corresponding ones in z
###### see what linespec does

nx=xg.size
ny=yg.size
nz=zg.size

x = x.flatten()
y = y.flatten()
z = z.flatten()

n = nx*ny*nz

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
####################Lmech = 10**z
######## get the values of Lmech from the database (see kkkkkkkkkk below)

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
if 1:
    f = file(outputDir+'results.out', 'w')
    for i in arange(n):
        f.write( format(ids[i]    , 's'  ) + ' ' + 
                 format(rho[i]    , '.3f') + ' ' + 
                 format(G0[i]     , '.3f') + ' ' + 
                 format(Lmech[i]  , '.3e') + ' ' + 
                 format(T[i]      , '.8f') + ' ' + 
                 format(mshErr[i] , 's'  ) + '\n')
    f.close()

print 'done'
###########################################
###########################################
###########################################kkkkkkkkkkkk
###########################################
###########################################
#---------------------------------------------------------------------------

runDirPath    = '/home/mher/ism/runs/oneSided/testDynamicSingleMesh/'
databasePath  = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid4/'

gridsRes      = 10
metallicity   = 1.0
plotRangenG0  = [[0,6],[0,6]]
#-----------------chemical network parameters------------------------
rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
removeManual  = ['13CH3']

nWorker = 3  # number of proccesses
pdr     = interface.pdrInterface( number_of_workers = nWorker, redirection='none') 

outputDir   = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/'
metallicity = 1.0   # in terms of solar metallicity
 
#---------------------amuse parameters------------------------------------------------
pdr.set_outputDir                  (outputDir + 'meshes/');
pdr.set_species_fName              ("/home/mher/ism/speciesInfo/species.inp");
pdr.set_underUbundant_fName        ("/home/mher/ism/speciesInfo/underabundant.inp");
pdr.set_rate99_fName               ("/home/mher/ism/speciesInfo/rate99.inp");
pdr.set_selfSheilding_CO_fName     ("/home/mher/ism/speciesInfo/self_shielding_CO.inp");
pdr.set_rotationalCooling_baseName ("/home/mher/ism/speciesInfo/rotationalcooling/rotcool");
pdr.set_vibrationalCooling_baseName("/home/mher/ism/speciesInfo/vibrationalcooling/vibcool");
pdr.set_database_fName             ("/home/mher/ism/database/database2.dat");
pdr.set_zeta                       (5.0e-17);
pdr.set_S_depletion                (200.0);
pdr.set_TTol                       (1e-3);
pdr.set_CTol                       (1e-3);
pdr.set_metalicity                 (metallicity);
pdr.set_AvMax                      (20.0);
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
arxv = meshArxv( metallicity = metallicity )
arxv.readDb( databasePath )
arxv.checkIntegrity()
print 'time reading %f' % (time() - t0)

#------------------------------------------------------------------
# read and setting up the chemical network used in the 
t0 = time()
# settin up the orignial netowrk
net = chemicalNetwork(rxnFile, baseSpecs, UMISTVER = 'umist99')
# reading the species to be removed from a file
net.removeSpecies( underAbunFile = underAbunFile )
net.removeSpecies( species = removeManual )
# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = specNumFile)
arxv.setChemicalNetwork(net) # assiginig the chemical network to the archive
#-------------------------------------------------------------------
# plotting stuff

#pyl.ioff()
#arxv.plotGrid(gridsRes, lgammaMechSec, radexParms)

quantity       = ['therm', 'heating']
plotRange_nG0  = [[0,6],[0,6]]
slabIdx        = 0
res            = [100,100]
lgammaMechSec  = -30.0
log            = True

f   = arxv.construct3DInterpolationFunction(quantity = quantity, slabIdx = slabIdx)
grd = arxv.computeInterpolated2DGrid(quantity = quantity, 
                                     slabIdx  = slabIdx, 
                                     ranges   = plotRange_nG0,
                                     res      = res, 
                                     zSec     = -30, 
                                     fInterp  = f)

#pyl.ion()

arxv.showGrid(quantity = quantity, 
              slabIdx  = slabIdx, 
              ranges   = plotRange_nG0,
              res      = res, 
              zSec     = -30,
              log      = log)

lnModel  = 5.0
lG0Model = 3.0
x = np.array([lnModel])
y = np.array([lG0Model])
z = np.array([-30])
data = np.array([x,y,z]).T

lgSurf = f(data) 
print lgSurf 

"""
rho   = lnModel
G0    = lG0Model
Lmech = lgSurf

# initializing variables and running the models
pdr.initialize()
ids,err = pdr.add_mesh(rho, 1000.0, G0, Lmech)

pdr.calc_equilibrium()
mshErr,err = pdr.get_errorFlags(ids)
T,err = pdr.get_temperature(ids)
"""
print 'done'
###########################################
###########################################
###########################################
###########################################