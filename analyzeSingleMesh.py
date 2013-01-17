from numpy import *
from time import *
import sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
from meshUtils import *

home = '/home/mher'
#------------------mesh data file-------------------------------------
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000-0.1.7'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/meshes/mesh.dat-id-000000'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat'
#meshFname      = '/data1/mher/ism/runs/oneSided/testOneSidedPDRGrid/meshes/mesh.dat-id-000000'
#meshFname      = '/data1/mher/ism/runs/oneSided/runSingleMesh/meshes/mesh.dat-10000percent'
#meshFname     = '/home/mher/ism/runs/oneSided/tests/M2-30/mesh.dat-adaptive-0.1.10'
#meshFname     = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/meshes/mesh.dat-id-000119'
meshFname     = '/home/mher/ism/runs/oneSided/testOneSidedPDRGrid/meshes/mesh.dat-id-000000'
metallicity   = 0.1

#-----------------chemical network parameters------------------------
chemParms  = {
              'rxnFile'       : home + '/ism/code/ismcpak/data/rate99Fixed.inp',
              'specNumFile'   : home + '/ism/code/ismcpak/data/species.inp',
              'underAbunFile' : home + '/ism/code/ismcpak/data/underabundant.inp',
              'removeManual'  : ['13CH3'],
              'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
              'umistVer'      : 'umist99',
              }

#importing the module which holds the definitions of the base species
baseSpecies = __import__(chemParms['baseSpecies'])
baseSpecs = baseSpecies.baseSpecies()

#------------------------------------------------------------------
net = chemicalNetwork.chemicalNetwork(chemParms['rxnFile'], 
                                      baseSpecs, 
                                      UMISTVER = 'umist99')

# reading the species to be removed from a file
net.removeSpecies( underAbunFile = chemParms['underAbunFile'] )
net.removeSpecies( species = chemParms['removeManual'] )
# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = chemParms['specNumFile'])

# setting up the PDR slab
m = mesh(meshFname, net, metallicity)
TMean1, nColls, NCOLVG1, = m.getRadexParameters('CO', -1.0)
print TMean1, nColls, NCOLVG1

m = mesh( meshFname, net, metallicity )
m.plot()

pyl.show()
print 'done'

# 378.37764026 1182.21211321 8.05680609553e+18