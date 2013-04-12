# in this test file, we load a database of PDR meshes and extract information
# at conditions determined by the SPH simulations
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

home = '/home/mher'

#path to sph snapshot (for the time being the .npz file
#containing the n,G0,gm, x,y,z of the gas particles
fPathDataSPH = home + '/ism/runs/galaxies/test1/test/xyz.npz'

#PDR database file
pdrDatabaseDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-1.0/'
metallicity        = 1.0
#-----------------chemical network parameters------------------------
chemParms  = {
              'rxnFile'       : home + '/ism/code/ismcpak/data/rate99Fixed.inp',
              'specNumFile'   : home + '/ism/code/ismcpak/data/species.inp',
              'underAbunFile' : home + '/ism/code/ismcpak/data/underabundant.inp',
              'removeManual'  : ['13CH3'],
              'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
              'umistVer'      : 'umist99',
              }
############################################################################
############################################################################
########################importing stuff#####################################
from numpy import *
import time, sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab
from meshUtils import *
from scipy      import interpolate
############################################################################
#################reading and setting up the pdr database####################
############################################################################

# reading the archive
print 'setting up the archive'
t0 = time()
arxvPDR = meshArxv(dirPath = pdrDatabaseDirPath )
arxvPDR.readDb(check = True)
print 'time reading data %f' % (time() - t0)

# read and setting up the chemical network used in the 
t0 = time()
arxvPDR.setupChemistry( parms = chemParms)
print 'time setting up the chemistry %f' % (time() - t0)

arxvPDR.set_metallicity(metallicity)
arxvPDR.set_default_attributes()
arxvPDR.set_grid_axes_quantity_values()
############################################################################
##################loading the sph simulation data###########################
############################################################################

f = np.load(fPathDataSPH)
lx, ly, lz = f['lx'], f['ly'], f['lz']
px, py, pz = f['px'], f['py'], f['pz']
m, tKin, h2frac = f['m'], f['tKin'], f['h2frac']

############################################################################
########### making an interpolation function for a certain quantitiy########
############################################################################

#getting some quantity from the meshes
q = np.zeros( arxvPDR.nMeshes)

for pdrMeshIdx, m in enumerate(arxvPDR.meshes):
    specIdx = arxvPDR.chemNet.species['e-'].num
    #q[pdrMeshIdx] = m['state']['abun'][specIdx][0]
    q[pdrMeshIdx] = m['state']['gasT'][0]

data = np.array([arxvPDR.grid_x, arxvPDR.grid_y, arxvPDR.grid_z]).T  #3D coordinates
        
f = interpolate.LinearNDInterpolator(data, q) # getting the interpolation function

#data = np.array([ [nGas], [G0], [gMech] ]).T

#print f(data)

#result=gridinterpolator("gasT", dens, G0, Gmech)

