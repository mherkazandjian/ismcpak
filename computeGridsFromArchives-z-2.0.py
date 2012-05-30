from numpy import *
from time import *
import sys
import pylab as pyl

from chemicalNetwork import *
from mesh import *
from meshUtils import *
from enumSpecies import *

# runs
#   uniformSweep2-z-2/
#   uniformSweep2-highRes-z-1.0
#   uniformSweep2-z-0.5
#   
#   uniformSweep2-z-1.0-no-mech
#   uniformSweep2-z-2-no-mech
#   uniformSweep2-z-0.5-no-mech
# 
#---------------------------Archive parameters-----------------------
#runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'

def runGrid(runDirPath, metallicity, lgammaMechSec, RadexSpecStr, RadexTransition ):
    #runDirPath    = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/'
    gridsRes      = 10
    #lgammaMechSec = -30.0
    #metallicity   = 2.0
    
    radexParms    = { 'radexPath'         : '/home/mher/ism/code/radex/Radex/bin/radex',  
                      'molDataDirPath'    : '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                      'specStr'           : RadexSpecStr,
    #                  'xH2_Min'          : 2*0.0000000001
                      'xH2_Min'           : -1.0,
                      'collisionPartners' : ['H2','H+','H'],
    #                  'collisionPartners' : ['H2','H','H+','e-']
    #                  'collisionPartners' : ['H2']
    #                  'collisionPartners' : ['H2','H+','e-','H']
                      'plotTransitionInGrid' : RadexTransition
                    }
    #-----------------chemical network parameters------------------------
    rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
    specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
    underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
    removeManual  = ['13CH3']
    
    
    # elements and basic species from which all the other species are made
    baseSpec = [  specie('CRPHOT', specType = -1, charge=0 , init=1),
                  specie('PHOTON', specType = -1, charge=0 , init=1),
                  specie('CRP'   , specType = -1, charge=0 , init=1),
                  specie('PAH'   , specType = 0 , charge=0 , init=1),
                  specie('H2V'   , specType = 0 , charge=0 , init=1, comp = [ ['H',2] ]),
                  specie('13C'   , specType = 0 , charge=0 , init=1),
                  specie('Na'    , specType = 0 , charge=0 , init=1),
                  specie('Mg'    , specType = 0 , charge=0 , init=1),
                  specie('Si'    , specType = 0 , charge=0 , init=1),
                  specie('Cl'    , specType = 0 , charge=0 , init=1),
                  specie('Fe'    , specType = 0 , charge=0 , init=1),
                  specie('He'    , specType = 0 , charge=0 , init=1),
                  specie('H'     , specType = 0 , charge=0 , init=1),
                  specie('M'     , specType = -1, charge=0 , init=1),
                  specie('C'     , specType = 0 , charge=0 , init=1),
                  specie('N'     , specType = 0 , charge=0 , init=1),
                  specie('O'     , specType = 0 , charge=0 , init=1),
                  specie('P'     , specType = 0 , charge=0 , init=1),
                  specie('S'     , specType = 0 , charge=0 , init=1),
                  specie('e-'    , specType = 0 , charge=-1, init=1) ]
    #------------------------------------------------------------------
    # reading the archive
    print 'setting up the archive'
    t0 = time()
    arxv = meshArxv( metallicity = metallicity )
    arxv.readDb( runDirPath )
    arxv.checkIntegrity()
    print 'time reading %f' % (time() - t0)
    #------------------------------------------------------------------
    # read and setting up the chemical network used in the 
    t0 = time()
    # settin up the orignial netowrk
    net = chemicalNetwork(rxnFile, baseSpec, UMISTVER = 'umist99')
    # reading the species to be removed from a file
    net.removeSpecies( underAbunFile = underAbunFile )
    net.removeSpecies( species = removeManual )
    # reading the species number and their corresponding indies and abundances from ascii files
    net.assignNumbersToSpecies(fileName = specNumFile)
    arxv.setChemicalNetwork(net) # assiginig the chemical network to the archive
    #-------------------------------------------------------------------
    # plotting stuff
    
    arxv.saveGridsToFiles(gridsRes, lgammaMechSec, radexParms)
    print 'done'
#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

species = (
              ('CN',[0,1,2,5,9,15]), 
          )
#----------------------------------------------------------------------------
z = 2.0
archive = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2-no-mech/'
for (specStr, transitions) in species:
    for transition in transitions:
        print specStr, transition
        runGrid(archive, z, -30.0, specStr, transition )
archive = '/home/mher/ism/runs/oneSided/uniformSweep2-z-2/'
for gm in [-24.0, -23.0, -22.0, -21.0, -20.0, -19.0, -18.0, -17.0, -16.0]:
    for (specStr, transitions) in species:
        for transition in transitions:
            print specStr, transition
            runGrid(archive, z, gm, specStr, transition )
