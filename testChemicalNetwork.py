from numpy import *
from time import *
from chemicalNetwork import *
from ctypes import *

T        = 800.0
Av       = 20.0
nDens    = 10**5.5
G0       = 10**5
findReact = ''
findProd  = 'HCN'
specAbunFname = '../code/M4-GM-16-Av19.0.out'

zeta     = 5e-17
albedo   = 0.6
rxnFile  = 'data/RATE06.txt'
abunSpecFname = '../speciesInfo/speciesNames.inp'


# elements and basic species from which all the other species are made
baseSpec = [  specie('CRPHOT', specType = -1, charge=0 , init=1),
              specie('PHOTON', specType = -1, charge=0 , init=1),
              specie('CRP'   , specType = -1, charge=0 , init=1),
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
              specie('F'     , specType = 0 , charge=0 , init=1),
              specie('P'     , specType = 0 , charge=0 , init=1),
              specie('S'     , specType = 0 , charge=0 , init=1),
              specie('e-'    , specType = 0 , charge=-1, init=1) ]

t0 = time()
net = chemicalNetwork(rxnFile, baseSpec)

#net.printReactions([1,2,3])
#net.printReactions(format = "status id hash type rxn abg trng acc ref" )
#net.printSpecies()

a={}
print a
a = {'x': baseSpec[0] }
print a
a['y'] = 5
print a
a['y'] = 8
print a

print net.specDict

print 'z' in a
print 'x' in a

print net.specDict
"""""
net.setAbundancesFromFile(abunSpecFname, specAbunFname )

net.setCloudParms(T, zeta, Av, albedo, nDens, G0)
net.computeReactionConstants()
net.computeRates()

#net.printReactions(ids, format = "id type rxn cst rate" )

print 'time elapsed : ', time() - t0, 'sec'

ids = net.filterReactions(findReact, findProd)
idsSorted = net.sortRxnsDecreasingRates(ids)

print '################################################################'
net.printReactions(idsSorted, format = "id type rxn cst rate trng" )
"""