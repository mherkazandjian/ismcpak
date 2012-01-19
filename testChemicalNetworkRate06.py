from numpy import *
from time import *
from chemicalNetwork import *
import sys

T        = 800.0
Av       = 20.0
nDens    = 10**5.5
G0       = 10**5
findReact = ''
findProd  = 'HCN'
specAbunFname = 'data/abun.out'

zeta     = 5e-17
albedo   = 0.6
rxnFile  = 'data/RATE06.txt'
abunSpecFname = 'data/speciesNumAndName.inp'


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
net = chemicalNetwork(rxnFile, baseSpec, UMISTVER = 'umist06')


# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = 'data/speciesNumAndName.inp')
net.setAbundances( specAbunFname  )
#net.printSpecies()

#net.computeReactionConstants()
#net.printReactions(format = 'rxn')
#net.printReactions(format = 'rxnNumeric')
#net.species['H'].num = -991
#net.printReactions(format = 'rxnNumeric')
#print 'done'

#net.printReactions(format = 'rxn')

#print id(net.)

print 'filter reactions with certain products and reactants'
inds = net.filterReactions(withReacts = ['PHOTON'], withProds = ['H', 'e-'])
net.printReactions(inds, format = 'type rxn')

print '-----------------------------------------------------------------------------' 
print 'filter reactions with certain reactants only'
inds = net.filterReactions(withReacts = ['HNC'])
net.printReactions(inds[0:5], format = 'type rxn')
net.species['HCN'].show()

print '-----------------------------------------------------------------------------' 
print 'changing the index of a specie'
net.printReactions([10], format = 'type rxn')
net.printReactions([10], format = 'type rxnNumeric')
net.species['H2O'].num = 51
net.printReactions([10], format = 'type rxn')
net.printReactions([10], format = 'type rxnNumeric')


print 'changing the abundance of H throught the net.abun numpy array'
idxH = net.species['H'].getNum()
print net.species['H'].getAbun()
print idxH, net.species['H'].getAbun(), net.abun[ idxH ]
print 'changing the abundance of H throught species.setAbun method'
net.abun[ idxH ][0] = 5.0  
print idxH, net.species['H'].getAbun(), net.abun[ idxH ]
net.species['H'].setAbun(55.0)
print idxH, net.species['H'].getAbun(), net.abun[ idxH ]

"""""
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