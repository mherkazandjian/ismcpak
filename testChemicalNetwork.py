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
specAbunFname = 'data/M4-GM-16-Av19.0.out'

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
net = chemicalNetwork(rxnFile, baseSpec)

# filter reactions with certain products and reactants
inds = net.filterReactions(withReacts = ['PHOTON'], withProds = ['H', 'e-'])
net.printReactions(inds, format = 'type rxn')

print '-----------------------------------------------------------------------------' 
# filter reactions with certain reactants only
inds = net.filterReactions(withReacts = ['HNC'])
net.printReactions(inds[0:5], format = 'type rxn')

print '-----------------------------------------------------------------------------' 
# changing the index of a specie
net.printReactions([10], format = 'type rxn')
net.printReactions([10], format = 'type rxnNumeric')
net.species['H2O'].num = 51
net.printReactions([10], format = 'type rxn')
net.printReactions([10], format = 'type rxnNumeric')

net.abun = np.zeros(50)

net.species['H'].num  = 0
net.species['H'].abun = -44.0
print net.species['H'].abun

net.abun[0] = 55.0
net.species['H'].abun = net.abun[net.species['H'].num]
print net.species['H'].abun
 
#net.computeReactionConstants()
#net.printReactions(format = 'rxn')
#net.printReactions(format = 'rxnNumeric')
#net.species['H'].num = -991
#net.printReactions(format = 'rxnNumeric')
#print 'done'

#net.printReactions(format = 'rxn')

net.assignNumbersToSpecies(fileName = 'data/speciesNumAndName.inp')

asdasd
rxn = net.reactions[0]
print rxn.str
rxn.display(format='rxnNumeric')
print rxn.reactants
print rxn.products

asdasda


#net.setAbundancesFromFile(abunSpecFname, specAbunFname )

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