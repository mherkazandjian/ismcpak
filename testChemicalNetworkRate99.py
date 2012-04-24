from numpy import *
from time import *
from chemicalNetwork import *
import sys

rxnFile       = '/home/mher/ism/code/ismcpak/data/rate99Fixed.inp'
specNumFile   = '/home/mher/ism/code/ismcpak/data/species.inp'
specAbunFile  = '/home/mher/ism/code/ismcpak/data/abun.out'
underAbunFile = '/home/mher/ism/code/ismcpak/data/underabundant.inp'
csvFile       = '/home/mher/ism/code/ismcpak/data/network.csv'

T        = 800.0
Av       = 20.0
nDens    = 10**5.5
G0       = 10**5
findReact = ''
findProd  = 'HCN'

zeta     = 5e-17
albedo   = 0.6


# elements and basic species from which all the other species are made
baseSpec = [  specie('CRPHOT', specType = -1, charge=0 , init=1),
              specie('PHOTON', specType = -1, charge=0 , init=1),
              specie('CRP'   , specType = -1, charge=0 , init=1),
              specie('PAH'   , specType = 0 , charge=0 , init=1),
              specie('H2V'   , specType = 0 , charge=0 , init=1),
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

t0 = time()
# settin up the orignial netowrk
net = chemicalNetwork(rxnFile, baseSpec, UMISTVER = 'umist99')
# setting the cloud parameters
net.setCloudParms(T, zeta, Av, albedo, nDens, G0)
# dumping the Gephi file for visualizing the network
net.writeNetworkGephiCSV(csvFile)        
# finding indentical reactions and pinting their IDs
sets = net.findIdenticalReactions()

for set in sets:
    for ind in set:
        print net.reactions[ind].str
    print '-----------------------------'
    
# set some species as inactive in the chemical network and assign new
# indicies from a file
#----------------------------------------------------------

# reading the species to be removed from a file
fObj =  open(underAbunFile, 'r')
setAsInactive = []
for line in fObj:
    line = line.split()
                
    specStr = line[1]
    setAsInactive.append(strip(specStr))
fObj.close() 

net.removeSpecies( setAsInactive )

# reading the species number and their corresponding indies and abundances from ascii files
net.assignNumbersToSpecies(fileName = specNumFile)

# in the original species.inp 13CH3 is missing, so we remove it from the network
# and re-assign numbers
net.removeSpecies( ['13CH3'] )
net.assignNumbersToSpecies(fileName = specNumFile)

# set the abundances
net.setAbundances( specAbunFile  )

print 'filter reactions with certain products and reactants'
inds = net.filterReactions(withReacts = ['PHOTON'], withProds = ['H', 'e-'])
net.printReactions(inds, format = 'type rxn')

print '-----------------------------------------------------------------------------' 
print 'filter reactions with certain reactants only'
inds = net.filterReactions(withReacts = ['HNC'])
net.printReactions(inds[0:5], format = 'type rxn')
print 'specie detailed info:'
net.species['HCN'].show()

#net.species['HCN'].hasBaseSpecies(['H'])
x =[ [spec[0], spec[1].hasBaseSpecies(['N', 'H'])] for spec in net.species.items()] 
for itm in x:
    if itm[1] == True:
        print itm 
print x

"""
print 'MESSING WITH THE NETWORK MANUALLY, JUST TO CHECK IF SPECIES ARE MAPPED CORRECTLY'
print '-----------------------------------------------------------------------------' 
print 'changing the index of a specie'
inds = net.filterReactions(withReacts = ['H2O'], withProds = ['e-'])
net.printReactions(inds, format = 'type rxn')
net.printReactions(inds, format = 'type rxnNumeric')
net.species['H2O'].num = 51
print 'printing now the numeric reactions wit hthe new number for H2O'
net.printReactions(inds, format = 'type rxnNumeric')

print '-----------------------------------------------------------------------------'
print 'changing the abundance of H throught the net.abun numpy array'
idxH = net.species['H'].getNum()
net.species['H'].show()
print idxH, net.abun[ idxH ]
print 'changing the abundance of H throught species.setAbun method'
net.abun[ idxH ][0] = 5.0
net.species['H'].show()
print idxH, net.abun[ idxH ]
net.species['H'].setAbun(55.0)
net.species['H'].show()
print idxH, net.abun[ idxH ]
"""

#net.computeReactionConstants()
"""
net.computeRates()

#net.printReactions(ids, format = "id type rxn cst rate" )

print 'time elapsed : ', time() - t0, 'sec'

ids = net.filterReactions(findReact, findProd)
idsSorted = net.sortRxnsDecreasingRates(ids)

print '################################################################'
net.printReactions(idsSorted, format = "id type rxn cst rate trng" )
"""