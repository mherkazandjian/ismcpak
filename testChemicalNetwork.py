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

rxnFile  = 'data/rate99Fixed.inp'
version  = 'umist99'

#rxnFile  = 'data/RATE06.txt'
#version  = 'umist06'
abunSpecFname = 'data/speciesNumAndName.inp'


baseSpecies = __import__('baseSpeciesDefault')
baseSpecs = baseSpecies.baseSpecies()


t0 = time()
net = chemicalNetwork(rxnFile, baseSpecs, UMISTVER = version)

# filter reactions with certain products and reactants
inds = net.filterReactions(withReacts = ['PHOTON'], withProds = ['H', 'e-'],
                           show = True, fmt = 'type rxn')

print '-----------------------------------------------------------------------------' 
# filter reactions with certain reactants only
inds = net.filterReactions(withReacts = ['HNC'])
net.printReactions(inds[0:5], fmt = 'type rxn')

print '-----------------------------------------------------------------------------'
# finding some reactions involving H2 and H2O and changing the index of H2O
inds = net.filterReactions( withReacts=['H2'], withProds=['H2O'])
print 'reactions involving H2 and H2O'
net.printReactions(inds, fmt = 'type rxn')
print 'the numerical representation of the same reactions'
net.printReactions(inds, fmt = 'type rxnNumeric')

originalIndex = net.species['H2O'].num 
print 'changing the index of H2O from %d to 999' % net.species['H2O'].num 
net.species['H2O'].num = 999
print 'the numerical representation of the same reactions (with updated index of H2O)'
net.printReactions(inds, fmt = 'type rxnNumeric')
net.species['H2O'].num = originalIndex # changing the index to what it was
print '-----------------------------------------------------------------------------'

#checking if the mapping of the abundances of net.species[:]._abun[0] to net._abun[:] has been
#done correctly by changing the abunance through net.abun[:] and then through net.species[:]._abun[0]
print 'setting the abundance of H to -44 through the net.species object'
net.species['H'].abun(-44.0)
print 'index H = %d, abun H = %e (from net.abun[])' % (net.species['H'].num, net.abun[net.species['H'].num])

print 'setting the abundance of H to -88 through the net.abun[]'
net.abun[net.species['H'].num] = -88.0
print "index H = %d, abun H = %e (from net.species['H'])" % (net.species['H'].num, net.abun[net.species['H'].num])

asdasd 
#net.computeReactionConstants()
#net.printReactions(fmt = 'rxn')
#net.printReactions(fmt = 'rxnNumeric')
#net.species['H'].num = -991
#net.printReactions(fmt = 'rxnNumeric')
#print 'done'

#net.printReactions(fmt = 'rxn')

net.assignNumbersToSpecies(fileName = 'data/speciesNumAndName.inp')

asdasd
rxn = net.reactions[0]
print rxn.str
rxn.display(fmt='rxnNumeric')
print rxn.reactants
print rxn.products

asdasda


#net.setAbundancesFromFile(abunSpecFname, specAbunFname )

"""""

net.setCloudParms(T, zeta, Av, albedo, nDens, G0)
net.computeReactionConstants()
net.computeRates()

#net.printReactions(ids, fmt = "id type rxn cst rate" )

print 'time elapsed : ', time() - t0, 'sec'

ids = net.filterReactions(findReact, findProd)
idsSorted = net.sortRxnsDecreasingRates(ids)

print '################################################################'
net.printReactions(idsSorted, fmt = "id type rxn cst rate trng" )
"""