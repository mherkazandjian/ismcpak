import numpy
from chemicalNetwork import chemicalNetwork
import sys

#--------state of the gas--------------
T        = 800.0
Av       = 20.0
nDens    = 10**5.5
G0       = 10**5
zeta     = 5e-17
albedo   = 0.6
#--------------------------------------
rxnFile       = 'data/rate99Fixed.inp'
version       = 'umist99'

#rxnFile  = 'data/RATE06.txt'
#version  = 'umist06'

speciesNums      = 'data/species.inp'
specAbunFname    = 'data/abun.out'            #file containing the abundances in the same ordfer of speciesNumAndName.inp 
underAbunFile    = 'data/underabundant.inp'   #rxns with those species are scratched (moved to the removed list)
manualRemoveList = ['13CH3']                  #rxns with those species are scratched (moved to the removed list)
#------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

baseSpecies = __import__('baseSpeciesDefault')
baseSpecs = baseSpecies.baseSpecies()

net = chemicalNetwork(rxnFile, baseSpecs, UMISTVER = version)

# filter reactions with certain products and reactants
print '\n\nreactions with PHOTO as reactant and H and e- as products'
IDs = net.filter_reactions(withReacts = ['PHOTON'], withProds = ['H', 'e-'],
                           show = True, fmt = 'type rxn')

#finding indentical reactions and printing them
sets = net.find_identical_reactions()
print '\n\nthe following are the identical reactions in the rxn file'
for ids in sets:
    print 'these reactions have the same hash code : ', ids
    net.print_reactions(IDs=ids, fmt='id hash type rxn abg trng')
    print '-------------------------------------------------------'

print '-----------------------------------------------------------------------------' 
# filter reactions with certain reactants only
IDs = net.filter_reactions(withReacts = ['HNC'])
net.print_reactions(IDs[0:5], fmt = 'type rxn')

print '-----------------------------------------------------------------------------'
# finding some reactions involving H2 and H2O and changing the index of H2O
IDs = net.filter_reactions( withReacts=['H2'], withProds=['H2O'])
print 'reactions involving H2 and H2O'
net.print_reactions(IDs, fmt = 'type rxn')
print 'the numerical representation of the same reactions'
net.print_reactions(IDs, fmt = 'type rxnNumeric')

originalIndex = net.species['H2O'].num 
print 'changing the index of H2O from %d to 999' % net.species['H2O'].num 
net.species['H2O'].num = 999
print 'the numerical representation of the same reactions (with updated index of H2O)'
net.print_reactions(IDs, fmt = 'type rxnNumeric')
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

#assigning new numbers to the species from a file listing the species numbers and strings
print '------------------------------------------------------------------------------' 
net.assign_numbers_to_species(fileName = speciesNums)
# reading the species to be removed from a file
net.remove_species( underAbunFile = underAbunFile)
#removing a specie manually from the network
net.remove_species( species = manualRemoveList )
#re-assign numbers to species(this is not neccessary, just to check if 
#there are species which do not have a number and that remove_species
#does what is supposed to do 
net.assign_numbers_to_species(fileName = speciesNums) 

#setting the abundances from a file (after this, all the species in the 
#network (being used [not the one which was read since we modified it
#by removing some species]) should have a float value for the abundances
#excpet PHOTON, CRP, CRPHOT, M 
net.set_abundances(fromFile = specAbunFname) 

#printing all the reactions involving a CRP (with rates and rxn constants)
print 'reactions of type CP with decreasing reaction rates'
idsCP = net.filter_reactions(withType='CP', show=True, fmt='id type rxn abg')

print 'done'