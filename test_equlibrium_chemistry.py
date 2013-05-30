import numpy
from chemicalNetwork import chemicalNetwork
import sys

##
##
#--------state of the gas--------------
state = {'T'          : 400.0,
        'Av'          : 20.0,
        'nDens'       : 10**5.5,
        'G0'          : 10**5,
        'zeta'        : 5e-17,
        'albedo'      : 0.6,
        'Tdust'       : 300.0,
        'metallicity' : 1.0,
        'PHI_PAH'     : 0.5,
        'beta_CO'     : 0.8,
        'beta_13CO'   : 0.8,
        'beta_H2'     : 0.8}
#--------------------------------------
rxnFile       = 'data/rate99Fixed.inp'
version       = 'umist99'

speciesNums      = 'data/species.inp'
specAbunFname    = 'data/abun.out'            #file containing the abundances in the same ordfer of speciesNumAndName.inp 
underAbunFile    = 'data/underabundant.inp'   #rxns with those species are scratched (moved to the removed list)
manualRemoveList = ['13CH3']                  #rxns with those species are scratched (moved to the removed list)
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

baseSpecies = __import__('baseSpeciesDefault')
baseSpecs = baseSpecies.baseSpecies()

net = chemicalNetwork(rxnFile, baseSpecs, UMISTVER = version)

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

print '------------------------------------------------------------------------------'

#scratching out reactions  (check if 3707 is removed in the original PDR code)
net.pop_reactions( [799, 3706, 3707, 4076, 4412, 4413] )
net.pop_reactions( net.filter_reactions(withType='TR') )

#changing temperature ranges of some of the reactions
net.setattr_rxn(5   , 'Tu', 700.0 ); net.setattr_rxn(6   , 'Tl', 700.0);
net.setattr_rxn(240 , 'Tu', 4150.0); net.setattr_rxn(241 , 'Tl', 4150.0);
net.setattr_rxn(240 , 'Tu', 4150.0); net.setattr_rxn(241 , 'Tl', 4150.0);

#setting the lower bounds of the temperature of some reactions
# [[ID1, Tlb1], [ID2, Tlb2]...[]]
ID_Tlb = [ [63 , 10.0 ], [64, 10.0  ], [80, 300.0 ],  [131, 100.0], 
           [196, 300.0], [210, 200.0], [234, 200.0],  [236, 200.0],
           [273, 200.0], [275, 300.0], [298, 200.0],  [337, 250.0],
           [343, 200.0], [348, 500.0], [354, 80.0 ],  [360, 200.0],
           [363, 300.0], [366, 200.0], [374, 25.0 ],  [379, 150.0],
           [387, 300.0], [393, 13.0 ], [400, 300.0],  [414, 200.0],
           [425, 200.0], [432, 200.0], [4278,300.0]
         ]

for ID, Tlb in ID_Tlb:
    print 'setting a lower bound to the T to be used in computing the rxn cst'
    net.setattr_rxn(ID, 'Tlb', Tlb)
    print '       ',
    net.print_reactions(ID, fmt='type rxn trng Tlb')

#in looking at chemcial_balance.c we see that all the reactions which have no
#Tl and Tu set, the rates are valid for all temperatures, so we set those to
#the min and max of UMIST Tl=10 and Tu=41000
for rxn in net.reactions:
    if rxn.Tl == -1 and rxn.Tu == -1:
        rxn.Tl = 10.0
        rxn.Tu = 41000.0
        print 'set trange manually ', 
        rxn.display(fmt='rxn trng')

#setting the functions which compute the reactions constants
print 'setting the functions which will compute the reaction constants'
net.set_rxn_cst_functions()

#checking if all the reactions have the functions set
found = 0
for rxn in net.reactions:
    if rxn.compute_rxn_cst_func == None:
        found = 0
        rxn.display(fmt='id rxn')
        
if found > 0:
    raise ValueError('one or more reaction does not have a rxn cst function set...see above..')
else:
    print 'reaction constant functions for all the reactions are set'
    
#now that the ranges of some of the reactions are set/modified, we can merge them
#and set the functions which compute the constants give a temperature 
net.merge_identical_reactions()
net.set_all_rxn_bounds()

#setting the state of the gass
net.set_environment_state(**state)

##################DONE FULLY SETTING UP THE NETWORK###############################
####NOW WE PROCEED BY COMPUTING THE JACOBIAN AND THE EQUILIBRIUM SOLUTION#########

net.compute_rxn_constants()
net.compute_rxn_rates()

#computing the rate of change of the abundance of each active species. Each entry with index
#'i' corresponds to a species in net.species whose number is species[SPECSTR].num
nSpecsActive = net.nSpecsActive
dxdt  = numpy.zeros(nSpecsActive, 'f8')                 #the rate of change of the abundance of the active species
jac   = numpy.zeros((nSpecsActive,nSpecsActive), 'f8')  #the jacobian

#for each reaction in the network, see what are the reactants and the products and 
#accroding to increment/decrement the rate of change in dxdt 
for i, rxn in enumerate(net.reactions):
    
    #rate of the reaction in the temperature range
    rate = rxn.rxn_in_trng.rate
    
    for specStr in rxn.reactants + rxn.products:
        
        spec = net.species[specStr]
        
        if spec.type >= 0:
            if specStr in rxn.reactants:
                dxdt[spec.num] -= rate
            else:
                dxdt[spec.num] += rate
    

#checking for the conservation of the number of 'baseSpecies' (which are active)
for specStr, spec in net.species.items():
    
    for comp in spec.comp:
        
        compStr   = comp[0]
        compCount = comp[1]
        
        if compStr == 'PAH':
            net.species[specStr].show()
        
        #index corresponding to the 'abundance' of the base specie in net.base_abun array
        idx = net.baseSpecies[compStr].num 
        
        net.base_abun[idx] += compCount * spec.abun()

for specStr, spec in net.baseSpecies.items():
    print specStr, spec.abun()

print 'done'