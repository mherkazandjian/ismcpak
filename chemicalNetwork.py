from string import *
from numpy import *
import StringIO
import re
from specie import *
from reaction import *
from time import *
from collections import *

"""
     self.species   # tuple containing the unique species objects in the network 
     self.__init__(fileName=None, baseSpecies=None)
     self.readDatabase(fileName)
     self.printFile()
     self.printReactions()
     self.setRxnAttributes(lineStr)
     self.parseReactions()
     self.analyzeNetwork(baseSpecies)
     self.updateSpeciesInReactions()
     self.setup(databaseFname, baseSpecies):
     # returns a counter object containing the enum of the species
     enumList = self.speciesEnum()  
    
     TODO:
         write a method which takes as input the species file used and filters
         out all the reaction which do not contain these species and re-numbers
         the species in the number representation of the chemical network
"""

### write a method which reads the abundances  and sets the abundance of each specie object
### this can be implemented by appending a new abudance column to species.inp 

# class definition for the chemical network
# ----------------------------------------
class chemicalNetwork(specie, reaction):

    def __init__(self, fileName=None, baseSpecies=None):
        self.nReactions=0              # number of reactions in the network
        self.reactions=[]              # reaction object list in the netowrk
        self.species=[]                # specie object list in the network
        self.fileName = fileName       # path of the ascii reaction file of the network
        self.baseSpecies = baseSpecies # object list of the base species
        self.nDens = None              # hydrogen gas ambient density
        self.T = None                  # gas temperature
        self.zeta = None               # cosmic ray ionization rate
        self.albedo = None
        self.Av = None
        self.G0 = None
        self.specDict = {}  # replace the self.species with this
        if (self.fileName != None) and (self.baseSpecies != None):
            self.setup(fileName, baseSpecies)

    # set the temp, zeta, albed, Av
    def setGasTemperature(self, T): self.T = T
    def setCrRate(self, zeta)     : self.zeta = zeta
    def setAlbed(self, albedo)    : self.albedo = albedo
    def setAv(self, Av)           : self.Av = Av
    def setG0(self, G0)           : self.G0 = G0
    def setDens(self, nDens)      : self.nDens = nDens
    
    def setCloudParms(self, T, zeta, Av, albedo, nDens, G0):
        self.setGasTemperature(T)
        self.setCrRate(zeta)
        self.setAlbed(albedo)
        self.setAv(Av)
        self.setG0(G0)
        self.setDens(nDens)
        
    # read and parses the UMIST 2006 reaction file without the header and assigns the 
    # variable fileStr and counts the number of reactions in the file 
    # sets : self.fileStr
    #        self.nRxn
    def readNetworkFile(self, fileName):
        file = open(fileName, 'r')
        print 'Opened chemical reaction network file : ' + fileName

        nRxn = 0
        fileStr = ''           # string that will contain the whole file
        for line in file:
            if line[0]=='#':
                continue
            fileStr += line
            nRxn = nRxn + 1
        
        self.fileStr = fileStr
        self.nRxn    = nRxn
        print '     Read ', self.nRxn, 'reactions'

    # prints the reaction file without the header
    def printFile(self):
        print self.fileStr

    # prints the reaction file without the header
    # rxnIds : an integer list of the reaction indecies to be printed
    # format : a string for the format of the reactions (see reaction object doc)
    def printReactions(self, rxnIds = None, format = None ):

        if rxnIds == None: # print all the reactions
            print 'total number of reactions = ', len(self.reactions)
            for rxn in self.reactions:
                rxn.display( format = format )
        else: # print reactions with specified ids
            for rxnId in rxnIds:
                for rxn in self.reactions:
                    if rxnId == rxn.id:
                        rxn.display(  format = format )
                        break

    # print all the species
    def printSpecies(self):
        for spec in self.species:
            print spec.str + ' : ',
            for com in spec.com:
                print com,
            print " charge : %+d " % spec.charge,
            print " abun = %e" % spec.abun
            print "",

    # returns the index of the specie which matches the input string
    def getSpecieIndex(self, specStr):
        i = 0
        for spec in self.species:
            if spec.str == specStr:
                return i
            i = i + 1
        print 'Specie ' + specStr + ' not found'
        return -1

    # set the abundances of the species from an ascii input file 
    # column 1 : index, column 2 : specie string, column 3 : abundance relative to H
    def setAbundancesFromFile(self, speciesFileName, abunFileName):
        
        print 'Setting the abundances of the species:'
        # reading the input species whose abundance is to be assigned
        specFile = open(speciesFileName, 'r')
        specsStrRead = []
        for line in specFile:
            lineParsed = line.split(' ')
            specStr = strip(lineParsed[1])
            specsStrRead.append(specStr)

        # reading the abundances of the species whose abundance is to be assigned
        abunFile = open(abunFileName, 'r')
        abunRead = []
        for line in abunFile:
            abun = float64(line)
            abunRead.append( abun )

        # assigining the abundances of the species to the specie objects in the network
        i = 0
        for i in arange(len(specsStrRead)):
            indx = self.getSpecieIndex(specsStrRead[i])
            if indx != -1:
                #print indx, specsStrRead[i], self.species[indx].str
                self.species[indx].abun = abunRead[i]
            i = i + 1


    # pasrse the UMIST reaction line and return the reaction object
    def setRxnAttributes(self, lineStr):
        lineParsed = lineStr.split(',')
        rxn = reaction()
        rxn.setAllFromRxnStrArr( lineParsed )
        return rxn
    
    # define rxn object list for all the reaction from the ascii database lines
    def parseReactions(self):

        print 'Parsing reactions...',

        for line in (self.fileStr).splitlines():
            rxn = self.setRxnAttributes(line)    # create a new  rxn object from reaction line
            self.reactions.append(rxn)           # append the new rxn to the reactions
            
        print 'complete'

    # count all the species and define them
    def analyzeNetwork(self, baseSpecies):

        print 'Analyzing the checmical network'
        # append the basic species
        uniqueSpecStr = []
        for baseSpec in baseSpecies:
            uniqueSpecStr.append(baseSpec.str)
            self.species.append(baseSpec)
            self.specDict[baseSpec.str] = baseSpec
            
        # getting the unique species from all the reactions 
        for rxn in self.reactions:

            # collecting all the species involved in a certain reaction
            specsStrAllThisRxn = [ rxn.R1.str, rxn.R2.str, rxn.R3.str,
                                   rxn.P1.str, rxn.P2.str, rxn.P3.str, rxn.P4.str ]

            # checking if the species in the reaction are already in the unique set or not
            for specStr in specsStrAllThisRxn:
                
                if len(specStr) != 0 and (specStr not in self.specDict):
                    newSpec = specie(specStr, specType=1)   # making the specie obj
                    newSpec.getComponents(baseSpecies)      # parsing the specie
                    self.specDict [specStr] = newSpec
                    
                if len(specStr) != 0:
                    if uniqueSpecStr.count(specStr) == 0:
                        uniqueSpecStr.append(specStr) # adding the unique spec to the tmp str buff
                        newSpec = specie(specStr, specType=1)   # making the specie obj
                        newSpec.getComponents(baseSpecies)  # parsing the specie
                        self.species.append(newSpec)  # appending the new spec to the spec list

        for i in arange(len(self.species)):
            self.species[i].num = i

        print '     Found ', len(self.species), ' unique species in the network'


    # construct the numeric representation of the species
    def updateSpeciesInReactions(self):

        print 'Updating the species in the reaction and setting the numbers....',
        
        # collect all the species strings into one string list
        allSpecStrList = []
        for spec in self.species:
            allSpecStrList.append(spec.str)

        # returns the specie object corresponding to the specie string
        def replaceSpecie( spec ):
            if len(spec.str) != 0:
                specIdx = allSpecStrList.index(spec.str)
                return self.species[specIdx]
            else:
                return spec

        # so far, the specie objects in the reactions have only the specie string
        # here we replance the species with their full objects with the parsed species
        for i in arange(len(self.reactions)):
            rxn = self.reactions[i]

            rxn.R1 = replaceSpecie(rxn.R1)
            if rxn.R1.str != '':
                rxn.reactants.append(rxn.R1)

            rxn.R2 = replaceSpecie(rxn.R2)
            if rxn.R2.str != '':
                rxn.reactants.append(rxn.R2)

            rxn.R3 = replaceSpecie(rxn.R3)
            if rxn.R3.str != '':
                rxn.reactants.append(rxn.R3)
            
            rxn.P1 = replaceSpecie(rxn.P1)
            if rxn.P1.str != '':
                rxn.products.append(rxn.P1)

            rxn.P2 = replaceSpecie(rxn.P2)
            if rxn.P2.str != '':
                rxn.products.append(rxn.P2)

            rxn.P3 = replaceSpecie(rxn.P3)
            if rxn.P3.str != '':
                rxn.products.append(rxn.P3)

            rxn.P4 = replaceSpecie(rxn.P4)
            if rxn.P4.str != '':
                rxn.products.append(rxn.P4)
            
            self.reactions[i] = rxn
        print 'compelte'
    # each reaction has 7 species involved max. In the UMIST2006 there are ~450 species. Each
    # specie can be represented by 9 bits ( max 512 unsigned integer). So each reaction can
    # have a tag composed of 9*7=63 bits, which is unique as long as the species involved in
    # the reactions are different.
    def setReactionHashcodes(self):
        
        base=uint64(512)
        allHashCodes = []

        for rxn in self.reactions:

            if rxn.R1.str != '':
                rxn.hash += (base**0)*uint64(rxn.R1.num)
            if rxn.R2.str != '':
                rxn.hash += (base**1)*uint64(rxn.R2.num)
            if rxn.R3.str != '':
                rxn.hash += (base**2)*uint64(rxn.R3.num)

            if rxn.P1.str != '':
                rxn.hash += (base**3)*uint64(rxn.P1.num)
            if rxn.P2.str != '':
                rxn.hash += (base**4)*uint64(rxn.P2.num)
            if rxn.P3.str != '':
                rxn.hash += (base**5)*uint64(rxn.P3.num)
            if rxn.P4.str != '':
                rxn.hash += (base**6)*uint64(rxn.P4.num)
    
            allHashCodes.append(rxn.hash)

#        # looking for identical reactions
#        for j in arange(len(self.reactions)):
#            rxn = self.reactions[j]
#            found=0
#            if allHashCodes[j] == -1: # if already been matched skip it
#                continue
#            for i in arange(len(allHashCodes[j+1:-1])):
#                if (allHashCodes[i] == rxn.hash) and (rxn.id != self.reactions[i].id):
#                    self.reactions[i].display(full=1)
#                    allHashCodes[i]=-1  # replacing the hashcode by -1 to indicate it has been matched
#                    found=1
#            if found == 1:
#                rxn.display(full=1)
#                print '--------------------------------------------------'


    # setup all the chemical network
    def setup(self, networkFname, baseSpecies):
        self.readNetworkFile(networkFname) # read the database into a buffer
        self.parseReactions()              # parse the reaction lines from the database
        self.analyzeNetwork(baseSpecies)   # fileter the unique species and parse them in their components 
        self.updateSpeciesInReactions()    # replacing the spcie objects in the reactions with the parsed species
        self.setReactionHashcodes()        # compute the hashcodes of the reactions

    #  write this method
    #        self.getDuplicateReaction()       # matched the hashcodes to determine which reactions are duplicated

        
    # method that returns the indecies of the reactions containing the input specie
    def filterReactions(self, reactSpecStr, productSpecStr):

        rxnIndsFound = []
        
        for rxn in self.reactions:

            inReacts = 0
            inProds  = 0

            # checking if the specie we are looking for is in the reactants
            if len(reactSpecStr) != 0 :
                for spec in rxn.reactants:
                    if reactSpecStr == spec.str:
                        inReacts = 1
                        break
            else:
                inReacts = 1

            # checking if the specie we are looking for is in the products
            if len(productSpecStr) != 0:
                for spec in rxn.products:
                    if productSpecStr == spec.str:
                        inProds = 1
                        break
            else:
                inProds = 1

            if inReacts and inProds:
                rxnIndsFound.append( rxn.id )

        return rxnIndsFound


    # compute the reaction constants
    def computeReactionConstants(self):
        
        for rxn in self.reactions:
            
            # compute the constant for reactions involving a photon
            if rxn.reactants[1].str == 'PHOTON':
                rxn.cst = rxn.alpha * exp( - rxn.gamma * self.Av )
            else:
                # compute the constant for reactions involving CRP
                if rxn.reactants[1].str == 'CRP':
                    rxn.cst = rxn.alpha
                else:
                    # compute the constant for reactions involving CRPHOTON
                    if rxn.reactants[1].str == 'CRPHOTON':
                        rxn.cst = rxn.alpha * ((self.T / 300.0)**rxn.beta) * rxn.gamma / ( 1.0 - self.albedo )
                    else:
                        # compute the constant for two body reactions
                        rxn.cst = rxn.alpha * ((self.T / 300.0)**rxn.beta) * exp( - rxn.gamma / self.T )


    # compute the reaction constants and rates
    def computeRates(self):

        for rxn in self.reactions:
            
            rxn.rate = rxn.cst
            
            # compute the rate for reactions involving a photon
            if rxn.reactants[1].str == 'PHOTON':
                if rxn.reactants[0].abun != None:
                    rxn.rate = rxn.rate * (self.nDens * rxn.reactants[0].abun)
                else:
                    rxn.rate = None
            else:
                # compute the rate for reactions involving CRP
                if rxn.reactants[1].str == 'CRP':
                    if rxn.reactants[0].abun != None:
                        rxn.rate = rxn.rate * (self.nDens * rxn.reactants[0].abun)
                    else:
                        rxn.rate = None
                else:
                    # compute the constant for reactions involving CRPHOTON
                    if rxn.reactants[1].str == 'CRPHOTON':
                        if reactants[0].abun != None:
                            rxn.rate = rxn.rate * (self.nDens * rxn.reactants[0].abun)
                        else:
                            rxn.rate = None
                    else:
                        # compute the constant for two body reactions
                        for spec in rxn.reactants:
                            if spec.abun != None:
                                rxn.rate = rxn.rate*(self.nDens * spec.abun)
                            else:
                                rxn.rate = None
                                break

    # sort the reactions with decreasing absolute rates
    def sortRxnsDecreasingRates(self, rxnIds):
        
        rates = []
        ids   = []

        for rxnId in rxnIds:
            for rxn in self.reactions:
                if rxnId == rxn.id and rxn.rate != None:
                    rates.append( rxn.rate )
                    ids.append( rxn.id )

        def sortIndices( data ):
            return sorted( range(len(data)), key = data.__getitem__, reverse=True)

        inds = sortIndices( rates )

        idsSorted = []
        for i in inds:
            idsSorted.append( ids[i] )

        return idsSorted

    def speciesEnum(self): #kkkkkkkkkkkk
        cnt = Counter()
        for spec in self.species:
            print spec.str
#            cnt = self.speciesEnum()  
