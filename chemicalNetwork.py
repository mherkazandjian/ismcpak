from string import *
import numpy as np

from specie import *
from reaction import *

"""
     self.species   # tuple containing the unique species objects in the network 
     self.__init__(fileName=None, baseSpecies=None)
     self.readDatabase(fileName)
     self.printFile()
     self.printReactions()
     self.setRxnAttributes(lineStr)
     self.parseReactions()
     self.getUniqueSpecies(baseSpecies)
     self.updateSpeciesInReactions()
     self.setup(databaseFname, baseSpecies):
     self.changeSpeciesNums( fileName )
     
     TODO:
         write a method which takes as input the species file used and filters
         out all the reaction which do not contain these species and re-numbers
         the species in the number representation of the chemical network
"""

### write a method which reads the abundances  and sets the abundance of each specie object
### this can be implemented by appending a new abudance column to species.inp 
### ;;;; as a check, dump the remainders of the parsed strings to a file...it should
#        mainly be black other than commas...etc...otherwise, there are stuff not taken
#        into account, this can be done by erasing stuff from self.fileStr after the whole
#        original network has been worked out

# class definition for the chemical network
# ----------------------------------------
class chemicalNetwork(specie, reaction):

    def __init__(self, fileName=None, baseSpecies=None):
        self.nReactions=None             # number of reactions in the network
        self.reactions=[]              # reaction object list in the netowrk
        self.nSpecs = None
        self.species={}                # specie object list in the network
        self.fileName = fileName       # path of the ascii reaction file of the network
        self.baseSpecies = baseSpecies # object list of the base species
        self.nDens = None              # hydrogen gas ambient density
        self.T = None                  # gas temperature
        self.zeta = None               # cosmic ray ionization rate
        self.albedo = None
        self.Av = None
        self.G0 = None
        self.abun = None
        
        if (self.fileName != None) and (self.baseSpecies != None):
            self.setup(fileName, baseSpecies)
        
    # setup all the chemical network
    def setup(self, networkFname, baseSpecies):
        self.readNetworkFile(networkFname) # read the database into a buffer
        self.parseReactions()              # parse the reaction lines from the database
        self.getUniqueSpecies(baseSpecies) # fileter the unique species and parse them in their components
        self.assignNumbersToSpecies()      # assingns the index of the specie in the dictionary
        self.setReactionHashcodes()        # compute the hashcodes of the reactions
        
    # read and parses the UMIST 2006 reaction file without the header and assigns the 
    # variable fileStr and counts the number of reactions in the file 
    # sets : self.fileStr
    #        self.nRxn
    def readNetworkFile(self, fileName):
        fObj = open(fileName, 'r')
        print 'Opened chemical reaction network file : ' + fileName

        nRxn = 0
        fileStr = ''           # string that will contain the whole file
        for line in fObj:
            if line[0]=='#':
                continue
            fileStr += line
            nRxn = nRxn + 1
        
        self.fileStr = fileStr
        self.nRxn    = nRxn
        
        fObj.close()
        print '     Read ', self.nRxn, 'reactions'

    # define rxn object list for all the reaction from the ascii database lines
    # here the reactant and product objects are not set
    def parseReactions(self):

        print 'Parsing reactions...',

        for line in (self.fileStr).splitlines():
            rxn = self.setRxnAttributes(line)    # create a new  rxn object from reaction line
            self.reactions.append(rxn)           # append the new rxn to the reactions
            
        print 'complete'

    # construct the the speices dict for all the species in the network
    #    self.species
    # also here all the unique species in the network are put in that dictionary
    # and the components of the species are also extracted
    def getUniqueSpecies(self, baseSpecies):

        print 'Analyzing the checmical network.....',
        # put the base species in the dictionary
        for baseSpec in baseSpecies:
            baseSpec.getComponents(baseSpecies)
            self.species[baseSpec.str] = baseSpec
                
        # getting the unique species from all the reactions 
        for rxn in self.reactions:
            for spec in rxn.species.values():
                if spec.str not in self.species: 
                    spec.type = 1   # setting the object as a non-base specie
                    spec.getComponents(baseSpecies)      # parsing the specie
                    self.species [spec.str] = spec
                                        
        self.nSpecs = len(self.species)
        print 'Found ', self.nSpecs, ' unique species in the network'
        
        # assign the the species object in the reactions to pointer to those
        # in self.species. This saves memory and time, since only
        # modifying self.species would modify all the corresponding
        # species all the reactions whereever they occure 
        self.updateSpeciesInReactions()
        
        # defining the array where the abundances of the species are stored
        # and each entry is mapped to self.species[].abun, such that
        # self.species['XXX'].abun == self.abun[ self.species['XXX'].num ] 
        self.abun = np.zeros( len(self.species) )


    # setting the number_tag of the species objects in the reaction objects
    # such that it would be possible to construct the numeric representation of
    # the reactions, i.e instead of having : 
    #           HCN  +  H   ->  HNC  +   H
    # we get    '56' +  '9' ->   27  +  '9'
    # assuming HCN, HNC and H have indicies 56, 27 and 9 respectively
    # WARNING : This should be called everytime an entry in self.species
    #           is replaced (i.e it addrress changes)
    def updateSpeciesInReactions(self):

        print 'constructing numeric represntation of the network....',
        
        # setting self.reactions[i].species[:].num 
        for i in np.arange( len(self.reactions) ):
            for specStr in (self.reactions[i].reactants + self.reactions[i].products):
                self.reactions[i].species[specStr] = self.species[specStr]

        print 'compelte'

    # set all the numbers of the species self.species[].num to None
    def setSpeciesNumsToNone(self):
        
        for specStr in self.species:      
            self.species[specStr].num = None
        
    # assing the self.species[:].num to their index in the dictionary
    # and update these indecis in the chemical network as well
    # WARNING : it might happen that by mistake some species number is 
    #           not specified in the input file containing the new number,
    #           so it will be None. Or that the file contains a species 
    #           that is not in the chemical network, these will also have
    #           None as numbers
    def assignNumbersToSpecies(self, fileName = None):

        self.setSpeciesNumsToNone()
                
        if fileName == None:
            print 'setting species number in the order they appear in self.species'
            # assign the number 
            i = 0
            for specStr in self.species:
                self.species[specStr].num = i
                i += 1
        else:
            # it might be a good idea to put this in the mesh database thing or mesh object
            # as a new dictionary and pass the dictionary object to this instead of a file
            print 'setting species number in the order they appear in the file: ' + fileName
            
            fObj =  open(fileName, 'r')
            for line in fObj:
                line = line.split()
                
                specStr = line[1]
                numStr  = line[0]
                if specStr in self.species:
                    self.species[ specStr ].num = int32(numStr)
                else:
                    print '%-12s from the file is not a species in the current network' % specStr 
#                print line[0], line[1]
        
            fObj.close() 
        
        for spec in self.species.values():
            if spec.num == None:
                print 'index of %-12s is missing in the file, it is kept as None' % spec.str
             
        self.mapSpeciesAbundances()
        
        print 'assigned species number and mapped the abundances'
        
    # map self.species[].abun to self.abun
    def mapSpeciesAbundances(self):
        
        for specStr in self.species:
            self.species[specStr].abun = self.abun[ self.species[specStr].num ] 

    def setAbundances(self, fromFile = None, fromMesh = None):
        
        print 'Setting the abundances of the species:'
        # reading the input species whose abundance is to be assigned
        specFile = open(speciesFileName, 'r')
        specsStrRead = []
        for line in specFile:
            lineParsed = line.split(' ')
            specStr = strip(lineParsed[1])
            specsStrRead.append(specStr)
            
        print specsStrRead
        aasdasd
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
                self.species[indx].abun = abunRead[i] ########
            i = i + 1

"""
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
    CONTINUE IMPLEMENTING SETTING THE ABUNDANCES FROM A MESH AND FROM AN ASCII FILE!!!
"""
    # set the abundances of the species from an ascii input file
    # containg the species name and the corresponding abundance 
    # column 1 : index (not used), column 2 : specie string, column 3 : abundance relative to H
    # ;;;;;;; continue working on this.......
    # pass as a paramter a dictionary specifying the specie and its abundance ex.. 
    #     spec1 : 1e-5
    #     spec2 : 0.034
    #           .
    #           .
    #          etc
    def setAbundancesFromMesh(self, mesh):
        
        print 'Setting the abundances of the species:'
        # reading the input species whose abundance is to be assigned
        specFile = open(speciesFileName, 'r')
        specsStrRead = []
        for line in specFile:
            lineParsed = line.split(' ')
            specStr = strip(lineParsed[1])
            specsStrRead.append(specStr)
            
        print specsStrRead
        aasdasd
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
                self.species[indx].abun = abunRead[i] ########
            i = i + 1
 
    # each reaction has 7 species max involved (UMIST). In the UMIST2006 there are ~450
    # species. Each specie can be represented by a 9 bit integer ( we use an unsigned integer)
    # which can have a max value of 512. So each reaction can have a tag composed of 9*7=63 
    # bits, which is unique as one of the involved in the reactions is. For example, the 
    # reaction :
    #       (0)   (1)    (2)        (3)   (4)    (5)    (6)
    #       221 + 286          ---> 249 + 369
    # wpuld have a hash code :
    #       (512^0)*221 + (512^1)*286 + (512^3)*249 + (512^4)*369   
    def setReactionHashcodes(self):

        print 'setting reaction hashcodes....',        
        base=uint64(512)
        allHashCodes = []

        for rxn in self.reactions:
            
            i = 0
            for specStr in rxn.reactants:
                rxn.hash += (base**i)*uint64( rxn.species[specStr].num ) 
                i += 1

            i = 0
            for specStr in rxn.products:
                rxn.hash += ( base**(3 + i) )*uint64( rxn.species[specStr].num )
                i += 1
    
            allHashCodes.append(rxn.hash)

        print 'complete'
 
#    def findIdenticalReactions(self):
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

    #  write this method
    #        self.getDuplicateReaction()       # matched the hashcodes to determine which reactions are duplicated


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
        for spec in self.species.values():
            spec.show()
            
    # returns the index of the specie which matches the input string
    def getSpecieIndex(self, specStr):
        if specStr in self.species: #i = 0
            return self.species[specStr].num
        else:
            print 'Specie ' + specStr + ' not found'
            return None


    # pasrse the UMIST reaction line and return the reaction object
    def setRxnAttributes(self, lineStr):
        lineParsed = lineStr.split(',')
        rxn = reaction()
        rxn.setAllFromRxnStrArr( lineParsed )
        return rxn
        
    # method that returns the indecies of the reactions containing the input specie
    # withReacts : a list of strings holding the reactants to be used as a filter
    # withProds  : a list of strings holding the reactants to be used as a filter
    # ex :
    # inds = filterReactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'])
    def filterReactions(self, withReacts = None, withProds = None):

        rxnIndsFound = []
        
        for rxn in self.reactions:

            inReacts = True
            inProds  = True

            # checking if the specie we are looking for is in the reactants
            if withReacts != None:
                for specStr in withReacts:
                    if specStr in rxn.reactants:
                        inReacts = inReacts and True
                    else:
                        inReacts = inReacts and False

            # checking if the specie we are looking for is in the reactants
            if withProds != None:
                for specStr in withProds:
                    if specStr in rxn.products:
                        inProds = inProds and True
                    else:
                        inProds = inProds and False
                

            if inReacts and inProds:
                rxnIndsFound.append( rxn.id )

        return rxnIndsFound


    # compute the reaction constants
    def computeReactionConstants(self):

        # it might be a good idea to put a check if
        # temperature,A_v are set...         
        for rxn in self.reactions:
            
            # compute the constant for reactions involving a photon
            if 'PHOTON' in rxn.reactants:
                rxn.cst = rxn.alpha * exp( - rxn.gamma * self.Av )
            else:
                # compute the constant for reactions involving CRP
                if 'CRP' in rxn.reactants:
                    rxn.cst = rxn.alpha
                else:
                    # compute the constant for reactions involving CRPHOTON
                    if 'CRPHOTON' in rxn.reactants:
                        rxn.cst = rxn.alpha * ((self.T / 300.0)**rxn.beta) * rxn.gamma / ( 1.0 - self.albedo )
                    else:
                        # compute the constant for two body reactions
                        rxn.cst = rxn.alpha * ((self.T / 300.0)**rxn.beta) * exp( - rxn.gamma / self.T )


    # compute the reaction constants and rates assuming all the abundances are set
    # adapted to new variables but need to check it!!!!
    def computeRates(self):

        for rxn in self.reactions:
            
            rxn.rate = rxn.cst
            # use a dictionary for the types of reactions instead 
            # of if/else statments
            
            reactants = rxn.reactants
            species   = rxn.species
            
            # compute the rate for reactions involving a photon
            if 'PHOTON' in rxn.reactants:
                rxn.rate *= (self.nDens * species[reactants[0]].abun)
            else:
                # compute the rate for reactions involving CRP
                if 'CRP' in rxn.reactants:
                    rxn.rate *= self.nDens * species[reactants[0]].abun
                else:
                    # compute the constant for reactions involving CRPHOTON
                    if 'CRPHOTON' in rxn.reactants:
                            rxn.rate *= self.nDens * species[reactants[0]].abun
                    else:
                        # compute the constant for two body reactions
                        for specStr in rxn.reactants:
                            rxn.rate *= self.nDens * species[specStr].abun

    # sort the reactions with decreasing absolute rates
    # adapted to new variables but need to check it!!!!
    # adapted to new variables but need to check it!!!!
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

    # change the number(index) 
    #     self.species[:].num
    #  of the species according to those in the file fileName which contains
    #        specie_idx  speciesStr 
    #def changeSpeciesNums( fileName ):