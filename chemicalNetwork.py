from string import *
import numpy as np

from specie import *
from reaction import *

"""
   *if a specie 'XXX' is not in :
       self.speciesRemoved
    and self.species['XXX'].num = None
    then the reactions which have 'XXX' as a specie are still in the 
    network but they are not used in computing the reaction, like for
    example PHOTON, it does not have an abundance and we do not keep
    track of the number of photons, so it need not have a number with
    which its abudnance can be tracked in self.abun[]   
   
   * reaction hash codes are set based on the original numbers assigned to
     species and it doesnt depend on the order the reactnats or the 
     products appear in the reaction input file. When computing the
     hash code, the species are picked in the increasing order as 
     they appeat in the sorted string array (using the sorted function)
     
     self.species   # tuple containing the unique species objects in the network 
     self.__init__(fileName=None, baseSpecies=None)
     self.readDatabase(fileName)
     self.printFile()
     self.printReactions()
     self.setRxnAttributes(lineStr)
     self.parseReactions()
     self.getUniqueSpecies(baseSpecies)
     self.mapReactionSpeciesToSpecieObjectList()
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

    def __init__(self, fileName=None, baseSpecies=None, UMISTVER = None):
        self.nReactions=None           # number of reactions in the network
        self.reactions=[]              # reaction object list in the netowrk

        self.nReactionsRemoved=None  # number of reactions in the network
        self.reactionsRemoved=[]     # reaction object list in the netowrk
        
        self.nSpecs = None             # number of active species in the network
        self.species={}                # active specie object list in the network
        
        self.nSpecsRemoved = None    # number of not active species in the network
        self.speciesRemoved={}       # not active specie object list in the network
        
        self.fileName = fileName       # path of the ascii reaction file of the network
        self.baseSpecies = baseSpecies # object list of the base species
        self.nDens = None              # hydrogen gas ambient density
        self.T = None                  # gas temperature
        self.zeta = None               # cosmic ray ionization rate
        self.albedo = None
        self.Av = None
        self.G0 = None
        self.abun = None               # (nSpecies, 1) np array, it is by default initialized to -1
        
        self.pstr = '' # just for printintg purposes
        
        if UMISTVER != None:
            self.umistVer = UMISTVER
        else:
            self.umistVer = None
            
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
            nRxn += 1
        
        self.fileStr = fileStr
        self.nReactions = nRxn
        
        fObj.close()
        print '     Read ', self.nReactions, 'reactions'

    # define rxn object list for all the reaction from the ascii database lines
    # here the reactant and product objects are not set
    def parseReactions(self):

        print 'Parsing reactions...',

        if self.umistVer == None:
            str  = 'Error : cannot parse network reaction, unknown database\n'
            str += '                   please set database version.'                 
            raise NameError(str)
        
        print 'using '+ self.umistVer + ' format ',
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
        self.mapReactionSpeciesToSpecieObjectList()
        
        # defining the array where the abundances of the species are stored
        # and each entry is mapped to self.species[].abun, such that
        # self.species['XXX'].abun == self.abun[ self.species['XXX'].num ]
        # the abundances are intiialized to -1
        self.initializeDefaultAbunBuffer() 
        
        self.nReactionsRemoved = len(self.reactionsRemoved)
        self.nSpecsRemoved = len(self.speciesRemoved)
        

    # mapping the species objects in all the reactions to the species
    # object list in the network object (this object)...this is useful
    # since instead of having a copy of a species which might occure in
    # many reaction, it is just referenced to the species in the species list
    # it is also useful since changin an attirbute of a species in the list
    # would be visible to all the occurances in the network. such as...
    # setting the number_tag of the species objects in the reaction objects
    # such that it would be possible to construct the numeric representation of
    # the reactions, i.e instead of having : 
    #           HCN  +  H   ->  HNC  +   H
    # we get    '56' +  '9' ->   27  +  '9'
    # assuming HCN, HNC and H have indicies 56, 27 and 9 respectively
    # and this number would be visible to all the instances where these
    # species occure
    # WARNING : This should be called everytime an entry in self.species
    #           is replaced (i.e it addrress changes) for ex by call
    def mapReactionSpeciesToSpecieObjectList(self):

        print 'constructing numeric represntation of the network....',
        
        # setting self.reactions[i].species[:].num 
        for i in np.arange( len(self.reactions) ):
            for specStr in (self.reactions[i].reactants + self.reactions[i].products):
                self.reactions[i].species[specStr] = self.species[specStr]

        print 'compelte'

    # remove species from the network. The species strings which appear
    # in the list are moved to speciesRemoved list and the reactions with 
    # the removed species are moved to reactionRemoved.
    def removeSpecies(self, specStrList):
        
        print 'removing species from the network....'

        for specStr in specStrList:
            # moving the specie to be removed to self.speciesRemoved
            # and update the numbers of spceis and removed one
            self.speciesRemoved[ specStr ] = self.species.pop(specStr)
            self.nSpecs = len(self.species)
            self.nSpecsRemoved = len(self.speciesRemoved) 
      
            # moving the reactions containing the species to be removed to the
            # self.reactionsRemoved list and updating the numbers. Starting from
            # the last reaction and working backwards in the reaction list
            n = self.nReactions
            for i in np.arange(n):
                if specStr in self.reactions[n-i-1].species:
                    self.reactionsRemoved.append( self.reactions.pop(n-i-1) )
                    
            self.nReactions = len(self.reactions)
            self.nReactionsRemoved = len(self.reactionsRemoved)

        self.pstr = '     '
        self.initializeDefaultAbunBuffer()
        self.assignNumbersToSpecies()
        
        print 'complete'

            
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
            print self.pstr + 'setting species number in the order they appear in self.species'
            self.pstr = '     '
            # assign the number 
            i = 0
            for specStr in self.species:
                self.species[specStr].num = i
                i += 1
        else:
            # it might be a good idea to put this in the mesh database thing or mesh object
            # as a new dictionary and pass the dictionary object to this instead of a file
            print self.pstr + 'setting species number in the order they appear in the file: ' + fileName
            self.pstr = '     '
            
            fObj =  open(fileName, 'r')
            i = 0
            for line in fObj:
                line = line.split()
                
                specStr = line[1]
                if specStr in self.species:
                    self.species[ specStr ].num = int32(i)
                else:
                    print self.pstr + '%-12s from the file is not a species in the current network' % specStr 
                i += 1
            fObj.close() 
        
        print  self.pstr + 'finished assigning the new species numbers'
        
        # printintg the species which did not get a number assigned to
        i = 0
        for spec in self.species.values():
            if spec.num == None:
                print self.pstr + 'index of %-12s is missing in the file, it is kept as ' % spec.str, spec.num
                i += 1
                
        if i != 0:
            print self.pstr + 'A total of %d species will have no numbers' % i 
        
        self.mapSpeciesAbundances()
        self.pstr = ''        
        
    # map self.species[].abun to self.abun
    # WARNING : species which do not have a number, get their abundance set to None
    def mapSpeciesAbundances(self):
        
        for specStr in self.species:
            specNum = self.species[specStr].num
            if specNum != None:
                self.species[specStr].abun = self.abun[ specNum ]
            else:
                self.species[specStr].abun = None
        print self.pstr + 'mapped spcies abundances to the abundance array'
            
    # each reaction has 7 species max involved (UMIST). In the UMIST2006 there are ~450
    # species. Each specie can be represented by a 9 bit integer ( we use an unsigned integer)
    # which can have a max value of 512. So each reaction can have a tag composed of 9*7=63 
    # bits, which is unique as one of the involved in the reactions is. For example, the 
    # reaction :
    #       (0)   (1)    (2)        (3)   (4)    (5)    (6)
    #       286 + 221          ---> 249 + 369  + 56
    # wpuld have a hash code :
    #       (512^0)*221 + (512^1)*286 + (512^3)*56 + (512^4)*249 + (512^5)*369   
    def setReactionHashcodes(self):

        print 'setting reaction hashcodes....',        
        base=uint64(512)

        for rxn in self.reactions:
            
            i = 0
            for specStr in sorted(rxn.reactants):
                rxn.hash += np.uint64((base**i))*np.uint64( rxn.species[specStr].num ) 
                i += 1

            i = 0
            for specStr in sorted(rxn.products):
                rxn.hash += np.uint64( base**(3 + i) )*np.uint64( rxn.species[specStr].num )
                i += 1
            
        print 'complete'

    # finds identical reactions which have the same reactants and the same
    # products and returns a tuple holding the reaction indecis
    # ( (ind1,ind2), (ind1,ind2)...)
    def findIdenticalReactions(self):
        
        sets = ()
        n = self.nReactions
        hashes = np.zeros( n, dtype = np.uint64 )
        inds   = np.zeros( n, dtype = np.uint32 )
        
        # collect hashes numpy arrays        
        i = 0
        for rxn in self.reactions:
            hashes[i] = rxn.hash
            inds[i]   = i
            i += 1
        
        # indices of sorted hashes in increasing order
        indsSorted = np.argsort(hashes, kind='quicksort')
        hashes     = hashes[indsSorted]
        inds       = inds[indsSorted]  
        
        # looping over the sorted hashes array and storing into
        # tuple the indicies of the corresponding reactions with the same
        # hash
        i = 0
        while i < n - 1:
            hi = hashes[i]
            
            if hashes[i+1] == hi:
                indsSameHash = (inds[i],)
                # checking for the hashes of the reaction below it to see
                # if there is more than one wit hthe same hash
                j=1
                while hashes[i+j] == hi:
                    indsSameHash += (inds[i+j],)
                    j += 1 
                
                i += j
                
                sets += (indsSameHash,)
            else:
                i += 1
                
        return sets
    
    #  write this method
    #        self.getDuplicateReaction()       # matched the hashcodes to determine which reactions are duplicated

    # set the abundances of the species from an ascii input file
    # or from a numpy array. 
    # the file contains one abundance on each line. each value will be set to the
    # corresponding vlaue in self.abun at the index corresponding to the line number
    # (staring with index 0). There should be as much lines as there are species with
    # assigned numbers. Also the numpy array should have the same length as the species
    # with numbers. 
    # WARNING : since self.specie[].abun are mapped to self.abun, in chanfing
    #           all the values of the array the method self.copyAbundancesFromArray( array )
    #           should be used, which sets the values one by one, instead of
    #           replacing self.abun with another array when self.abun = newAbun
    #           which removes the mapping.
    def setAbundances(self, fromFile = None, fromArray = None):
        
        print 'Setting the abundances of the species from:',
        
        newAbunArr = None
        
        if fromFile != None:
            print 'file : %s' % fromFile            
            newAbunArr = np.fromfile(fromFile, dtype = np.float64, sep = " ")
            
        if fromArray != None:
            print 'file : %s' % fromArray            
            newAbunArr = fromArray
        
        self.copyAbundancesFromArray( newAbunArr )

    def copyAbundancesFromArray(self, array):
        n = len(array)
        self.abun[0:n, 0] = array

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
        
        if self.umistVer == 'umist06':
            lineParsed = lineStr.split(',')
        
        if self.umistVer == 'umist99':
            u06Str = self.convertRxnLineToU06FormatStr( lineStr )
            lineParsed = u06Str.split(',')
            
        rxn = reaction()
        rxn.setAllFromRxnStrArr( lineParsed )
        return rxn

    # converts a umist99 reaction string line to one as it would be formatted in the umist06 database  
    def convertRxnLineToU06FormatStr( self, line ):
        
        # function which returns the type of the reaction based on its ID (table 5 in umist99 papper)
        def getReactionType( rxnID ):
                        
            if rxnID >= 1 and rxnID <= 433:
                return 'NN'
            if rxnID >= 434 and rxnID <= 2606:
                return 'IN'
            if rxnID >= 2607 and rxnID <= 3144:
                return 'CE'
            if rxnID >= 3145 and rxnID <= 3175:
                return 'II'
            if rxnID >= 3176 and rxnID <= 3606:
                return 'DR'
            if rxnID >= 3607 and rxnID <= 3631:
                return 'RR'
            if rxnID >= 3632 and rxnID <= 3678:
                return 'AD'
            if rxnID >= 3679 and rxnID <= 3760:
                return 'RA'
            if rxnID >= 3761 and rxnID <= 3916:
                return 'PH'
            if rxnID >= 3917 and rxnID <= 3927:
                return 'CP'
            if rxnID >= 3928 and rxnID <= 4059:
                return 'CR'
            if rxnID >= 4060 and rxnID <= 4077:
                return 'CL'
            if rxnID >= 4078 and rxnID <= 4107:
                return 'TR'
            if rxnID >= 4108 and rxnID <= 4113:
                return '-'

            return ''
        
        pos = 0
        newLine = ''
        
        w = 4; 
        cmp = line[pos:(pos+4)].strip() # index
        newLine += cmp +','; pos += w
        
        rxnID = np.int32(cmp)
        cmp = getReactionType( rxnID )
        newLine += cmp +',';  

        w = 9;
        cmp = line[pos:(pos+10)].strip() #r1
        newLine += cmp +',';  pos += w
        w = 9;
        cmp = line[pos:(pos+w)].strip() #r2
        newLine += cmp +',';   pos += w
        w = 9;
        cmp = line[pos:(pos+w)].strip() #r3
        newLine += cmp +',';   pos += w
        w = 9;
        cmp = line[pos:(pos+w)].strip() #p1
        newLine += cmp +',';   pos += w
        w = 9;
        cmp = line[pos:(pos+w)].strip() #p3
        newLine += cmp +',';   pos += w
        w = 6;
        cmp = line[pos:(pos+w)].strip() #p4
        newLine += cmp +',';   pos += w
        w = 5;
        cmp = line[pos:(pos+w)].strip() #alpha
        newLine += cmp +',';   pos += w
        w = 8;
        cmp = line[pos:(pos+w)].strip() #beta
        newLine += cmp +',';   pos += w
        w = 8;
        cmp = line[pos:(pos+w)].strip() #beta
        newLine += cmp +',';   pos += w
        w = 10;
        cmp = line[pos:(pos+w)].strip() #gamma
        newLine += cmp +',';   pos += w
        w = 1;
        cmp = line[pos:(pos+w)].strip() #kind of data
        newLine += cmp +',';   pos += w
        w = 5;
        cmp = line[pos:(pos+w)].strip() #T_low
        newLine += cmp +',';   pos += w
        w = 5;
        cmp = line[pos:(pos+w)].strip() #T_high
        newLine += cmp +',';   pos += w
        w = 1;
        cmp = line[pos:(pos+w)].strip() #accuracy code
        newLine += cmp +',';   pos += w
        w = 4;
        cmp = line[pos:(pos+w)].strip() #refCode
        newLine += cmp +',';   pos += w

        #print '---------------'
        #print line
        #print newLine
        
        return newLine
    
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
        
    def initializeDefaultAbunBuffer(self):
        self.abun = np.zeros( (len(self.species), 1) ) - 1
        print self.pstr + 'Reset default abundances.'

    # set all the numbers of the species self.species[].num to None
    def setSpeciesNumsToNone(self):
        
        for specStr in self.species:      
            self.species[specStr].num = None
