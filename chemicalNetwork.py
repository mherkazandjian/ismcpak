import string
import numpy
from specie import specie
from reaction import reaction

class chemicalNetwork(specie, reaction):
    """If a specie 'XXX' is not in self.speciesRemoved and self.species['XXX'].num = None
       then the reactions which have 'XXX' as a specie are still in the network but they are
       not used in computing the reaction, like for example PHOTON, it does not have an 
       abundance and we do not keep track of the number of photons, so it need not have a
       number with which its abudnance can be tracked in self.abun[]   
       
       Reaction hash codes are set based on the original numbers assigned to species and it doesnt 
       depend on the order the reactnats or the products appear in the reaction input file. When 
       computing the hash code, the species are picked in the increasing order as they appeat in
       the sorted string array (using the sorted function).
        
       .. todo::  Write a method which merges reactions with the same hashcode into one reaction 
        where in different temperature ranges different temperatures are set.
        
       .. todo::  Write a method which checks if the reactions are balanced or not.
              
       .. todo::  Write a method which reads the abundances and sets the abundance of each specie object 
        this can be implemented by appending a new abudance column to species.inp.

       .. todo::  Write a method which takes as input the species file used and filters out all the
        reaction which do not contain these species and re-numbers the species in the number representation of the chemical network.
           
       .. todo::  As a check, dump the remainders of the parsed strings to a file...it should
        mainly be black other than commas...etc...otherwise, there are stuff not taken into account,
        this can be done by erasing stuff from self.fileStr after the whole original network has 
        been worked out.
    """
    def __init__(self, fileName = None, baseSpecies = None, UMISTVER = None):
        self.nReactions = None         #: number of reactions in the network.
        self.reactions = []            #: reaction object list in the netowrk.

        self.nReactionsRemoved = None  #: number of reactions removed from the originally read network.
        self.reactionsRemoved = []     #: a list of the reaction removed from the originally read netowrk.
        
        self.nSpecs = None             #: number of active species in the network.
        self.species = {}              #: a dictionary of the active specie object list in the network.
        
        self.nSpecsRemoved = None    #: number of non active species in the network.
        self.speciesRemoved = {}     #: a dictionary of the non active species object in the network.
        
        self.fileName = fileName       #: path of the ascii file from which the reactions of the network were read.
        self.baseSpecies = baseSpecies #: a list of the species of the base objects.
        self.nDens = None              #: density of the ambient hydrogen gas.
        self.T = None                  #: gas temperature.
        self.zeta = None               #: cosmic ray ionization rate.
        self.albedo = None             #: Albedo 
        self.Av = None                 #: visual extinction
        self.G0 = None                 #: FUV fux in terms of G0
        self.abun = None               #: a numpy array of shape (nSpecies, 1), it is by default initialized to -1.
        
        self.__pstr = '' # just for printintg purposes
        
        if UMISTVER != None:
            self.umistVer = UMISTVER  
        else:
            self.umistVer = None
            """string of the version of UMIST used. The two supported values are 'umist06'
               and 'umist99'.
            """
            
        if (self.fileName != None) and (self.baseSpecies != None):
            self.setup(fileName, baseSpecies)
        
    def setup(self, networkFname, baseSpecies):
        """Setup all the chemical network."""

        self.readNetworkFile(networkFname)   # read the database into a buffer
        self.parseReactions()                # parse the reaction lines from the database
        self.getUniqueSpecies(baseSpecies)   # fileter the unique species and parse them in their components
        self.assignNumbersToSpecies()        # assingns the index of the specie in the dictionary
        self.setReactionHashcodes()          # compute the hashcodes of the reactions
        self.updateReactionsTypes()          # make sure reaction types are correct
        self.checkReactionsTypes()           # make sure reaction types are correct
         
    def readNetworkFile(self, fileName):
        """Read and parses the UMIST 2006 or UMIST 99 reaction file (without the header) and assigns 
           the variable fileStr and counts the number of reactions in the file sets : self.fileStr
           self.nRxn.
        """

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

    def parseReactions(self):
        """define rxn object list for all the reaction from the ascii database lines
           here the reactant and product objects are not set.
        """

        print 'Parsing reactions...',

        if self.umistVer == None:
            strng  = 'Error : cannot parse network reaction, unknown database\n'
            strng += '                   please set database version.'                 
            raise NameError(strng)
        
        print 'using '+ self.umistVer + ' format...',
        for line in (self.fileStr).splitlines():
            rxn = self.setRxnAttributes(line)    # create a new  rxn object from reaction line
            self.reactions.append(rxn)           # append the new rxn to the reactions
            
        print 'complete'

    def getUniqueSpecies(self, baseSpecies):
        """Construct the the speices dict for all the species in the network self.species
           also here all the unique species in the network are put in that dictionary and
           the components of the species are also extracted
        """
        print 'Analyzing the checmical network.....',
        # put the base species in the dictionary
        for baseSpec in baseSpecies:
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
        # and each entry is mapped to self.species[]._abun, such that
        # self.species['XXX']._abun == self.abun[ self.species['XXX'].num ]
        # the abundances are intiialized to -1
        self.initializeDefaultAbunBuffer() 
        
        self.nReactionsRemoved = len(self.reactionsRemoved)
        self.nSpecsRemoved = len(self.speciesRemoved)
        

    def mapReactionSpeciesToSpecieObjectList(self):
        """Mapping the species objects in all the reactions to the species object list in
           the network object (this object)...this is useful since instead of having a copy 
           of a species which might occure in many reaction, it is just referenced to the
           species in the species list it is also useful since changin an attirbute of a 
           species in the list would be visible to all the occurances in the network. such
           as... setting the number_tag of the species objects in the reaction objects such
           that it would be possible to construct the numeric representation of the reactions,
           i.e instead of having :
            
                     HCN  +  H   ->  HNC  +   H
                     
           we get    
           
                     '56' +  '9' ->   27  +  '9'
                     
           assuming HCN, HNC and H have indicies 56, 27 and 9 respectively
           and this number would be visible to all the instances where these
           species occure.
           
           .. warning::  This should be called everytime an entry in self.species is replaced 
            (i.e it addrress changes) for ex by call.
           
           .. warning:: It is highly recommended to use the specie.abun(x) method to set the abun
            of a species. Since it is easy to get confused and set the abundance via specie._abun = x
            which would break the mapping with net.abun. The way to change the abundance without
            breaking the mapping would be to use specie.abun[0] = x, since specie._abun is an ndarray
            of shape (1,) and setting the value via specie._abun = x would assign a new object instead
            of a new value.
        """
        print 'constructing numeric represntation of the network....',
        
        # setting self.reactions[i].species[:].num 
        for i in numpy.arange( len(self.reactions) ):
            for specStr in (self.reactions[i].reactants + self.reactions[i].products):
                self.reactions[i].species[specStr] = self.species[specStr]

        print 'compelte'

    def removeSpecies(self, species = None, underAbunFile = None):
        """Remove species from the network. The species strings which appear in the list are
           moved to speciesRemoved list and the reactions with the removed species are moved
           to reactionRemoved.
        """
        print 'removing species from the network....', 
        
        if species == None and underAbunFile != None:
            print 'using the list from the file %s' % underAbunFile
            # reading the species to be removed from a file
            fObj =  open(underAbunFile, 'r')
            setAsInactive = []
            for line in fObj:
                line = line.split()
                
                specStr = line[1]
                setAsInactive.append(string.strip(specStr))
            fObj.close() 

            print 'species to be removed ', setAsInactive
            species = setAsInactive
            
        if species != None:
            specStrList = species
        
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
                for i in numpy.arange(n):
                    if specStr in self.reactions[n-i-1].species:
                        self.reactionsRemoved.append( self.reactions.pop(n-i-1) )
                    
                self.nReactions = len(self.reactions)
                self.nReactionsRemoved = len(self.reactionsRemoved)

        self.__pstr = '     '
        self.initializeDefaultAbunBuffer()
        self.assignNumbersToSpecies()
        
        print 'complete'
            
    def assignNumbersToSpecies(self, fileName = None):
        """Assing the self.species[:].num to their index in the dictionary and update
           these indecis in the chemical network as well.
           
           .. warning::  it might happen that by mistake some species number is not specified
             in the input file containing the new number, so it will be None. Or that the
             file contains a species that is not in the chemical network, these will also
             have None as numbers.
             
        """
        self.setSpeciesNumsToNone()
                
        if fileName == None:
            print self.__pstr + 'setting species number in the order they appear in self.species'
            self.__pstr = '     '
            # assign the number 
            i = 0
            for specStr in self.species:
                self.species[specStr].num = i
                i += 1
        else:
            # it might be a good idea to put this in the mesh database thing or mesh object
            # as a new dictionary and pass the dictionary object to this instead of a file
            print self.__pstr + 'setting species number in the order they appear in the file: ' + fileName
            self.__pstr = '     '
            
            fObj =  open(fileName, 'r')
            i = 0
            for line in fObj:
                line = line.split()
                
                specStr = line[1]
                if specStr in self.species:
                    self.species[ specStr ].num = numpy.int32(i)
                else:
                    print self.__pstr + '%-12s from the file is not a species in the current network' % specStr 
                i += 1
            fObj.close() 
        
        print  self.__pstr + 'finished assigning the new species numbers'
        
        # printintg the species which did not get a number assigned to
        i = 0
        for spec in self.species.values():
            if spec.num == None:
                print self.__pstr + 'index of %-12s is missing in the file, it is kept as ' % spec.str, spec.num
                i += 1
                
        if i != 0:
            print self.__pstr + 'A total of %d species will have no numbers' % i 
        
        self.mapSpeciesAbundances()
        self.__pstr = ''        
        
    def mapSpeciesAbundances(self):
        """Map self.species[]._abun to self.abun, so that changing the abundace of the element
           'X' in self.species[] would change the value in self.abun[self.specie['X'].num] and 
           vise versa.
        
           .. warning::  species which do not have a number, get their abundance set to None.
        
        """

        for specStr in self.species:
            specNum = self.species[specStr].num
            if specNum != None:
                self.species[specStr]._abun = self.abun[ specNum ]
            else:
                self.species[specStr]._abun = None
        print self.__pstr + 'mapped spcies abundances to the abundance array'
        
    def setReactionHashcodes(self):
        """Each reaction has 7 species max involved (UMIST). In the UMIST2006 there are ~450
           species. Each specie can be represented by a 9 bit integer ( we use an unsigned integer)
           which can have a max value of 512. So each reaction can have a tag composed of 9*7=63 
           bits, which is unique as one of the involved in the reactions is. For example, the 
           reaction :
           
                 (0)   (1)    (2)        (3)   (4)    (5)    (6)
                 286 + 221          ---> 249 + 369  + 56
                 
           would have a hash code :
           
                 (512^0)*221 + (512^1)*286 + (512^3)*56 + (512^4)*249 + (512^5)*369
                    
        """

        print 'setting reaction hashcodes....',        
        base=numpy.uint64(512)

        for rxn in self.reactions:
            
            i = 0
            for specStr in sorted(rxn.reactants):
                rxn.hash += numpy.uint64((base**i))*numpy.uint64( rxn.species[specStr].num ) 
                i += 1

            i = 0
            for specStr in sorted(rxn.products):
                rxn.hash += numpy.uint64( base**(3 + i) )*numpy.uint64( rxn.species[specStr].num )
                i += 1
            
        print 'complete'

    def updateReactionsTypes(self): 
        """Checks and updates the types of the reactions."""
        
        # maybe it is better if the dictonary definig the reaction types
        # based on what is in the reactants here instead of in the reaction class
        for rxn in self.reactions:
            rxn.updateType()
        print 'updated reaction types'
        
    def checkReactionsTypes(self): 
        """Checks and updates the types of the reactions."""
        
        for rxn in self.reactions:
            if rxn.type == '':
                strng = 'reaction types, reaction :\n%s\n has no type set' % rxn.str
                raise NameError(strng)

        print 'reaction check completed...all passed.'
                    
    def findIdenticalReactions(self):
        """Finds identical reactions which have the same reactants and the same products
           and returns a tuple holding the reaction indecis ( (ind1,ind2), (ind1,ind2)...).
        """
        
        sets = ()
        n = self.nReactions
        hashes = numpy.zeros( n, dtype = numpy.uint64 )
        inds   = numpy.zeros( n, dtype = numpy.uint32 )
        
        # collect hashes numpy arrays        
        i = 0
        for rxn in self.reactions:
            hashes[i] = rxn.hash
            inds[i]   = i
            i += 1
        
        # indices of sorted hashes in increasing order
        indsSorted = numpy.argsort(hashes, kind='quicksort')
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
    
    def setAbundances(self, fromFile = None, fromArray = None):
        """Set the abundances of the species from an ascii input file or from a numpy array. 
           The file contains one abundance on each line. each value will be set to the corresponding
           vlaue in self.abun at the index corresponding to the line number (staring with index 0). 
           There should be as much lines as there are species with assigned numbers. Also the numpy
           array should have the same length as the species with numbers. 
           
           .. warning::  since self.specie[]._abun are mapped to self.abun, in chanfing all the values
             of the array the method self.copyAbundancesFromArray( array ) should be used, which
             sets the values one by one, instead of replacing self.abun with another array when
             self.abun = newAbun which removes the mapping.
        """
        
        print 'Setting the abundances of the species from:',
        
        newAbunArr = None
        
        if fromFile != None:
            print 'file : %s' % fromFile            
            newAbunArr = numpy.fromfile(fromFile, dtype = numpy.float64, sep = " ")
            
        if fromArray != None:
            print 'file : %s' % fromArray            
            newAbunArr = fromArray
        
        self.copyAbundancesFromArray( newAbunArr )

    def copyAbundancesFromArray(self, array):
        n = len(array)
        self.abun[0:n, 0] = array

    def printFile(self):
        """prints the reaction file without the header."""

        print self.fileStr

    def printReactions(self, rxnIds = None, fmt = None ):
        """Prints the reaction file without the header. 
        
        :param int rxnIds: an integer list of the reaction indecies to be printed.
        
        :param string fmt: a string for the format of the reactions (see reaction object doc).
        
        """

        if rxnIds == None: # print all the reactions
            print 'total number of reactions = ', len(self.reactions)
            for rxn in self.reactions:
                rxn.display( fmt = format )
        else: # print reactions with specified ids
            for rxnId in rxnIds:
                for rxn in self.reactions:
                    if rxnId == rxn.ID:
                        rxn.display(  fmt = fmt )
                        break

    def printSpecies(self, removed = None):
        """print all the species"""
        
        if removed == None:
            specs = self.species.values()
        else:          
            specs = self.speciesRemoved.values()
            
        for spec in specs: 
            spec.show()
            
    def getSpecieIndex(self, specStr):
        """Returns the index of the specie which matches the input string."""

        if specStr in self.species: #i = 0
            return self.species[specStr].num
        else:
            print 'Specie ' + specStr + ' not found'
            return None

    def setRxnAttributes(self, lineStr):
        """Pasrse the UMIST reaction line and return the reaction object."""

        if self.umistVer == 'umist06':
            lineParsed = lineStr.split(',')
        
        if self.umistVer == 'umist99':
            u06Str = self.convertRxnLineToU06FormatStr(lineStr)
            lineParsed = u06Str.split(',')
            
        rxn = reaction()
        rxn.setAllFromRxnStrArr( lineParsed )
        return rxn

    def convertRxnLineToU06FormatStr( self, line ):
        """Converts a umist99 reaction string line to one as it would be formatted in the umist06 database."""  
        
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
        cpmnt = line[pos:(pos+4)].strip() # index
        newLine += cpmnt +','; pos += w
        
        rxnID = numpy.int32(cpmnt)
        cpmnt = getReactionType( rxnID )
        newLine += cpmnt +',';  

        w = 9;
        cpmnt = line[pos:(pos+10)].strip() #r1
        newLine += cpmnt +',';  pos += w
        w = 9;
        cpmnt = line[pos:(pos+w)].strip() #r2
        newLine += cpmnt +',';   pos += w
        w = 9;
        cpmnt = line[pos:(pos+w)].strip() #r3
        newLine += cpmnt +',';   pos += w
        w = 9;
        cpmnt = line[pos:(pos+w)].strip() #p1
        newLine += cpmnt +',';   pos += w
        w = 9;
        cpmnt = line[pos:(pos+w)].strip() #p3
        newLine += cpmnt +',';   pos += w
        w = 6;
        cpmnt = line[pos:(pos+w)].strip() #p4
        newLine += cpmnt +',';   pos += w
        w = 5;
        cpmnt = line[pos:(pos+w)].strip() #alpha
        newLine += cpmnt +',';   pos += w
        w = 8;
        cpmnt = line[pos:(pos+w)].strip() #beta
        newLine += cpmnt +',';   pos += w
        w = 8;
        cpmnt = line[pos:(pos+w)].strip() #beta
        newLine += cpmnt +',';   pos += w
        w = 10;
        cpmnt = line[pos:(pos+w)].strip() #gamma
        newLine += cpmnt +',';   pos += w
        w = 1;
        cpmnt = line[pos:(pos+w)].strip() #kind of data
        newLine += cpmnt +',';   pos += w
        w = 5;
        cpmnt = line[pos:(pos+w)].strip() #T_low
        newLine += cpmnt +',';   pos += w
        w = 5;
        cpmnt = line[pos:(pos+w)].strip() #T_high
        newLine += cpmnt +',';   pos += w
        w = 1;
        cpmnt = line[pos:(pos+w)].strip() #accuracy code
        newLine += cpmnt +',';   pos += w
        w = 4;
        cpmnt = line[pos:(pos+w)].strip() #refCode
        newLine += cpmnt +',';   pos += w

        return newLine
    
    def filterReactions(self, withReacts = None, withProds = None, show = None, **kwargs):
        """Method that returns the indicies of the reactions containing certain reactants
           or products or both. Either one of those two keywords should be provided.
           
           :param withReacts: a string or a list of strings holding the reactants to be used as a filter.
           
           :param withProds: a string or a list of strings holding the products to be used as a filter.
           
           :param bool show: if this is True, the filtered reactions are printed. If it is False (or None). 
              only the ids are returned. By default *kwargs* is passed to the function which prints the reactions,
              so the format of the printed reactions can be determined by passing the format
              argument which internallty is passed to self.printReactions.
           
           for example:
           
           .. code-block:: python

               #filter with prodcuts only
               inds = filterReactions(withProds = 'H2')
           
               #filter with reactants only
               inds = filterReactions(withReacts = ['H','CN-'])
           
               #filter with some reactants and products
               inds = filterReactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'])

               #filter with some reactants and products
               inds = filterReactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'], show = True)

               #filter with some reactants and products
               inds = filterReactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'], 
                                      show = True, format = 'type rxn')
               
        """
        
        rxnIndsFound = []
        
        #checking for the type of the input parameter of the reactants
        if hasattr(withReacts, '__iter__'):
            lookForReacts = withReacts
        else:
            lookForReacts = [withReacts] 

        #checking for the type of the input parameter of the products
        if hasattr(withProds, '__iter__'):
            lookForProds = withProds
        else:
            lookForProds = [withProds] 
        
        for rxn in self.reactions:

            inReacts = True
            inProds  = True

            # checking if the specie we are looking for is in the reactants
            if withReacts != None:
                for specStr in lookForReacts:
                    if specStr in rxn.reactants:
                        inReacts = inReacts and True
                    else:
                        inReacts = inReacts and False

            # checking if the specie we are looking for is in the reactants
            if withProds != None:
                for specStr in lookForProds:
                    if specStr in rxn.products:
                        inProds = inProds and True
                    else:
                        inProds = inProds and False
                

            if inReacts and inProds:
                rxnIndsFound.append( rxn.ID )

        if show == True:
                self.printReactions(rxnIndsFound, **kwargs)
                
        return rxnIndsFound

    def computeReactionConstants(self):
        """Compute the reaction constants from the state of the gas in the environment.

           .. todo::  CONTINUE FROM HERE!!!!!!!! DO IT THE SAME WAY IT IS DONE IN THE PDR CODE THINK OF A MORE ELEGANT WAY AFTERWARDS 
        """
        
        def computePhotoRxnCnst(rxn):
            return rxn.alpha * numpy.exp( - rxn.gamma * self.Av )
        def computeCrpRxnCnst(rxn):
            return rxn.alpha
        def computeRrpPhotonRxnCnst(rxn):
            return rxn.alpha * ((self.T / 300.0)**rxn.beta) * rxn.gamma / ( 1.0 - self.albedo )
        def computeTwoBodyRxnRate(rxn):
            return rxn.alpha * ((self.T / 300.0)**rxn.beta) * numpy.exp( - rxn.gamma / self.T ) 

        reactionTypes = { 'PHOTON'   : computePhotoRxnCnst,
                          'CRP'      : computeCrpRxnCnst,
                          'CRPHOT'   : computeRrpPhotonRxnCnst,
                          'twoBody'  : computeTwoBodyRxnRate }
        
        # it might be a good idea to put a check if
        # temperature,A_v are set...         
        for rxn in self.reactions:
            id = rxn.ID
            print id
            
    def computeRates(self):
        """Compute the reaction constants and rates assuming all the abundances are set
           adapted to new variables but need to check it!!!!
        """

        for rxn in self.reactions:
            
            rxn.rate = rxn.cst
            # use a dictionary for the types of reactions instead 
            # of if/else statments
            
            reactants = rxn.reactants
            species   = rxn.species
            
            # compute the rate for reactions involving a photon
            if 'PHOTON' in rxn.reactants:
                rxn.rate *= (self.nDens * species[reactants[0]].abun())
            else:
                # compute the rate for reactions involving CRP
                if 'CRP' in rxn.reactants:
                    rxn.rate *= self.nDens * species[reactants[0]].abun()
                else:
                    # compute the constant for reactions involving CRPHOTON
                    if 'CRPHOTON' in rxn.reactants:
                            rxn.rate *= self.nDens * species[reactants[0]].abun()
                    else:
                        # compute the constant for two body reactions
                        for specStr in rxn.reactants:
                            rxn.rate *= self.nDens * species[specStr].abun()

    def sortRxnsDecreasingRates(self, rxnIds):
        """Sort the reactions with decreasing absolute rates adapted to new variables but need to check it!!!!
           adapted to new variables but need to check it!!!!
        """
        rates = []
        ids   = []

        for rxnId in rxnIds:
            for rxn in self.reactions:
                if rxnId == rxn.ID and rxn.rate != None:
                    rates.append( rxn.rate )
                    ids.append( rxn.ID )

        def sortIndices( data ):
            return sorted( range(len(data)), key = data.__getitem__, reverse=True)

        inds = sortIndices( rates )

        idsSorted = []
        for i in inds:
            idsSorted.append( ids[i] )

        return idsSorted

    def setCloudParms(self, T, zeta, Av, albedo, nDens, G0):
        self.setGasTemperature(T)
        self.setCrRate(zeta)
        self.setAlbed(albedo)
        self.setAv(Av)
        self.setG0(G0)
        self.setDens(nDens)
        
    def initializeDefaultAbunBuffer(self):
        self.abun = numpy.zeros( (len(self.species), 1) ) - 1
        print self.__pstr + 'Set abundances to default values of -1.'

    # set all the numbers of the species self.species[].num to None
    def setSpeciesNumsToNone(self):
        
        for specStr in self.species:      
            self.species[specStr].num = None
            
    def writeNetworkGephiCSV(self, fName):
        """Writes a CSV file where each in each reaction the all the specie/product 
           pairs combination is treated as a node->edge pairs. For example :
           the nodes->edge pairs in the reaction :
           
                C6H7+ +   e-   ----> C6H2 + H2   + H2  + H   (in ascii format  )
                 354  +  144          82  + 369  + 369 + 221 (in numeric format)
                 
           node, edge (list written to the CSV file
           
              - 354,82
              - 354,369
              - 354,369
              - 354,221
              - 144,82
              - 144,369
              - 144,369
              - 144,221
              
           :NOTE: species which have their number as None are not included, such as 
             PHOTON, CRP...
        """
        
        fObj = open(fName, "w")
        
        fObj.write('source;target\n')
        
        for rxn in self.reactions:
            for reactStr in rxn.reactants:
                reactNum = self.species[reactStr].num
                for prodStr in rxn.products:
                    prodNum  = self.species[prodStr].num
                    if reactNum == None or prodNum == None:
                        continue
                    else:
                        outStr = '%s;%s\n' % (reactStr, prodStr)
                    fObj.write(outStr)
                    
        fObj.close()
        print 'Wrote the Gephi readable CSV file.'
    
    # set the temp, zeta, albed, Av
    def setGasTemperature(self, T): self.T = T
    def setCrRate(self, zeta)     : self.zeta = zeta
    def setAlbed(self, albedo)    : self.albedo = albedo
    def setAv(self, Av)           : self.Av = Av
    def setG0(self, G0)           : self.G0 = G0
    def setDens(self, nDens)      : self.nDens = nDens
