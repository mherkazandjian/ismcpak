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
       
       The reaction types supported are those of UMIST2006. They are only used when computing the
       reaction constants and rates. When a UMIST99 reaction network file is provided, they are 
       assigned types bases on the reaction ID. 
       
       .. warning:: The user might need to check the method self.convert_rxn_line_from_U99_to_U06_format_str to
        make sure it is done correctly in case the reaction file is customised.
        
        The types of reaction and their abbreviations are listed in Table-7 of the UMIST2006 paper.   
        
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

        self.read_network_file(networkFname)   # read the database into a buffer
        self.parse_reactions()                # parse the reaction lines from the database
        self.get_unique_species(baseSpecies)   # fileter the unique species and parse them in their components
        self.assign_numbers_to_species()        # assingns the index of the specie in the dictionary
        self.set_reaction_hash_codes()          # compute the hashcodes of the reactions
        self.update_reactions_types()          # make sure reaction types are correct
        self.check_reactions_types()           # make sure reaction types are correct
         
    def read_network_file(self, fileName):
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

    def parse_reactions(self):
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
            rxn = self.set_rxn_attributes(line)    # create a new  rxn object from reaction line
            self.reactions.append(rxn)           # append the new rxn to the reactions
            
        print 'complete'

    def get_unique_species(self, baseSpecies):
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
                    spec.get_components(baseSpecies)      # parsing the specie
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

    def remove_species(self, species = None, underAbunFile = None):
        """Remove species from the network. All reactions containing those species are moved to 
           the self.reactionsRemoved list. Also those species are moved from self.species to 
           self.removedSpecies list.
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
        self.assign_numbers_to_species()
        
        print 'complete'
            
    def assign_numbers_to_species(self, fileName = None):
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
        
        # getting the list species which did not get a number assigned to
        specsWithNoNum = []
        for spec in self.species.values():
            if spec.num == None:
                specsWithNoNum.append(spec.str)
                
        #printing the species with no numbers (10 per line)
        if len(specsWithNoNum) != 0:
            print self.__pstr + 'A total of %d species will have no numbers' % len(specsWithNoNum)
            
            print self.__pstr,
            
            for i, specStr in enumerate(specsWithNoNum):
                print ' %s' % specStr,
                if ((i+1) % 10) == 0:
                    print '\n%s' % self.__pstr,
            
            print 
                     
        self.map_species_abundances()
        self.__pstr = ''        
        
        
    def map_species_abundances(self):
        """Map self.species[]._abun to self.abun, so that changing the abundace of the element
           'X' in self.species[] would change the value in self.abun[self.specie['X'].num] and 
           vise versa.
        
           .. warning::  species which do not have a number, get their abundance set to None.
             i.e they are not mapped. 
        
        """

        for specStr in self.species:
            specNum = self.species[specStr].num
            if specNum != None:
                self.species[specStr]._abun = self.abun[ specNum ]
            else:
                self.species[specStr]._abun = None
        print self.__pstr + 'mapped spcies abundances to the abundance array'
        
    def set_reaction_hash_codes(self):
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

    def update_reactions_types(self): 
        """Checks and updates the types of the reactions. Reactions whose types are not
           within those of UMIST2006 and are custom types for example, get their types
           set based on their reactants. This is done by calling the method 
           reaction.update_type() for each reaction."""
        
        # maybe it is better if the dictonary definig the reaction types
        # based on what is in the reactants here instead of in the reaction class
        for rxn in self.reactions:
            rxn.update_type()
        print 'updated reaction types'
        
    def check_reactions_types(self): 
        """Checks the types of the reactions and raises an error if a reaction without 
           a type is detected. This should be called before trying to compute constants."""
        
        for rxn in self.reactions:
            if rxn.type == '':
                strng = 'reaction types, reaction :\n%s\n has no type set' % rxn.str
                raise NameError(strng)
        
        print 'reaction check completed...all passed.'
                    
    def find_identical_reactions(self):
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
    
    def set_abundances(self, fromFile = None, fromArray = None):
        """Set the abundances of the species from an ascii input file or from a numpy array.
        
           :param string fromFile: The path of the file which contains one abundance which will be
             assigned to the species. The file should have one abundance on each line. Each value 
             will be set to the corresponding vlaue in self.abun[] at the index corresponding to 
             the line number (staring with index 0). There should be as much lines as there are 
             species with assigned numbers. 
             
           :param numpy.ndarray fromArray: A numpy array holding the abundances of the speceies. The 
             values of this array will be copies element-wise to self.abun[]. They must have the same
             length. 
           
           .. warning::  since self.specie[]._abun are mapped to self.abun, in chanfing all the values
             of the array the method self.copy_abundances_from_array( array ) should be used, which
             sets the values one by one, instead of replacing self.abun with another array when
             self.abun = newAbun which removes the mapping.
        """
        
        print 'Setting the abundances of the species from:',
        
        newAbunArr = None
        
        if fromFile != None:
            print 'file : %s' % fromFile            
            newAbunArr = numpy.fromfile(fromFile, dtype = numpy.float64, sep = " ")
        
        if fromArray != None:
            print 'array (id) : %d' % id(fromArray)            
            newAbunArr = fromArray
        
        self.copy_abundances_from_array( newAbunArr.flatten() )

    def copy_abundances_from_array(self, array):
        n = len(array)
        self.abun[0:n, 0] = array[:]

    def printFile(self):
        """prints the reaction file without the header."""

        print self.fileStr

    def print_reactions(self, rxnIDs = None, fmt = None):
        """Prints the reaction file without the header. NOTE: this should not be
           confused with reaction index (i.e the position in the self.reactions[] list). 
        
        :param int rxnIDs: an integer or an integer list (or numpy array) of reaction IDs to be printed.
        
        :param string fmt: a string for the format of the reactions (see reaction object doc). IF
           this is not passed, a default format of 'type rxn' is used. 
        
        """
        
        if fmt == None:
            fmt = 'type rxn'
            
        if rxnIDs == None: # print all the reactions
            print 'total number of reactions = ', len(self.reactions)
            for rxn in self.reactions:
                rxn.display(fmt = fmt)
        else: # print reactions with specified ids
                        
            #if rxnIDs is just an interger, putting it into a list
            if hasattr(rxnIDs, '__iter__') == False:
                rxnIDs = numpy.array([rxnIDs])
            
            for rxnID in rxnIDs:
                
                found = False

                for rxn in self.reactions:
                    if rxnID == rxn.ID:
                        found = True 
                        rxn.display(fmt = fmt)
                        break
            
                if found == False:
                    raise ValueError('No reaction found with ID = %d' %rxnID)
                    
    def print_species(self, removed = None):
        """print all the species"""
        
        if removed == None:
            specs = self.species.values()
        else:          
            specs = self.speciesRemoved.values()
            
        for spec in specs: 
            spec.show()
            
    def get_specie_index(self, specStr):
        """Returns the index of the specie which matches the input string."""

        if specStr in self.species: #i = 0
            return self.species[specStr].num
        else:
            print 'Specie ' + specStr + ' not found'
            return None

    def set_rxn_attributes(self, lineStr):
        """Pasrse the UMIST reaction line and return the reaction object."""

        if self.umistVer == 'umist06':
            rxnStr06 = lineStr
        
        if self.umistVer == 'umist99':
            rxnStr06 = self.convert_rxn_line_from_U99_to_U06_format_str(lineStr)
        
        rxn = reaction()
        rxn.set_all_attr_from_rxn_str(rxnStr06)
        return rxn

    def convert_rxn_line_from_U99_to_U06_format_str( self, line ):
        """Converts a umist99 reaction string line to one as it would be formatted in the umist06
           database. Reactions withs IDs >= 4108 and <= 4113 are assigned the type SU short for  
           Sundries.
        """  
        
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
                return 'SU'

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
    
    def filter_reactions(self, withReacts = None, withProds = None, withType = None, 
                               show = None, sort = None, **kwargs):
        """Method that returns the IDs (as a numpy int32 array) of the reactions containing certain reactants
           or products or both. Either one of those two keywords should be provided. The
           matching is done using &&, so all the entries must be present in the reaction.
           
           :param withReacts: a string or a list of strings holding the reactants to be used as a filter.
           
           :param withProds: a string or a list of strings holding the products to be used as a filter.
           
           :param bool show: if this is True, the filtered reactions are printed. If it is False (or None). 
              only the ids are returned. By default *kwargs* is passed to the function which prints the reactions,
              so the format of the printed reactions can be determined by passing the format
              argument which internallty is passed to self.print_reactions.
           
           :param string withType: the reaction type against which filtering will be done.           
           
           :param bool sorted: if this is set to True, the reactions are printed in decreasnig 
             reaction rate order.

           for example:
           
           .. code-block:: python

               #filter with prodcuts only
               inds = filter_reactions(withProds = 'H2')
           
               #filter with reactants only
               inds = filter_reactions(withReacts = ['H','CN-'])
           
               #filter with some reactants and products
               inds = filter_reactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'])

               #filter with some reactants and products
               inds = filter_reactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'], show = True)

               #filter with some reactants and products
               inds = filter_reactions(withReacts = ['H','CN-'], withProds = ['HCN','e-'], 
                                      show = True, format = 'type rxn')
               
        """
        
        rxnIDsFound = []
        
        #checking if withType is of the correct type
        if hasattr(withType, '__iter__'):
            raise ValueError('withType can not be a list or an iterable object, it must be a string.')
        
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
            
            #checking for the type
            if withType != None:
                if rxn.type == withType:
                    hasType = True
                else:
                    hasType = False
            else:
                hasType = True
                
            #appending the reaction ID to a list
            if inReacts and inProds and hasType:
                rxnIDsFound.append( rxn.ID )

        if show == True:
            
            if sort == True:
                rxnIDsFound = self.sort_rxns_decreasing_rates(rxnIDsFound)
    
            self.print_reactions(rxnIDsFound, **kwargs)
                
        return numpy.array(rxnIDsFound)

    def compute_rxn_constants(self):
        """Compute the reaction constants from the state of the gas in the environment. The
           equation used in computing the rate is based on the type of the reaction (reaction.type).
           In case a custom reaction is added make sure to modify this method to compute the 
           reaction constant correctly.
           
           .. note:: The attributes self.T, self.Av and self.albedo must be set before call this method.
           
           .. todo:: eventually this should be done in the reaction calss, where the state of the 
            environment is passed as an argument which is used to compute the constant.
            
           .. todo:: document the actual equations (in latex) for the reaction constants. As in, reaction
             type and the equation for it.

        """

        #definitng the functions which compute the reaction constants. The general
        #signature of those functions should be as follows:
        #    compute_bla_bla_bla_cnst(rxn, state)
        #where state is the dict which holds the parameters used to compute the rate
        #in the functions.
        def compute_photo_rxn_cnst(rxn, state):
            return rxn.alpha * numpy.exp( - rxn.gamma * state['Av'] )
        def compute_CP_rxn_cnst(rxn, state):
            return rxn.alpha
        def compute_CR_rxn_cnst(rxn, state):
            return rxn.alpha * ((state['T'] / 300.0)**rxn.beta) * (rxn.gamma / ( 1.0 - state['albedo'] )) * (state['zeta']/1.36e-17) 
        def compute_nbody_rxn_cnst(rxn, state):
            return rxn.alpha * ((state['T'] / 300.0)**rxn.beta) * numpy.exp( - rxn.gamma / state['T'] ) 

        if self.T == None or self.Av == None or self.albedo == None:
            raise ValueError('values of attributes of T, Av and albedo are not set.')
                
        #a dictonary which sets the function for computing the reaction
        #constant based on the reaction type (see table-7 in umist06 paper
        #and table-5 in umist99)
        compute_cst = {
                       #two or three body reactions
                       'IN': compute_nbody_rxn_cnst, #umist06
                       'NN': compute_nbody_rxn_cnst, #umist06
                       'CE': compute_nbody_rxn_cnst, #umist06
                       'II': compute_nbody_rxn_cnst, #umist06
                       'DR': compute_nbody_rxn_cnst, #umist06
                       'RR': compute_nbody_rxn_cnst, #umist06
                       'AD': compute_nbody_rxn_cnst, #umist06
                       'RA': compute_nbody_rxn_cnst, #umist06
                       'CD': compute_nbody_rxn_cnst, #umist06
                       'CI': compute_nbody_rxn_cnst, #umist06
                       'IM': compute_nbody_rxn_cnst, #umist06
                       'CL': compute_nbody_rxn_cnst, #umist06
                       'TR': compute_nbody_rxn_cnst, #umist99
                       'SU': compute_nbody_rxn_cnst, #umist99 (Sundries)
                       'NB': compute_nbody_rxn_cnst, #general type (NBODY reactions)
                       #photoprocess
                       'PH': compute_photo_rxn_cnst, #umist06
                       #cosmic ray photon (CRP)
                       'CP': compute_CP_rxn_cnst,    #umist06
                       #cosmic ray photon (CRPHOT)
                       'CR': compute_CR_rxn_cnst,    #umist06
                       } 

        #the state of the environment (gas)
        state = {'Av':self.Av, 'T':self.T, 'albedo':self.albedo, 'zeta':self.zeta}
                
        # it might be a good idea to put a check if
        # temperature,A_v are set...         
        for rxn in self.reactions:
            rxn.cst = compute_cst[rxn.type](rxn, state)
            
    def compute_rxn_rates(self):
        """Compute the reaction rates assuming all the abundances are set and the reaction
           constants are already computed. 
           
           ..warning:: The rates are computed by multiplying the reaction constant with the
            DENSITY of all the reactants.  The reactants which have an abundance of None are
            skipped (this is useful for the reactions which have CRP, PHOTON, M and CRPHOT
            as reactants.  So in case custom reactions are added, there is a caveat where
            if a reaction includes any of the 4 mentioned 'species' in addition to a thrid
            reactant the rate might be computed incorrectly. For example, a reaction like : 
            
                    CO + PHOTON + XYZ ->  ....
                    
            would have the rate computed as : rxn.cst * density(CO) * density(XYZ)
            which might not be what we want to do.
        """

        for rxn in self.reactions:

            reactants = rxn.reactants
            species   = rxn.species

            
            #mutiplyuing the reaction constant with the density of the first reactant
            #this is always the case for all the reactions including the PH, CP and CR ones. 
            rate = rxn.cst * (self.nDens * species[reactants[0]].abun())

            #if the reaction is not a PH or CP or CR, the other reactants also must
            #be taken into account in computing the rate          
            # compute the rate for reactions involving a photon
            for specStr in rxn.reactants:
                abun = species[specStr].abun()
                 
                if abun != None:
                    rate *= (self.nDens * abun)

            if rate == None:
                raise ValueError('reaction %s has a rate equal to None. Something has gone wrong while computing its rate' % rxn.str)

            rxn.rate = rate

    def sort_rxns_decreasing_rates(self, rxnIDs):
        """Sort the reactions with decreasing absolute rates.
        """
        
        rates = []
        ids   = []

        for rxnId in rxnIDs:
            for rxn in self.reactions:
                if rxnId == rxn.ID and rxn.rate != None:
                    rates.append( rxn.rate )
                    ids.append( rxn.ID )
        
        #sorting the reaction with decreasing rates
        rates = numpy.abs(numpy.array(rates))
        ids   = numpy.array(ids)
        
        inds  = numpy.argsort(rates)
        idsSorted = ids[inds]
        
        return idsSorted[::-1]
    
    def get_unique_rxn_types(self):
        """Returns a dict holding the type and the counts of the available
           reactions types in the network"""
        
        rxnTypes = {}
        
        for rxn in self.reactions:
            
            if rxn.type in rxnTypes: 
                rxnTypes[rxn.type] += 1
            else:
                rxnTypes[rxn.type] = 1
            
        return rxnTypes
        
    def set_environment_state(self, T = None, zeta = None, Av = None, 
                                    albedo = None, nDens = None, G0 = None):
        """This method is used to set the properties of the environment, mainly needed
           to compute the reaction constants and rates.
        """ 
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
