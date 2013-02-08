import numpy
from specie import specie

# class definition for a single reaction
# --------------------------------------
class reaction():
    """This class parses a reaction line read (as a string) into its components 
       (reactants, prodcuts, temperature ranges, constants...etc..). It also
       provides a method to print the reaction in a pretty format. 
    
       .. todo:: write a method which checks if the reaction is balanced or not
    """
    def __init__(self):
        self.str = ''      #: The string representation of the reaction
        self.status = -1   #: NA 
        self.ID   = -1     #: ID of the reaction, read from the reaction file
        self.hash = numpy.uint64(0) #: a number generated from the products and the reactants (useful to check for identical reactions efficiently)
        #--------------
        self.type = ''  #: umist (or customly defined ) string indicating the type of the reaction.
        #--------------
        self.reactants  = []  #: a list of strings of the reactants
        self.nReactants = None  #: number of reactants
        
        #--------------
        self.products  = [] #: a list of strings of the products
        self.nProducts = None #: number of the products
        
        self.species = {} #: a dictionary whose keys are the reactants and the products.
        self.nSpec = None #: number of species involved in the reaction. 
        
        #--------------
        self.alpha = 0.0 #: the first constant in the reaction line.
        self.beta  = 0.0 #: the second constant in the reaction line.
        self.gamma = 0.0 #: the third constant in the reaction line.
        self.cst   = None #: the reaction constant computed from alpha, beta, gamma (and maybe other parameters)
        self.rate  = None #: the reaction rate computed from self.cst and the abundances (and maybe other parameters) 
        self.src   = '' #: a string indicating the literature source of the reaction.
        self.Tl = -1.0 #: lower temperature range where the reaction is valid
        self.Tu = -1.0 #: upper temperature range where the reaction is valid
        #--------------
        self.accuracy = ''  #: umist error code of accuracy of the reaction constatn
        self.naccuracy = '' #: numeric error of the reaction constant
        #--------------
        self.refCode = '' #: ??? (not used anywhere)
        self.complements = None
        """list of other reactions objects which have the same hash codes, but are 
           have reaction constants defined over other (non-overlapping) complementary
           temperature ranges. If this is None, then the reaction does not have any
           complements.
        """
        
        self.all_rxns = None
        "same as self.complements with self appended to it"

        self.bounds = {'below':None, 'above':None} 
        """dictonary which helps decide which reaction to use when the the temperature
           at which the reaction constant will be computed is not withing the range 
           defined in the reaction file. (See set_bound_rxns)
           
           .. code-block:: python
           
               self.bouds['below']['T']    #the minimum temperature where the reaction is defined
               self.bouds['below']['rxn']  #the reaction which has the constants defined for this minimum temperature
               
               self.bouds['above']['T']    #the maximum temperature where the reaction is defined
               self.bouds['above']['rxn']  #the reaction which has the constants defined for this maximum temperature
   
        """  
        
        self.Tlb = None;
        """This is the temperature value used as a minimum for computing the reaction constant.
           i.e if a temperature lower than this is passed to compute the reaction constant, this
           self.Tlb is used instead. This is usually use by the function which computes the 
           constant self.compute_cst_func
        """
        
        self.compute_rxn_cst_func = None
        """The function which computes the reaction constant. It should take as first  
           argument a dictonary which has all the parameters it needs to compute 
           the reaction constant.
        """
        
        self.rxn_in_trng = None #: the reaction which is chosen to compute the reaction rate
         
    def set_all_attr_from_rxn_str(self, rxnStr):
        """set all the attributes from the string array of the reaction line"""

        self.str = rxnStr
        
        rxnStr = rxnStr.split(',')
        
        self.active()
        self.setId  (numpy.int32(rxnStr[0]))
        self.setType(rxnStr[1])
        
        self.set_reactants_from_strings([rxnStr[2], rxnStr[3], rxnStr[4]])
        self.set_products_from_strings([rxnStr[5], rxnStr[6], rxnStr[7], rxnStr[8]])
        self.setAlpha(numpy.float64(rxnStr[9]))
        self.setBeta(numpy.float64(rxnStr[10]))
        self.setGamma(numpy.float64(rxnStr[11]))
        self.setSrc(rxnStr[12] )
        # be more careful in setting the upper and lower temperature warnings
        Tl = numpy.float64(rxnStr[13]) if len(rxnStr[13]) > 0 else -1 
        self.setTl(Tl)
        Tu = numpy.float64(rxnStr[14]) if len(rxnStr[14]) > 0 else -1 
        self.setTu(Tu)
        self.setAccuracy(rxnStr[15])
        self.setRefCode(rxnStr[16])

    def update_type(self):
        """if CUSTOM reactions are added to the reaction file, sometimes the type 
           of the reaction is missing. Here, the type of the reaction is set based
           on what is in the reactants. It is set to CP is it is a CRP reaction,
           CR is it is a CRPHOT reaction and to PH if it is a PHOTO reaction.
           otherwise, it is set to NB indicating a two or three body reaction
        """  
        
        if len(self.type) == 0:
            
            typeDict = { 'CRP'    : 'CP',
                         'CRPHOT' : 'CR',
                         'PHOTON' : 'PH' }
            defualtType = 'NB'
            
            s1 = set( typeDict.keys() )
            inter = s1.intersection( set(self.reactants))
            
            if len(inter) == 1:
                reactStr = inter.pop()
                newType = typeDict[reactStr]
                #print 'setting the type of the reaction to UMIST type %s bec on of the reactants is reactants %s' % (newType, reactStr)
                self.type = newType
            else:
                if len(inter) >= 2:
                    #str = 'Error : reaction cannot have two matching entries from typeDict'
                    raise NameError(str)
                else:
                    if self.type == '':
                        #print 'setting reaction type to default'
                        self.type = defualtType

    def set_rxn_cst_func(self, functionDict):
        """method that computes the reaction constant.
        
           .. todo:: SETS THE FUNCTION WHICH COMPUTES THE REACTION CONSTANT FROM A FUNCTION NAME DICTIONARY
           .. warning:: not implmented yet. 
        """  
        pass

    def compute_reaction_constant(self, state):
        """method that computes the reaction constant using the constants from the the 
           appropriate temperature range.  The computed value is set to self.cst
           and to the reaction which was used to compute the constant.
        
           :param dict state: a dictionary holding all the parameters needed by the
            functions which compute the constants to do their job. 
            
        """
        
        rxn = self.get_rxn_in_trange(state['T'])
        cst = rxn.compute_rxn_cst_func(rxn, state)
        
        rxn.cst = cst
        self.rxn_in_trng = rxn
    
    def compute_reaction_rate(self, nDens):
        """computes the reaction rate from the reaction constant and the
           abundances of the reactants. The rate is set to self.rxn_in_trng.rate
           and not to self.rate
           
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
        
        reactants = self.reactants
        species   = self.species

        #mutiplyuing the reaction constant with the density of the first reactant                                                                                                                                                        
        #this is always the case for all the reactions including the PH, CP and CR ones.                                                                                                                                                 
        rate = self.rxn_in_trng.cst * (nDens * species[reactants[0]].abun())

        #if the reaction is not a PH or CP or CR, the other reactants also must                                                                                                                                                          
        #be taken into account in computing the rate                                                                                                                                                                                     
        # compute the rate for reactions involving a photon                                                                                                                                                                              
        for specStr in reactants:
            abun = species[specStr].abun()

            if abun != None:
                rate *= (nDens * abun)

        if rate == None:
            raise ValueError('reaction %s has a rate equal to None. Something has gone wrong while computing its rate' % rxn.str)

        self.rxn_in_trng.rate = rate
    
        
    def active(self):
        """sets status to 1 indication the reaction IS being used"""
        self.status=1
        
    def passive(self):
        """sets status to 0 indication the reaction IS NOT being used"""
        self.status=0
        
    def setId(self, ID):
        """sets ID and type"""
        self.ID = ID
        
    def setType(self, typ):
        """Set the type"""
        self.type=typ
        
    def set_reactants_from_strings(self, reactantsStrArr):
        """set the reactants"""
        for specStr in reactantsStrArr:
            if len(specStr) != 0:
                self.species[ specStr] =  specie(specStr)
                self.reactants.append( specStr )
        
        self.nReactants = len(self.reactants)
            
    def set_products_from_strings(self, productsStrArr):
        """set the products"""

        for specStr in productsStrArr:
            if len(specStr) != 0:
                self.species[ specStr] =  specie(specStr)
                self.products.append( specStr ) 
            
        self.nProducts = len(self.products)
        self.nSpecs = self.nReactants + self.nProducts
    
    def update_specie_in_reaction(self, specStr, specObj):
        self.species[specStr] = specObj
        
    def setAlpha(self, alpha):
        """set alpha"""
        self.alpha=alpha
    def setBeta(self, beta):
        """set beta"""
        self.beta=beta
    def setGamma(self, gamma):
        """set gamma"""
        self.gamma=gamma
    def setSrc(self, src):
        """set src"""
        self.src=src
    def setTl(self, Tl):
        """set Tl"""
        self.Tl = Tl
    def setTu(self, Tu):
        """set Tu"""
        self.Tu = Tu
    def setAccuracy(self, acc):
        """set the accuracu and the refcode"""
        self.accuracy = acc
    def setRefCode(self, refCode):
        """set the reference code"""
        self.refCode = refCode
        
    def has_ID(self, ID_check):
        """Returns True if the reaction ID is the same as the input ID. and if the ID
           matches one of the complement reactions (if any) in self.complements.
           
           .. todo:: see if using rxn.has_ID( ID_check) for the cmplements below
             is a good idea. Snice it might be recursive if self.complements are linked
             to the current reactions!!.
        """
        
        IDs = [self.ID]
        
        if self.has_complements():
            for rxn in self.complements:
                IDs.append(rxn.ID)
                    
        if ID_check in IDs:
            return True
        else:
            return False
    
    def set_bound_rxns(self):
        """In this method, we setup a dictionary of complementaty reactions. They
           keys of the dictonary are 'below' and 'above'. Those point to the
           reaction objects which are used when the temperature is above or 
           below the gloab ranges where the reactions are defined: For example, if 
           a reaction with ranges [30,300] for rxn1 and [300, 1000] for rxn2, 
           and [1000,30000] for rxn3. Then 'above' would be rxn3 and below would
           be rxn1.  Such that the rates of those are used whenever the temperature
           range is not withing 30 and 30000.
        """
        
        below = {'T':self.Tl, 'rxn':self}
        above = {'T':self.Tu, 'rxn':self}

        #collecting the reactions into a list (including self and the complements)
        if self.complements != None:
            rxns = self.all_rxns
        else:
            rxns = [self]
        
        #finding and setting the min and max and the corresponding reactions   
        for rxn in rxns:
            if rxn.Tl < below['T']:
                below['rxn'] = rxn
            if rxn.Tu > above['T']:
                above['rxn'] = rxn

        self.bounds['below'] = below
        self.bounds['above'] = above
            
    def get_rxn_in_trange(self, T):
        """Returns the reaction among self or self.complements (if any) where the
           temperature range of the reaction is within the value (or range) passed.
        """

        #if the temperature is out of bounds
        if T < self.bounds['below']['T']:
            return self.bounds['below']['rxn']
        if T > self.bounds['above']['T']:
            return self.bounds['above']['rxn']
        
        #collecting the reactions into a list (including self and the complements)
        if self.complements != None:
            rxns = self.all_rxns
        else:
            rxns = [self]
        
        #checking if the temperature ranges match one
        #of those of the reaction and it complements
        #(Note that the Tu comparison is not inclusive)      
        for rxn in rxns:
            if T >= rxn.Tl and T < rxn.Tu:
                return rxn
        
        #the code should not get here
        raise ValueError("""can determin what reaction to use to compute the reaction 
                            constant. This might be because there are gaps between the
                            complements (defining the temperature ranges""")
    
    def get_rxn_with_ID(self, ID):
        """returns the reaction which mateched ID. A value error is raise if it is not found.
           Also a vlaue error is raise if self.complement == None. It is also assumed that
        """
        
        if self.complements == None:
            raise ValueError("This method should be called only if self.complements is NOT None.")
        else:
            if self.ID == ID:
                return self
            else:
                for rxn in self.complements:
                    if rxn.ID == ID:
                        return rxn
                    
        raise ValueError("""The code should not get here. If it does, it means that the ID %d  
                            which is being looked for is not in the composite reacation.""" % ID)
    
    def has_complements(self):
        """Returns True if self.complements is not None. Returns False if it is None"""
        if self.complements == None:
            return False
        else:
            return True
                
    def display(self, fmt = None ):
        """display a reaction based on the input format requested.
        
           .. todo:: document the parameter fmt
           
        """

        # functions which print the compoments of a reaction
        def printStatus () : print "%d"   % self.status,
        def printId     () : print "%04d" % self.ID,
        def printHash   () : print "%20d" % self.hash,
        def printType   () : print "%2s"  % self.type,
        def printReacts () :
            # collecting the strings to be printed into a tuple and filling the 
            # reactants which do not exist by empty strings
            cmpnts = ()
            fmtStr = ""
            for i in [0,1,2]:
                if i < self.nReactants:
                    cmpnts  += (self.reactants[i],)
                    fmtStr  += "%-10s  "   
                else:
                    cmpnts += ('',)
                    fmtStr += "%-10s  " 
            print fmtStr % cmpnts,
        def printProds  () :
            # collecting the strings to be printed into a tuple and filling the 
            # species which do not exist by empty strings
            cmpnts = ()
            fmtStr = ""

            for i in [0,1,2,3]:
                if i < self.nProducts:
                    cmpnts += (self.products[i],)
                    fmtStr += "%-10s  "   
                else:
                    cmpnts += ('',) 
                    fmtStr += "%-10s  "   
            print fmtStr % cmpnts,
        def printRxn    () : 
            printReacts()
            print ' --> ',
            printProds(),
        def printReactsNumeric () :
            # collecting the numeric reprentations of the species to be printed into 
            # a tuple and filling the spcies which do not exist by empty strings
            cmpnts = ()
            fmtStr = ""
            for i in [0,1,2]:
                if i < self.nReactants:
                    
                    n = self.species[ self.reactants[i] ].num 
                    
                    if n != None:
                        cmpnts += ( self.species[ self.reactants[i] ].num,)
                        fmtStr += "%-10d  "
                    else:
                        cmpnts += ( self.species[ self.reactants[i] ].str,)
                        fmtStr += "%-10s  "
                else:
                    cmpnts += ('',)
                    fmtStr += "%-10s  "
                    
            print fmtStr % cmpnts,
        def printProdsNumeric  () :
            # collecting the numeric reprentations of the species to be printed into 
            # a tuple and filling the spcies which do not exist by empty strings
            cmpnts = ()
            fmtStr = ""
            for i in [0,1,2,3]:
                if i < self.nProducts:
                    
                    n = self.species[ self.products[i] ].num 
                    if n != None:
                        cmpnts += (self.species[ self.products[i] ].num,)
                        fmtStr += "%-10d"    
                    else:
                        cmpnts += (self.species[ self.products[i] ].str,)
                        fmtStr += "%-10s  "   
                else:
                    cmpnts += ('',)
                    fmtStr += "%-10s  "   
            print fmtStr % cmpnts,
        def printRxnNumeric    () :
            printReactsNumeric()
            print ' --> ',
            printProdsNumeric(),

        def printABG    () : print "|%+-5.2e %+-5.2e %+-5.2e|" % (self.alpha, self.beta, self.gamma), 
        def printTrng   () : print "|%-5d %-5d|" % (self.Tl, self.Tu), 
        def printAcc    () : print " %s " % ( self.accuracy ),
        def printRef    () : print " %s " % ( self.refCode ),
        def printCst    () :
            if self.cst == None:
                print "    NA   ",  
            else:
                print " %+-15.8e" % ( self.cst ),
        def printRate   () :
            if self.rate != None:
                print " %+-15.8e" % ( self.rate ),
            else:
                print "        NA       ",
        def printAbun   () : 

            print ''
            for spec in self.reactants:
                if spec._abun != None:
                    print "%10s : %+-15.8e" % ( spec.str, spec._abun )
                else:
                    print "     NA     "
            for spec in self.products:
                if spec._abun != None:
                    print "%10s : %+-15.8e" % ( spec.str, spec._abun )
                else:
                    print "     NA     "
                    
        def printStr():
            print self.str

        def printTlb():
            if self.Tlb != None:
                print " %.4e " % self.Tlb,
            else:
                print " None ",
            
        action = {
                 "status"     : printStatus,
                 "id"         : printId,
                 "hash"       : printHash,
                 "type"       : printType,
                 "reacts"     : printReacts,
                 "prods"      : printProds,
                 "rxn"        : printRxn,
                 "rxnNumeric" : printRxnNumeric,
                 "abg"        : printABG,
                 "trng"       : printTrng,
                 "acc"        : printAcc,
                 "ref"        : printRef,
                 "cst"        : printCst, 
                 "rate"       : printRate,
                 "abun"       : printAbun,
                 "str"        : printStr,
                 "Tlb"        : printTlb,
            }
        
        if self.complements != None:
            print '---------------------------------------------------------------------' 
        
        if fmt != None:    
            fmtStrSplt = fmt.split()
            for fmtComp in fmtStrSplt:
                action[fmtComp.strip()]()
            print 
        else:
            print self.str
        
        #printing also the complement reactions
        if self.complements != None:
            for rxn in self.complements:
                rxn.display(fmt = fmt)
            print '---------------------------------------------------------------------' 

        
        