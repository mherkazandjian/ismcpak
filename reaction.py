import string
import numpy
from specie import specie

# class definition for a single reaction
# --------------------------------------
class reaction():
    """asdasdadasda
    
       .. todo:: write a method which checks if the reaction is balanced or not
    """
    def __init__(self):
        self.str = ''  #: The string representation of the reaction
        self.status = -1
        self.ID   = -1
        self.hash = numpy.uint64(0)
        #--------------
        self.type = ''  # umist string for the reaction
        self.ntype = '' # numeric type/code of reaction
        #--------------
        self.reactants  = []  # a string list of the reactants
        self.nReactants = None
        
        #--------------
        self.products  = [] # a string list of the products
        self.nProducts = None
        
        self.species = {} # a dictionary of the species
        self.nSpec = None
        
        #--------------
        self.alpha = 0.0
        self.beta  = 0.0
        self.gamma = 0.0
        self.cst   = None
        self.rate  = None
        self.src   = ''
        self.Tl = -1.0
        self.Tu = -1.0
        #--------------
        self.accuracy = ''  # umist error code of accuracy of the reaction constatn
        self.naccuracy = '' # numeric error of the reaction constant
        #--------------
        self.refCode = ''

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

    def get_reaction_constant(self, parameters):
        """method that computes the reaction constant
        
           .. todo:: RETURNS THE REACTION CONSTANT BASED ON THE TEMPERATURE AND THE REST OF THE PARAMETERS.
           .. warning:: not implemented yet 
        """  
        pass
    
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
            }
        
        if fmt != None:    
            fmtStrSplt = fmt.split(' ')
            for fmtComp in fmtStrSplt:
                action[fmtComp.strip()]()
            print 
        else:
            print self.str
        
