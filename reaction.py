from string import *
import numpy as np
import StringIO
import re
from specie import *

"""
     methods : self.__init__()
               self.active()
               self.passive()
               self.setAllFromRxnStrArr(rxnStr)
               self.display(full=None)
               self.setId(id)
               self.setType(type)
               self.setSpeciesFromStrings( )
               self.updateSpecies( [R1, R2, R3, P1, P2, P3, P4] )
               self.setReactantsFromString(R1, R2, R3)
               self.setProductsFromString(P1, P2, P3, P4)
               self.setAlpha(alpha)
               self.setBeta(beta)
               self.setGamma(gamma)
               self.setSrc(src)
               self.setTl(Tl)
               self.setTu(Tu)
               self.setAccuracy(acc)
               self.setRefCode(refCode)
"""

# class definition for a single reaction
# --------------------------------------
class reaction():
    def __init__(self):
        self.str = ''
        self.status = -1
        self.id   = -1
        self.hash = int64(0)
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

    # set all the attributes from the string array of the reaction line
    def setAllFromRxnStrArr(self, rxnStr):
        self.str = rxnStr
        self.active()
        self.setId  ( int32(rxnStr[0]) )
        self.setType( rxnStr[1]        )
        
        self.setReactantsFromStrings( [ rxnStr[2], rxnStr[3], rxnStr[4] ] )
        self.setProductsFromStrings(  [ rxnStr[5], rxnStr[6], rxnStr[7], rxnStr[8] ])
        self.setAlpha( float64(rxnStr[9])  )
        self.setBeta(  float64(rxnStr[10]) )
        self.setGamma( float64(rxnStr[11])  )
        self.setSrc( rxnStr[12] )
        self.setTl( float64(rxnStr[13]) )
        self.setTu( float64(rxnStr[14]) )
        self.setAccuracy( rxnStr[15] )
        self.setRefCode( rxnStr[16] )
    
    # sets status to 1 indication the reaction IS being used
    def active(self):
        self.status=1
    # sets status to 0 indication the reaction IS NOT being used
    def passive(self):
        self.status=0
    # sets id and type
    def setId(self, id):
        self.id=id
    def setType(self, type):
        self.type=type
        
    # set the reactants
    def setReactantsFromStrings(self, reactantsStrArr):
        
        for specStr in reactantsStrArr:
            if len(specStr) != 0:
                self.species[ specStr] =  specie(specStr)
                self.reactants.append( specStr )
        
        self.nReactants = len(self.reactants)
            

    # set the products
    def setProductsFromStrings(self, productsStrArr):

        for specStr in productsStrArr:
            if len(specStr) != 0:
                self.species[ specStr] =  specie(specStr)
                self.products.append( specStr ) 
            
        self.nProducts = len(self.products)
        self.nSpecs = self.nReactants + self.nProducts
    
    def updateSpecieInReaction(self, specStr, specObj):
        self.species[specStr] = specObj
        
    # set alpha, beta, gamma
    def setAlpha(self, alpha):
        self.alpha=alpha
    def setBeta(self, beta):
        self.beta=beta
    def setGamma(self, gamma):
        self.gamma=gamma
    # set src
    def setSrc(self, src):
        self.src=src
    # set Tl, Tu
    def setTl(self, Tl):
        self.Tl = Tl
    def setTu(self, Tu):
        self.Tu = Tu
    # set the accuracu and the refcode
    def setAccuracy(self, acc):
        self.accuracy = acc
    def setRefCode(self, refCode):
        self.refCode = refCode
 
    def display(self, format = None ):

        # methods which print the compoments of a reaction
        def printStatus () : print "%d"   % self.status,
        def printId     () : print "%04d" % self.id,
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
            for i in [0,1,2,4]:
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
                    cmpnts += ( self.species[ self.reactants[i] ].num,)
                    fmtStr += "%-10d  "    
                else:
                    cmpnts += ('',)
                    fmtStr += "%-10s  "   
            print fmtStr % cmpnts,
        def printProdsNumeric  () :
            # collecting the numeric reprentations of the species to be printed into 
            # a tuple and filling the spcies which do not exist by empty strings
            cmpnts = ()
            fmtStr = ""
            for i in [0,1,2,4]:
                if i < self.nProducts:
                    cmpnts += (self.species[ self.products[i] ].num,)
                    fmtStr += "%-10d"    
                else:
                    cmpnts += ('',)
                    fmtStr += "%-10s  "   
            print fmtStr % cmpnts
        def printRxnNumeric    () :
            printReactsNumeric()
            print ' --> ',
            printProdsNumeric(),

        def printABG    () : print "|%+-5.2e %+-5.2e %+-5.2e|" % (self.alpha, self.beta, self.gamma), 
        def printTrng   () : print "|%-5d %-5d|" % (self.Tl, self.Tu), 
        def printAcc    () : print " %s " % ( self.accuracy ),
        def printRef    () : print " %s " % ( self.refCode ),
        def printCst    () : print " %+-15.8e" % ( self.cst ),
        def printRate   () :
            if self.rate != None:
                print " %+-15.8e" % ( self.rate ),
            else:
                print "        NA       ",
        def printAbun   () : 

            print ''
            for spec in self.reactants:
                if spec.abun != None:
                    print "%10s : %+-15.8e" % ( spec.str, spec.abun )
                else:
                    print "     NA     "
            for spec in self.products:
                if spec.abun != None:
                    print "%10s : %+-15.8e" % ( spec.str, spec.abun )
                else:
                    print "     NA     "


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
            "abun"       : printAbun}


        if format != None:
            fmtStrSplt = format.split(' ')
            for fmtComp in fmtStrSplt:
                action.get( strip(fmtComp) )()
            print 
        else:
            print 'NOT IMPLEMENTED YET, IMPLENT PRINTINTG THE WHOLE UMIST LINE'
        
