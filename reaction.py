from string import *
from numpy import *
import StringIO
import sys
import subprocess
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
               self.setReactants(R1, R2, R3)
                 self.setR1(R1)
                 self.setR2(R2)
                 self.setR3(R3)
               self.setReactants(P1, P2, P3, P4)
                 self.setP1(P1)
                 self.setP2(P2)
                 self.setP3(P3)
                 self.setP4(P4)
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
        self.R1   = None        
        self.R2   = None
        self.R3   = None
        self.reactants = []
        #--------------
        self.P1   = None
        self.P2   = None
        self.P3   = None
        self.P4   = None
        self.products = []
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
    # set the reactant strings
    def setR1(self, R1):
        self.R1=R1
    def setR2(self, R2):
        self.R2=R2
    def setR3(self, R3):
        self.R3=R3
    def setReactants(self, R1, R2, R3):
        self.R1=specie(R1)
        self.R2=specie(R2)
        self.R3=specie(R3)
    # set the product strings
    def setP1(self, P1):
        self.P1=P1
    def setP2(self, P2):
        self.P2=P2
    def setP3(self, P3):
        self.P3=P3
    def setP4(self, P4):
        self.P4=P4
    def setProducts(self, P1, P2, P3, P4):
        self.P1=specie(P1)
        self.P2=specie(P2)
        self.P3=specie(P3)
        self.P4=specie(P4)
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

    # set all the attributes from the string array of the reaction line
    def setAllFromRxnStrArr(self, rxnStr):
        self.active()
        self.setId  ( int32(rxnStr[0]) )
        self.setType( rxnStr[1]        )
        self.setReactants( rxnStr[2], rxnStr[3], rxnStr[4] )
        self.setProducts( rxnStr[5], rxnStr[6], rxnStr[7], rxnStr[8] )
        self.setAlpha( float64(rxnStr[9])  )
        self.setBeta(  float64(rxnStr[10]) )
        self.setGamma( float64(rxnStr[11])  )
        self.setSrc( rxnStr[12] )
        self.setTl( float64(rxnStr[13]) )
        self.setTu( float64(rxnStr[14]) )
        self.setAccuracy( rxnStr[15] )
        self.setRefCode( rxnStr[16] )

    def display(self, format = None ):

        # methods which print the compoments of a reaction
        def printStatus () : print "%d"   % self.status,
        def printId     () : print "%04d" % self.id,
        def printHash   () : print "%20d" % self.hash,
        def printType   () : print "%2s"  % self.type,
        def printReacts () : print "%-10s  %-10s  %-10s" % (self.R1.str, self.R2.str, self.R3.str),
        def printProds  () : print "%-10s  %-10s  %-10s  %-10s" % (self.P1.str, self.P2.str, self.P3.str, self.P4.str),
        def printRxn    () : 
            printReacts()
            print ' --> ',
            printProds(),
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
            "status" : printStatus,
            "id"     : printId,
            "hash"   : printHash,
            "type"   : printType,
            "reacts" : printReacts,
            "prods"  : printProds,
            "rxn"    : printRxn,
            "abg"    : printABG,
            "trng"   : printTrng,
            "acc"    : printAcc,
            "ref"    : printRef,
            "cst"    : printCst, 
            "rate"   : printRate,
            "abun"   : printAbun}


        if format != None:
            fmtStrSplt = format.split(' ')
            for fmtComp in fmtStrSplt:
                action.get( strip(fmtComp) )()
            print 
        else:
            print 'NOT IMPLEMENTED YET, IMPLENT PRINTINTG THE WHOLE UMIST LINE'
        
