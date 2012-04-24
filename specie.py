from string import *
from numpy import *
import re

""" methods :  self.__init__(specStr, specType=None, charge=None, init=None)
               self.show()
               self.getComponents(baseSpec)
               self.get Abun(), setAbun( x ), getNum()
"""

# class definition for a single specie
# --------------------------------------
class specie():
    def __init__(self, specStr, specType=None, charge=None, init=None, comp=None):
        self.type = specType   # -2 => ignore specie
                               # -1 => do not conisder it as a specie
                               #  0 => basic specie
                               # 1 => composite specie
        self.str  = specStr    # string representation of the specie
        self.num  = None         # numeric representation of the specie

        if comp == None:
            if specType == 0:               # a list holding the indicies of the basic species making 
                self.comp  = [[specStr,1]]   # it up and the number number of each sub-specie
            else:
                self.comp = []
        else:
            self.comp = comp

        self.charge = charge
        self.abun = None     # abundance of the specie relative to total H nuclei
        self.init = None     # flag which indicates if the specie is inetialized or not
        self.active = None

    # method that prints the string representation of the specie
    def show(self):

        if self.num != None:
            print 'num = %04d ' % self.num,
            print '%-12s charge = %+-d  type = %+-d abun = %+-e' % (self.str, self.charge, self.type, self.getAbun() ),
        else:
            print 'num = NA   ',   
            print '%-12s charge = NA  type = %+-d abun = NA' % (self.str, self.type),

        print self.comp

    # method that parses the components of species by counting the elements in each specie
    # make sure baseSpec has the species with String in a decreasing order in length
    # (except for the charge for example e-), i.e in this order
    #     longest string names ( CRPHOTON, PHOTON, CRP)
    #     regular elements (Na, Cl,...)
    #     elements with a single letter Upper case (H, F, C...)
    #     species denoted with lower case letters ( e-)
    def getComponents(self, baseSpec):

        specStr = self.str

        for bs in baseSpec:
            # work with one base specie at a time
            nFound = int32(0)
            while True:
                regex = '(' + bs.str + ')' + '([0-9]{0,2})'
                m = re.search( regex, specStr)
                if m:
                    # updating the number of elements found
                    if len(m.string[m.start(2):m.end(2)])==0:
                        nFound = nFound + 1
                    else:
                        nFound = nFound + int32(m.string[m.start(2):m.end(2)])
                        
                    # removing the found element from the string
                    specStr = re.sub( regex, ' ', specStr, count=1)
                else:
                    break

            if nFound >= 1:
                self.comp.append( [bs.str, nFound] ) 
        
        # by now, all the elements should be extracted and counted and only
        # the charge char (if any) should be in the string
        specStr = strip(specStr)
        if specStr == '':
            self.charge = 0
        elif specStr == '+':
            self.charge = 1
        elif specStr == '-':
            self.charge = -1

        # succcessfully processes
        self.init = 1
        
    # returns True if the base species in specList ['XX','YY'] are in self.comp
    #     ex : print net.species['HCN'].hasComponents(['H','N','C'])
    #     ex : print net.species['13C'].hasComponents(['C'])
    def hasBaseSpecies(self, specStrList ):
        
        hasAllComponents = True
        
        cmpntsStrList = [ strng for strng, n in self.comp ]
        
        for specStr in specStrList:
            if specStr not in cmpntsStrList:
                hasAllComponents &= False  
                
        return hasAllComponents
        
    def getAbun(self):
        if self.abun != None:
            return self.abun[0]
        else:
            return None
        
    def setAbun(self, abun):
        self.abun[0] = abun
    
    def getNum(self):
        return self.num
