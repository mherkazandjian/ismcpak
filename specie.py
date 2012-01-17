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
    def __init__(self, specStr, specType=None, charge=None, init=None):
        self.type = specType   # -2 => ignore specie
                               # -1 => do not conisder it as a specie
                               #  0 => basic specie
                               # 1 => composite specie
        self.str  = specStr    # string representation of the specie
        self.num  = None         # numeric representation of the specie

        if specType == 0:               # a list holding the indicies of the basic species making 
            self.comp  = [[specStr,1]]   # it up and the number number of each sub-specie
        else:
            self.comp = []

        self.charge = charge
        self.abun = None     # abundance of the specie relative to total H nuclei
        self.init = None     # flag which indicates if the specie is inetialized or not
        self.active = None

    # method that prints the string representation of the specie
    def show(self):

        if self.num != None:
            print '%04d ' % self.num,
        else:
            print 'None ',  
        
        print '%-12s charge = %+-d  type = %+-d abun = %+-e' % (self.str, self.charge, self.type, self.getAbun() )
        """
        #printStr +=  
        #
        #    numPrnt = float(nan)
        #else:
        #    numPrnt = self.num
            
        #print self.num, self.str, self.charge, self.type, self.abun
        print '%03d %-12s charge = %+-d  type = %d abun = %e' % (self.num, self.str, self.charge, self.type, self.abun[0]) 
    #, self.comp
    # method that parses the components of species by counting the elements in each specie
    # make sure baseSpec has the species with String in a decreasing order in length
    # (except for the charge for example e-), i.e in this order
    #     longest string names ( CRPHOTON, PHOTON, CRP)
    #     regular elements (Na, Cl,...)
    #     elements with a single letter Upper case (H, F, C...)
    #     species denoted with lower case letters ( e-)
    """
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
        
    def getAbun(self):
        return self.abun[0]
    def setAbun(self, abun):
        self.abun[0] = abun
    
    def getNum(self):
        return self.num
