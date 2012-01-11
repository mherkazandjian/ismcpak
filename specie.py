from string import *
from numpy import *
import re

""" methods :  self.__init__(specStr, specType=None, charge=None, init=None)
               self.show()
               self.getComponents(baseSpec)

"""

# class definition for a single specie
# --------------------------------------
class specie():
    def __init__(self, specStr, specType=None, charge=None, init=None):
        self.type = specType   # , -1 => do not conisder it as a specie, 0 => basic specie, 1 => composite specie
        self.str  = specStr    # string representation of the specie
        self.num  = -1         # numeric representation of the specie

        if specType == 0:               # a list holding the indicies of the basic species making 
            self.comp  = [[specStr,1]]   # it up and the number number of each sub-specie
        else:
            self.comp = []

        self.charge = charge
        self.abun = None     # abundance of the specie relative to total H nuclei
        self.init = None     # flag which indicates if the specie is inetialized or not

    # method that prints the string representation of the specie
    def show(self):
        print self.str, ' : ', self.comp, ' charge', self.charge, 'type', self.type

    # method that parses the components of species by counting the elements in each specie
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
