from string import *
from numpy import *
import re

# class definition for a single specie
# --------------------------------------
class specie():
    """A class for defining species.
    """
    def __init__(self, specStr, specType=None, charge=None, init=None, comp=None):
        self.type = specType
        """   
             + -2 => ignore specie
             + -1 => do not conisder it as a specie
             + 0 => basic specie
             + 1 => composite specie
             
        """
        
        self.str  = specStr    #: string representation of the specie
        self.num  = None       #: numeric representation of the specie

        if comp == None:
            if specType == 0:               #: a list holding the indicies of the basic species making 
                self.comp  = [[specStr,1]]  #: it up and the number number of each sub-specie
            else:
                self.comp = []
        else:
            self.comp = comp

        self.charge = charge
        self._abun = None
        """
            Abundance of the specie relative to total H nuclei. It is either None or an ndarray of shape (1,)
        """
        self.init = None     #: flag which indicates if the specie is inetialized or not
        self.active = None

    def show(self):
        """method that prints the string representation of the specie."""

        chargeStr = '%+-d' % self.charge if self.charge != None else 'NA'
        numStr    = '%04d' % self.num if self.num != None else '%-4s' % 'NA'
        typeStr   = '%+-d' % self.type if self.type != None else 'NA'
        abunStr   = '%+-.2e' % self.abun() if self.abun() != None else '%-9s' % 'NA'
        
        print 'num = %s %-10s charge = %s type = %s abun = %s elements = ' % (numStr, self.str, chargeStr, typeStr, abunStr),
        print self.comp
        
    def get_components(self, baseSpec):
        """Method that parses the components of species by counting the elements in each specie
           make sure baseSpec has the species with String in a decreasing order in length
           (except for the charge for example e-), i.e in this order longest string names ( CRPHOTON, PHOTON, CRP)
           regular elements (Na, Cl,...) elements with a single letter Upper case (H, F, C...)
           species denoted with lower case letters ( e-).
        """

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
        
    def hasBaseSpecies(self, specStrList ):
        """Returns True if the base species in specList ['XX','YY'] are in self.comp
           ex : 
           
              print net.species['HCN'].hasComponents(['H','N','C']) 
              print net.species['13C'].hasComponents(['C'])
              
        """
        
        hasAllComponents = True
        
        cmpntsStrList = [ strng for strng, n in self.comp ]
        
        for specStr in specStrList:
            if specStr not in cmpntsStrList:
                hasAllComponents &= False  
                
        return hasAllComponents
        
    def abun(self, *args):
        """This works both ways to set and to get the value of self._abun. When passed
           without an argment it returns self._abun. When passed with an argument it
           sets the value of the argument to self._abun[0].
           
           .. code-block:: python
           
               x = spec.abun()  # return spec._abun[0]
               spec.abun(0.53)  # set spec._abun[0] to 0.53
               
           .. warning:: It is highly recommended to use the specie.abun(x) method to set the abun
              of a species. Since it is easy to get confused and set the abundance via specie._abun = x
              which would break the mapping with net.abun. The way to change the abundance without
              breaking the mapping would be to use specie.abun[0] = x, since specie._abun is an ndarray
              of shape (1,) and setting the value via specie._abun = x would assign a new object instead
              of a new value.

            .. todo:: see how to include the documentation of _abun with sphinx
        """

        if len(args) == 0:        
            if self._abun != None:
                return self._abun[0]
            else:
                return None
        else:
            self._abun[0] = args[0]
        
    def get_num(self):
        return self.num