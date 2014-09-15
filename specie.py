from string import *
from numpy import *
import re

# class definition for a single specie
# --------------------------------------
class specie():
    """A class for defining species. If the species is a baseSpecie (specType=0), the all the keywords
    (except comp) are mandatory.
    """
    def __init__(self, specStr, specType=None, charge=None, init=None, comp=None, mass = None):
        
        self.type = None
        """   
             + -2 => ignore specie
             + -1 => do not conisder it as a specie
             + 0 => basic specie
             + 1 => composite specie
             
        """

        #checking keywords and arguments        
        if specStr  == None: raise ValueError('specStr is a mandatory argument.')
        
        #if it is a base sepcie, need to provide the charge, the mass
        if specType == 0:
            
            if charge  == None: raise ValueError('charge is a mandatory keyword when specType = 0')
            if init    == None: raise ValueError('init is a mandatory keyword when specType = 0')
            if mass    == None: raise ValueError('mass is a mandatory keyword when specType = 0')
            
            self.type = 0
            
        self.charge = charge #: the charge of the species (integer) 
        self._abun  = None   #: Abundance of the specie relative to total H nuclei. It is either None or an ndarray of shape (1,)
        self.init   = None   #: flag which indicates if the specie is initialized or not
        
        self.str     = specStr #: string representation of the specie
        self.num     = None    #: the index of the specie (a numeric representation of the species mapping it to an array index for example)
        self.comp    = None    
        """"a list of [specStr, N, baseSpecObj] for each compoenent making up the specie
        where N is the number of times the base species occures (basically baseSpecObj.str = specStr)
        """
        self.mass = mass
        
        #setting up the components list
        if comp == None:
            
            if specType == 0: #if no components have been passed and the specie is a base species, the components are set to [self.specStr, 1, self]         
                self.comp  = [[specStr, 1, self]]  #otherwise, to [] to be parsed maybe later.
            else:
                self.comp = []
        else:
            
            if len(comp[0]) != 3:
                raise ValueError('input component must be a list of lists eahc of length 3 [specStr, count, object|string]')
            else:
                
                if comp[0][2] == 'self':
                    comp[0][2] = self
                    
            self.comp = comp

        #setting the mass of the specie
        if mass != None:
            self.mass = float64(mass)

                        
    def show(self):
        """method that prints the string representation of the specie."""

        chargeStr = '%+-d' % self.charge if self.charge != None else 'NA'
        numStr    = '%04d' % self.num if self.num != None else '%-4s' % 'NA'
        typeStr   = '%+-d' % self.type if self.type != None else 'NA'
        abunStr   = '%+-.2e' % self.abun() if self.abun() != None else '%-9s' % 'NA'
        mass      = '%-10.3f' % self.mass if self.mass != None else '%-4s' % 'NA'
        
        print 'num = %s %-10s charge = %s type = %s abun = %s mass = %s amu | elements = ' % (numStr, self.str, chargeStr, typeStr, abunStr, mass),
        print self.comp
        
    def get_components(self, baseSpecs, get_mass = None):
        """Method that parses the components of species by counting the elements in each specie
           make sure baseSpecs has the species with String in a decreasing order in length
           (except for the charge for example e-), i.e in this order longest string names ( CRPHOTON, PHOTON, CRP)
           regular elements (Na, Cl,...) elements with a single letter Upper case (H, F, C...)
           species denoted with lower case letters ( e-).
        """

        specStr = self.str  #specStr gets a copy (this is not assignment by reference)
        
        self.comp = []
        
        for baseSpecStr, baseSpec in baseSpecs.items():
            
            # work with one base specie at a time
            nFound = int32(0)
            while True:
                
                regex = '(' + baseSpecStr + ')' + '([0-9]{0,2})' #look for AAANN, i.e text followed by numbers
                m = re.search(regex, specStr) #if not found, m = None

                if m:
                    # updating the number of elements found (by looking at the number after the species string)
                    if len(m.string[m.start(2):m.end(2)])==0: #component string not followed by a number => there is one compoenet of it, for example HCN and we are looking for H
                        nFound += 1 
                    else:
                        nFound += int32(m.string[m.start(2):m.end(2)]) #component string followed by a number => there is more than one compoenet of it, for example the '3' in H3CN and we are looking for H
                        
                    # removing the found element from the string
                    specStr = re.sub(regex, ' ', specStr, count=1)
                else:
                    break

            if nFound >= 1:
                self.comp.append( [baseSpecStr, nFound, baseSpec] ) 
        
        # by now, all the elements should be extracted and counted and only
        # the charge char (if any) should be in the string
        specStr = strip(specStr)
        if specStr == '':
            self.charge = 0
        elif specStr == '+':
            self.charge = 1
            specStr = specStr.replace('+','')
        elif specStr == '-':
            self.charge = -1
            specStr = specStr.replace('-','')

        #by now, specStr should be a blank string, otherwise, it contains components that
        #are not in baseSpecies (or something went wrong)
        specStr = specStr.strip() 
        if len(specStr) != 0:
            raise ValueError("""After parsing self.str = %s the following part '%s' was
             not traced back to a base species. Please make sure that all the base 
             species are provided.""" % (self.str, specStr))        

        if get_mass == True:
            self.compute_mass()
            
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
    
    def compute_mass(self):
        """computes the weight of the speicies. Nothing is done for base species because the mass
        is a requirement when the species is instantiated.
        """
        
        if self.type == 0:
            pass
        else:
            
            #computing the mass of the non base specie
            if len(self.comp) == 0:
                raise ValueError('cannot compute masses from, compoenents must be computed before calling this method') 
            
            mass = float64(0.0)
            
            for compStr, count, baseSpec in self.comp:
                mass += baseSpec.mass * count
                
            self.mass = mass
            
    def copy(self):
        """returns a new object whose attribute values are a copy of self"""
        
        #setting up a new object which will hold the copied attributes
        specNew = specie(self.str)
        
        #copying basic attributes
        specNew.charge = self.charge        
        specNew.init  = self.init
        specNew.num = self.num
        specNew.mass = self.mass
        specNew._abun = self._abun
        specNew.type = self.type
        
        #copying the components attribute (if any)
        if len(self.comp) != None:
            
            compNew = []
            
            for comp in self.comp:
                compNew.append( [comp[0], comp[1], comp[2]] )
                
            specNew.comp = compNew
            
        return specNew