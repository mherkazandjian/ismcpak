from string import *
from numpy import *
import StringIO
import re
from collections import namedtuple

'''
Created on 21 Dec 2011

this module is used to read and dump the species names into C or IDL
enum data structures.
  
  
     arguments : 
        - inName                : the text file containing the species and their indicies

    language binding : the language to which the binding will be generated
                                   IDL code :
                                     print, s._HCO_P
                                   C,   (not implemented yet)
                                   PYTHON 
                                    If this is set to None, the a class containing
                                    the species as attributes is returned, where
                                    the attribute values are the indicies in the 
                                    input file. A species ( say HCO+ ) index is obtained
                                    as : 
                                       s = enumSpecies(inName)
                                       print s.xHCO_P
    NOTES : first three lines are ignored
@author: mher
'''

class enumSpecies():
    
    def __init__(self, inName, offsetFromZero=None):
        self.inName  = inName     # the text file containing the species and their indicies
        self.buffer  = None        
        
        if offsetFromZero == None:
            self.offsetFromZero = 0
        else:
            self.offsetFromZero = offsetFromZero
        
        self.read()


    # method that read the species info file and stores ints info into tuples
    def read(self):
        
        fileIn  = open(self.inName, 'r')

        buff = []
        i = 0
        for line in fileIn:
            i = i + 1
            if i <= 3:
                continue
            else :
                lineParsed = line.split(' ')
                specName = lineParsed[1]
                specIdx  = int32(lineParsed[0])-self.offsetFromZero
                
                buff.append((specName, specIdx))
           
        self.buffer = buff 
    
    def genVarCodeIDL_or_C_Py(self, specName):
            
        m = re.search('[\+,\-]', specName)
        if m:
            charge = specName[m.start(0):m.end(0)]  
            specName = specName[0:-1]
            if charge == '+' :
                specName = specName + '_P'
            if charge == '-' :
                specName = specName + '_N'
        return specName
               
    def writeIDL_FILE(self, outName, offset=None):
        
        if offset == None:
            offset = 0

        fileOut = open(outName , 'w')

        fileOut.write( "function species \n" )
        fileOut.write( "species = { $ \n" )

        for item in self.buffer:
            varName = self.genVarCodeIDL_or_C_Py(item[0])
            fileOut.write( "_%s : %i ,$\n" % (varName, item[1] + offset) )

        fileOut.write( " _FOO : -1 }\n " )

        fileOut.write( "return, species \n" )
        fileOut.write( "end\n" )

        fileOut.close()

    def genPythonClass(self, offset=None):

        if offset == None:
            offset = 0
          
        attrNames = []
        for item in self.buffer:
            varName = 'x'+self.genVarCodeIDL_or_C_Py(item[0])
            attrNames.append(varName)
        buffTmp = ' '.join(attrNames)
        
        cls = namedtuple('specEnum', attrNames)

        for item in self.buffer:
            varName = 'x'+self.genVarCodeIDL_or_C_Py(item[0])
            setattr(cls, varName, item[1] - offset)

        return cls