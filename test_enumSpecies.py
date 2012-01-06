from enumSpecies import *


inName  = '/home/mher/ism/speciesInfo/species.inp'
outName = '/home/mher/outName.pro'

#reading the file
e = enumSpecies(inName)
 
# generating a file which when called returns an IDL structure
# whose tags are the species
e.writeIDL_FILE(outName)

# generating an object whose attributes are the species whose 
# values are the indidices in the file
s = e.genPythonClass(offset=1)
