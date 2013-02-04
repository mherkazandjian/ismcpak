from leidenLambda import molData 
import numpy as np
import re
from ismUtils import planckOccupation as ng

#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath, specie = 'CO')

# getting the data from a particular file
#CO = reader.get_specie(specStr = 'CO', inPath = 'xpol_new.dat')

for spec in reader.speciesInfo:
    print '----------------------------------------------'
    for partner in spec['transColl']['partnersList']:   
        print '          %s' % partner 
