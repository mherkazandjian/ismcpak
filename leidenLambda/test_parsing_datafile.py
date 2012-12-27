from leidenLambda import molData 
import numpy as np
import re
from ismUtils import planckOccupation as ng

#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath)

#selecting the one holding the info for p-NH3
"""
for specDict in reader.speciesInfo:
    if 'p-NH3' in specDict['specStr']:
        pNH3 = specDict
        print 'found p-NH3 info, filePath : %s' % pNH3['path']
        break
"""