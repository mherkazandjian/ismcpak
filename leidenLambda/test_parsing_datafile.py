from leidenLambda import molData 
import numpy as np
import re
from ismUtils import planckOccupation as ng

#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath, specie = 'CO')

# getting the data from a particular file
# data = reader.get_specie(specStr = 'CO', inPath = 'xpol_new.dat')
# print data[0]['path']  #file path of the matched specie
# specInfo = data[0] 
# molData.critical_density(specInfo =specInfo, upper=1, lower=0, T_kin = 20.0, collider = 'p-H2')

for spec in reader.speciesInfo:
    print '----------------------------------------------'
    for partner in spec['transColl']['partnersList']:   
        print '          %s' % partner 

#### add a class called specie with methods 
#       ciritical density
#       transition info give upper /lowe
#       in/out rate transition 
#        ...

#plot 
#print CO_13['transRad'][4][3], CO['transRad'][4][3]
