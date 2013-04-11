from leidenLambda import molData 
import numpy
import os
HOME = os.environ['HOME']
#------------------------------------------------

#reading the whole database of line info of species from LAMBDA
lambdaPath = HOME + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
reader = molData.reader(dirPath = lambdaPath, species = 'HCO+')

# getting the data from a particular file
# data = reader.get_specie(specStr = 'CO', inPath = 'xpol_new.dat')
# print data[0]['path']  #file path of the matched specie
# specInfo = data[0] 
# molData.critical_density(specInfo =specInfo, upper=1, lower=0, T_kin = 20.0, collider = 'p-H2')

#### add a class called specie with methods 
#       ciritical density
#       transition info give upper /lowe
#       in/out rate transition 
#        ...

#plot 
#print CO_13['transRad'][4][3], CO['transRad'][4][3]


#reader = molData.reader(dirPath = lambdaPath, species = ['CO', '13CO', 'HCO+', 
#                                                         'HCN', 'HNC', 'CS', 'CN'])

#hcop = reader.get_specie(specStr='HCO+', inPath='hco+.dat')
#cs   = reader.get_specie(specStr='CS', inPath='cs.dat')
#cn   = reader.get_specie(specStr='CN', inPath='cn.dat')
#co   = reader.get_specie(specStr='CO', inPath='co.dat')


kb_cgs = numpy.float64(1.38e-16) #boltzmann constant in cgs

"""
idx = 0
specInfo = reader.get_specie(specStr='CO', inPath='co.dat')
transRad = specInfo['transRad'][idx]

print transRad['u'], transRad['l'] 
print transRad['nu']
"""