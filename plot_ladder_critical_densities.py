#this is a test file for leidenLambda.py.
import sys, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from leidenLambda import molData 
import numpy
import re
from ismUtils import planckOccupation as ng

#reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'

#species to be loaded
species = {
           'CO'    :  {'partner': 'o-H2', 'obj':None, 'nc':None},
           'HCN'   :  {'partner': 'H2'  , 'obj':None, 'nc':None},
           'HNC'   :  {'partner': 'H2'  , 'obj':None, 'nc':None},
           'HCO+'  :  {'partner': 'H2'  , 'obj':None, 'nc':None},
           'CS'    :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           '13C16O':  {'partner': 'p-H2', 'obj':None, 'nc':None},
           'CN'    :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           'SiO'   :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           'OH'    :  {'partner': 'p-H2', 'obj':None, 'nc':None},
#           'SO2'   :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           'SO'    :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           'N2H+'  :  {'partner': 'H2'  , 'obj':None, 'nc':None},
#           'HCS+'  :  {'partner': 'H2'  , 'obj':None, 'nc':None},
           }

#transitions indicies for which the critical densities will be computed
transInds = numpy.arange(0, 15)

#kinetic temperature at which the critical densities will be computed
T = 100.0

#--------------------------------------------------------------------------------------------
#loading the species
reader = molData.reader(dirPath = lambdaPath, species = species.keys())

for specStr in species:
    species[specStr]['obj'] = reader.get_specie(specStr=specStr) 
    species[specStr]['nc']  = numpy.zeros(transInds.size) 

#getting the critical densities
for u, i in enumerate(transInds):
    for specStr in species:  
        species[specStr]['nc'][i] = species[specStr]['obj'].critical_density(trans=i, T_kin=T, collider=species[specStr]['partner'])

plts = []
lbls = []
#plotting the ncrical densities
for specStr in species:
    plt, = pylab.semilogy(transInds, species[specStr]['nc'])
    lbl = '%s_%s' % (specStr , species[specStr]['partner'])
    print lbl
    plts.append(plt)
    lbls.append(lbl)
    
pylab.xticks(transInds)
pylab.legend(plts, lbls)

pylab.show()