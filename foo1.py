# the models of paper1 (table-1)
#---------------------------------------------------------------------------------------------------
from numpy import *
from numpy.random import *
import pylab as pyl

from amuse.community.pdr import interface
from mesh import *
from chemicalNetwork import *
from enumSpecies import *
from ismUtils import *
import time

dx   = 1.0      # log10 density
xMin = 5.0
xMax = 6.1

dy   = 1.0     # log10 G0
yMin = 5.0  
yMax = 6.1

# factor of surface heating to be added as mechanical heating
z = [0.01, 0.05]

# generating the parameter space 
xg = arange(xMin, xMax, dx)
yg = arange(yMin, yMax, dy)
zg = arange(len(z))
x, y, z = mgrid[ xMin:(xMax+1e-10):dx, yMin:(yMax+1e-10):dy, 0:len(z):1]

print 'x'
print x
print 'y'
print y
print 'z'
print z
print x.flatten()
print y.flatten()
print z.flatten()
print '---------'
print z[:][:][:]
#print z.flatten()
