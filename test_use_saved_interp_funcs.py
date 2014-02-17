"""
<keywords>
interpolation functions, scipy, ismcpak
</keywords>
<description>
demonstration of using the saved interpolation functions
</description>
"""

import time
import os
import numpy
from numpy import random, array
import pylab
import meshUtils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/sph-db-z-1.0-low-res')

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=False)


F = arxvPDR.load_interp_func(info={'source':'radex'},
                             line='13CO1-0', 
                             quantity='fluxKkms')
 

logn  = 3.0
logG0 = 3.0
logGM = -50.0
Av    = 1.0

data = array([logn, logG0, logGM, Av]).reshape(4,1).T
v = F.get( data )


## testing the perforamce by interpoalting at the same point many times
nPts = 1000
t0 = time.time()
for i in numpy.arange(nPts):
    v = F.get( data )    
dt = time.time() - t0
print 'interpolating single points at a rate of %.2e points per second' % (nPts / dt)


## testing the perforamce by interpoalting at differnet locations in the 
## parameter space
nTrials = 50
nPtsPerTrial = 1
logn_rng = [0.0, 3.0]
logG0_rng = [3.0, 5.0]
logGM_rng = [-50.0, -20.0]
Av_rng = [1.0, 5.0]

t0 = time.time()
for i in numpy.arange(nTrials):
    
    logn  = random.rand(nPtsPerTrial)*\
               (logn_rng[1] - logn_rng[0]) + logn_rng[0]
    logG0 = random.rand(nPtsPerTrial)*\
               (logG0_rng[1] - logG0_rng[0]) + logG0_rng[0]
    logGM = random.rand(nPtsPerTrial)*\
               (logGM_rng[1] - logGM_rng[0]) + logGM_rng[0]
    Av    = random.rand(nPtsPerTrial)*\
               (Av_rng[1] - Av_rng[0]) + Av_rng[0]

    data = array([logn, logG0, logGM, Av]).reshape(4,nPtsPerTrial).T
    print data
    v = F.get( data )
    
dt = time.time() - t0
print 'interpolating points at a rate of ' +\
        '%.2e points per second' % ( (nPtsPerTrial * nTrials) / dt)


## plotting the output for an interpolated set of points while varying
## the input coordinates in one dimension

nPts = 20
logn = numpy.linspace(0.0, 6.0, nPts)
logG0 = numpy.ones(nPts)*3.0
logGM = numpy.ones(nPts)*(-25.0)
Av = numpy.ones(nPts)*(1.0)

# the quantity that will be plotted on the x-axis
x = logn

data = numpy.vstack((logn, logG0, logGM, Av)).T

t0 = time.time()
v = F.get( data )    
dt = time.time() - t0

print 'interpolating points at a rate of ' +\
        '%.2e points per second' % ( nPts / dt)

pylab.plot(x, v)
pylab.plot(x, v, 'o')
pylab.xlim(x.min(), x.max())
