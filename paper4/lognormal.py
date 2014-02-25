"""
<keywords>
lognormal, sample, distribution, PDF
</keywords>
<description>
play around with lognormal distributions
</description>
"""

from mylib.utils import histogram
import numpy
from numpy import log, log10, exp, sqrt, pi

import pylab

nPts = 1000000

#mean, disp = 10.0**-0.5, 10.0**2.1
#mean, disp = 10.0**2.2, 10.0**1.4   # mean and dispersion


#mu = 1.0 #10.0**-0.5
#sigma = 3.1 #10.0**2.1

#adasdad
#---------------------------------------------

################
# wada 2001 data
################
l10xdata = [0.0 ,  1.0,  2.0, 3.0,  4.0 ,  5.0]
l10ydata = [-1.5, -1.8, -2.3, -2.8, -3.4, -4.4]
################

var = disp**2

mu = log( (mean**2) / sqrt(var + mean**2) )
sigma = sqrt( log(1.0 + var / mean**2) )

print 'mu, sigma = ', mu, sigma

def lognorm_pdf(mu, sigma, x):
     
    fac = 1.0 / (sqrt(2 * pi) * x * sigma )
     
    return fac * exp( -0.5 * ( (log(x) - mu)/sigma )**2)
        
#x_sample = numpy.random.lognormal(mu, sigma, (1, nPts))

## re-constructing the distribution function from the sampled points
## in log10(x) to recover the log-normal distribution
#lx = log(x_sample)
#hist = histogram.hist_nd(lx, nbins=1000, loc=True)

## number of points in unit inteval (d(log(x))
#lxbin = hist.f.cntrd
#f = hist.f / nPts / hist.szBins

## the theoretical distribution
x_exact = 10.0**numpy.linspace(-20, 20, 1000)
pdf_exact = lognorm_pdf(mu, sigma, x_exact)

## plotting
#pylab.plot(lxbin, log(f),'o-')

#
pylab.plot(log10(x_exact), log10(pdf_exact), 'r--')

#plotting the wada data
#pylab.plot(l10xdata, l10ydata, '+')

pylab.xlim(-3, 7)
pylab.ylim(-7, -1)
pylab.show()