from scipy import stats
from numpy import *
from pylab import *

no = 10000.0
Ns = 500000.0

rv = stats.norm(loc=0.314, scale=2.16)

############
r_func = lambda n_o, N_s, n_n: ((1.0/float(N_s))*(n_o / n_n))**(1.0/3.0)
############
n_new = array([])

i = 0
while n_new.size < Ns:
    
    sample = exp(rv.rvs(100000))

    n_new = hstack((n_new, sample[ sample > no]))
    
    print i, n_new.size
    i += 1
    
print 'sample size = ', n_new.size
print 'keeping only %d of them' % Ns

print 'density of the original particle is ', no
print 'number of sampled points picked from the big sample = ', n_new.size
#area_new = ((1.0/float(Ns))*(no / n_new))**(2.0/3.0)

r_new = ((1.0/float(Ns))*(no / n_new))**(1.0/3.0)
area_new = (r_new**2.0) * rv.pdf(log(n_new))
#area_new = (((rv.pdf(log(n_new)))*(no / n_new))**(2.0/3.0))

print area_new.sum()

#######################











