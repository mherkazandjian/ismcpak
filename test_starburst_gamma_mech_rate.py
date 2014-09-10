import ismUtils
import numpy
from scipy import integrate

'''
#test1 : if y = sqrt(x), the integral from 1 to 4 is 14/3
x = numpy.linspace(1.0, 4.0, 10000.0)
y = numpy.sqrt(x)
v_integrand = integrate.simps(y, x)
err = 1.0 - v_integrand / (14.0/3.0)
print 'error = %e' % err
'''

print ismUtils.startburst_gamma_mech_rate(SFR=10.0, n_PDR=15.0, d_PDR=0.1, d_SB=100.0, E_SN=1e51, eta=0.1)