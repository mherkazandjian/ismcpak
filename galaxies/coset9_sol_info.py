import fi_utils
import numpy
 
# returns the probability density given the density in cm-3 (not in log scale)
# mu = 0.31 = log(1.36 cm^-3)     mu(in log10) = log10(1.36) = 0.134
# sigma = 2.16 = log(13.6)        sigma(in log10) = log10(13.6) = 1.13
nPDF_sph = fi_utils.density_distribution('lognormal', 
                                         mu=numpy.exp(0.31), 
                                         sigma=numpy.exp(2.16))

# returns the probability density given the density in cm-3 (not in log scale)
nPDF_sph_test = fi_utils.density_distribution('lognormal', mu=numpy.exp(0.31), sigma=numpy.exp(6.16))

# returns the radius of a cloud (in kpc) given the density in cm-3 
r_sph_kpc = lambda n: (30.0/1000.0)*n**(-1.0/3.0)