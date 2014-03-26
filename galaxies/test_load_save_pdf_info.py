from galaxies import fi_utils
from pylab import *
from scipy import *

gas = fi_utils.gas_set(10)

## loading the sampling info from the saved file
filename = '/home/mher/ism/runs/galaxies/coset2run4/coset-9-sol-ext-100/firun/weights_func.000004.npz'
info = gas.load_weights_function(filename)        


x, y = info['x_fit'], info['y_fit']

plot(log(x), y)

nPDF = fi_utils.density_distribution('lognormal', mu=10.0**0.13, sigma=10.0**0.93)
xs = 10.0**linspace(-3, 6, 100)
ys = nPDF.pdf(xs)
plot(log(xs), ys, 'r-')


nPDF = fi_utils.density_distribution('lognormal', mu=10.0**0.13, sigma=10.0**1.8)
xs = 10.0**linspace(-3, 6, 100)
ys = nPDF.pdf(xs)
plot(log(xs), ys, 'r:')


nPDF = fi_utils.density_distribution('lognormal', mu=10.0**0.13, sigma=10.0**2.7)
xs = 10.0**linspace(-3, 6, 100)
ys = nPDF.pdf(xs)
plot(log(xs), ys, 'r--')


show()