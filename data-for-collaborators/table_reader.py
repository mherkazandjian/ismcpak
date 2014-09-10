#description: sample reader of the table of PDR models and their emissions
#
#to get the emissions for a model with n = 1e3, G0 = 1e3, Z=1, Av=10
#           
#   model_data = data[ (data['n'] == 1e3) & (data['G0'] == 1e3) & (data['Av'] == 5.0) & (data['Z']==1) & (data['alpha'] <= 0.0) ]
#   
#   CO10 = model_data['CO_1-0']
#   HCN10 = model_data['HCN_1-0']
#   #and so on...
#
#-------------------------------------------------------------------------------------------
import numpy
import pylab

fileName = 'test.txt'


descr = [
         ('idx'    , 'i'),
         ('n'      , 'f'),
         ('G0'     , 'f'),
         ('gm'     , 'f'),
         ('alpha'  , 'f'),
         ('Z'      , 'f'),
         ('Av'     , 'f'),
         ('CO_1-0'  , 'f'),
         ('HNC_1-0' , 'f'),
         ('HCN_1-0' , 'f'),
         ('HCO+_1-0', 'f'),
         ('CS_2-1'  , 'f')
        ]
 
data = numpy.loadtxt(fileName, dtype=descr)

subset = data[ (data['n'] == 1e3) & (data['G0'] == 1e3) & (data['Av'] == 5.0) ]
x, y = subset['gm'], subset['CS_2-1'] / subset['CO_1-0'] 
inds = numpy.argsort(x)
plt1, = pylab.loglog( x[inds], y[inds], 'r' )


subset = data[ (data['n'] == 1e3) & (data['G0'] == 1e3) & (data['Av'] == 10.0) ]
x, y = subset['gm'], subset['CS_2-1'] / subset['CO_1-0'] 
inds = numpy.argsort(x)
plt2, = pylab.loglog( x[inds], y[inds], 'b' )

subset = data[ (data['n'] == 1e3) & (data['G0'] == 1e3) & (data['Av'] == 30.0) ]
x, y = subset['gm'], subset['CS_2-1'] / subset['CO_1-0'] 
inds = numpy.argsort(x)
plt3, = pylab.loglog( x[inds], y[inds], 'k' )

pylab.legend([plt1, plt2, plt3], ['Av = 5', 'Av = 10', 'Av = 30'])
pylab.xlabel(r'mechanical heating rate (erg/cm^3/s)')
pylab.ylabel(r'line ratio (CS(2-1)/CO(1-0)')
pylab.show()

