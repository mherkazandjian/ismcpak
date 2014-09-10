from despotic import cloud
import numpy
import pylab

specStr = 'CO'

mycloud = cloud()
mycloud.nH = 1.0e5          #gas density
mycloud.colDen = 2.0e22     #cloud column density
mycloud.sigmaNT = 4.0e5     #non-thermal velocity despersion 

mycloud.comp.xoH2 = 0.1     #ortho-H2 composition, xoH2 molecule per H nucleus
mycloud.comp.xpH2 = 0.4     #para-H2 composition, xpH2 molecule per H nucleus

mycloud.Tg = None           #cloud gas kinetic temperature
mycloud.Td = 0.0            #cloud dust temperature

mycloud.addEmitter(specStr, 1.0e-7)  #abudnace of the emitting species per H nucleus

mycloud.Tg = 40.0 
lines1 = mycloud.lineLum(specStr)

mycloud.Tg = 400.0 
lines2 = mycloud.lineLum(specStr)


u1 = [l['upper'] for l in lines1]
inten1 = numpy.array([l['intIntensity'] for l in lines1]) #intensity after subtracting the CMB contributuin

u2 = [l['upper'] for l in lines2]
inten2 = numpy.array([l['intIntensity'] for l in lines2])

pylab.plot(u1, inten1, '-o')
pylab.plot(u2, inten2, '-x')
pylab.xlabel('upper')
pylab.ylabel(r'erg.cm$^{-2}$.s$^{-1}$')
pylab.show()

"""" plotting the solutions of despotic using radex plotting utilities"""
"""
r = radex('','')
r.set_attributes_from_despotic(specStr, mycloud, lines1)
r.plotModel()
"""