import numpy as np
import pylab as pyl

# NOTE : THIS IS NOT RUNNING IN ECLIPSE...but it is running from a regular 
# python shell

x = np.load('/home/mher/ism/docs/paper02/lineData/CO-1-0-2.0-30.npy')
print x

im11 = pyl.imshow( x.transpose(), origin = 'bottom', extent=[ 0, 6, 0, 6], 
                   aspect='auto', interpolation='nearest')

cbar11 = pyl.colorbar(im11, ax=pyl.gca(), orientation = 'horizontal')
pyl.contour( x.transpose(), extent=[ 0, 6, 0, 6] )
pyl.show()
print 'ok'