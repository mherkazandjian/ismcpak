import numpy as np
import pylab as pyl

# NOTE : THIS IS NOT RUNNING IN ECLIPSE...but it is running from a regular 
# python shell

x = np.load('/home/mher/ism/docs/paper02/lineData/CO-1-0-2.0-30.npy')
x = np.log10(x)
print x

mn = -12 #np.min(x)
mx = -4  #np.max(x)
ticks_at = [mn, -7, mx]
print 'min = ', mn, 'max = ', mx

gridLeft = 0.1
gridBotm = 0.1
gridSz   = 0.5

cbarLeft = 0.1
cbarBotm = 0.8
cBarLengt = 0.8
cBarThick = 0.1

fig = pyl.figure(1, figsize = (6,6) )
# plotting the color coded grid
grdAxs = fig.add_axes([gridLeft, gridBotm, gridSz, gridSz])
im11 = grdAxs.imshow( x.transpose(), origin = 'bottom', extent=[ 0, 6, 0, 6], 
                   aspect='auto', interpolation='nearest',
                   vmin = mn, vmax = mx )

# creating and plotting the colorbar
cbarAxs = fig.add_axes([cbarLeft, cbarBotm, cBarLengt, cBarThick])
"""
for tick in cbarAxs.xaxis.get_major_ticks():
    tick.label1On = False
    tick.label2On = False
for tick in cbarAxs.yaxis.get_major_ticks():
    tick.label1On = False
    tick.label2On = False
"""

cbar = pyl.colorbar(im11, cax = cbarAxs, ax = grdAxs, orientation = 'horizontal',
                    format = '%d', ticks = ticks_at)

#cbar.ax.set_yticklabels([str(mn), '0', str(mx)])

# plotting the contour lines
##grdAxs.contour( x.transpose(), extent=[ 0, 6, 0, 6] )

pyl.show()

fig.savefig('/home/mher/1.eps')
print 'ok'