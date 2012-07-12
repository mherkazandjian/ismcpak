import numpy as np
import pylab as pyl

# NOTE : THIS IS NOT RUNNING IN ECLIPSE...but it is running from a regular 
# python shell

baseDir = '/home/mher/ism/docs/paper02/lineData' 
specStr = 'HCN'
upperStr = '1'
lowerStr = '0'
metallicity = 2
gMech    = -30

fName = ('%s/%s-%s-%s-%.1f%d.npy') % (baseDir,
                                      specStr,
                                      upperStr,
                                      lowerStr,
                                      metallicity,
                                      gMech)
#-----------------------------------------------------------------------
x = np.load( fName )
print x
x = np.log10(x)

mn = -12 #np.min(x)
mx = -1  #np.max(x)
ticks_at = [mn, -10, -8, -6, -4, -2, mx]
print 'min = ', mn, 'max = ', mx

gridLeft = 0.15
gridBotm = 0.1
gridSz   = 0.7

cbarLeft = 0.1
cbarBotm = 0.88
cBarLengt = 0.8
cBarThick = 0.05

fig = pyl.figure(1, figsize = (6,6) )
# plotting the color coded grid
grdAxs = fig.add_axes([gridLeft, gridBotm, gridSz, gridSz])
im11 = grdAxs.imshow( x.transpose(), origin = 'bottom', extent=[ 0, 6, 0, 6], 
                   aspect='auto', interpolation='nearest',
                   vmin = mn, vmax = mx )
grdAxs.set_xlabel("$\log_{10} [ n / cm^{-3} ]$")
grdAxs.set_ylabel("$\log_{10} G_0$")



# creating and plotting the colorbar
cbarAxs = fig.add_axes([cbarLeft, cbarBotm, cBarLengt, cBarThick])
cbarAxs.text(0.3, 1.2, "$\log_{10}$ [ CO (%s-%s) / $erg\,cm^{-3}\,s^{-1}$ ] " %(upperStr, lowerStr))

cbar = pyl.colorbar(im11, cax = cbarAxs, ax = grdAxs, orientation = 'horizontal',
                    format = '%d', ticks = ticks_at)

# plotting the contour lines
##grdAxs.contour( x.transpose(), extent=[ 0, 6, 0, 6] )

pyl.show()

fig.savefig( fName.replace('.npy', '.eps') )
print 'ok'