import numpy as np
import pylab as pyl
from removeAxesLabels import *

# NOTE : THIS IS NOT RUNNING IN ECLIPSE...but it is running from a regular 
# python shell


nx = 3
ny = 1

specStr     = 'HNC'
metallicity =  2.0
gMech       = -22

baseDir = '/home/mher/ism/docs/paper02/lineData'

title = "$\log_{10}$ [ %s line intensities / $erg\,cm^{-3}\,s^{-1}$ ] " % specStr 

filesNfo = (  
             (baseDir, specStr, '01', '00', metallicity, gMech), # [0][0]
             (baseDir, specStr, '02', '01', metallicity, gMech), # [0][1]
             (baseDir, specStr, '04', '03', metallicity, gMech), # [0][2]
           )
         
gridsLeft  = 0.06
gridsBotm  = 0.14
gridsRight = 0.98
gridsTop   = 0.80
gridsWspsc = 0.05
gridsHspsc = 0.05

cbarLeft = 0.2
cbarBotm = 0.88
cBarLengt = 0.6
cBarThick = 0.03


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

outFname = '/home/mher/ism/docs/paper02/figs/%s-%d-%.1f.eps' % (specStr, gMech, metallicity)

fNameFormat = ('%s/%s-%s-%s-%.1f%d.npy')

fig, axsGrids = pyl.subplots( ny, nx, sharex = True, sharey = True, figsize = (8, 3.5) )
pyl.subplots_adjust( left = gridsLeft , bottom = gridsBotm, right = gridsRight, top = gridsTop, 
                     wspace = gridsWspsc, hspace = gridsHspsc)

#-----------------------------------------------------------------------

mn = -12 #np.min(x)
mx = -1  #np.max(x)
ticks_at = [mn, -10, -8, -6, -4, -2, mx]
print 'min = ', mn, 'max = ', mx

for i, nfo in enumerate(filesNfo):        
        print i, nfo
    
        fName = fNameFormat % nfo
        print fName
        
        x = np.load( fName )
        #print x
        x = np.log10(x)

        print '-------------------'
        axs = axsGrids[i]
        axs.text(4.5, 5, ('%s-%s') % (nfo[2],nfo[3]) )

        im  = axs.imshow( x.transpose(), origin = 'bottom', extent=[ 0, 6, 0, 6], 
                             aspect='auto', interpolation='nearest',
                             vmin = mn, vmax = mx )

# creating and plotting the colorbar

cbarAxs = fig.add_axes([cbarLeft, cbarBotm, cBarLengt, cBarThick])
cbar = pyl.colorbar(im, cax = cbarAxs, ax = axs, orientation = 'horizontal',
                    format = '%d', ticks = ticks_at)
#removeAllLabels(cbarAxs)
cbarAxs.text(0.2, 1.4, title)

# final touches on the plots
for axs in axsGrids:
    axs.set_xlabel("$\log_{10} [ n / cm^{-3} ]$")
axsGrids[0].set_ylabel("$\log_{10} G_0$")

removeSubplotLabels( nx, ny, axsGrids)

pyl.show()

fig.savefig( outFname )
print 'ok'
