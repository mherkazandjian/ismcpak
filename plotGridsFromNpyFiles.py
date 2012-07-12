import numpy as np
import pylab as pyl
from removeAxesLabels import *

# NOTE : THIS IS NOT RUNNING IN ECLIPSE...but it is running from a regular 
# python shell


nx = 3
ny = 2

outFname = '/home/mher/ism/docs/paper02/figs/CO-noMech.eps'
baseDir = '/home/mher/ism/docs/paper02/lineData'

title = "$\log_{10}$ [ CO line intensities / $erg\,cm^{-3}\,s^{-1}$ ] " 

filesNfo = (
           # row 1 (from top)
            (  
             (baseDir, 'CO', '1', '0', 2, -30), # [0][0]
             (baseDir, 'CO', '2', '1', 2, -30), # [0][1]
             (baseDir, 'CO', '3', '2', 2, -30), # [0][2]
            ),
           # row 2 (from top)
            (  
             (baseDir, 'CO', '6' , '5' , 2, -30), # [1][0]
             (baseDir, 'CO', '10', '9' , 2, -30), # [1][2]
             (baseDir, 'CO', '16', '15', 2, -30), # [1][3]
            )
           )

gridsLeft  = 0.06
gridsBotm  = 0.1
gridsRight = 0.98
gridsTop   = 0.85
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

fNameFormat = ('%s/%s-%s-%s-%.1f%d.npy')

fig, axsGrids = pyl.subplots( ny, nx, sharex = True, sharey = True, figsize = (8, 6) )
pyl.subplots_adjust( left = gridsLeft , bottom = gridsBotm, right = gridsRight, top = gridsTop, 
                     wspace = gridsWspsc, hspace = gridsHspsc)

baseDir = '/home/mher/ism/docs/paper02/lineData'
specStr = 'CO'
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

mn = -12 #np.min(x)
mx = -1  #np.max(x)
ticks_at = [mn, -10, -8, -6, -4, -2, mx]
print 'min = ', mn, 'max = ', mx

for i, nfoTpl in enumerate(filesNfo):
    for j, nfo in enumerate(nfoTpl):
        print i, j, nfo
        
        fName = fNameFormat % nfo
        print fName
        
        x = np.load( fName )
        #print x
        x = np.log10(x)

        print '-------------------'
        axs = axsGrids[i,j]
        axs.text(4.5, 5, ('%s-%s') % (nfo[2],nfo[3]) )

        im  = axs.imshow( x.transpose(), origin = 'bottom', extent=[ 0, 6, 0, 6], 
                             aspect='auto', interpolation='nearest',
                             vmin = mn, vmax = mx )


# creating and plotting the colorbar

cbarAxs = fig.add_axes([cbarLeft, cbarBotm, cBarLengt, cBarThick])
cbar = pyl.colorbar(im, cax = cbarAxs, ax = axs, orientation = 'horizontal',
                    format = '%d', ticks = ticks_at)
removeAllLabels(cbarAxs)
cbarAxs.text(0.3, 1.4, title)

# final touches on the plots
for axs in axsGrids[-1,:]:
    axs.set_xlabel("$\log_{10} [ n / cm^{-3} ]$")
for axs in axsGrids[:,0]:
    axs.set_ylabel("$\log_{10} G_0$")

removeSubplotLabels( nx, ny, axsGrids)

pyl.show()

fig.savefig( outFname )
print 'ok'
