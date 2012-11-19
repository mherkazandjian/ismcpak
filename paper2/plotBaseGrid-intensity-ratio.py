# plot a grid of line ratios into one panel with colorbar, contour lines (values are user specified)
#   grid1/grid2
#------------------------------------------------------------------------------------
import numpy as np
import pickle, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
import matplotlib.cm as cm
from fetchGridData import fetchRadexGrid

#-----------------------------------
#grid 1
specStr1       = 'HNC'
transition1    = '01-00'
#grid 2
specStr2       = 'HCN'
transition2    = '1-0'
dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
#-----------------------------------
cLevels       = [-2, -1, 0, 1, 2]
cbarTicks     = np.arange(-2,2,0.1)

dirname1       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr1
parmsFile1     = dirname1 + 'parms.out'
fileInfoFile1  = dirname1 + 'filesInfo.out'

dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
parmsFile2     = dirname2 + 'parms.out'
fileInfoFile2  = dirname2 + 'filesInfo.out'

colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/lineRatio-%s-%s-%s-%s-base.eps' % (specStr1, transition1, specStr2, transition2)
#====================================================================================== 

fig    = pyl.figure(0, figsize = (6, 6) )
ax1    = fig.add_axes( [0.25, 0.2 , 0.55, 0.55])
axCbar = fig.add_axes( [0.2, 0.8, 0.65, 0.02])

#==========================getting grid1 data and information======================================
grd1 = fetchRadexGrid( dirname  = dirname1, specStr = specStr1, 
                      gmechSec = 1e-10, transition = transition1)
#==========================getting grid2 data and information======================================
grd2 = fetchRadexGrid(dirname  = dirname2, specStr = specStr2, 
                      gmechSec = 1e-10, transition = transition2)

parms = pickle.load(open(parmsFile1))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

#computing the log of the ratios
grd = np.log10(10.0**grd1 / 10.0**grd2)

im = ax1.imshow(grd, extent = rangesLst, origin='lower', cmap = colormap)

CS = ax1.contour(grd, cLevels, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')
pyl.clabel(CS, fmt = '%.1f' )
cbar = pyl.colorbar(im, cax = axCbar, orientation = 'horizontal')

ax1.set_xlabel('$\log_{10} [ n_{gas} / (cm^{-3}) ] $')
ax1.set_ylabel('$\log_{10} [ G_0 ] $')

axCbar.set_title('$\log_{10}[ %s (%s) / %s (%s) ]$' % (specStr1, transition1, specStr2, transition2) )
cbar.set_ticks( cbarTicks )

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pyl.show()

print 'min-max log10 intensities = ', np.nanmin(grd), np.nanmax(grd) 
print 'done'
