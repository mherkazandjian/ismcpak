# plot a grid into one panel with colorbar, contour lines (values are user specified)
#------------------------------------------------------------------------------------
import numpy as np
import pickle, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
import matplotlib.cm as cm
from fetchGridData import fetchRadexGrid

specStr       = '13CO'
transition    = '6-5'
cLevels       = [-15, -10, -8, -7, -6, -6.5, -5.5, -5.2, -5, -4, -3, -2]
cbarTicks     = np.arange(-18,-1, 2)
dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr
parmsFile     = dirname + 'parms.out'
fileInfoFile  = dirname + 'filesInfo.out'
colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/%s-%s-base.eps' % (specStr,transition)
#====================================================================================== 

fig    = pyl.figure(0, figsize = (6, 6) )
ax1    = fig.add_axes( [0.25, 0.2 , 0.55, 0.55])
axCbar = fig.add_axes( [0.2, 0.8, 0.65, 0.02])

# getting the grid data
grd = fetchRadexGrid( dirname  = dirname, specStr = specStr, 
                      gmechSec = 1e-10, transition = transition, verbose = True)

#getting the range over which the data was generated
parms = pickle.load(open(parmsFile))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

im = ax1.imshow(grd, extent = rangesLst, origin='lower', cmap = colormap)

CS = ax1.contour(grd, cLevels, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')
pyl.clabel(CS, fmt = '%.1f' )
cbar = pyl.colorbar(im, cax = axCbar, orientation = 'horizontal')

ax1.set_xlabel('$\log_{10} [ n_{gas} / (cm^{-3}) ] $')
ax1.set_ylabel('$\log_{10} [ G_0 ] $')

axCbar.set_title('$\log_{10}[ %s (%s) erg.cm^{-2} s^{-1} ]$' % (specStr, transition) )
cbar.set_ticks( cbarTicks )

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pyl.show()

print 'min-max log10 intensities = ', np.nanmin(grd), np.nanmax(grd) 
print 'done'
