# plot a grid from a pdr quantity (see mesh data format) into one panel with
# colorbar, contour lines (values are user specified)
#------------------------------------------------------------------------------------
import numpy as np
import pickle, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
import matplotlib.cm as cm

subDirName    = 'surfaceHeating'
plotTitle     = 'surfaceHeating' 
quantity      = ['therm','heating']
slabIdx       = 0
cLevels       = np.arange(20) - 30
cbarTicks     = np.arange(-6,-4,0.1)
dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % subDirName
parmsFile     = dirname + 'parms.out'
fileInfoFile  = dirname + 'filesInfo.out'
colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/%s-%d-base.eps' % (subDirName, slabIdx)
#====================================================================================== 

fig    = pyl.figure(0, figsize = (6, 6) )
ax1    = fig.add_axes( [0.25, 0.2 , 0.55, 0.55])
axCbar = fig.add_axes( [0.2, 0.8, 0.65, 0.02])

#loading the pickle file holding the parameters used to generate the data (a dict)
parms = pickle.load(open(parmsFile))
#loading the picke file holding the information of the file (a list of dicts)
filesInfo = pickle.load(open(fileInfoFile))

print '-------------------------------------------'
print '   zSec       transition    filename '
print '-------------------------------------------'
for fileInfo in filesInfo:
    print '%1.3e       %s         %s' % (10.0**fileInfo['zSec'], fileInfo['slabIdx'], fileInfo['filename'])

ranges = parms['plotRanges']

#looking for the reference grid (zero gmech) for this transition
for fileInfo in filesInfo:
    if (np.fabs(10.0**fileInfo['zSec'] - 1e-10) < 1e-14):
        print 10.0**fileInfo['zSec'], fileInfo['slabIdx'], fileInfo['filename']
        fname = fileInfo['filename']        
grd = np.loadtxt(fname, dtype = np.float64)

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

axCbar.set_title(plotTitle)
cbar.set_ticks( cbarTicks )

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pyl.show()

print 'min-max log10 intensities = ', np.nanmin(grd), np.nanmax(grd) 
print 'done'
