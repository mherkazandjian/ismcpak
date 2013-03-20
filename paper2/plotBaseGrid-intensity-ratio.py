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
from mylib.utils.misc import scale as scale
from fetchGridData import fetchRadexGrid

#-----------------------------------
#grid 1
specStr1       = 'O'
transition1    = '1-0'
#grid 2
specStr2       = 'C'
transition2    = '1-0'
dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
#-----------------------------------
log_v_range   = [-2, 2] # log10 of the range of the values, also that of the cbar
log_cbarTicks = np.arange(log_v_range[0], log_v_range[1] + 0.01, 0.5)
log_cLevels   = np.log10([0.1, 0.5, 1, 2, 4, 6, 8, 10, 50])
#log_cLevels   = np.arange(log_v_range[0], log_v_range[1] + 0.01, 0.5)
showLogLabels = False # if False the bar ticks are displayes as 10^x and the contour
                      # show the actual values  

dirname1       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr1
parmsFile1     = dirname1 + 'parms.out'
fileInfoFile1  = dirname1 + 'filesInfo.out'

dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
parmsFile2     = dirname2 + 'parms.out'
fileInfoFile2  = dirname2 + 'filesInfo.out'

Av_max         = 10.0
colormap      = cm.jet
#imageSavePath = '/home/mher/ism/docs/paper02/src/figs/lineRatio-%s-%s-%s-%s-base.eps' % (specStr1, transition1, specStr2, transition2)
imageSavePath = '/home/mher/foo.eps'
#======================================================================================
 
width  = 6    #figure width (non normalized) 
height = 6.8  #figure heigh (non normalized)
as_rat = width/height #aspect ratio of the figure

ax_xs  = 0.09 #axses x start (normalized)
ax_ys  = 0.09 #axses y start (normalized)
ax_sz  = 0.90 #axses size (normalized)

cbar_xs = ax_xs        #colorbar x start
cbar_ys = ax_sz + 0.03 #colorbar x start
cbar_sc = 0.9          #scale of the width of the cbar (relative to the width of ax)
cbar_w  = 0.02         #width of the cbar (normalized)

fig    = pyl.figure(0, figsize = (width, height) )
fig.set_facecolor('white')
ax1    = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])
axCbar = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_sz, cbar_ys, cbar_sc*ax_sz - (0.5*(1.0-cbar_sc))*ax_sz, cbar_w])

#==========================getting grid1 data and information======================================
grd1, filename1 = fetchRadexGrid(dirname  = dirname1, specStr = specStr1, Av_max = Av_max,
                                 gmechSec = 1e-10, transition = transition1)
#==========================getting grid2 data and information======================================
grd2, filename2 = fetchRadexGrid(dirname  = dirname2, specStr = specStr2,  Av_max = Av_max,
                                 gmechSec = 1e-10, transition = transition2)

parms = pickle.load(open(parmsFile1))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

#computing the log of the ratios
grd = np.log10(10.0**grd1 / 10.0**grd2)

#setting the values which are outside the allowed ranges to the max 
#allowed ranges
indsBelow = np.where( grd < log_v_range[0])
grd[indsBelow] = log_v_range[0]
indsAbove = np.where( grd > log_v_range[1])
grd[indsAbove] = log_v_range[1]

grdScaled = scale(grd, 0.0, 1.0, log_v_range[0], log_v_range[1]) 

im = ax1.imshow(grdScaled, extent = rangesLst, origin='lower', cmap = colormap,
                vmin = 0.0, vmax = 1.0, norm = None)

#-----------------plotting the colorbar and setting up the colorbar-----------------
cbarData = np.linspace(0, 1, 500)
cbarv = []
for i in np.arange(50):
    cbarv.append( cbarData.copy() )
cbarv = np.array(cbarv)
im = axCbar.imshow(cbarv, aspect = 'auto', vmin= 0, vmax = 1, cmap = colormap,
                   extent = [log_v_range[0], log_v_range[1], 0, 1])

axCbar.axes.get_yaxis().set_ticks([]) #removing y ticklabels

if showLogLabels == False:
    cbarLabeslStrs = []
    for tickv in log_cbarTicks:
        cbarLabeslStrs.append('$10^{%.1f}$' % tickv)
        print '$10^{%.1f}$' % tickv
    axCbar.set_xticklabels( cbarLabeslStrs )
    
    titleStr = '$%s (%s)/%s (%s)$' % (specStr1,transition1,specStr2,transition2)
else: 
    titleStr = '$\log_{10}[%s (%s)/%s (%s)]$' % (specStr1,transition1,specStr2,transition2)

axCbar.set_title( titleStr)
    
#----------------------done plotting the colorbar-----------------------------------

if showLogLabels == True:
    CS = ax1.contour(grd, np.array(log_cLevels), 
                     extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                     origin='lower', 
                     colors = 'black')
else:
    CS = ax1.contour(10**grd, 10**log_cLevels, 
                     extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                     origin='lower', 
                     colors = 'black')    
pyl.clabel(CS, fmt = '%.1f' )

ax1.set_xlabel('$\log_{10} [ n_{gas} / (cm^{-3}) ] $')
ax1.set_ylabel('$\log_{10} [ G_0 ] $')

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pyl.show()

print 'min-max log10 intensities = ', np.nanmin(grd), np.nanmax(grd) 
print 'done'
