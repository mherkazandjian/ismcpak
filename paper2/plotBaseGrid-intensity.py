# plot a grid into one panel with colorbar, contour lines (values are user specified)
#------------------------------------------------------------------------------------
import numpy
import pickle, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import matplotlib.cm as cm
from fetchGridData import fetchRadexGrid

specStr       = 'HNC'
Av_max        = 10.0    
transition    = '04-03'
#transition    = '2_1.5-1_1.5'
#plotTitle     = '$\log_{10}$[ %s (%s) erg.cm^{-2} s^{-1} ]$' % (specStr, transition)
#plotTitle     = r'$\log_{10}$([CII] 158$\mu m$)' 
#plotTitle     = r'$\log_{10}$([OI] 63$\mu m$)' 
#plotTitle     = r'$\log_{10}$(%s(%s))' % (specStr,transition)
#plotTitle     = r'$\log_{10}$(HCO$^+$(%s))' % (transition)
#plotTitle     = r'$\log_{10}$(%s($2_{3/2}$-$1_{3/2}$))' % (specStr,)
plotTitle     = r'$\log_{10}$(%s(4-3))' % (specStr,)
cLevels       = [-10, -9, -8, -8.5, -7, -7.5, -6.5, -6, -5, -4, -3, -2]
cbarTicks     = numpy.hstack([numpy.arange(-18,-1, 2), []])
#cLevels       = numpy.arange(-18,-1, 1)
#cbarTicks     = numpy.arange(-18,-1, 2)
#------------------------------------------
#clip          = None
clip          = [-15.0, 0.0]
removeNans    = True
saveCleaned   = True
#-----------------------------------------
dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr
parmsFile     = dirname + 'parms.out'
fileInfoFile  = dirname + 'filesInfo.out'
colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/%s-%s-base.eps' % (specStr,transition)
#imageSavePath = '/home/mher/foo.eps'
#====================================================================================== 

width  = 3    #figure width (non normalized) 
height = 3.4  #figure heigh (non normalized)
as_rat = width/height #aspect ratio of the figure

ax_xs  = 0.17 #axses x start (normalized)
ax_ys  = 0.17 #axses y start (normalized)
ax_sz  = 0.75 #axses size (normalized)

cbar_xs = ax_xs        #colorbar x start
cbar_ys = ax_sz + 0.13 #colorbar x start
cbar_sc = 0.99         #scale of the width of the cbar (relative to the width of ax)
cbar_w  = 0.02         #width of the cbar (normalized)

fig    = pylab.figure(0, figsize = (width, height) )
fig.set_facecolor('white')
ax1    = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])
axCbar = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_sz, cbar_ys, cbar_sc*ax_sz - (0.5*(1.0-cbar_sc))*ax_sz, cbar_w])

# getting the grid data
grd, filename = fetchRadexGrid( dirname  = dirname, specStr = specStr, 
                                gmechSec = 1e-10, transition = transition,  
                                Av_max = Av_max, verbose = True)

#---------printing the min and max (excluding nans and infs)
indsNotInf = numpy.where(numpy.invert(numpy.isinf(grd)))
print 'min-max log10 intensities = ', numpy.nanmin(grd[indsNotInf]), numpy.nanmax(grd[indsNotInf]) 
#------------------------------------------------------------

#getting the range over which the data was generated
parms = pickle.load(open(parmsFile))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

#---------------discarding and cliping the value ranges---------------------------
xrange, yrange = ranges[0], ranges[1]
grdShape = grd.shape
grdPts = numpy.meshgrid(numpy.linspace(0,grdShape[0]-1,grdShape[0]),
                        numpy.linspace(0,grdShape[1]-1,grdShape[1]))

if removeNans:
    cond = True
    cond = numpy.bitwise_and(numpy.isnan(grd), cond)
    cond = numpy.bitwise_and(grdPts[0] < 4.0*(grdShape[0]/xrange[1]), cond)
    cond = numpy.bitwise_or(grd >= -1.0, cond)
    indsNan = numpy.where(cond)
    grd[indsNan] = -18.0

if clip != None:
    grd = grd.clip(clip[0], clip[1])
#------------done discarding and cliping the value ranges---------------------------

im = ax1.imshow(grd, extent = rangesLst, origin='lower', cmap = colormap)

CS = ax1.contour(grd, cLevels, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')
pylab.clabel(CS, fmt = '%.1f' )
cbar = pylab.colorbar(im, cax = axCbar, orientation = 'horizontal')

ax1.set_xlabel('$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='large')
ax1.set_ylabel('$\log_{10}$ [G$_0$]', size='large')

axCbar.set_title(plotTitle)
cbar.set_ticks( cbarTicks )

fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pylab.show()

if saveCleaned:
    newFname = filename+'-cleaned'
    numpy.savetxt(newFname, grd, fmt = '%1.4e')
    print 'saved the cleaned file into:\n   %s' % newFname
    
print 'done'
