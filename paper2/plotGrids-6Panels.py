# plot a grid into one panel with colorbar, contour lines (values are user specified)
#------------------------------------------------------------------------------------
import numpy
import pickle, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import matplotlib.cm as cm
from mylib.utils.misc import scale as scale
from mylib.utils.misc import rm_noise

from fetchGridData import fetchRadexGrid

specStr       = 'CO'

transition    = '1-0'
#transition    = '2_1.5-1_1.5'

plotTitle     = r'$\log_{10}$[$R$(%s(%s))]' % (specStr,transition)
#plotTitle     = r'$\log_{10}$[$R$(HCO$^+$(%s))]' % (transition)
#plotTitle     = r'$\log_{10}$(%s($2_{3/2}$-$1_{3/2}$))' % (specStr,)
#plotTitle     = r'$\log_{10}$[$R$(%s(4-3))]' % (specStr,)

Av_max        = 10.0
relGmech      = [[1e-3, 1e-2, 5e-2], [0.1, 0.5, 1.0 ] ]
v_range       = [-2, 2] # range of the values, also that of the cbar
cLevels       = [0]
cbarTicks     = numpy.arange(v_range[0], v_range[1], 1)
dirname       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr
parmsFile     = dirname + 'parms.out'
fileInfoFile  = dirname + 'filesInfo.out'
colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/%s-%s-gMech.eps' % (specStr,transition)
#imageSavePath = '/home/mher/foo.eps'
#=====================================================================================
ny = len(relGmech)
nx = len(relGmech[0])

# the indicies on the plot (x,y), x = horizontal pos, y = vertical position
panelInds = ( 
             ((0,0),(0,1),(0,2)),
             ((1,0),(1,1),(1,2)),
            ) 

#--------------setting up the axes---------------------
width  = 5    #figure width (non normalized) 
height = 4.6  #figure heigh (non normalized)
as_rat = width/height #aspect ratio of the figure

fig, panelAxs = pylab.subplots(ny, nx, sharex = True, sharey = True, figsize = (width,height))
fig.set_facecolor('white')

pylab.subplots_adjust(left = 0.1, bottom = 0.15, right = 0.98, top = 0.8,
                     wspace = 0.0, hspace = 0.0)
axCbar = fig.add_axes( [0.15, 0.87, 0.80, 0.04])
#------------done setting up axes----------------------

#loading the pickle file holding the parameters used to generate the data (a dict)
parms = pickle.load(open(parmsFile))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

# getting the reference grid data
refGrd, filenameRef = fetchRadexGrid(dirname  = dirname, specStr = specStr, Av_max = Av_max,
                                     gmechSec = 1e-10, transition = transition, verbose = True)

print '###############################################################'

# setting the grids to be displayed with the desired mechanical heating percentages         
dataGrds = [[None,None,None   ], [None,None,None] ]

mn = -1.0 #initializing the mn, mx to defaults before finding them
mx =  1.0

#---------------------------------------------------------------------------------
# looping over the indicies in the image to be plotted
# and gettting the data for each panel
#---------------------------------------------------------------------------------
for r in panelInds: # looping over the rows (r)
    for inds in r:  # looping over each entry in r
        
        i, j = inds
        print inds
        rGmech = relGmech[i][j] # relative gMech of the grid to be displayed
        
        print 'indicies = ', i,j
        print 'relative gmech = ', rGmech

        grd, filename = fetchRadexGrid( dirname  = dirname, specStr = specStr, Av_max = Av_max, 
                                        gmechSec = rGmech, transition = transition, verbose = False)

        #taking the ratios of the base grid with the current one                
        data = numpy.log10(10.0**grd / 10.0**refGrd)
        dataGrds[i][j] = data
                
        ax = panelAxs[i][j]
        #im = ax.imshow(dataGrds[i][j], extent = rangesLst, origin='lower')
        ax.text(1,5,'%.2f%%' % (100*rGmech), size='large')
                
        indsFinit = numpy.where(numpy.fabs(data) != numpy.inf)
        dataTmp = data[indsFinit]
        currMin = numpy.nanmin(dataTmp)
        currMax = numpy.nanmax(dataTmp)
        print 'min = %e, max = %e' % (currMin, currMax)

        mn = numpy.min([mn, currMin])
        mx = numpy.max([mx, currMax])
        
        print '##########################'
print 'global min = %f, global max = %f' % (mn,mx)

#-----------------------done getting the data for the panels------------------------

#-----------------plotting the colorbar and setting up the colorbar-----------------
cbarData = numpy.linspace(0, 1, 500)
cbarv = []
for i in numpy.arange(50):
    cbarv.append( cbarData.copy() )
cbarv = numpy.array(cbarv)
im = axCbar.imshow(cbarv, aspect = 'auto', vmin= 0, vmax = 1, cmap = colormap,
                   extent = [v_range[0], v_range[1], 0, 1])
axCbar.axes.get_yaxis().set_ticks([])
axCbar.set_title(plotTitle, size='large')
#----------------------done plotting the colorbar-----------------------------------

#-----------------------------------------------------------------------------------
# now the we have the data for each panel
#looping over the indicies of the figure to be plotted (the panels in the plot)
# and showin the panel in the image
for r in panelInds:
    for inds in r:
        
        i, j = inds
        print inds
        rGmech = relGmech[i][j]
        grd = dataGrds[i][j]
        
        print 'indicies = ', i,j
        print 'relative gmech = ', rGmech        

        #setting the values which are outside the allowed ranges to the max 
        #allowed ranges
        indsBelow = numpy.where( grd < v_range[0])
        grd[indsBelow] = v_range[0]
        indsAbove = numpy.where( grd > v_range[1])
        grd[indsAbove] = v_range[1]

        grdScaled = scale(grd, 0.0, 1.0, v_range[0], v_range[1]) 
        ax = panelAxs[i][j]
                
        im = ax.imshow(grdScaled, extent = rangesLst, origin='lower', cmap = colormap,
                       vmin = 0.0, vmax = 1.0, norm = None)
        
        cLevelsScaled = scale(cLevels, 0.0, 1.0, v_range[0], v_range[1])
        CS = ax.contour(grdScaled, cLevelsScaled, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')
        #pylab.clabel(CS, fmt = '%.1f' )

    print '-----------------------------'

fig.text(0.35, 0.05, '$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='large')
fig.text(0.02, 0.55, '$\log_{10}$ [G$_0$]', size='large', rotation = 'vertical')

fig.savefig(imageSavePath)
print 'saved image to file : %s' % imageSavePath

pylab.show()
print 'done'
