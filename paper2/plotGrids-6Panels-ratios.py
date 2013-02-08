# plot a grid of line ratios into one panel with colorbar, contour lines (values are user specified)
#   grid1/grid2
# for different values of relative Gmech.
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
specStr1       = 'CO'
transition1    = '6-5'
#grid 2
specStr2       = 'CO'
transition2    = '6-5'
dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
#-----------------------------------

relGmech      = [[1e-3, 1e-2,5e-2], [0.1, 0.5, 1.0 ] ]
log_v_range   = [-2, 2] # range of the values, also that of the cbar
#cLevels       = np.arange(log_v_range[0], log_v_range[1] + 0.1, 0.5)
cLevels       = [0, 1, 2, 3, 3.5]
cbarTicks     = np.arange(log_v_range[0], log_v_range[1], 1)
showLogLabels = False
verbose       = True

dirname1       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr1
parmsFile1     = dirname1 + 'parms.out'
fileInfoFile1  = dirname1 + 'filesInfo.out'

dirname2       = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/%s/' % specStr2
parmsFile2     = dirname2 + 'parms.out'
fileInfoFile2  = dirname2 + 'filesInfo.out'

colormap      = cm.jet
imageSavePath = '/home/mher/ism/docs/paper02/src/figs/lineRatio-%s-%s-%s-%s-gMech.eps' % (specStr1, transition1, specStr2, transition2)
#=====================================================================================

ny = len(relGmech)
nx = len(relGmech[0])

# the indicies on the plot (x,y), x = horizontal pos, y = vertical position
panelInds = ( 
             ((0,0),(0,1),(0,2)),
             ((1,0),(1,1),(1,2)),
            ) 

#--------------setting up the axes---------------------
fig, panelAxs = pyl.subplots(ny, nx, sharex = True, sharey = True, figsize = (6,6))
pyl.subplots_adjust( left = 0.15, bottom = 0.15, right = 0.95, top = 0.7,
                     wspace = 0.0, hspace = 0.0)
axCbar     = fig.add_axes( [0.2, 0.75, 0.7, 0.02])
#------------done setting up axes----------------------

parms = pickle.load(open(parmsFile1))
ranges = parms['plotRanges']
rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

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
        
        
        grd1 = fetchRadexGrid(dirname  = dirname1, specStr = specStr1, 
                              gmechSec = rGmech, transition = transition1, verbose = verbose)
        grd2 = fetchRadexGrid(dirname  = dirname2, specStr = specStr2, 
                              gmechSec = rGmech, transition = transition2, verbose = verbose)

        #taking the ratios of the base grid with the current one                
        data = np.log10(10.0**grd1 / 10.0**grd2)
        dataGrds[i][j] = data
                
        ax = panelAxs[i][j]
        #im = ax.imshow(dataGrds[i][j], extent = rangesLst, origin='lower')
        ax.text(1,5,'%.2f%%' % (100*rGmech))
                
        indsFinit = np.where(np.fabs(data) != np.inf)
        dataTmp = data[indsFinit]
        currMin = np.nanmin(dataTmp)
        currMax = np.nanmax(dataTmp)
        print 'min = %e, max = %e' % (currMin, currMax)

        mn = np.min([mn, currMin])
        mx = np.max([mx, currMax])
                
        print '-----------------------------'
print 'global min = %f, global max = %f' % (mn,mx)
#-----------------------done getting the data for the panels------------------------

#-----------------plotting the colorbar and setting up the colorbar-----------------
cbarData = np.linspace(0, 1, 500)
cbarv = []
for i in np.arange(50):
    cbarv.append( cbarData.copy() )
cbarv = np.array(cbarv)
im = axCbar.imshow(cbarv, aspect = 'auto', vmin= 0, vmax = 1, cmap = colormap,
                   extent = [log_v_range[0], log_v_range[1], 0, 1])

axCbar.axes.get_yaxis().set_ticks([])
axCbar.set_title('$\log_{10}[%s (%s)/%s (%s)]$' % (specStr1,transition1,specStr2,transition2))
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
        indsBelow = np.where( grd < log_v_range[0])
        grd[indsBelow] = log_v_range[0]
        indsAbove = np.where( grd > log_v_range[1])
        grd[indsAbove] = log_v_range[1]

        grdScaled = scale(grd, 0.0, 1.0, log_v_range[0], log_v_range[1]) 
        ax = panelAxs[i][j] 
        im = ax.imshow(grdScaled, extent = rangesLst, origin='lower', cmap = colormap,
                       vmin = 0.0, vmax = 1.0, norm = None)
        
        CS = ax.contour(grd, cLevels, 
                extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                origin='lower', 
                colors = 'black')
        pyl.clabel(CS, fmt = '%.1f' )

    print '-----------------------------'

fig.text(0.45, 0.08, '$\log_{10} [ n_{gas} / (cm^{-3}) ] $')
fig.text(0.08, 0.45, '$\log_{10} [ G_0 ] $', rotation = 'vertical')

fig.savefig(imageSavePath)
print 'saved image to file : %s' % imageSavePath
pyl.show()
print 'done'
