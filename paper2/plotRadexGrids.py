import numpy as np
import pickle
import pylab as pyl

dirname      = '/home/mher/ism/runs/oneSided/dynamicMeshTest1/analysis/CO/'
parmsFile    = dirname + 'parms.out-Wed-Nov--7-00:38:15-2012'
fileInfoFile = dirname + 'filesInfo.out-Wed-Nov--7-00:38:15-2012'
transition   = '1-0'

#loading the pickle file holding the parameters used to generate the data (a dict)
parms = pickle.load(open(parmsFile))
#loading the picke file holding the information of the file (a list of dicts)
filesInfo = pickle.load(open(fileInfoFile))

print '-------------------------------------------'
print '   zSec       transition    filename '
print '-------------------------------------------'
for fileInfo in filesInfo:
    print '%1.3e       %s         %s' % (10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename'])

ranges = parms['plotRanges']

#looking for the reference grid (zero gmech) for this transition
for fileInfo in filesInfo:
    if (np.fabs(10.0**fileInfo['zSec'] - 1e-10) < 1e-14) and fileInfo['transition'] == transition:
        print 10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename']
        fname = fileInfo['filename']        
refGrd = np.loadtxt(fname, dtype = np.float64)

rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

fig = pyl.figure(0, figsize = (12,12) )

ax1 = fig.add_axes( [0.05, 0.05, 0.35, 0.35])
im = ax1.imshow(refGrd, extent = rangesLst, origin='lower')
pyl.colorbar(im, shrink = 0.8, extend = 'both')


ax00 = fig.add_axes( [0.05      , 0.4, 0.2, 0.2])
ax01 = fig.add_axes( [0.05 + 0.2, 0.4, 0.2, 0.2])
ax02 = fig.add_axes( [0.05 + 0.4, 0.4, 0.2, 0.2])
ax10 = fig.add_axes( [0.05      , 0.6, 0.2, 0.2])
ax11 = fig.add_axes( [0.05 + 0.2, 0.6, 0.2, 0.2])
ax12 = fig.add_axes( [0.05 + 0.4, 0.6, 0.2, 0.2])
relAxs = [ [ax00, ax01, ax02], [ax10, ax11, ax12] ]



print '###############################################################'

# the indicies on the plot (x,y), x = horizontal pos, y = vertical position
allInds = ( 
           ((0,0),(0,1),(0,2)),
           ((1,0),(1,1),(1,2)),
          ) 

# setting the grids to be displayed with the desired mechanical heating percentages         
relGmech = [[1e-3, 1e-2,5e-2], [0.1, 0.5, 1.0 ] ]
dataGrds = [[None,None,None   ], [None,None,None] ]

mn = -1.0
mx =  1.0

# looping over the indicies in the image to be plotted
# and gettting the data for each panel
for r in allInds:
    for inds in r:
        
        i, j = inds
        print inds
        rGmech = relGmech[i][j] # relative gMech of the grid to be displayed
        
        print 'indicies = ', i,j
        print 'relative gmech = ', rGmech
        
        for fileInfo in filesInfo:
            
            #finding the grid to be plotted in the poistion i,j in the figure 
            if np.fabs(10.0**fileInfo['zSec'] - rGmech) < 1e-14 and fileInfo['transition'] == transition:
                print 10.0**fileInfo['zSec'], fileInfo['transition'], fileInfo['filename']
                fname = fileInfo['filename']  #<----this is the name of the file holding the data
                grd = np.loadtxt(fname, dtype = np.float64)
                
                #taking the ratios of the base grid with the current one                
                data = np.log10(np.fabs(1.0 - 10.0**grd / 10.0**refGrd))
                dataGrds[i][j] = data
                
                ax = relAxs[i][j]
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

print '###########################################################################'
# now the we have the data for each panel
#looping over the indicies of the figure to be plotted (the panels in the plot)
# and showin the panel in the image
for r in allInds:
    for inds in r:
        
        i, j = inds
        print inds
        rGmech = relGmech[i][j]
        grd = dataGrds[i][j]
        
        print 'indicies = ', i,j
        print 'relative gmech = ', rGmech        
        
        ax = relAxs[i][j] 
        im = ax.imshow(grd, extent = rangesLst, origin='lower')
                
    print '-----------------------------'

print 'global min = %f, global max = %f' % (mn,mx)
axCbarDataAx = fig.add_axes( [0.3,0.9,0.4,0.02])
axCbar2 = fig.add_axes( [0.3,0.95,0.4,0.02])

cbarData = np.linspace(mn,mx,500)
cbarv = []
for i in np.arange(50):
    cbarv.append( cbarData.copy() )
cbarv = np.array(cbarv)
im = axCbarDataAx.imshow(cbarv)

pyl.colorbar(im, cax = axCbar2, ax = axCbarDataAx, 
             orientation = 'horizontal', 
             norm = False)

pyl.show()
        
print 'done'
