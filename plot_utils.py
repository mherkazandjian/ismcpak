import numpy
from numpy import log10
from scipy import interpolate
import pylab
import os
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator 
import matplotlib

from mylib.utils.misc import scale as scale
import lineDict
import meshUtils
import sys
import collections
import mylib.utils.removeAxesLabels
from leidenLambda import molData

try:    
    from PyQt4 import QtGui, QtCore
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
    #import matplotlib.pyplot as plt
except:
    print 'can not use radex in gui mode..failed to import PyQt4, or matplotlib.backends.backend_qt4agg'

def plot_intensity_grid(grd, line=None, ranges=None, cmap=None, cLevels=None, 
                             cbarTicks=None, imageSaveDir=None, removeNans=None, 
                             clip=None, clip_max_n=None, 
                             f_interp_dim=None, **kwargs):
    '''plots and saves an intensity grid''' 
    
    width  = 3    #figure width (non normalized) 
    height = 3.4  #figure heigh (non normalized)
    as_rat = width/height #aspect ratio of the figure
    
    ax_xs  = 0.19 #axses x start (normalized)
    ax_ys  = 0.19 #axses y start (normalized)
    ax_sz  = 0.75 #axses size (normalized)
    
    cbar_xs = ax_xs        #colorbar x start
    cbar_ys = ax_sz + 0.15 #colorbar x start
    cbar_sc = 0.99         #scale of the width of the cbar (relative to the width of ax)
    cbar_w  = 0.02         #width of the cbar (normalized)
    
    fig    = pylab.figure(figsize = (width, height) )
    fig.set_facecolor('white')
    ax1    = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])
    axCbar = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_sz, cbar_ys, cbar_sc*ax_sz - (0.5*(1.0-cbar_sc))*ax_sz, cbar_w])
    
    #---------printing the min and max (excluding nans and infs)
    #indsNotInf = numpy.where(numpy.invert(numpy.isinf(grd)))
    #print 'min-max log10 intensities = ', numpy.nanmin(grd[indsNotInf]), numpy.nanmax(grd[indsNotInf]) 
    #------------------------------------------------------------
    
    rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

    #---------------discarding and cliping the value ranges---------------------------
    xrange, yrange = ranges[0], ranges[1]
    grdShape = grd.shape
    grdPts = numpy.meshgrid(numpy.linspace(0, grdShape[0]-1, grdShape[0]),
                            numpy.linspace(0, grdShape[1]-1, grdShape[1]))
    grdPts[0] = grdPts[0].T
    grdPts[1] = grdPts[1].T
    
    if removeNans != None and removeNans==True:
        cond = True
        
        cond = numpy.bitwise_and(numpy.isnan(grd), cond)   #mask of points which are nans (true => point is nan)
        
        if clip_max_n != None:
            cond = numpy.bitwise_and(grdPts[0] < clip_max_n*(grdShape[0]/xrange[1]), cond) #mask of points which are nans for n <
         
        cond = numpy.bitwise_or(grd >= -1.0, cond)
        indsNan = numpy.where(cond)
        grd[indsNan] = -18.0

    if clip != None:
        grd = grd.clip(clip[0], clip[1])
    #------------done discarding and cliping the value ranges---------------------------

    im = ax1.imshow(grd.T, extent = rangesLst, origin='lower', cmap = cmap)
    
    CS = ax1.contour(grd.T, cLevels, 
                    extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                    origin='lower', 
                    colors = 'black')
    
    pylab.clabel(CS, fmt = '%.1f' )
    cbar = pylab.colorbar(im, cax = axCbar, orientation = 'horizontal')
    
    ax1.set_xlabel('$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='large')
    ax1.set_ylabel('$\log_{10}$ [G$_0$]', size='large')
    
    plotTitle = r'$\log_{10}$(%s / erg cm$^{-2}$ s$^{-1}$)' %  lineDict.lines[line['code']]['latex']
    axCbar.set_title(plotTitle, size='small')
    cbar.set_ticks( cbarTicks )
    
    fig.canvas.set_window_title(line['code'])
    
    imageSavePath = os.path.join(imageSaveDir,line['code']+'-base.eps') 
    fig.savefig(imageSavePath)
    print 'wrote image : %s' % imageSavePath
    pylab.show()
    
    print 'done'
    
def plot_relative_change_intensity_grid(arxvPDR, line=None, ranges=None, relGmech_ref=None,
                                        relGmech=None, v_range=None, cLevels=None, 
                                        cmap=None, imageSaveDir=None, removeNans=None,
                                        clip=None, clip_max_n=None, Av_use=None, res=None, 
                                        f_interp_dim=None, **kwargs):
    '''plots 6 grids showing the relative change in the intesitiy as a function of the relative mechanical
    heating in the grids
    '''

    nx, ny = 3, 2
    
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

    rangesLst = (ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]) 

    ref_grd = arxvPDR.get_emission_grid_from_databases(line=line, Av_use=Av_use, ranges=ranges, 
                                                       res=res, z_sec=relGmech_ref, f_interp_dim=f_interp_dim)

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
    
            grd = arxvPDR.get_emission_grid_from_databases(line=line, Av_use=Av_use, ranges=ranges, 
                                                           res=res, z_sec=numpy.log10(rGmech), f_interp_dim=f_interp_dim)
    
            #taking the ratios of the base grid with the current one                
            data = numpy.log10(10.0**grd / 10.0**ref_grd)
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
    im = axCbar.imshow(cbarv, aspect = 'auto', vmin= 0, vmax = 1, cmap = cmap,
                       extent = [v_range[0], v_range[1], 0, 1])
    axCbar.axes.get_yaxis().set_ticks([])

    plotTitle = r'$\log_{10}$[R(%s)]' %  lineDict.lines[line['code']]['latex']
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
                    
            im = ax.imshow(grdScaled.T, extent = rangesLst, origin='lower', cmap = cmap,
                           vmin = 0.0, vmax = 1.0, norm = None)
            
            cLevelsScaled = scale(cLevels, 0.0, 1.0, v_range[0], v_range[1])
            
            CS = ax.contour(grdScaled.T, cLevelsScaled, 
                    extent=(ranges[0][0], ranges[0][1], ranges[1][0], ranges[1][1]), 
                    origin='lower', 
                    colors = 'black')
            #pylab.clabel(CS, fmt = '%.1f' )
    
        print '-----------------------------'
    
    fig.text(0.35, 0.05, '$\log_{10}$ [n$_{gas}$ / cm$^{-3}$]', size='large')
    fig.text(0.02, 0.55, '$\log_{10}$ [G$_0$]', size='large', rotation = 'vertical')
    
    imageSavePath = os.path.join(imageSaveDir, line['code']+'-gMech.eps') 
    fig.savefig(imageSavePath)
    
    fig.canvas.set_window_title(line['code'])
    
    pylab.show()
    print 'done'
    


class line_ratio(object):
    
    def __init__(self, 
                 arxvPDR,     #: an PDR mesh database object
                 line1={      #: a dict holding the info of the first line (should be defined in the lineDict module
                        'code' : 'O63',
                        'type' : 'pdr'
                        },  
                 line2={      #: a dict holding the info of the second line (should be defined in the lineDict module
                        'code' : 'C+158',
                        'type' : 'pdr'
                       }, 
                 ranges=[[0.0, 6.0], [0.0, 6.0]],  #ranges in n and G0
                 cmap={#: a dict holding the colormap, its range and whether the colors correpond to log quantities or not
                       'obj'     : matplotlib.cm.jet,
                       'v_range' : [-2, 2],
                       'log'     : True,
                      },
                 cbar={#: a dict holding the info about the ticks of the colorbar and their format
                       'ticks'  : numpy.linspace(-2, 2, 5),
                       'format' : '%.1f',  # None, '10^x', '%e' (or any format)
                      },
                 c_levels={#: a dict holding the info about the countour lines their format
                           'values' : numpy.log10([0.01, 0.03,  0.1,  0.3, 1,  3,  10,    30,   100]),
                           'format' : '%.2f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          }, 
                 Av_use=10.0, #: the Av of the grid
                 res=[100,100],  #: the resolution of the interpolated image
                 z_sec=-10.0, #: the section in mechanical heating (in the default we assume it is a db with relative gmech with alpha=1e-10
                 image_save_dir=None, 
                 removeNans=None, 
                 log_title=None,
                 clip=None, 
                 clip_max_n=None, 
                 f_interp_dim=None, 
                 only_init=None, 
                 only_get_data =None,
                 lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                 **kwargs):

        #setting the parameters and keywords as attributes
        self.arxvPDR = arxvPDR
        self.line1 = line1
        self.line2 = line2
        self.Av_use = Av_use
        self.z_sec = z_sec
        self.ranges = ranges
        self.cmap = cmap
        self.cbar = cbar
        self.c_levels = c_levels
        self.res = res
        self.clip = clip
        self.clip_max = clip_max_n
        self.f_interp_dim = f_interp_dim
        self.log_title = log_title
        self.image_save_dir = image_save_dir
        self.lambdaPath = lambdaPath
        
        #attributes related to the gui and where things will be plotted
        self.figure  = None
        self.axs     = None
        self.cbar_ax = None
        
        self.grd = None
        self.grd_scaled = None
        
        self.im = None
        self.cbar_im = None
        self.contours = None
        self.nc1_plt = None
        self.nc2_plt = None
                 
        if only_get_data != None and only_get_data == True:
            self.grd, self.grd_scaled = self.get_grid()
        else: 
            #defining and setting the values of the above attributes
            self.figure, self.axs, self.cbar_ax = self.init_plot()
                
            #plotting the grid with the accessories (using dummy data)
            self.setup_plot()
        
            if only_init != None and only_init == True:
                pass
            else:    
                #updating the plot
                self.update_plot()
                pylab.show()
    
    def get_critical_densities(self, T_kin=100.0):
        '''returns a list of the critical densities of the lines in self.line1, self.line2 for H2'''
        
        line1 = lineDict.lines[self.line1['code']]
        line2 = lineDict.lines[self.line2['code']]
        
        spec1Str = line1['specStr']
        spec2Str = line2['specStr']
        
        #an instance of the radex object (will be used to recover the species info paths, 
        #and hence critical densities of lines         
        reader = molData.reader(dirPath = self.lambdaPath, species = [spec1Str, spec2Str])
 
        nc1 = reader.get_line_mean_critical_density_H2(self.line1['code'], T_kin)
        nc2 = reader.get_line_mean_critical_density_H2(self.line2['code'], T_kin)
        
        return numpy.log10(nc1), numpy.log10(nc2)        
        
    def init_plot(self):
        '''initialized the plotting reigon, the colorbar and the axes'''
        
        width  = 3    #figure width (non normalized) 
        height = 3.4  #figure heigh (non normalized)
        as_rat = width/height #aspect ratio of the figure
        
        ax_xs  = 0.17 #axses x start (normalized)
        ax_ys  = 0.17 #axses y start (normalized)
        ax_sz  = 0.75 #axses size (normalized)
        
        cbar_xs = ax_xs        #colorbar x start
        cbar_ys = ax_sz + 0.15 #colorbar x start
        cbar_sc = 0.99         #scale of the width of the cbar (relative to the width of ax)
        cbar_w  = 0.02         #width of the cbar (normalized)
        
        fig = pylab.figure(figsize = (width, height) )
        axs = fig.add_axes([ax_xs, ax_ys*as_rat, ax_sz, ax_sz*as_rat])
        cbar_ax = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_sz, cbar_ys, cbar_sc*ax_sz - (0.5*(1.0-cbar_sc))*ax_sz, cbar_w])

        #setting the labels of the axes
        axs.set_xlabel(r'$\log_{10}$ [n$_{gas}$]')    
        axs.set_ylabel(r'$\log_{10}$ [G$_0$]')
        
        self.figure = fig
        self.axs = axs
        self.cbar_ax = cbar_ax
        
        self.figure.set_facecolor('white')
                
        return fig, axs, cbar_ax
        
    def get_grid(self):
        '''getting grid1 data and information'''
        
        grd1 = self.arxvPDR.get_emission_grid_from_databases(line=self.line1, Av_use=self.Av_use, ranges=self.ranges, res=self.res, z_sec=self.z_sec, f_interp_dim=self.f_interp_dim)
        grd2 = self.arxvPDR.get_emission_grid_from_databases(line=self.line2, Av_use=self.Av_use, ranges=self.ranges, res=self.res, z_sec=self.z_sec, f_interp_dim=self.f_interp_dim)
    
        #computing the log of the ratios
        if self. cmap['log'] == True:
            grd = numpy.log10(10.0**grd1 / 10.0**grd2)
        else:
            grd = 10.0**grd1 / 10.0**grd2

        v_range = self.cmap['v_range']
        
        #setting the values which are outside the allowed ranges to the max 
        #allowed ranges
        indsBelow = numpy.where( grd < v_range[0])
        grd[indsBelow] = v_range[0]
        indsAbove = numpy.where( grd > v_range[1])
        grd[indsAbove] = v_range[1]
        
        grd_scaled = scale(grd, 0.0, 1.0, v_range[0], v_range[1]) 

        return grd, grd_scaled
        
    def setup_plot(self):
        '''sets up the plot with dummy data to intialized the objects used for the plotting'''
        
        rangesLst = (self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1])
    
        #dummy data
        dummy_data = numpy.zeros(numpy.array(self.res, 'i')) 
        
        self.im = self.axs.imshow(dummy_data, extent=rangesLst, origin='lower',     
                                  cmap=self.cmap['obj'], vmin=0.0, vmax=1.0, norm=None)

        #-----------------plotting the colorbar and setting up the colorbar--------------------------------------
        #dummy cbar data
        dummy_cbar_data = numpy.zeros(numpy.array((50, 500), 'i')) 
        
        self.cbar_im = self.cbar_ax.imshow(dummy_cbar_data, aspect='auto', vmin=0, vmax=1, 
                                           cmap=self.cmap['obj'], 
                                           extent = [self.cmap['v_range'][0], self.cmap['v_range'][1], 0.0, 1.0])

        self.cbar_ax.axes.get_yaxis().set_ticks([]) #removing y ticklabels
        #----------------------done plotting the colorbar------------------------------------------------------
        
        
        self.contours = self.axs.contour(dummy_data, numpy.array(self.c_levels['values']), 
                                         extent=(self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1]),    
                                         origin='lower', colors = 'black')
        
        #plotting dummy vertical line which will eventually hold the plot 
        #locations of the critical densities of each line
        #self.nc1_plt, = self.axs.plot([0, 0], [0, 6], 'b-.', lw=2)
        #self.nc2_plt, = self.axs.plot([0, 0], [0, 6], 'r-.', lw=2)
        
    def update_plot(self):
        '''sets the actual data to the plot'''

        self.grd, self.grd_scaled = self.get_grid() 
            
        self.figure.canvas.set_window_title(self.line1['code']+'/'+self.line2['code'])
        
        self.im.set_data(self.grd_scaled.T)
        
        self.update_colorbar()
        
        self.update_contours()
        
        #self.update_positions_of_critical_density_lines()
        
        self.save_figure()
        
    def update_colorbar(self):

        cbarData = numpy.linspace(0, 1, 500)
        cbarv = []
        for i in numpy.arange(50):
            cbarv.append( cbarData.copy() )
        cbarv = numpy.array(cbarv)
        
        self.cbar_im.set_data(cbarv)

        self.cbar_ax.set_xticks( self.cbar['ticks'] )
        
        #setting the labels of the colorbar    
        cbarLabeslStrs = []
        
        for tickv in self.cbar['ticks']:
            if self.cbar['format'] == '10^x':
                cbarLabeslStrs.append(r'$10^{%.1f}$' % tickv)
            elif self.cbar['format'] == None:
                cbarLabeslStrs.append(tickv)
            else:
                cbarLabeslStrs.append(self.cbar['format'] % tickv)

        #setting the title of the colorbar
        titleStr = r'%s/%s' % (lineDict.lines[self.line1['code']]['latex'], 
                               lineDict.lines[self.line2['code']]['latex'])
        
        #showing the critical densities of each line involved
        
        if self.log_title != None and self.log_title==True:
            titleStr = r'$\log_{10}$(%s)' % titleStr
        
        self.cbar_ax.set_title( titleStr)
        self.cbar_ax.set_xticklabels( cbarLabeslStrs )

    def update_contours(self):

        clabel_font_size = 10

        #deleting the countour lines (if they are found)
        if 'collections' in dir(self.contours):
            for c in self.contours.collections:
                paths = c.get_paths()
                del paths[:]
            #deleting the labels
            for txt in self.contours.labelTexts:
                txt.set_text('')

        self.contours = self.axs.contour(self.grd.T, numpy.array(self.c_levels['values']), 
                                         extent=(self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1]),    
                                         origin='lower', colors = 'black', linestyle='-')

        if self.c_levels['format']=='10^x':
            pylab.clabel(self.contours, fmt = ticker.LogFormatterMathtext())
        elif self.c_levels['format']==None:
            pylab.clabel(self.contours)
        elif self.c_levels['format']=='strs':  
            fmt={}
            for l,s in zip(self.contours.levels, self.c_levels['strs']):
                fmt[l]= s
                
            pylab.clabel(self.contours, fmt=fmt, inline=True, fontsize=clabel_font_size)
    
        else:
            
            if 'pow10' in self.c_levels and self.c_levels['pow10'] == True:
                
                #getting the contour labels as strings with the format provided by the user
                c_level_strs = []
                for v in self.c_levels['values']:

                    if numpy.int32(v) <= -3:
                        fstr = '%.0e' 
                    if numpy.int32(v) == -2:
                        fstr = '%.2f' 
                    if numpy.int32(v) == -1:
                        fstr = '%.1f' 
                    if numpy.int32(v) == 0:
                        fstr = '%.1f' 
                    if numpy.int32(v) == 1:
                        fstr = '%.0f' 
                    if numpy.int32(v) == 2:
                        fstr = '%.0f' 
                    if numpy.int32(v) >= 3:
                        fstr = '%.0e' 
                
                    c_level_strs.append(fstr % 10.0**v)
  
                fmt={}
                for l, s in zip(self.contours.levels, c_level_strs):
                    fmt[l]= s
                    
                pylab.clabel(self.contours, fmt=fmt, inline=True, fontsize=clabel_font_size)
  
            else:
                pylab.clabel(self.contours, fmt = self.c_levels['format'], fontsize=clabel_font_size)
        
        
        #getting the critical densities for the lines
    
    def update_positions_of_critical_density_lines(self):
        '''updates the locations of the lines marking the critical densities'''
        nc1, nc2 = self.get_critical_densities()
        
        self.nc1_plt.set_xdata([nc1, nc1])
        self.nc1_plt.set_ydata([self.ranges[1][0], self.ranges[1][1]])

        self.nc2_plt.set_xdata([nc2, nc2])
        self.nc2_plt.set_ydata([self.ranges[1][0], self.ranges[1][1]])
    
    def save_figure(self):
        
        if self.image_save_dir != None:
            image_save_path = os.path.join(self.image_save_dir, self.line1['code']+'-'+self.line2['code']+'.eps') 
            self.figure.savefig(image_save_path)
            
    def get_value_at_coords(self, x, y):
        '''retuns the value of the grid and the provided x and y values'''
        
        #getting the indicies
        x_ind = scale(x, 0, self.grd.shape[0], self.ranges[0][0], self.ranges[0][1])
        y_ind = scale(y, 0, self.grd.shape[1], self.ranges[0][0], self.ranges[0][1])

        return self.grd[x_ind, y_ind]


class line_ratio_app(QtGui.QWidget):
    '''A QT GUI application for displaying line ratio grids from a PDR grid
    
    .. todo:: add buttons, checkboxes..etc.. for the other line_ratio object
    which are not present in the current gui. 
    
    .. todo:: add an option to plot contours by pressing on the figure itself.
    '''
    
    def __init__(self, parent=None, arxvPDR=None, initial_parms={}, **kwargs):
        
        super(line_ratio_app, self).__init__(parent)
        
        #intitializing the line_ratio object where the maps will be plotted
        self.line_ratio = line_ratio(arxvPDR, only_init=True, **initial_parms)
        
        self.initUI()
        
    def initUI(self):
        
        #-----------------------------------------------------------------------
        # a figure instance to plot on
        self.figure = self.line_ratio.figure

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        self.canvas.mpl_connect('button_press_event', self.on_button_down_in_axes)
        #-----------------------------------------------------------------------
    
        self.init_input_fields()
        
        # Just some button connected to `plot` method
        self.button = QtGui.QPushButton('Plot')
        self.button.clicked.connect(self.plot)
        
        self.set_layout()
                
        self.setGeometry(300, 300, 1000, 700)
        self.setWindowTitle('line ratio explorer')    
        self.show()

    def init_input_fields(self):
        '''defining the input fields'''
        
        self.lbl_line1 = QtGui.QLabel('line1')
        default_text = self.line_ratio.line1['code']
        self.qle_line1 = QtGui.QLineEdit(default_text, self)
        self.qle_line1.editingFinished.connect(self.line1_editing_finished)
            
        self.lbl_type1 = QtGui.QLabel('type1')
        default_text = self.line_ratio.line1['type']
        self.qle_type1 = QtGui.QLineEdit(default_text, self)
        self.qle_type1.editingFinished.connect(self.type1_editing_finished)

        self.lbl_line2 = QtGui.QLabel('line2')
        default_text = self.line_ratio.line2['code']
        self.qle_line2 = QtGui.QLineEdit(default_text, self)
        self.qle_line2.editingFinished.connect(self.line2_editing_finished)

        self.lbl_type2 = QtGui.QLabel('type2')
        default_text = self.line_ratio.line2['type']
        self.qle_type2 = QtGui.QLineEdit(default_text, self)
        self.qle_type2.editingFinished.connect(self.type2_editing_finished)

        self.lbl_Av = QtGui.QLabel('Av')
        default_text = '%.2f' % self.line_ratio.Av_use
        self.qle_Av = QtGui.QLineEdit(default_text, self)
        self.qle_Av.editingFinished.connect(self.Av_editing_finished)

        self.lbl_g_mech = QtGui.QLabel('g_mech')
        default_text = '%.2e' % self.line_ratio.z_sec
        self.qle_g_mech = QtGui.QLineEdit(default_text, self)
        self.qle_g_mech.editingFinished.connect(self.g_mech_editing_finished)

        self.lbl_contours = QtGui.QLabel('contours')
        default_text = ''
        self.qle_contours = QtGui.QLineEdit(default_text, self)
        self.qle_contours.editingFinished.connect(self.set_contours_editing_finished)

    def set_layout(self):
        '''laying out the gui'''
        
        grid = QtGui.QGridLayout()
        grid.setSpacing(10)
        
        grid.addWidget(self.toolbar, 1, 0, 1, 8) #location 1,0 on the grid spanning 1 row and 3 columns

        #--------------------------------------------
        grid.addWidget(self.lbl_line1, 1, 8)
        grid.addWidget(self.qle_line1, 1, 9, 1, 6)

        grid.addWidget(self.lbl_type1, 2, 8)
        grid.addWidget(self.qle_type1, 2, 9, 1, 6)

        grid.addWidget(self.lbl_line2, 3, 8)
        grid.addWidget(self.qle_line2, 3, 9, 1, 6)

        grid.addWidget(self.lbl_type2, 4, 8)
        grid.addWidget(self.qle_type2, 4, 9, 1, 6)

        grid.addWidget(self.lbl_Av, 5, 8)
        grid.addWidget(self.qle_Av, 5, 9, 1, 6)
        
        grid.addWidget(self.lbl_g_mech, 6, 8)
        grid.addWidget(self.qle_g_mech, 6, 9, 1, 6)

        grid.addWidget(self.lbl_contours, 7, 8)
        grid.addWidget(self.qle_contours, 7, 9, 1, 6)
        #---------------------------------------------
        
        grid.addWidget(self.button, 15, 14, 1, 1)
        
        grid.addWidget(self.canvas, 2, 0, 8, 8)
        self.setLayout(grid) 

    def line1_editing_finished(self):
        if self.qle_line1.isModified():
            self.line_ratio.line1['code'] = str(self.qle_line1.text())
        self.qle_line1.setModified(False)

    def type1_editing_finished(self):
        if self.qle_type1.isModified():
            self.line_ratio.line1['type'] = str(self.qle_type1.text())
        self.qle_type1.setModified(False)

    def line2_editing_finished(self):
        if self.qle_line2.isModified():
            self.line_ratio.line2['code'] = str(self.qle_line2.text())
        self.qle_line2.setModified(False)

    def type2_editing_finished(self):
        if self.qle_type2.isModified():
            self.line_ratio.line2['type'] = str(self.qle_type2.text())
        self.qle_type2.setModified(False)

    def Av_editing_finished(self):
        if self.qle_Av.isModified():
            value = numpy.float64(eval(str(self.qle_Av.text())))
            self.line_ratio.Av_use = value
        self.qle_type2.setModified(False)

    def g_mech_editing_finished(self):
        if self.qle_g_mech.isModified():
            value = numpy.float64(eval(str(self.qle_g_mech.text())))            
            self.line_ratio.z_sec = value
        self.qle_g_mech.setModified(False)

    def set_contours_editing_finished(self):
        if self.qle_contours.isModified():
            value = numpy.array(eval('[' + str(self.qle_contours.text()) + ']'), 'f8')            
            self.line_ratio.c_levels['values'] = value
        self.qle_contours.setModified(False)
    
    def on_button_down_in_axes(self, event):
        
        # get the x and y coords, flip y from top to bottom
        xd, yd = event.xdata, event.ydata

        if self.line_ratio.axs.contains(event)[0]:

            #left click (print the value at that point and add the contour line)            
            v = self.line_ratio.get_value_at_coords(xd, yd)
            
            if self.line_ratio.cmap['log'] == True:
                vlog = v
                v = 10.0**vlog
                self.line_ratio.c_levels['values'] = numpy.hstack((vlog, self.line_ratio.c_levels['values'])) 
            else:
                vlog = numpy.log10(v)
                v = v
                self.line_ratio.c_levels['values'] = numpy.hstack((v, self.line_ratio.c_levels['values']))

            #updating the contours of the contour plot                
            self.line_ratio.update_contours()
            self.canvas.draw()
              
            print 'line ratio at log10(n)=%.1f and log10(y)=%.1f is %f ( = 10^{%.2f} ) ' % (xd, yd, v, vlog)        
            
            
            #setting the contour values to the countour line edit object
            text = ''
            for value in self.line_ratio.c_levels['values']:
                text += '%.2f, ' % value
            self.qle_contours.setText(text)
 
    def plot(self):
        print 'updating the plot'
        self.line_ratio.update_plot()
        self.canvas.draw()

class bar_plot_line_ratios(object):
    
    def __init__(
                 self, 
                 arxvPDR = None,
                 ranges = None,
                 names = None,
                 names_pos = None,
                 grid_coords = None,
                 gm_v = numpy.array([0.1, 1.0, 5.0, 10.0, 50.0])/100.0,
                 colors = [            'k', 'g', 'b', 'c' , 'y',   'r'],
                 line_ratios = None,
                 Av_max = None,
                 bar_width = 0.1,
                 image_save_dir=None,
                 plot_title = None,
                ):
        
        self.arxvPDR = arxvPDR
        self.gm_v   =  gm_v
        self.colors =  colors
        self.line_ratios = line_ratios
        self.Av_max = Av_max
        self.bar_width = bar_width
        self.grid_coords = grid_coords
        self.image_save_dir = image_save_dir
        self.plot_title = plot_title
        self.names_pos = names_pos
        
        self.n_subplots = len(names)
        
        fig, axs = pylab.subplots(self.n_subplots, 1, sharex = True, sharey = False, figsize = (6,12))
        pylab.subplots_adjust(left = 0.15, bottom = 0.1, right = 0.98, top = 0.9,
                             wspace = 0.0, hspace = 0.03) 
        
        #####################################################################################
        for i in numpy.arange(self.n_subplots):
            
            pylab.subplot(self.n_subplots, 1, self.n_subplots-i)
            info = self.plot_ratios_bars(arxvPDR, ranges[i], names[i], log_n = grid_coords[i][0], log_G0 = grid_coords[i][1])

            #setting the tick locations
            pylab.gca().set_xticks(numpy.arange(len(self.line_ratios)))
            
            if i == 0: 
                #setting the tick labels 
                pylab.gca().set_xticklabels(info['ratios'].keys(), rotation = 45, fontsize = 10)
            else:
                mylib.utils.removeAxesLabels.removeAll_xLabels(pylab.gca())
            
            #removing the last tick on the y axis for each subplot
            yticks = pylab.gca().axes.yaxis.get_major_ticks()
            yticks[-1].label1On = False

        legend = pylab.legend(info['rects'], info['strings'], 
                             bbox_to_anchor = (-0.0, 1.1, 1.0, .102), loc = 3,  
                             ncol=5, mode = 'expand', borderaxespad=0.0,
                             title = r"$\alpha$", fontsize=10)

        fig.canvas.set_window_title(self.plot_title)

        pylab.show()

        if self.image_save_dir != None:
            image_save_path = os.path.join(self.image_save_dir, self.plot_title + '.eps') 
            fig.savefig(image_save_path)

        
    
    def get_intensities_and_ratios(self, idx, pdrMeshObj):
        '''gets the intensities and the ratios for a mesh with the specified index in the pdr arxv. A
        tuple is returned where the first item is a dict object of lines from lineDict involved in the ratios 
        with an extra key 'flux' holding the flux of each line. The second item is a dict of the line ratios'''
        
        def get_lines_info():
            '''Gets the line info from the line ratio codes'''
    
            lines = dict()
            
            for line_ratio in self.line_ratios:
                
                #getting the species string from the line ratio fron the precdefined line codes
                line1_code, line2_code = line_ratio.split('/')
                line1 = lineDict.lines[line1_code]
                line2 = lineDict.lines[line2_code]
                lines[line1_code] = line1
                lines[line2_code] = line2
    
            return lines
        
        def read_radex_dbs_of_involves_lines(lines):
            '''reads the corresponding radex DBs if neccessary'''
                    
            #reading the dbs of the radex lines only
            for line in lines.values():
                if line['type'] == 'radex-lvg':
                    self.arxvPDR.readDbsRadex(species = [line['specStr']], Av = self.Av_max) 
        #
    
        def computeRatio(fluxes, rStr):
            xStr, yStr = rStr.split('/')
            ratio = fluxes[xStr]/fluxes[yStr]
            return ratio
        #
        
        def get_fluxes(lines):
            '''gets the flux of each line from the databases'''
            
            for line in lines.values():
                
                #flux of a radex line from a radex database
                if line['type'] == 'radex-lvg':
                    transitions = self.arxvPDR.radexDbs['%.2f' % self.Av_max][line['specStr']]['meshes'][idx]
                    if transitions == None:
                        line['flux'] = None
                    else:
                        line['flux'] = transitions[line['radexIdx']]['fluxcgs'] 
                
                #flux computed from a pdr mesh
                if line['type'] == 'pdr':
                    line['flux'] = pdrMeshObj.compute_integrated_quantity(line['ismcpak'], Av_range = [0.0, self.Av_max])
        #
        
        def get_line_ratios_values(lines):
            '''returns the values of the specified line ratios'''
            ratios = collections.OrderedDict()
            
            for line_ratio in self.line_ratios:
                
                #getting the species string from the line ratio fron the precdefined line codes
                line1_code, line2_code = line_ratio.split('/')
                
                flux1 = lines[line1_code]['flux']
                flux2 = lines[line2_code]['flux']
                
                if flux1 == None or flux2 == None:
                    ratios[line_ratio] = None
                else: 
                    ratios[line_ratio] = flux1 / flux2
                
            return ratios
        #
        
        lines = get_lines_info()
        
        read_radex_dbs_of_involves_lines(lines)
        
        get_fluxes(lines)
    
        line_ratio_values = get_line_ratios_values(lines)    
            
        return (lines, line_ratio_values)

    def plot_ratios_bars(self, arxvPDR, ylim, modelName, log_n = None, log_G0 = None):
        
        barWidth = self.bar_width
        gm_v = self.gm_v
        
        pylab.ylim( ylim )
        pylab.xlim( [0, len(self.line_ratios)] )
        
        allRects = []
        legendStrs = []
        
        for i, gm in enumerate(gm_v):
            
            print 'generating plots for model %s n, G0, gm = ' % modelName, log_n, log_G0, numpy.log10(gm)
            idx        = arxvPDR.get_mesh_index(x=log_n, y=log_G0, z=numpy.log10(gm) )
            pdrMeshObj = arxvPDR.get_mesh_data(x=log_n , y=log_G0, z=numpy.log10(gm) )
    
            lines, line_ratio_values = self.get_intensities_and_ratios(idx, pdrMeshObj)

            inds   = numpy.arange(len(line_ratio_values))
            values = numpy.log10(numpy.array(line_ratio_values.values(),'f8'))
    
            rect = pylab.bar(inds + i*barWidth, values, width = barWidth, bottom = 0, color = self.colors[i])
            allRects.append(rect)
            legendStrs.append(gm)
            print '      success'
        
        #pylab.yscale('log')
        pylab.text(self.names_pos[0], ylim[0] + (ylim[1] - ylim[0])*self.names_pos[1], modelName)        
        pylab.ylabel(r"$log_{10}$[line ratio]")
        pylab.grid(True)
        
        return {'rects' : allRects, 'strings' : legendStrs, 'ratios' : line_ratio_values}
    
class line_ratio_grid_of_grid():
    
    def __init__(self,
                 arxvPDR,       #: an PDR mesh database object
                 ratios=None,   #: a dict holding the ratios to be plotted (y-axis of the grid of the grids). for ex, ratios = ('CO2-1/CO1-0', 'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0')
                 zsecs=None,   #: a dict holding the section ins mechanical heating to be displayed (as the x-axis of the grid of grids). For ex, zsecs= (1e-10, 0.01, 0.05, 0.1, 0.25)
                 zsecType='alpha', #: a string specifying what is on the zaxis. choices are  'alpha' or 'CR'  
                 ranges=[[0.0, 6.0], [0.0, 6.0]],  #ranges in n and G0 for each grid
                 cmap={#: a dict holding the colormap, its range and whether the colors correpond to log quantities or not
                       'obj'     : matplotlib.cm.jet,
                       'v_range' : [-2, 2],
                       'log'     : True,
                      },
                 cbar={#: a dict holding the info about the ticks of the colorbar and their format
                       'ticks'  : numpy.linspace(-2, 2, 5),
                       'format' : '%.1f',  # None, '10^x', '%e' (or any format)
                      },
                 c_levels={#: a dict holding the info about the countour lines their format
                           'values' : numpy.log10([0.01, 0.03,  0.1,  0.3, 1,  3,  10,    30,   100]),
                           'format' : '%.1f',  # None, '10^x', '%e' (or any format), 'strs'
                           #'strs'   : ['0.01',  '0.1', '1', '10', '30', '100'],
                          }, 
                 Av_use=10.0, #: the Av of the grid
                 em_unit = 'cgs',
                 res=[100,100],  #: the resolution of the interpolated image
                 image_save_dir=None, 
                 removeNans=None, 
                 clip=None, 
                 clip_max_n=None, 
                 f_interp_dim=None, 
                 only_init=None, 
                 fname=None,                 
                 **kwargs):

        #setting the parameters and keywords as attributes
        self.arxvPDR = arxvPDR
        self.ratios = ratios
        self.zsecs = zsecs
        self.zsecType = zsecType
        self.Av_use = Av_use
        self.ranges = ranges
        self.em_unit = em_unit
        self.cmap = cmap
        self.cbar = cbar
        self.c_levels = c_levels
        self.res = res
        self.clip = clip
        self.clip_max = clip_max_n
        self.f_interp_dim = f_interp_dim
        self.image_save_dir = image_save_dir
        self.fname = fname
        
        #attributes related to the gui and where things will be plotted
        self.nx = None
        self.ny = None
        self.figure  = None
        self.axs     = None
        self.cbar_ax = None
        
        self.line_ratio_objs = None #:holds the grids of grids of the line ratios
        self.obs_ratios = None #: holds the observed ratios (set using set_obs_ratios method) 
        
        self.im = None
        self.cbar_im = None
        self.contours = None
        
        #defining and setting the values of the above attributes
        self.figure, self.axs, self.cbar_ax = self.init_plot()

        self.line_ratio_objs = numpy.ndarray(self.axs.shape, dtype='object') 
        
        #getting the data for all the subplots
        self.get_data_all()

        #plot the data in the suplots
        self.plot_grids()

        if self.fname != None and self.image_save_dir != None:
            imageSavePath = os.path.join(self.image_save_dir, self.fname) 
            self.figure.savefig(imageSavePath)
            print 'wrote image : %s' % imageSavePath
        
    def init_plot(self):
        '''initialized the plotting reigon, the colorbar and the axes'''

        self.nx = len(self.zsecs)
        self.ny = len(self.ratios)
        
        width  = 6    #figure width (non normalized) 
        height = width*(float(self.ny)/float(self.nx)) #figure heigh (non normalized)
        
        ax_xs = 0.12 #axses x start (normalized)
        ax_ys = 0.12 #axses y start (normalized)
        ax_xe = 0.72  #axses x end   (normalized)
        ax_ye = 0.8  #axses y end   (normalized)
        
        cbar_xs = ax_xs        #colorbar x start
        cbar_ys = ax_ye + 0.1  #colorbar y start
        cbar_sc = 0.8          #scale of the width of the cbar (relative to the width of ax)
        cbar_w  = 0.04         #width of the cbar (normalized)

        ax_xsz = ax_xe - ax_xs #width of the x-axes
        ax_ysz = ax_ye - ax_ys #height of the y-axes
        subplt_xsz = ax_xsz/self.nx #width of the subplots (all the subplots have this width)
        subplt_ysz = ax_ysz/self.ny #height of the subplots (all the subplots have this height)

        self.ax_xs  = ax_xs
        self.ax_ys  = ax_ys
        self.ax_xe  = ax_xe
        self.ax_ye  = ax_ye
        self.ax_xsz = ax_xsz
        self.ax_ysz = ax_ysz 
        self.subplt_xsz = subplt_xsz
        self.subplt_ysz = subplt_ysz
        
        fig, axs = pylab.subplots(self.ny, self.nx, sharex = False, sharey = False, figsize=(width, height),
                      subplot_kw = {'xlim':[self.ranges[0][0], self.ranges[0][1] ],  #keywords passed to figure.add_subplot() 
                                    'ylim':[self.ranges[1][0], self.ranges[1][1] ],
                                    'aspect':'equal',
                                    'adjustable':'datalim',
                                    }
                      )
        
        #plotting the axes labels and some other labels
        pylab.figtext(ax_xsz/2.0 + ax_xs,            0.02    , r'$\log_{10}$ [$n_{gas}$]', rotation=0 , horizontalalignment='center', verticalalignment='bottom')
        pylab.figtext(0.01              , ax_ysz/2.0 + ax_xs , r'$\log_{10}$ [$G_0$]'    , rotation=90, horizontalalignment='left', verticalalignment='center')
        pylab.figtext(ax_xe - 0.1,  0.9, r'A$_V$' + '\n %.2f' % self.Av_use)                                            
                                                   
        pylab.subplots_adjust(left=ax_xs, bottom=ax_ys, right=ax_xe, top=ax_ye, wspace=0.0, hspace=0.0)
        
        #plottuing the colorbar
        def plot_colorbar():
            
            dummy_cbar_data = numpy.zeros(numpy.array((50, 500), 'i'))
             
            cbar_ax = fig.add_axes([cbar_xs + (0.5*(1.0-cbar_sc))*ax_xsz, cbar_ys, cbar_sc*ax_xsz - (0.5*(1.0-cbar_sc))*ax_xsz, cbar_w])
            cbar_ax.set_title(r'$\log_{10}$[line ratio]')
            cbar_im = cbar_ax.imshow(dummy_cbar_data, aspect='auto', vmin=0, vmax=1, 
                                     cmap=self.cmap['obj'], 
                                     extent = [self.cmap['v_range'][0], self.cmap['v_range'][1], 0.0, 1.0])
            
            cbar_ax.axes.get_yaxis().set_ticks([]) #removing y ticklabels
    
            cbarData = numpy.linspace(0, 1, 500)
            cbarv = []
            for i in numpy.arange(50):
                cbarv.append( cbarData.copy() )
            cbarv = numpy.array(cbarv)
            
            cbar_im.set_data(cbarv)
    
            cbar_ax.set_xticks( self.cbar['ticks'] )
            
            #emphasising the contour line styles on the colorbar
            for i, v in enumerate(self.c_levels['values']):
                if i % 2 == 0:
                    cbar_ax.plot([v,v], [0,1], linestyle='solid', linewidth=2, color='k')
                else:
                    cbar_ax.plot([v,v], [0,1], linestyle='dashed', linewidth=2, color='k')
            
            #setting the labels of the colorbar    
            cbarLabeslStrs = []
            
            for tickv in self.cbar['ticks']:
                if self.cbar['format'] == '10^x':
                    cbarLabeslStrs.append(r'$10^{%.1f}$' % tickv)
                elif self.cbar['format'] == None:
                    cbarLabeslStrs.append(tickv)
                else:
                    cbarLabeslStrs.append(self.cbar['format'] % tickv)
            
            cbar_ax.set_xticklabels( cbarLabeslStrs )
        #
        
        cbar_ax = plot_colorbar()
        
        self.figure = fig
        self.axs = axs
        self.cbar_ax = cbar_ax
        
        self.figure.set_facecolor('white')
        
        return fig, axs, cbar_ax
    
    def get_data_all(self):
        '''gets and sets all the data for all the grids'''
        
        for r, row in enumerate(self.axs):
            for c, sub_ax in enumerate(row):
                
                ratio = self.ratios[r]
                z_sec = numpy.log10(self.zsecs[c])
                code1, code2 = ratio.split('/') 
                lr = line_ratio(self.arxvPDR, 
                                line1={'code':code1, 'type':'radex-lvg', 'em_unit':self.em_unit},
                                line2={'code':code2, 'type':'radex-lvg', 'em_unit':self.em_unit},
                                z_sec=z_sec,
                                Av_use=self.Av_use,
                                res=self.res,
                                cmap=self.cmap,
                                ranges=self.ranges,
                                only_get_data=True,
                                )
                
                self.line_ratio_objs[r][c] = lr
            #
        #
        
        mylib.utils.removeAxesLabels.removeSubplotLabels(self.axs, keep_first_last_labels=True)
        
    def plot_grids(self):
        '''plot the images in the subplots'''

        rangesLst = (self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1])

        for r, row in enumerate(self.axs):
            for c, sub_ax in enumerate(row):
                
                grd = self.line_ratio_objs[r][c].grd
                grd_scaled = self.line_ratio_objs[r][c].grd_scaled
                
                if rangesLst[0] == 0.0 and rangesLst[2] == 0.0:
                    sub_ax.set_xticks([1,3,5])
                    sub_ax.set_yticks([1,3,5])

                if rangesLst[0] == 3.0 and rangesLst[2] == 3.0:                
                    sub_ax.set_xticks([3, 4, 5])
                    sub_ax.set_yticks([3, 4, 5])
                
                #plotting the image
                sub_ax.imshow(grd_scaled.T, extent=rangesLst, origin='lower',    
                              cmap=self.cmap['obj'], vmin=0.0, vmax=1.0, norm=None)
                
                #plotting the contour lines
                clines = sub_ax.contour(grd.T, self.c_levels['values'], 
                                        extent=(self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1]), 
                                        origin='lower', 
                                        colors = 'black')
                
                #making every other contour line dashed
                for cline in clines.collections:
                    pylab.setp(cline, linestyle='solid', color='k')

                #making every other contour line dashed
                for cline in clines.collections[1::2]:
                    pylab.setp(cline, linestyle='dashed')
                     
                #displaying the title at the top of the subplots in the top row 
                if r == 0:
                    
                    if self.zsecType == 'alpha':
                        title_str = r'  $\alpha$' + '\n' + '%.2f' % self.zsecs[c] #alpha string
                    if self.zsecType == 'zeta':
                        title_str = r'  $\zeta$' + '\n' + '%.1e' % self.zsecs[c] #alpha string
                        
                    sub_ax.set_title(title_str, 
                                 size='small', rotation=0, 
                                 horizontalalignment='center', 
                                 verticalalignment='bottom')

                #printing the label for each row (line1/line2) 
                if c+1==self.nx:
                    text_x = self.ax_xe + 0.05
                    text_y = self.ax_ye - self.subplt_ysz*(r+0.5)
                    pylab.figtext(text_x, text_y, 
                                  self.ratios[r],
                                  rotation='horizontal', 
                                  size='small',
                                  horizontalalignment='left', verticalalignment='center'
                                 )
    

    def set_obs_ratios(self, obs_ratios):
        '''sets the obseravations object to the attribute and checks if all the line ratios 
         grids are available in self.line_ratio_objs
        ''' 
        
        #array holding the boolean as a check if the obs_ratios are in the 
        #self.line_ratio_objs
        check = numpy.zeros(len(obs_ratios),'b') 
        
        for r, row in enumerate(self.axs):
            for c, sub_ax in enumerate(row):
                lr = self.line_ratio_objs[r][c]
                
                for i, ratio in enumerate(obs_ratios):
                    if ratio == lr.line1['code'] + '/' + lr.line2['code']:
                        check[i] = True
                        break
        
        if numpy.prod(check) == False:
            strng = 'line ratios below are not in the grid of grids: \n   '
            inds = numpy.where(check == False)
            for i, ind in enumerate(inds):
                strng += (obs_ratios.keys()[ind] + ' ')
            
            raise ValueError(strng)
        else:
            print 'all input ratio grids are available' 
        
        #all ok, setting the observerd ratios as an attribute
        self.obs_ratios = obs_ratios
        
    def fit_observations_visual(self, obs_ratios=None):
        
        if obs_ratios != None:
            self.set_obs_ratios(obs_ratios)
        else:
            obs_ratios=self.obs_ratios

        #number of columns
        n_col = self.axs.shape[1]
        
        #making the truth grids
        truth_grids = numpy.ndarray((n_col), 'object')
        
        #intializing the truth grids
        for i in numpy.arange(n_col): truth_grids[i] = numpy.ones((self.res[0],self.res[1]), 'b') 

        #getting the truth grids (product of the grids which satisfy the observed 
        #ranges in the grids)            
        for r, row in enumerate(self.axs):
            for c, sub_ax in enumerate(row):
                
                lr = self.line_ratio_objs[r][c]
                
                #getting the line ratio from the 
                for i, ratio in enumerate(obs_ratios):
                    
                    if ratio == '%s/%s' % (lr.line1['code'],lr.line2['code']):
                        v, e = obs_ratios[ratio]['v'], obs_ratios[ratio]['e']
                        truth_grids[c] *= ((lr.grd >= numpy.log10(v-e)) *  (lr.grd <= numpy.log10(v+e)))
            #
        #

        fig, panelAxs = pylab.subplots(1, n_col, sharex = True, sharey = True, figsize = (6, 1.5),
                                       subplot_kw = {'xlim':[self.ranges[0][0], self.ranges[0][1] ],  #keywords passed to figure.add_subplot()     
                                                     'ylim':[self.ranges[1][0], self.ranges[1][1] ],
                                                     'aspect':'equal',
                                                     'adjustable':'datalim',
                                                     'autoscale_on' : False,                                                     
                                                     }
                                       )

        pylab.subplots_adjust(left=0.2, bottom=0.25, right=0.8, top=0.7, wspace=0.0, hspace=0.0)        
        fig.set_facecolor('white')

        rangesLst = (self.ranges[0][0], self.ranges[0][1], self.ranges[1][0], self.ranges[1][1])
        
        for i in numpy.arange(n_col):
            
            panelAxs[i].imshow(truth_grids[i].T, origin='lower', extent = rangesLst)

            if rangesLst[0] == 0.0 and rangesLst[2] == 0.0:
                panelAxs[i].set_xticks([1,3,5])
                panelAxs[i].set_yticks([1,3,5])

            if rangesLst[0] == 3.0 and rangesLst[2] == 3.0:                
                panelAxs[i].set_xticks([3, 4, 5])
                panelAxs[i].set_yticks([3, 4, 5])
                    
        #making the string to be printed at the top of the plot
        '''
        strng = ''
        strng += 'Av = %.1f' % self.Av_use
        for r in obs_ratios:
            strng += ' %s' % r
        pylab.figtext(0.1, 0.9, strng)
        '''
        pylab.show()
        
        return fig, panelAxs