# purpose : - produce bar plots for the line emissions for different lines, column
#             column densities...etc..
#           - here focusing on CO lines
# keywords: plot, line ratio, intensity, molecular, CO, fine structure
#--------------------------------------------------------------------------------

import numpy
import sys, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import meshUtils
import mesh
import mylib.utils.removeAxesLabels as axisUtils 
import collections

#########################################parameters##########################################################
home = '/home/mher'

metallicity = 2.0
Av_max      = 29.5

spec1Str = 'HNC'
spec2Str = 'HCN'

#imageSavePath = '/home/mher/ism/docs/paper02/src/figs/bar-plots-lineRatios-%s-%s-z-%.1f.eps' % (spec1Str, spec2Str, metallicity)
imageSavePath = '/home/mher/foo.eps'

parms = {
         #path to the database files
         'dirPath'       : home + '/ism/runs/oneSided/singleModels-z-%.1f/' % metallicity,
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
         'plotRanges'    : [[0,6],[0,6  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[0,6],[0,6],[-51, -15]],  # uniform gmech
         
         'plotGrids'     : False,
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : False,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  
                             'path'                 : home + '/ism/code/radex/Radex/bin/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : 'CO',
                             'freqRange'            : [0, 50000],
                             #'xH2_Min'              : 2*0.0000000001
                             'xH2_Min'              : -1.0,
                             'collisionPartners'    : ['H2','H+','H','e-','He'],
                             'tBack'                : 2.73,
                             'lineWidth'            : 1.0,
                             'verbose'              : True, 
                             'maxDisplayTranistion' : 20,
                             ###----------extra convergence params-----------------------
                             'checkOutputIntegrity' : False,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 1e-2,
                             #'popDensSumTol'        : 10,
                             'changeFracTrial'      : 0.001,
                             'nMaxTrial'            : 100,
                            },
        }
##############################################setting up and reading the data###############################################################

# reading the archive
arxv = meshUtils.meshArxv(readDb = True, **parms) 

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech = parms['relativeGmech']) 

#reading all the available precomputed radex databases
arxv.readDbsRadex(species = ['CO','13CO','HCN','HNC','HCO+','CS','CN'], Av = Av_max)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_intensities_and_ratios(idx):
    
    def computeRatio(fluxes, rStr):
        xStr, yStr = rStr.split('/')
        ratio = fluxes[xStr]/fluxes[yStr]
        return ratio
        
    #get the intensities from a model
    flux = collections.OrderedDict()

    specStr = spec1Str
    transitions = arxv.radexDbs['%.2f' % Av_max][specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:    
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[2]['fluxcgs'] 

    specStr = spec2Str
    transitions = arxv.radexDbs['%.2f' % Av_max][specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[2]['fluxcgs'] 

    ratios = collections.OrderedDict()

    rStr = spec1Str + '(4-3)/' + spec1Str + '(1-0)'; ratios[rStr] = computeRatio(flux, rStr)    
    rStr = spec2Str + '(4-3)/' + spec2Str + '(1-0)'; ratios[rStr] = computeRatio(flux, rStr)
    rStr = spec1Str + '(1-0)/' + spec2Str + '(1-0)'; ratios[rStr] = computeRatio(flux, rStr)
    rStr = spec1Str + '(4-3)/' + spec2Str + '(4-3)'; ratios[rStr] = computeRatio(flux, rStr)
    rStr = spec1Str + '(4-3)/' + spec2Str + '(1-0)'; ratios[rStr] = computeRatio(flux, rStr)
    
    return (flux, ratios)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_ratios_bars(arxv, ylim, modelName, log_n = None, log_G0 = None):
    global barWidth, colors, gm_v

    pylab.ylim( ylim )

    allRects = []
    legendStrs = []
    
    for i, gm in enumerate(gm_v):
        
        print 'generating plots for n, G0, gm = ', log_n, log_G0, numpy.log10(gm)
        idx = arxv.get_mesh_index(x = log_n, y = log_G0, z = numpy.log10(gm) )
        fluxes, ratios = get_intensities_and_ratios(idx)
        
        if ratios != None:
            inds   = numpy.arange(len(ratios))
            values = numpy.log10(ratios.values())

            rect = pylab.bar(inds + i*barWidth, values, width = barWidth, bottom = 0, color = colors[i])
            allRects.append(rect)
            legendStrs.append(gm)
            print '      success'
        else:
            print '      failed <-----------------'
            pass
    
    #pylab.yscale('log')
    pylab.text(0.8, ylim[1]*0.8, modelName)        
    pylab.ylabel(r"$log_{10}$[line ratio]")
    pylab.grid(True)
    
    return {'rects' : allRects, 'strings' : legendStrs, 'ratios' : ratios}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig, axs = pylab.subplots(6, 1, sharex = True, sharey = False, figsize = (6,12))
pylab.subplots_adjust(left = 0.15, bottom = 0.1, right = 0.98, top = 0.9,
                     wspace = 0.0, hspace = 0.0)
pylab.subplot(616) 
barWidth = 0.1

gm_v   = numpy.array([0.1, 1.0, 5.0, 10.0, 50.0])/100.0
colors = [            'k', 'g', 'b', 'c' , 'y',   'r']

#####################################################################################
info = plot_ratios_bars(arxv, [-1.0, 1.0], 'MA1', log_n = 1.0, log_G0 = 2.0)
pylab.xticks(range(len(info['ratios'].keys())))
pylab.gca().set_xticklabels(info['ratios'].keys(), rotation = 45, fontsize = 10)

pylab.subplot(615)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-0.6, 1.0], 'MA2', log_n = 2.0, log_G0 = 1.0)


pylab.subplot(614)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-1.0, 0.5], 'M1', log_n = 3.0, log_G0 = 3.0)

pylab.subplot(613)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-1.0, 0.5], 'M2', log_n = 3.0, log_G0 = 5.0)

pylab.subplot(612)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-0.5, 1.0], 'M3', log_n = 5.5, log_G0 = 3.0)

pylab.subplot(611)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-1.0, 1.0], 'M4', log_n = 5.5, log_G0 = 5.0)

legen = pylab.legend(info['rects'], info['strings'], 
                   bbox_to_anchor = (-0.1, 1.1, 1.1, .102), loc = 3, 
                   ncol=5, mode = 'expand', borderaxespad=0.0,
                   title = r"$\alpha$")

pylab.show()

fig.savefig(imageSavePath)
print 'saved image to path %s', imageSavePath
