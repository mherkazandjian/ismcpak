# purpose : - produce bar plots for the line emissions for different lines, column
#             column densities...etc..
#           - for multiple species to be plotted on the same bar plot
# keywords: plot, line ratio, intensity, molecular, CO, fine structure
#--------------------------------------------------------------------------------

import numpy
import sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab
import meshUtils
import mesh
import mylib.utils.removeAxesLabels as axisUtils 
import collections

#-----------
#set the spcie and the lines to be extracted from the data in the sections that
#look line this one.. in the sections below
"""
    specStr = 'HCO+'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:    
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[3]['fluxcgs'] 
"""
#set the line ratios to be displayed in the section which looks like this below as well
"""
    rStr = 'HCO+(4-3)/HCO+(1-0)'; ratios[rStr] = computeRatio(flux, rStr)    
    rStr = 'CN(3-2)/CN(1-0)';     ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CS(4-3)/CS(1-0)';     ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CS(1-0)/HCO+(1-0)';   ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CS(4-3)/HCO+(4-3)';   ratios[rStr] = computeRatio(flux, rStr)
"""
#########################################parameters##########################################################
home = '/home/mher'

metallicity = 2.0

imageSavePath = '/home/mher/ism/docs/paper02/src/figs/bar-plots-lineRatios-misc1-z-%.1f.eps' % (metallicity)
#imageSavePath = '/home/mher/foo.eps'


parms = {
         #path to the database files
         'dirPath'       : home + '/ism/runs/oneSided/singleModels-z-%.1f/' % metallicity,
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
         'plotRanges'    : [[0,6],[0,6  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[0,6],[0,6],[-51, -15]],  # uniform gmech
         
         'plotGrids'     : False,
         'radex'         : { 'use'                  : True,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
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
arxv.readDbsRadex(species = ['CO','13CO','HCN','HNC','HCO+','CS','CN'])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_intensities_and_ratios(idx):
    
    def computeRatio(fluxes, rStr):
        xStr, yStr = rStr.split('/')
        ratio = fluxes[xStr]/fluxes[yStr]
        return ratio
        
    #get the intensities from a model
    flux = collections.OrderedDict()

    specStr = 'HCO+'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:    
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[3]['fluxcgs'] 

    specStr = 'HCN'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[2]['fluxcgs'] 

    specStr = 'CO'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(2-1)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(3-2)']  = transitions[2]['fluxcgs'] 

    specStr = '13CO'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(2-1)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(3-2)']  = transitions[2]['fluxcgs'] 

    specStr = 'CN'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(3-2)']  = transitions[1]['fluxcgs'] 

    specStr = 'CS'
    transitions = arxv.radexDbs[specStr]['meshes'][idx]
    if transitions == None:
        return (None, None)
    else:
        flux[specStr + '(1-0)']  = transitions[1]['fluxcgs'] 
        flux[specStr + '(4-3)']  = transitions[1]['fluxcgs'] 



    ratios = collections.OrderedDict()

    rStr = 'HCO+(1-0)/CO(1-0)';   ratios[rStr] = computeRatio(flux, rStr)    
    rStr = 'HCO+(1-0)/13CO(1-0)'; ratios[rStr] = computeRatio(flux, rStr)    
    rStr = 'HCN(1-0)/CO(1-0)';    ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'HCN(1-0)/HCO+(1-0)';  ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'HCN(4-3)/HCO+(4-3)';  ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CN(1-0)/HCN(1-0)';    ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CN(3-2)/HCN(1-0)';    ratios[rStr] = computeRatio(flux, rStr)
    rStr = 'CS(1-0)/CO(1-0)';     ratios[rStr] = computeRatio(flux, rStr)
    
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

#fig, axs = pylab.figure(figsize = (6,12))
fig, axs = pylab.subplots(4, 1, sharex = True, sharey = False, figsize = (8,10))
barWidth = 0.1

gm_v   = numpy.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])/100.0
colors = [            'k', 'g', 'b', 'c' , 'y',   'm',  'r']

#####################################################################################

pylab.subplot(414) 
info = plot_ratios_bars(arxv, [-2.0, 1.0], 'M1', log_n = 3.0, log_G0 = 3.0)
pylab.xticks(range(len(info['ratios'].keys())))
pylab.gca().set_xticklabels(info['ratios'].keys(), rotation = 45, fontsize = 8)

pylab.subplot(413)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-3.0, 1.5], 'M2', log_n = 3.0, log_G0 = 5.0)

pylab.subplot(412)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-2.0, 2.0], 'M3', log_n = 5.5, log_G0 = 3.0)

pylab.subplot(411)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-2.0, 2.0], 'M4', log_n = 5.5, log_G0 = 5.0)

legen = pylab.legend(info['rects'], info['strings'], 
                   bbox_to_anchor = (-0.1, 1.1, 1.2, .102), loc = 3,  
                   ncol=7, mode = 'expand', borderaxespad=0.0,
                   title = r"$\alpha$")
pylab.show()

fig.savefig(imageSavePath)
print 'saved image to path %s', imageSavePath
