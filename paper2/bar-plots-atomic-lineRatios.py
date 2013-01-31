# purpose : - produce bar plots for the line emissions for different lines, column
#             column densities...etc..
#           - here focusing on atomic lines
# keywords: plot, line ratio, intensity, atomic, fine structure
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

#########################################parameters##########################################################
home = '/home/mher'

metallicity = 0.5

imageSavePath = '/home/mher/ism/docs/paper02/src/figs/bar-plots-lineRatios-atomic-z-%.1f.eps' % metallicity
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

#setting up a mesh object for temporary use (utility stuff)
pdrMeshObj = mesh.mesh(chemNet = arxv.chemNet, metallicity = arxv.metallicity) 

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech = parms['relativeGmech']) 

#reading all the available precomputed radex databases
###arxv.readDbsRadex(species = ['CO'])

############################extracting what we want and plotting####################################
###arxv.use_radexDb('CO')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_intensities_and_ratios(pdrMeshObj):
    #get the intensities from a model
    flux = collections.OrderedDict()
    
    quantity = ['fineStructureCoolingComponents','O','rate','1-0'] #OI 63um
    flux['O63'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C','rate','1-0'] # CI 609um
    flux['C609'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C','rate','2-1'] # CI 369um
    flux['C369'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C+','rate','1-0'] # CII 158um
    flux['C+158'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    
    ratios = collections.OrderedDict()
    
    ratios['O63/C369'] = flux['O63']/flux['C369']
    ratios['O63/C609'] = flux['O63']/flux['C609']
    ratios['O63/C+158'] = flux['O63']/flux['C+158']
    ratios['C+158/C369'] = flux['C+158']/flux['C369']
    ratios['C+158/C609'] = flux['C+158']/flux['C609']
    ratios['C369/C609'] = flux['C369']/flux['C609']

    return (flux, ratios)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_ratios_bars(arxv, ylim, modelName, log_n = None, log_G0 = None):
    global barWidth, colors, gm_v

    pylab.ylim( ylim )

    allRects = []
    legendStrs = []
    
    for i, gm in enumerate(gm_v):
        
        data = arxv.get_mesh_data(x = log_n, y = log_G0, z = numpy.log10(gm) )
        pdrMeshObj.setData( data )
        fluxes, ratios = get_intensities_and_ratios(pdrMeshObj)
        
        inds   = numpy.arange(len(ratios))
        values = numpy.log10(ratios.values())
        
        rect = pylab.bar(inds + i*barWidth, values, width = barWidth, bottom = 0, color = colors[i])
        allRects.append(rect)
        legendStrs.append(gm)
    
    #pylab.yscale('log')
    pylab.text(0.8, ylim[1]*0.8, modelName)        
    pylab.ylabel(r"$log_{10}$[line ratio]")
    pylab.grid(True)
    
    return {'rects' : allRects, 'strings' : legendStrs, 'ratios' : ratios}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#fig, axs = pylab.figure(figsize = (6,12))
fig, axs = pylab.subplots(6, 1, sharex = True, sharey = False, figsize = (6,12))
#ax = fig.add_subplot(616)
pylab.subplot(616) 
#ax.set_position([0.1, 0.2, 0.8, 0.7])
barWidth = 0.13

gm_v   = numpy.array([0.1, 1.0, 5.0, 10.0, 50.0])/100.0
colors = [                  'k',  'g', 'b', 'c', 'y',   'r']

#####################################################################################
info = plot_ratios_bars(arxv, [-1.0, 1.0], 'MA1', log_n = 1.0, log_G0 = 1.0)
pylab.gca().set_xticklabels(info['ratios'].keys(), rotation = 45, fontsize = 10)

pylab.subplot(615)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-1.0, 2.0], 'MA2', log_n = 2.0, log_G0 = 2.0)

pylab.subplot(614)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [-0.5, 2.5], 'M1', log_n = 3.0, log_G0 = 3.0)

pylab.subplot(613)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [0.0, 3.0], 'M2', log_n = 3.0, log_G0 = 5.0)

pylab.subplot(612)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [0.0, 4.0], 'M3', log_n = 5.5, log_G0 = 3.0)

pylab.subplot(611)
axisUtils.removeAll_xLabels(pylab.gca())
info = plot_ratios_bars(arxv, [0.0, 5.0], 'M4', log_n = 5.5, log_G0 = 5.0)


legen = pylab.legend(info['rects'], info['strings'], 
                   bbox_to_anchor = (-0.1, 1.1, 1.2, .102), loc = 3,  
                   ncol=5, mode = 'expand', borderaxespad=0.0,
                   title = r"$\alpha$")


pylab.show()

fig.savefig(imageSavePath)

print 'done'