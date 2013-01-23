# purpose : - produce bar plots for the line emissions for different lines, column
#             column densities...etc..
#           - here focusing on atomic lines
# keywords: plot, line ratio, intensity, atomic, fine structure
#--------------------------------------------------------------------------------

from numpy import *
from time import *
import sys, os
if 'particle3' in os.uname():
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
from meshUtils import *
import mylib.utils.removeAxesLabels as axisUtils 

#########################################parameters##########################################################
home = '/home/mher'

metallicity = 2.0

#imageSavePath = '/home/mher/ism/docs/paper02/src/figs/bar-plots-lineRatios-atomic-z-%.1f.eps' % metallicity
imageSavePath = '/home/mher/foo.eps'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/singleModels-z-%.1f/' % metallicity,
         
         # reference database
         'runDirPath2'   : home + '/ism/runs/oneSided/surfaceGrid-z-%.1f-high-res-no-gmech/' % metallicity,
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'plotRanges'    : [[0,6],[0,6  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[0,6],[0,6],[-51, -15]],  # uniform gmech
         'metallicity'    : metallicity,

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
                             #'collisionPartners'    : ['H2','H','H+','e-'],
                             #'collisionPartners'    : ['H2'],
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
         #-----------------chemical network parameters------------------------
         'chemistry'     : {
                            'rxnFile'       : home + '/ism/code/ismcpak/data/rate99Fixed.inp',
                            'specNumFile'   : home + '/ism/code/ismcpak/data/species.inp',
                            'underAbunFile' : home + '/ism/code/ismcpak/data/underabundant.inp',
                            'removeManual'  : ['13CH3'],
                            'baseSpecies'   : 'baseSpeciesDefault', #name of the module holding the base species
                            'umistVer'      : 'umist99',
                           }
        }
##############################################setting up and reading the data###############################################################

# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv( **parms )
arxv.readDb( check = True)
print 'time reading data %f' % (time() - t0)

#setting up a mesh object for temporary use (utility stuff)
pdrMeshObj = mesh(chemNet = arxv.chemNet, metallicity = arxv.metallicity) 

# read and setting up the chemical network used in the 
t0 = time()
arxv.setupChemistry()
print 'time setting up the chemistry %f' % (time() - t0)

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech         = parms['relativeGmech'], 
                                   referenceDatabasePath = parms['runDirPath2'] )

#reading all the available precomputed radex databases
###arxv.readDbsRadex(species = ['CO'])

############################extracting what we want and plotting####################################
###arxv.use_radexDb('CO')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_intensities_and_ratios(pdrMeshObj):
    #get the intensities from a model
    flux = {}
    
    quantity = ['fineStructureCoolingComponents','O','rate','1-0'] #OI 63um
    flux['O63'] = (1.0/(2.0*np.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C','rate','1-0'] # CI 609um
    flux['C609'] = (1.0/(2.0*np.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C','rate','2-1'] # CI 369um
    flux['C369'] = (1.0/(2.0*np.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    quantity = ['fineStructureCoolingComponents','C+','rate','1-0'] # CII 158um
    flux['C+158'] = (1.0/(2.0*np.pi))*pdrMeshObj.compute_integrated_quantity(quantity)
    
    ratios = {}
    
    ratios['O63/C369'] = flux['O63']/flux['C369']
    ratios['O63/C609'] = flux['O63']/flux['C609']
    ratios['O63/C+158'] = flux['O63']/flux['C+158']
    ratios['C369/C609'] = flux['C369']/flux['C609']
    ratios['C+158/C369'] = flux['C+158']/flux['C369']
    ratios['C+158/C609'] = flux['C+158']/flux['C609']

    return (flux, ratios)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def plot_ratios_bars(arxv, ylim, modelName, log_n = None, log_G0 = None):
    global barWidth, colors, gm_v

    pyl.ylim( ylim )

    allRects = []
    legendStrs = []
    
    for i, gm in enumerate(gm_v):
        
        #very low density models 
        data = arxv.get_mesh_data(x = log_n, y = log_G0, z = np.log10(gm) )
        pdrMeshObj.setData( data )
        fluxes, ratios = get_intensities_and_ratios(pdrMeshObj)
        
        inds   = np.arange(len(ratios))
        values = np.log10(ratios.values())
        
        rect = pyl.bar(inds + i*barWidth, values, width = barWidth, bottom = 0, color = colors[i])
        allRects.append(rect)
        legendStrs.append(gm)
    
    #pyl.yscale('log')
    pyl.text(0.8, ylim[1]*0.8, modelName)        
    pyl.ylabel(r"$log_{10}$[line ratio]")
    pyl.grid(True)
    
    return {'rects' : allRects, 'strings' : legendStrs, 'ratios' : ratios}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#fig, axs = pyl.figure(figsize = (6,12))
fig, axs = pyl.subplots(6, 1, sharex = True, sharey = False, figsize = (6,12))
#ax = fig.add_subplot(616)
pyl.subplot(616) 
#ax.set_position([0.1, 0.2, 0.8, 0.7])
barWidth = 0.13

gm_v   = np.array([0.1, 1.0, 5.0, 10.0, 50.0])/100.0
colors = [                  'k',  'g', 'b', 'c', 'y',   'r']

#####################################################################################
info = plot_ratios_bars(arxv, [-0.5, 1.0], 'MA1', log_n = 1.0, log_G0 = 1.0)
pyl.gca().set_xticklabels(info['ratios'].keys(), rotation = 45, fontsize = 10)

pyl.subplot(615)
axisUtils.removeAll_xLabels(pyl.gca())
info = plot_ratios_bars(arxv, [-0.5, 2.0], 'MA2', log_n = 2.0, log_G0 = 2.0)

pyl.subplot(614)
axisUtils.removeAll_xLabels(pyl.gca())
info = plot_ratios_bars(arxv, [0.0, 3.0], 'M1', log_n = 3.0, log_G0 = 3.0)

pyl.subplot(613)
axisUtils.removeAll_xLabels(pyl.gca())
info = plot_ratios_bars(arxv, [0.0, 4.0], 'M2', log_n = 3.0, log_G0 = 5.0)

pyl.subplot(612)
axisUtils.removeAll_xLabels(pyl.gca())
info = plot_ratios_bars(arxv, [0.0, 5.0], 'M3', log_n = 5.5, log_G0 = 3.0)

pyl.subplot(611)
axisUtils.removeAll_xLabels(pyl.gca())
info = plot_ratios_bars(arxv, [0.0, 5.0], 'M4', log_n = 5.5, log_G0 = 5.0)


legen = pyl.legend(info['rects'], info['strings'], 
                   bbox_to_anchor = (-0.1, 1.1, 1.2, .102), loc = 3,  
                   ncol=5, mode = 'expand', borderaxespad=0.0,
                   title = r"$\alpha$")


pyl.show()

fig.savefig(imageSavePath)

print 'done'