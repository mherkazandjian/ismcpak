from numpy import *
from time import *
import sys, os
if os.uname()[1] in ['ruineraa', 'particle3']:
    import matplotlib
    matplotlib.use('Qt4Agg')
import pylab as pyl
from meshUtils import *
from mesh import mesh

#########################################parameters##########################################################
home = '/home/mher'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
        
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plot'          : False, 
         'showGrids'     : False,
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
        }
#############################################################################################################

# reading the archive
print 'setting up the archive'
t0 = time()
arxv = meshArxv(readDb = True, **parms)
print 'time reading data %f' % (time() - t0)

# setting the x,y,z quantities to be used for ploting
arxv.set_grid_axes_quantity_values(relativeGmech = parms['relativeGmech']) 

def get_line_vs_Av(specStr, transition_indx, log10nDens, log10G0, gmFrac, Avs):
    
    parms['radex']['specStr'] = specStr

    #getting the pdr model
    arxv.get_mesh_data(x = log10nDens, y = log10G0, z = np.log10(gmFrac))
    
    flux = []
    
    #getting the emissions for different Avs
    for Av in Avs:
        
        print 'Av = %f' % Av
        arxv.computeAndSetRadexCurves(compute_only=True, Av_range = [0.0, Av])
        flux.append(arxv.radexObj.transitions['fluxcgs'][transition_indx])

        
    return numpy.array(flux)
    


spec1Str    = 'HCN'
trans1_indx = 0

spec2Str    = 'HNC'
trans2_indx = 0

log10nDens  = 5.0
log10G0     = 5.0
gmFrac      = 0.01    
Avs         = numpy.linspace(2.0, 30.0, 15.0)

pylab.semilogy([], [])
pylab.ylim([0.001, 1.0])
pylab.hold(True)

plots = []
titles = []

#for gmFrac in [0.01, 0.1, 0.25, 0.5, 1.0]:
for gmFrac in [0.01, 0.05, 0.1]:
    flux1 = get_line_vs_Av(spec1Str, trans1_indx, log10nDens, log10G0, gmFrac, Avs)
    flux2 = get_line_vs_Av(spec2Str, trans2_indx, log10nDens, log10G0, gmFrac, Avs)
    
    p, = pylab.semilogy(Avs, flux2/flux1, '-')
    
    plots.append(p)
    titles.append( '%.1f' % (gmFrac*100) + '%')

pylab.legend(plots, titles)
pylab.show()


print 'done'