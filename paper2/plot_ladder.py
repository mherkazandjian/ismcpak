# plot the ladder of transitions (tested for CO data) emission and optical depth
#------------------------------------------------------------------------------------

import numpy
import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')

import pylab
import meshUtils
#########################################parameters##########################################################
home = '/home/mher'

specStr_PDR   = 'CO'
specStr_Radex = 'CO'

imageSavePath = '/home/mher/ism/docs/paper02/src/figs/CO-ladder.eps'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
                                  # False => 3rd dim is gMech 
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plotRanges'    : [[-1,7],[-1,7  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[-4,7],[-4,7],[-51, -15]],  # uniform gmech
         
         'plot'          : True, 
         'showGrids'     : True,
         'gridsInfo'     : { '00' : {#some quantity
                                    'show'     : True,
                                    #'quantity' : ['state', 'gasT'],
                                    'quantity' : ['therm', 'heating'],
                                    'slabIdx'  : 0,
                                    },
                             '01' : {# abundance 
                                    'show'     : True, 
                                    'quantity' : ['state', 'abun'],
                                    'slabIdx'  : 0,
                                    'specStr'  : specStr_PDR,
                                    },
                             '10' : { # column density
                                    'show'     : True,
                                    'maxAv'    : 30,
                                    'specStr'  : specStr_PDR,
                                    },
                             '11' : { # line intensitities
                                     'show'           : True,
                                     ##------------------comment those if radex parms is 'pdr' is selected below this------------                                    
                                     #'type'           : 'pdr', #if type = pdr, quantity should point to a valid destination in the dtype in arxv.meshes[i]
                                     #'quantity'      : ['fineStructureCoolingComponents','C','rate','1-0'], # for use with 'pdr'
                                     #'specStr'        : 'C',     # database to be restored/computed
                                     ##-----------comment those radex parms if 'pdr' is selected above this--------------
                                     'type'           : 'radex',
                                     'specStr'        : specStr_Radex,     # database to be restored/computed
                                     'transitionIndx' : 0,
                                     'quantity'       : 'fluxcgs',
                                     #----------------end radex parms---------------------------------------------------
                                     'showContours'   : True,
                                     'Av_max'         : 10.0,  #the maximum Av to be used  
                                    },
                           },
         'gridsRes'      : 100,
         'nThreads'      : 1,
         'meshPltAvRng'  : [0, 30.0], #plotting range as a function of Av
          
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : False,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  #range which will be used in extracting data needed by radex from the PDR models
                                                                    #(only relevent to constructing databases)
                             'path'                 : home + '/ism/code/radex/Radex/bin-gcc/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : specStr_Radex,
                             'freqRange'            : [0, 50000],
                             #'xH2_Min'              : 2*0.0000000001
                             'xH2_Min'              : -1.0,
                             #'collisionPartners'    : ['H2','H+','H','e-','He'],
                             #'collisionPartners'    : ['H2','H','H+','e-'],
                             'collisionPartners'    : ['H2'],
                             'use_pdr_gas_den_H2'   : True,   #<----------
                             'tBack'                : 2.73,
                             'lineWidth'            : 1.0,
                             'verbose'              : False,
                             'maxDisplayTranistion' : 20,
                             ###----------extra convergence params-----------------------
                             'checkOutputIntegrity' : False,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 1e-2,
                             #'popDensSumTol'        : 10,
                             'changeFracTrial'      : 0.01,
                             'nMaxTrial'            : 100,
                            },
        }
#############################################################################################################

# reading the archive
print 'setting up the archive'
t0 = time.time()
arxv = meshUtils.meshArxv(readDb = True, **parms)
print 'time reading data %f' % (time.time() - t0)

arxv.set_attributes(**parms)
arxv.set_default_attributes()

"""
if parms['plot']:
    # plotting stuff
    pylab.ioff()
    arxv.plotGrids()
    pylab.show()
"""

print 'done'

fig, axs = pylab.subplots(1, 2, figsize = (12, 6))

 
pylab.subplot(121)

pylab.hold(True)
indx1=arxv.get_mesh_index(x=3.0, y=3.0, z=-10)
transRadex = arxv.meshesRadex[indx1]
p1, = pylab.semilogy(numpy.int32(transRadex[:]['upper']), transRadex[:]['fluxcgs'], 'r')

indx2=arxv.get_mesh_index(x=3.0, y=3.0, z=-3)
transRadex = arxv.meshesRadex[indx2]
p2, = pylab.semilogy(numpy.int32(transRadex[:]['upper']), transRadex[:]['fluxcgs'], 'b')

indx3=arxv.get_mesh_index(x=3.0, y=3.0, z=-2)
transRadex = arxv.meshesRadex[indx3]
p3, = pylab.semilogy(numpy.int32(transRadex[:]['upper']), transRadex[:]['fluxcgs'], 'g')

indx4=arxv.get_mesh_index(x=3.0, y=3.0, z=-1)
transRadex = arxv.meshesRadex[indx4]
p4, = pylab.semilogy(numpy.int32(transRadex[:]['upper']), transRadex[:]['fluxcgs'], 'k')

pylab.hold(False)

pylab.show()

pylab.axis([0, 20, 1e-10, 1])
pylab.legend([p1, p2, p3, p4], ['0%','0.1%','1%','10%'])
pylab.xlabel('J')
pylab.ylabel('$\log_{10} [ flux / erg.cm^{-2}.s^{-1} ]$')


pylab.subplot(122)

pylab.hold(True)
indx1=arxv.get_mesh_index(x=3.0, y=3.0, z=-10)
transRadex = arxv.meshesRadex[indx1]
p1, = pylab.plot(numpy.int32(transRadex[:]['upper']), transRadex[:]['tau'], 'r')

indx2=arxv.get_mesh_index(x=3.0, y=3.0, z=-3)
transRadex = arxv.meshesRadex[indx2]
p2, = pylab.plot(numpy.int32(transRadex[:]['upper']), transRadex[:]['tau'], 'b')

indx3=arxv.get_mesh_index(x=3.0, y=3.0, z=-2)
transRadex = arxv.meshesRadex[indx3]
p3, = pylab.plot(numpy.int32(transRadex[:]['upper']), transRadex[:]['tau'], 'g')

indx4=arxv.get_mesh_index(x=3.0, y=3.0, z=-1)
transRadex = arxv.meshesRadex[indx4]
p4, = pylab.plot(numpy.int32(transRadex[:]['upper']), transRadex[:]['tau'], 'k')

pylab.hold(False)

pylab.axis([0, 20, 0, 500])
pylab.legend([p1, p2, p3, p4], ['0%','0.1%','1%','10%'])
pylab.xlabel('J')
pylab.ylabel('$tau$')


fig.savefig(imageSavePath)
print 'wrote image : %s' % imageSavePath
pylab.show()

print 'done'
