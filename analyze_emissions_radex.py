"""
see how emissions (or emission radiative properties) are effected as a function of different 
parameters using radex. For example, see how emission intensity changes with temperature, optical
depth, column density...etc..
"""

import radex
import numpy
import pylab

# path of the radex excutable
radexPath      = '/home/mher/ism/code/radex/Radex/bin-gcc/radex'  
molDataDirPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'

# parameters that will be passed to radex single partner
inFile = { 'specStr'                : 'HCN'      ,
           'freqRange'              : [0, 0]    ,
           'tKin'                   : 30.0     ,
           'collisionPartners'      : ['H2']    ,
           'nDensCollisionPartners' : [1e3]     ,
           'tBack'                  : 2.73      ,
           'molnDens'               : 1e14      ,
           'lineWidth'              : 1.0       ,
           'runAnother'             : 1         }

#transition info to be plotted
transIndx = 1

#-------x1 (quantity to be plotted wit curves on the same plot with different color--------
x1 = [10, 50, 200, 500]  #changing the kinetic temperature
qx1 = 'tKin' 

#x1 = [[1e3], [1e4], [1e5], [1e6]]   #changing the density of the colliding specie
#qx1 = 'nDensCollisionPartners' 

#-------quantity to be considered as the one on the x-axis---------------------------------
x2 = 10.0**numpy.linspace(12.0, 24.0, 20)
qx2 = 'molnDens'

plt_quantities = ['fluxcgs', 'fluxKkms', 'tau', 'Tex']

#-------------------------------------------------------------------------------------------
# creating and setup the radex process instance
radexObj = radex.radex(radexPath, molDataDirPath)

# setting put the parameters, running and parsing the output
radexObj.setInFile( inFile )


fig, axs = pylab.subplots(nrows=len(plt_quantities), ncols=1, sharex=True, figsize=(6,10))

plts, titles = [], [] #stuff for making the legend

#run radex for all the values x2 for each value of x1 and 
#plot the curves as a function of x2 for each x1.
for i, v1 in enumerate(x1):
    
    #array which will store the computed values
    y = numpy.zeros(x2.size, radexObj.generateTransitionDtype().descr)  

    inFile[qx1] = v1
    
    #plotting stuff varying the x axis variable
    for j, v2 in enumerate(x2):
        
        inFile[qx2] = v2
        
        #radexObj.run( checkInput = True, verbose = True)        #<<< running radex without error checks
        radexObj.run_mutiple_trials(1.0, 1e-2, 1e-2, 100, strict = True, verbose = False) #<<<< running radex with checks
        
        if radexObj.getStatus() &  radexObj.FLAGS['SUCCESS'] :
            radexObj.parseOutput()
            trans = radexObj.getTransition(transIndx)  # getting the info of the transition
    
            for q in plt_quantities:
                y[j][q] = trans[q]
            
            print 'radex converged %s = %e, %s = %e' % (qx1, v1, qx2, v2)
        
        else:
            if radexObj.getStatus() &  radexObj.FLAGS['ITERWARN']:
                print 'did not converge'
                
                print 'warnings'
                print '--------'
                print radexObj.warnings

                print 'radex did not converg %s = %e, %s = %e' % (qx1, v1, qx2, v2)
                
                y[j][q] = numpy.nan
        
    for k, q in enumerate(plt_quantities):
        plt, = axs[k].loglog(x2, y[q])
        
        if k == 0:
            plts.append(plt)
            titles.append(qx1 + ' = %e' % v1)
            
#setting up the labels
for k, q in enumerate(plt_quantities):
    axs[k].set_ylabel(q)
axs[-1].set_xlabel(qx2)
    
pylab.figtext(0.1, 0.95, 'tranistion = %d' % transIndx, size='xx-large')
fig.legend(plts, titles)

pylab.show()
