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

#specStr_PDR   = 'HCO+'
#specStr_Radex = 'HCO+'

imageSavePath = '/home/mher/foo.eps'

parms = {
         #path to the database files
         'dirPath'      : home + '/ism/runs/oneSided/dynamicMeshTest1/',
         'relativeGmech' : True,  # True  => 3rd dim is the gMech/gSurface(gMech=0)
         'min_gMech'     : 1e-50, # set the mimum value of gMech to be used in the ref arxive
         
         'plotRanges'    : [[-1,7],[-1,7  ],[-12, 6]],     # adaptive gMech 
         #'plotRanges'     : [[-3,7],[-3,7],[-51, -15]],  # uniform gmech
         
         'plot'          : False, 
         'showGrids'     : False,
         'gridsRes'      : 100,         
         'radex'         : { 'use'                  : True,
                             'loadAllDbs'           : True,
                             ###-----------radex database parms-----------------
                             'compute'              : False, #if true, runns radex on all meshes
                             'writeDb'              : False, #if true, writes the computed stuff to a db
                             'Av_range'             : [0.0, 10.0],  #range which will be used in extracting data needed by radex from the PDR models
                                                                    #(only relevent to constructing databases)
                             'path'                 : home + '/ism/code/radex/Radex/bin/radex',  
                             'molDataDirPath'       : home + '/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles',
                             'specStr'              : None,#specStr_Radex,
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
                             'checkOutputIntegrity' : True,  # if true, check the radex output (sometimes although it converges, the numbers do not make sense)                             
                             'popDensSumExpected'   : 1.0, 
                             'popDensSumTol'        : 1e-2,
                             #'popDensSumTol'        : 10,
                             'changeFracTrial'      : 0.01,
                             'nMaxTrial'            : 100,
                            },
        }
#############################################################################################################

# reading the archive
#--------------------
print 'setting up the archive'
t0 = time.time()
arxv = meshUtils.meshArxv(readDb = True, **parms)
print 'time reading data %f' % (time.time() - t0)

#----------
def ratio_grid_gm_Av(logn, logg0, spec1, trans1, spec2, trans2, arxv):
    """
    for each gm, getting the emissions as a function of Av
    for a model with a certain logn and G0 for two line with
    spec1(trans1)/spec2(trans2) and returning the ratio as a 
    2D array, gm_all, Av_all. For example:
    
    ratio, gms, Avs = ratio_grid_gm_Av(3.0, 3.0, 'CO', 1, 'CO', 0, arxv)
    
    The ratio is that from spec1(trans1)/spec2(trans2)
    
    ratio[i,j] is the ratio for gm = gm_all[i] and Av_all[j]
    """
    
    #getting the unique sections in gmech and Av
    gm_all = arxv.get_unique_grid_z_sections(1e-13)
    Av_all = numpy.array([float(AvStr) for AvStr in arxv.radexDbs.keys()])

    n_gm = gm_all.size
    n_Av = Av_all.size
    
    data1 = numpy.zeros((n_gm, n_Av), dtype=numpy.float64) 
    data2 = numpy.zeros((n_gm, n_Av), dtype=numpy.float64)
    
    for i, gm in enumerate(gm_all):        
        for j, Av in enumerate(Av_all):
            
            Av_this = Av_all[j]
            
            #index of the PDR mesh 
            indMin = arxv.get_mesh_index(x = logn, y = logg0, z = gm_all[i])
            
            arxv.use_radexDb(Av_this, spec1)
            if arxv.meshesRadex[indMin] != None:
                data1[i,j] = arxv.meshesRadex[indMin]['fluxcgs'][trans1]
            else:
                continue
    
            arxv.use_radexDb(Av_this, spec2)
            if arxv.meshesRadex[indMin] != None:
                data2[i,j] = arxv.meshesRadex[indMin]['fluxcgs'][trans2]
            else:
                continue
    
    #the ratio for all the selected models for the two lines
    return data1/data2, gm_all, Av_all
##

def plot_ratios(logn, logg0, ratio1, ratio2, arxv, xrng, yrng):
    r1, gms1, Avs1 = ratio_grid_gm_Av(logn, logg0, ratio1['spec1'], ratio1['trans1'], ratio1['spec2'], ratio1['trans2'], arxv)
    r2, gms2, Avs2 = ratio_grid_gm_Av(logn, logg0, ratio2['spec1'], ratio2['trans1'], ratio2['spec2'], ratio2['trans2'], arxv)
    
    fig = pylab.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    
    gm_plot = numpy.log10([1e-10, 0.01, 0.05, 0.1, 0.5, 1.0])
    colors  = ['k', 'g', 'b', 'c', 'y', 'r']
    titles  = []
    plts    = []
    
    pylab.hold(True)
     
    for i, gm in enumerate(gm_plot):
        
        ind_gm = numpy.where(numpy.fabs(1.0 - gm/gms1) < 1e-6)[0]
        
        if ind_gm.size != 0:
            
            print gms1[ind_gm]
            inds = numpy.argsort(Avs1)
            
            #plt1, = ax.plot(Avs1[inds], numpy.log10(r1[ind_gm, inds]), colors[i])
            #plt2, = ax.plot(Avs1[inds], numpy.log10(r2[ind_gm, inds]), colors[i]+'--')
            plt1, = ax.semilogy(Avs1[inds], r1[ind_gm, inds], colors[i])
            plt2, = ax.semilogy(Avs1[inds], r2[ind_gm, inds], colors[i]+'--')
            plts.append(plt1)
            titles.append('%.2f' % (10.0**gms1[ind_gm]))

    
    ax.legend(plts, titles)
    
    ax.set_ylim(yrng[0], yrng[1])
    ax.set_xlim(xrng[0], xrng[1])
    ax.set_xlabel(r'A$_V$')#, size='small')
    ax.set_ylabel('line ratio')#, size='small')

    #plotting the intersection lines
    """
    x1, y1 = (6.51, 0.254)
    ax.semilogy([5, x1], [y1, y1], 'r--')
    ax.semilogy([x1, x1], [y1, 0.0001], 'r--')
        

    x2, y2 = (23.1, 0.254)
    ax.semilogy([5, x2], [y2, y2], 'r--')
    ax.semilogy([x2, x2], [y1, 0.0001], 'r--')        

    ax.semilogy([x1], [y1], 'ko', markersize=12)
    ax.text(x1, 0.05, 'Av = %.1f' % x1, size=16)

    ax.semilogy([x2], [y2], 'ko', markersize=12)
    ax.text(x2, 0.05, 'Av = %.1f' % x2, size=16)
    """
    
    pylab.hold(False)
    
    pylab.show()

    fig.savefig(imageSavePath)
    
    return fig

fig = plot_ratios(3.0, 3.0, 
                  {'spec1':'CO'  , 'trans1':1, 'spec2':'CO', 'trans2':0},  
                  {'spec1':'13CO', 'trans1':0, 'spec2':'CO', 'trans2':0},
                  arxv,
                  [5.0,30.0],
                  [1e-2,100]
                  )
        

print 'done'