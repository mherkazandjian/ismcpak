import numpy
from numpy import log10

import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')
#matplotlib.use('PS')
from mylib.utils.misc import fetchNestedDtypeValue

import pylab
import meshUtils
import lineDict

#########################################parameters##########################################################
home = '/home/mher'

parameters = {
             'dirPath'    : home + '/ism/runs/tests/oneSidedGrid/',
             'range_logn'  : [0.0, 6.0],    #range in log n of models to be considered
             'range_logG0' : [0.0, 6.0],    #range in log G0 of models to be considered
             'range_alpha' : [0.0, 1],      #range in alpha of models to be considered (not in log10)
             'Av_values'   : [0.01, 0.1, 1.0, 2.0], #values in Av to be considered
             #the emissions to be extracted from the LVG models
             'emissions'   : (
                              ['lvg', 'radex' , 'CO'  , '1-0'],   
                              ['lvg', 'radex' , 'CO'  , '2-1'],        
                              ['lvg', 'radex' , 'CO'  , '3-2'],
                              ['lvg', 'radex' , 'CO'  , '4-3'],   
                             #
                              ['lvg', 'radex' , '13CO', '1-0'],   
                              ['lvg', 'radex' , '13CO', '2-1'],   
                              ['lvg', 'radex' , '13CO', '3-2'],   
                              ['lvg', 'radex' , '13CO', '4-3'],   
                              #
                              ['lvg', 'radex' , 'HCO+', '1-0'],   
                              ['lvg', 'radex' , 'HCO+', '2-1'],   
                              ['lvg', 'radex' , 'HCO+', '3-2'],   
                              ['lvg', 'radex' , 'HCO+', '4-3'],
                              #
                              ['lvg', 'radex' , 'HCN', '1-0'],   
                              ['lvg', 'radex' , 'HCN', '2-1'],   
                              ['lvg', 'radex' , 'HCN', '3-2'],   
                              ['lvg', 'radex' , 'HCN', '4-3'],
                              #
                              ['pdr', 'leiden', 'C'   , '609'],
                              ['pdr', 'leiden', 'C'   , '369'],
                             ),
             #the column densities to be written
             'col_dense'   : (
                              ['pdr', 'CO'],   
                              ['pdr', 'H2'],
                             ),
             #integrated quantities to be written to the file
             'int_quantities' : (
                                 ['total_heating', ['therm', 'heating']],
                                 ['heating_CR', ['heating', 'cr']],
                                ), 
             #quantities at a certain Av
             'quantities_at_AV' : (
                                   ['T_kin', ['state', 'gasT']],
                                  ),
            # custom quantities at a certain Av
            'custom_quantities_at_Av' : (
                                         ['T_radex'],
                                        ),
             'outFname'    : home + '/tmp/models-z-1.0.txt',
             }

# 
def main(parms):
    
    ## reading the archive
    print 'setting up the archive'
    t0 = time.time()
    arxv = meshUtils.meshArxv(dirPath=parms['dirPath'], readDb=True, relativeGmech=False, min_gMech=1e-50)
    print 'time reading data %f' % (time.time() - t0)

    ##loading all the radex model info
    for line_model, code, specStr, trans  in parms['emissions']:
        if line_model == 'lvg':            
            arxv.readDbsRadex(Av = parms['Av_values'], species = specStr)

    fObj = open(parms['outFname'], 'w')

    ## writing the column numbers
    fObj.write('#c1    c2         c3         c4         c5      c6      c7    ')

    for i, p in enumerate(parms['emissions'] + 
                          parms['col_dense'] + 
                          parms['int_quantities'] +
                          parms['quantities_at_AV'] +
                          parms['custom_quantities_at_Av']):
        fObj.write('c%-15d' % (i+7+1))
    fObj.write('\n')
    
    ## writing the column quantities
    fObj.write('#indx  n_gas      G0         gmech      alpha   Z       Av    ')
    for i, p in enumerate(parms['emissions']):
        fObj.write('%-15s ' % (p[2]+p[3]))
    for i, p in enumerate(parms['col_dense']):
        fObj.write('%-15s ' % ('N(%s) ' % (p[1])))
    for i, q in enumerate(parms['int_quantities']):
        fObj.write('%-15s ' % q[0])
    for i, q in enumerate(parms['quantities_at_AV']):
        fObj.write('%-15s ' % q[0])
    for i, q in enumerate(parms['custom_quantities_at_Av']):
        fObj.write('%-15s ' % q[0])
    fObj.write('\n')

    idx = 0 #an index written at the beginning of each line
    
    ## looping over all the meshes and extracting the info needed and writing them to the ascii file
    for i, msh in enumerate(arxv.meshes):
        
        hdr = msh['hdr']
        n, G0, gm = hdr['nGas'], hdr['G0'], hdr['gammaMech']
        alpha = 10.0**arxv.grid_z[i]
        
        if log10(n) == 5.0 and log10(G0) == 5.0:
            print  log10(n), log10(G0), gm
        
        ## selecting models with specified ranges        
        if log10(n) < parms['range_logn'][0] or log10(n) > parms['range_logn'][1]:
            continue
        if log10(G0) < parms['range_logG0'][0] or log10(G0) > parms['range_logG0'][1]:
            continue
        if alpha < parms['range_alpha'][0] or alpha > parms['range_alpha'][1]:
            continue
        
        ## for each model, looping over the specied values of Av
        for Av in parms['Av_values']:

            outStr = '' #string which will hold the output that will be written to the ascii data file
            
            #           c1      c2      c3    c4      c5       c6      c7
            outStr += '%06d %-+9.3e %-+9.3e %-+9.3e %-+5.4f %-+6.3f  %+05.1f ' % (idx, n, G0, gm, alpha, arxv.metallicity, Av)

            ## getting the emissions for the lines specieified in parms['emissions'] and writing them to each row in the output file
            for line_model, code, specStr, trans in parms['emissions']:
                
                ## getting the info of a line from the corresponding radex model
                if line_model == 'lvg' and code == 'radex':
                    
                    arxv.use_radexDb(Av, specStr)

                    if arxv.meshesRadex[i] != None:
                        
                        transIdx = lineDict.lines[specStr + trans]['radexIdx'] 
                        flux =  arxv.meshesRadex[i][transIdx]['fluxcgs']
                        outStr += '%-+15.3e ' % flux 
                    else:
                        outStr += '%-+15.3e ' % numpy.nan
                        
                ## getting info of a line from the PDR 
                if line_model == 'pdr' and code == 'leiden':
                    
                    ## using the current mesh as the temporary mesh in the mesh arxv object
                    arxv.mshTmp.setData(msh)
                    quantity = lineDict.lines[specStr + trans]['ismcpak']
                    flux =  arxv.mshTmp.compute_integrated_quantity(quantity, Av_range = [0.0, Av])
                    outStr += '%-+15.3e ' % flux 

            ## getting the column density of each specie in parms['col_dense'] and writing them to each row in the output file
            for model_type, specStr in parms['col_dense']:
                
                ## getting the info of a line from the corresponding radex model
                if model_type == 'pdr':
                    
                    ## using the current mesh as the temporary mesh in the mesh arxv object
                    arxv.mshTmp.setData(msh)
                    colDens = arxv.mshTmp.getColumnDensity(specsStrs = [specStr], maxAv = Av)
                    outStr += '%-+15.3e ' % colDens[0] 
            #
            
            ## getting the column density of each specie in parms['col_dense'] and writing them to each row in the output file
            for quantity_name, quantity in parms['int_quantities']:

                ## using the current mesh as the temporary mesh in the mesh arxv object
                arxv.mshTmp.setData(msh)
                q = arxv.mshTmp.compute_integrated_quantity(quantity, Av_range=[0.0, Av])
                outStr += '%-+15.3e ' % q
            #

            ## getting a certain quantity at a certain Av from the meshes and writing them to each row in the output file
            for quantity_name, quantity in parms['quantities_at_AV']:

                ## using the current mesh as the temporary mesh in the mesh arxv object
                q_vs_Av_this_mesh = fetchNestedDtypeValue(msh, quantity)
                Av_this_mesh = fetchNestedDtypeValue(msh, ['state', 'Av'])
                
                ## getting the index in the quantity array nearest to the desired Av values
                ind = numpy.argmin( numpy.fabs(Av_this_mesh - Av) )
                
                q = q_vs_Av_this_mesh[ind]
                
                outStr += '%-+15.3e ' % q
            #
            
            ## writing custom quantities computed at this Av
            for quantity in parms['custom_quantities_at_Av']:
                
                if quantity[0] == 'T_radex':
                    
                    arxv.mshTmp.setData(msh)

                    radex_parms = arxv.mshTmp.get_radex_parameters( speciesStr = 'CO', 
                                                                    threshold  = -1.0,
                                                                    Av_range   = [0.0, Av])
                    
                    q = radex_parms[0]

                    outStr += '%-+15.3e ' % q
                    

            outStr += '\n'
            fObj.write(outStr)
            
            idx += 1
        #
    #
      
    fObj.close()
    return arxv
    print 'done'
#


arxv = main(parameters)
