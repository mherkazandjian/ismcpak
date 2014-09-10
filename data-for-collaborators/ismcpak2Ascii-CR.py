import numpy
from numpy import log10

import time
import sys
import os
import matplotlib
matplotlib.use('Qt4Agg')
#matplotlib.use('PS')

import pylab
import meshUtils
import lineDict

#########################################parameters##########################################################
home = '/home/mher'

parameters = {'dirPath'    : home + '/ism/runs/oneSided/uniformGrid-z-1.0-no-gm-CR-sweep/',
             'range_logn'  : [0.0, 6.0],     #range in log n of models to be considered
             'range_logG0' : [0.0, 6.0],     #range in log G0 of models to be considered
             'range_CR'    : [5e-17, 5e-13], #range in CR of models to be considered (not in log10)
             'Av_values'   : [5.0, 10.0, 30.0], #values in Av to be considered

             #useful for variable CR runs
             'grid_qx'       : ['hdr','nGas'],
             'grid_qy'       : ['hdr','G0'],
             'grid_qz'       : ['from_meshes_info', 'parms', 3, 'CR_rate'],  # 3 indicates the 4th column in self.infoAll['parms']
             
             'emissions'   : (#the emissions to be extracted from the LVG models
                              ['lvg', 'radex' , 'CO'  , '1-0'],   
                              ['lvg', 'radex' , 'CO'  , '2-1'],   
                              ['lvg', 'radex' , 'CO'  , '3-2'],   
                              ['lvg', 'radex' , 'CO'  , '4-3'],   
                              ['lvg', 'radex' , 'CO'  , '5-4'],  
                              ['lvg', 'radex' , 'CO'  , '6-5'], 
                              ['lvg', 'radex' , 'CO'  , '7-6'],   
                              ['lvg', 'radex' , 'CO'  , '8-7'],   
                              ['lvg', 'radex' , 'CO'  , '9-8'],   
                              ['lvg', 'radex' , 'CO'  , '10-9'],  
                              ['lvg', 'radex' , 'CO'  , '11-10'], 
                              ['lvg', 'radex' , 'CO'  , '12-11'], 
                              ['lvg', 'radex' , 'CO'  , '13-12'],
                              ['lvg', 'radex' , 'CO'  , '14-13'], 
                              ['lvg', 'radex' , 'CO'  , '15-14'],
                              #
                              ['lvg', 'radex' , '13CO', '1-0'],   
                              ['lvg', 'radex' , '13CO', '2-1'],   
                              ['lvg', 'radex' , '13CO', '3-2'],   
                              ['lvg', 'radex' , '13CO', '4-3'],   
                              ['lvg', 'radex' , '13CO', '5-4'],  
                              ['lvg', 'radex' , '13CO', '6-5'], 
                              ['lvg', 'radex' , '13CO', '7-6'],   
                              #
                              ['lvg', 'radex' , 'HCO+', '1-0'],   
                              ['lvg', 'radex' , 'HCO+', '4-3'],   
                              ['lvg', 'radex' , 'HCO+', '7-6'],   
                              ['lvg', 'radex' , 'HCN' , '1-0'],
                              ['lvg', 'radex' , 'HCN' , '3-2'],
                              ['lvg', 'radex' , 'HCN' , '4-3'],
                              ['lvg', 'radex' , 'HNC' , '1-0'],
                              ['lvg', 'radex' , 'HNC' , '3-2'],
                              ['pdr', 'leiden', 'C'   , '609'],
                              ['pdr', 'leiden', 'C'   , '369'],
                             ),
             'col_dense'   : (#the column densities to be written
                              ['pdr', 'CO'],   
                              ['pdr', 'H2'],
                             ), 
             'outFname'    : home + '/tmp/test.txt',
             }

# 
def main(parms):
     
    # reading the archive
    print 'setting up the archive'
    t0 = time.time()
    arxv = meshUtils.meshArxv(readDb=True, **parameters)
    print 'time reading data %f' % (time.time() - t0)

    #loading all the radex model info
    for line_model, code, specStr, trans  in parameters['emissions']:
        if line_model == 'lvg':            
            arxv.readDbsRadex(Av = parms['Av_values'], species = specStr)

    fObj = open(parms['outFname'], 'w')

    #writing the header of the file
    fObj.write('#c1    c2         c3         c4         c5      c6      c7    ')    
    for i, p in enumerate(parameters['emissions'] + parameters['col_dense']):
        fObj.write('c%-10d' % (i+7+1))
    fObj.write('\n')
    
    fObj.write('#indx  n_gas      G0         gmech      CR   Z       Av    ')
    for i, p in enumerate(parameters['emissions']):
        fObj.write('%-10s ' % (p[2]+p[3]))
    for i, p in enumerate(parameters['col_dense']):
        fObj.write('N(%-7s) ' % (p[1]))
    fObj.write('\n')

    idx = 0 #an index written at the beginning of each line
    
    #looping over all the meshes and extracting the info needed and writing them to the ascii file
    for i, msh in enumerate(arxv.meshes):
        
        hdr = msh['hdr']
        n, G0, gm = hdr['nGas'], hdr['G0'], hdr['gammaMech']
        CR = 10.0**arxv.grid_z[i]
        
        #selecting models with specified ranges        
        if log10(n) < parms['range_logn'][0] or log10(n) > parms['range_logn'][1]:
            continue
        if log10(G0) < parms['range_logG0'][0] or log10(G0) > parms['range_logG0'][1]:
            continue
        if CR < parms['range_CR'][0] or CR > parms['range_CR'][1]:
            continue
        
        #for each model, looping over the specied values of Av
        for Av in parms['Av_values']:

            outStr = '' #string which will hold the output that will be written to the ascii data file
            
            #           c1      c2      c3    c4      c5       c6      c7
            outStr += '%06d %-+9.3e %-+9.3e %-+9.3e %-+5.4f %-+6.3f  %+05.1f ' % (idx, n, G0, gm, CR, arxv.metallicity, Av)
                        
            #getting the emissions for the lines specieified in parms['emissions'] and writing them to each row in the output file
            for line_model, code, specStr, trans in parms['emissions']:
                
                #getting the info of a line from the corresponding radex model
                if line_model == 'lvg' and code == 'radex':
                    
                    arxv.use_radexDb(Av, specStr)

                    if arxv.meshesRadex[i] != None:
                        
                        transIdx = lineDict.lines[specStr + trans]['radexIdx'] 
                        flux =  arxv.meshesRadex[i][transIdx]['fluxcgs']
                        outStr += '%+10.3e ' % flux #c8->12
                    else:
                        outStr += '%+10.3e ' % numpy.nan
                        
                #getting info of a line from the PDR 
                if line_model == 'pdr' and code == 'leiden':
                    #using the current mesh as the temporary mesh in the mesh arxv object
                    arxv.mshTmp.setData(msh)
                    quantity = lineDict.lines[specStr + trans]['ismcpak']
                    flux =  arxv.mshTmp.compute_integrated_quantity(quantity, Av_range = [0.0, Av])
                    outStr += '%+10.3e ' % flux #c13,c14

            #getting the column density of each specie in parms['col_dense'] and writing them to each row in the output file
            for model, specStr in parms['col_dense']:
                
                #getting the info of a line from the corresponding radex model
                if model == 'pdr':
                    
                    #using the current mesh as the temporary mesh in the mesh arxv object
                    arxv.mshTmp.setData(msh)
                    colDens = arxv.mshTmp.getColumnDensity(specsStrs = [specStr], maxAv = Av)
                    outStr += '%+10.3e ' % colDens[0] #c15,c16
            #
            
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
