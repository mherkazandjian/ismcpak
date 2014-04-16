"""
<keywords>
interpolation functions, scipy, ismcpak, radex, plot, interp, emission, database
</keywords>
<description>
plot stuff from saved interpolation functions
</description>
"""

import os
import numpy
import pylab
import meshUtils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/sph-db-z-1.0-low-res')

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=False)


F = arxvPDR.load_interp_func(info={'source':'radex'},
                             line='HCN1-0', 
                             quantity='fluxKkms')
 

## making a plot of the emission for fixed Av, G0 for varying gmech as a function of n 
if True:
    gms = [-50.0, -30.0, -24.0, -23.0, -22.5, -22.0,]
    Av  = 10.0 
    G0  = 2.0
    
    for gm in gms:
         
        logns  = numpy.linspace(0.0, 6.0, 100.0)
        logG0s = numpy.ones(logns.shape, 'f8')*G0
        logGMs = numpy.ones(logns.shape, 'f8')*gm
        Avs    = numpy.ones(logns.shape, 'f8')*Av
        
        data = numpy.array([logns, logG0s, logGMs, Avs]).T
        v = F.get( data )
        
        pylab.semilogy(logns, 10.0**v, '-', label=gm)
    
    pylab.legend(loc=0)
    pylab.ylabel('flux')
    pylab.show()
    print 'done'


if False:
    ## plotting sections of the DB
    #log10n  = numpy.linspace(-1.0, 6.0, 50)
    #log10G0 = numpy.linspace(-1.0, 6.0, 50)
    #Av      = 10.0
    #gmech   = -24.0
    
    log10n  = numpy.linspace(-1.0, 6.0, 50)
    log10G0 = numpy.linspace(-1.0, 6.0, 50)
    Av      = numpy.linspace(1.0, 30.0, 50)
    gmech   = numpy.linspace(-25.0, -20.0, 50)
    
    data = numpy.meshgrid(
                          log10n,   # log10(n)  
                          log10G0,  # log10(G0)
                          [gmech],  # log10(gmech)
                          Av,       # Av
                         )
    
    data = numpy.array([ dim.flatten() for dim in data], 'f8').T
    
    v = F.get(data).reshape(log10n.size, log10G0.size)
    
    pylab.figure()
    pylab.imshow(v, extent=[log10n.min(), log10n.max(), log10G0.min(), log10G0.max()], origin='lower')

    
pylab.show()
