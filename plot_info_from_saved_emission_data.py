"""
<keywords>
scipy, ismcpak, radex, plot, interp, emission, database
</keywords>
<description>
plot stuff from saved emission data database
</description>
"""

import os
import numpy
import pylab
import meshUtils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

v, data = arxvPDR.get_emission_from_all_radex_dbs_for_Av_range(line = 'HCN1-0', 
                                                               Avs  = 'all', 
                                                               quantity = 'fluxKkms')

n, g0, alpha, Av = data.T

## making a plot of the emission for fixed Av, G0 for varying gmech as a function of n 
if True:

    pylab.figure()
        
    alphas = numpy.log10([1e-10, 0.01, 0.05, 0.25, 0.75])
    AvSec  = 10.0
    G0Sec  = 5.0
    
    for alphaThis in alphas:
         
        inds = numpy.where( 
                           (g0 >= (G0Sec - 1e-6))*(g0 <= (G0Sec + 1e-6))*
                           (Av >= (AvSec - 1e-6))*(Av <= (AvSec + 1e-6))*
                           (alpha >= (alphaThis - 1e-6))*(alpha <= (alphaThis + 1e-6))
                          )
        
        lognPlt = n[inds]
        vPlt = v[inds]
        
        inds = numpy.argsort(lognPlt)
        lognPlt = lognPlt[inds]
        vPlt = vPlt[inds]
        
        plt, = pylab.semilogy(lognPlt, 10.0**vPlt, '-', label=10.0**alphaThis)
        pylab.semilogy(lognPlt, 10.0**vPlt, plt.get_color()+'o-')
    
    pylab.legend(loc=0)
    pylab.ylabel('flux')
    pylab.show()
    print 'done'
    
pylab.show()
