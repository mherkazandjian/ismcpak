# utility function useful for ISM stuff
import numpy as np

def Av2NH(Av, Z):
    """Computes the column density of hydrogen nuclei given an Av using the formula
    in paper1 (add ref)
    """
    return (Av*1.87e21)/Z
def NH2Av(NH, Z):
    """The inverse of :data:`Av2NH`"""
    return NH*(1.0/1.87e21)*Z

    
def AvToLength(Av, nGas, Z):
    """ convert the input visual extinction value to the length in cm (L) of the the NH 
        column (not the actual column density) this is the inverse of Eq-4 
        in paper01 where  :math:`N_H = L.n_{gas}` 
        
        :param Av: The Av (in mag units).
        :type Av: numpy.float64
        :param nGas: The gas density in :math:`cm^{-3}`.
        :type nGas: numpy.float64
        :param Z: The metllicity in terms of solar metallicity (:math:`Z_{\odot}`).
        :type Z: numpy.float64
    """
    return (Av*1.87e21)/(nGas*Z)

def LengthToAv(NH, nGas, Z):
    """ convert the column density to visual extinction value in Av. This is
        Eq-4 in paper01 implemented in fucntion :data:`LengthToAv`
        
    :param NH: The H column density (in :math:`cm^{-2}`).
    :type NH: numpy.float64
    :param nGas: The gas density in :math:`cm^{-3}`.
    :type nGas: numpy.float64
    :param Z: The metllicity in terms of solar metallicity (:math:`Z_{\odot}`).
    :type Z: numpy.float64
    """
    return NH*nGas*Z / 1.87e21

def getSlabThicknessFromAv(slabsAv, nGas, Z):
    """ given the Av of the sub-slabs of a plane parallel model, 
        it computes the thickness of each sub-slab of the descretized
        model. All the parameters are of the same type as :data:`AvToLength`
        except for :
        
        :param slabsAv: An array holding the Av (in mag units) for the staring positions of the sub-slabs.
        :type slabsAv: numpy.ndarry([], dtype = numpy.float64)
        :returns: The thickness of each slab in cm
        :rtype: numpy.ndarry([], dtype = numpy.float64)
     """
    dAvSlabs = slabsAv[1::] - slabsAv[0:-1]
    dxSlabs  = AvToLength( dAvSlabs, nGas, Z)
    return dxSlabs
    
def planckOccupation(h, nu, kb, T):
    """computes the planck function for input parameters.
       :keywords: beta, nu,  
    """
    x = h * nu / (kb * T)
    return 1.0 / (np.exp(x) - 1.0)

def ortho_para_abundance_at_eq(tkin, xH2):
    """compute the ortho to para abundance (relative to xH2) at thermodynamic equilibirum given a kinetic temperature.
    i.e here, they add up to """
    rop = 9.0*np.exp(-170.6/tkin)
    
    xoH2 = xH2 / (1.0 + 1.0 / rop)
    xpH2 = xH2 - xoH2
    
    return xoH2, xpH2
