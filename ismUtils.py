# utility function useful for ISM stuff

def AvToLength(Av, nGas, Z):
    """ convert the input visual extinction value to the length (L) of the the NH 
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
    
    