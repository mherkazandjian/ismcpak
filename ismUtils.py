# utility function useful for ISM stuff

# equation 4 in Kazandian et. al (from Bohlin et. al 1978)
def AvToLength(Av, nGas, Z):
    return (Av*1.87e21)/(nGas*Z)

# slabAv : (ndarray of length at least 2) 
#          visual extintion at the beginning of each slab
# nGas   : number density of hydrogen nuclei
# Z      : metallicity
def getSlabThicknessFromAv(slabsAv, nGas, Z):
    
    dAvSlabs = slabsAv[1::] - slabsAv[0:-1]
    dxSlabs  = AvToLength( dAvSlabs, nGas, Z)
    return dxSlabs
    
    