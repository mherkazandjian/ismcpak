import numpy as np
import pylab as pyl

class mesh( ):
    """
    This class read a PDR mesh and returns an object containing all the data of that mesh
    
    the class attributes are the same as the xxxxFormats defined below.
        
        m = mesh(path)
        m.hdr     contains the info in the header         ( dtype )
        m.data    contains temperatures, abundances....   ( dtype )
        m.cooling contains all the cooling info           ( dtype )
        m.heating contains all the heating info           ( dtype )
    """

    def __init__(self, fileName=None):
        
        if fileName != None:
            self.fName = fileName
        
            # reading on the header, constructing the dull mesh dtype
            # and re-reading the whole file into a mesh dtype
            dtHdr = np.dtype( self.headerFormat() )
            hdr   = np.fromfile(self.fName, dtype = dtHdr, count = 1 )
            self.hdr = hdr[0]
            # re-reading the whole file, header and the data
            dtMesh        = self.constructMeshDtype( hdr['nSpecs'], hdr['nSteps'])
            data          = np.fromfile( fileName, dtype = dtMesh, count = 1)
            self.data     = data[0]
        else:
            self.data = None
            self.hdr  = None
            
            #plotting stuff

        self.fig = None
        self.axs = None

    def meshFormat(self, nSpecs, nSteps):
        return [ 
                   ('hdr'    , np.dtype( self.headerFormat ()                ), 1),
                   ('state'  , np.dtype( self.stateFormat  ( nSpecs, nSteps) ), 1),
                   ('heating', np.dtype( self.heatingFormat( nSteps)         ), 1),
                   ('cooling', np.dtype( self.coolingFormat( nSteps)         ), 1),
                ]
    
    def headerFormat(self):
        return [ 
                  ('version'  , np.int32   ),
                  ('G0'       , np.float64 ),
                  ('nGas'     , np.float64 ),
                  ('gammaMech', np.float64 ),
                  ('nSteps'   , np.int32   ),
                  ('nSpecs'   , np.int32   ),
               ]
        
    def stateFormat(self, nSpecs, nSteps):
        return [
                  ('gasT'  , np.float64, (nSteps, 1) ),
                  ('dustT' , np.float64, (nSteps, 1) ),
                  ('Av'    , np.float64, (nSteps, 1) ),
                  ('abun'  , np.float64, (nSpecs, nSteps) ),
               ]

    def heatingFormat(self, nSteps):
        return [
                  ('photo'    , np.float64, (nSteps, 1) ),
                  ('cIon'     , np.float64, (nSteps, 1) ),
                  ('molHydro' , np.float64, (nSteps, 1) ),
                  ('H2pump'   , np.float64, (nSteps, 1) ),
                  ('ggColl'   , np.float64, (nSteps, 1) ),
                  ('visc'     , np.float64, (nSteps, 1) ),
                  ('cr'       , np.float64, (nSteps, 1) ),
               ]

    def coolingFormat(self, nSteps):
        return [
                  ('metaStable'    , np.float64, (nSteps, 1) ),
                  ('fineStructure' , np.float64, (nSteps, 1) ),
                  ('roVib'         , np.float64, (nSteps, 1) ),
                  ('recom'         , np.float64, (nSteps, 1) ),
                  ('lymanAlpha'    , np.float64, (nSteps, 1) ),
               ]

    def constructMeshDtype(self, nSpecs, nSteps):
        return np.dtype( self.meshFormat( nSpecs, nSteps ) ) 
    
    def getData(self):
        return self.data
    
    def saveData(self, fName):
        self.data.tofile( fName )

    def toBuffer(self, fObj):
        self.data.tofile( fObj )

    def setData(self, data):
        self.data = data
    
    def setFigure(self, figObj, axObj):
        self.fig = figObj
        self.axs = axObj
        
    def plot(self, chemNet):
        
        data = self.data
        spcs = chemNet.species
        
        if self.fig == None:
            self.fig, self.axs = plt.subplots(2, 2, sharex=True, sharey=False) 

            
        pyl.subplot(221)
        pyl.hold(False)
        pyl.semilogy(data['state']['Av'],  data['state']['gasT'] , 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['dustT'], 'b' )
        pyl.axis([0, 20, 1, 100000])
        pyl.xlabel('$A_V$')
        pyl.ylabel('$T(K)$')
        pyl.text(0.4, 1e3, '$T_{gas}$' , color='r')
        pyl.text(0.4, 1e4, '$T_{dust}$', color='b')
        pyl.text(0.8, 5e5,'$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (np.log10(data['hdr']['G0']), np.log10(data['hdr']['nGas']), np.log10(data['hdr']['gammaMech']) ) )

        pyl.subplot(222)
        pyl.hold(False)
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H+'].num], 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H'].num],   'g' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H2'].num],  'b' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['e-'].num], 'c' )
        pyl.axis([0, 20, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$H^+$' , color='r')
        pyl.text(0.4, 1e-9, '$H$'   , color='g')
        pyl.text(0.4, 1e-8, '$H_2$'  , color='b')
        pyl.text(0.4, 1e-7, '$e^-$'  , color='c')
        pyl.xlabel('$A_V$')
        pyl.ylabel('abun')

        pyl.subplot(223)
        pyl.hold(False)
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['C+'].num], 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['C'].num],   'g' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['CO'].num],  'b' )
        pyl.axis([0, 20, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$C^+$' , color='r')
        pyl.text(0.4, 1e-9, '$C$'   , color='g')
        pyl.text(0.4, 1e-8, '$CO$'  , color='b')
        pyl.xlabel('$A_V$')
        pyl.ylabel('abun')

        pyl.subplot(224)
        pyl.hold(False)
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['HCN'].num], 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['HNC'].num],   'g' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['HCO+'].num],  'b' )
        pyl.axis([0, 20, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$HCN$' , color='r')
        pyl.text(0.4, 1e-9, '$HNC$'   , color='g')
        pyl.text(0.4, 1e-8, '$HCO^+$'  , color='b')
        pyl.xlabel('$A_V$')
        pyl.ylabel('abun')
