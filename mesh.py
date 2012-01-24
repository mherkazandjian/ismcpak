import numpy as np
import pylab as pyl

class mesh( ):
    """
    This class read a PDR mesh and returns an object containing all the data of that mesh
    
    the class attributes are the same as the xxxxFormats defined below.
        
        m = mesh(path)
        m.hdr     contains the info in the header
        m.data    contains temperatures, abundances....
        m.cooling contains all the cooling info
        m.heating contains all the heating info
    """

    def __init__(self, fileName=None):
        
        if fileName != None:
            self.fName = fileName
        
            # reading on the header, constructing the dull mesh dtype
            # and re-reading the whole file into a mesh dtype
            dtHdr = np.dtype( self.headerFormat() )
            hdr   = np.fromfile(self.fName, dtype = dtHdr, count = 1 )
            # re-reading the whole file, header and the data
            dtMesh        = self.constructMeshDtype( hdr['nSpecs'], hdr['nSteps'])
            data          = np.fromfile( fileName, dtype = dtMesh, count = 1)
            self.data     = data[0]
            
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
        
    def plot(self, eSpcs):
        
        data = self.data
        
        if self.fig == None:
            self.fig, self.axs = plt.subplots(2, 2, sharex=True, sharey=False) 

        self.fig.text(0.3, 0.04, '$A_V$')
        self.fig.text(0.7, 0.04, '$A_V$')
        self.fig.text(0.04, 0.25, 'Abun', rotation = 90)
        self.fig.text(0.04, 0.75, 'T(K)', rotation = 90)
        self.fig.text(1.0 - 0.08, 0.25, 'Abun', rotation = 90)
        self.fig.text(1.0 - 0.08, 0.75, 'Abun', rotation = 90)


        self.axs[0,0].axis([0, 20, 1, 10000])
        self.axs[0,0].semilogy(data['state']['Av'],  data['state']['gasT'] , 'r' )
        self.axs[0,0].semilogy(data['state']['Av'],  data['state']['dustT'], 'b' )
        self.fig.text(0.4, 0.85, '$T_{gas}$' , color='r')
        self.fig.text(0.4, 0.82, '$T_{dust}$', color='b')

        self.axs[0,1].axis([0, 20, 1e-12, 2])
        self.axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH_P], 'r' )
        self.axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH],   'g' )
        self.axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH2],  'b' )
        self.axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xe_N], 'c' )
        self.fig.text(0.85, 0.84, '$H^+$' , color='r')
        self.fig.text(0.85, 0.81, '$H$'   , color='g')
        self.fig.text(0.85, 0.78, '$H_2$' , color='b')
        self.fig.text(0.85, 0.75, '$e^-$' , color='c')


        self.axs[0,1].axis([0, 20, 1e-12, 2])
        self.axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xC_P], 'r' )
        self.axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xC],   'g' )
        self.axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xCO],  'b' )
        self.fig.text(0.4, 0.44, '$C^+$' , color='r')
        self.fig.text(0.4, 0.41, '$C$'   , color='g')
        self.fig.text(0.4, 0.38, '$CO$'  , color='b')

        self.axs[0,1].axis([0, 20, 1e-12, 2])
        self.axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHCN], 'r' )
        self.axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHNC],   'g' )
        self.axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHCO_P],  'b' )
        self.fig.text(0.8, 0.44, '$HCN$'   , color='r')
        self.fig.text(0.8, 0.41, '$HNC$'   , color='g')
        self.fig.text(0.8, 0.38, '$HCO^+$' , color='b')
        
    # clears the plotted lines in all the subplots
    def clearLinesInSubPlots(self):
        
        if self.fig != None:
            for axs in self.axs.flat:
                n = len(axs.lines)
                i = 0
                while i < n:
                    axs.lines.pop(0)
                    i += 1
        