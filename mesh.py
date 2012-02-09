import numpy as np
import pylab as pyl

class mesh( ):
    """
    This class read a PDR mesh and returns an object containing all the data of that mesh
    
    the class attributes are the same as the xxxxFormats defined below.
        
        m = mesh(path)
        m.hdr     contains the info in the header         ( dtype )
        
        m.data    contains everything abt the mesh        ( dtype ) 
          data['hdr']     = m.hdr
          data['state']   = temperatures, abundances.... 
          data['cooling'] = cooling information
          data['heating'] = heating information
    """

    def __init__(self, fileName=None):
        
        if fileName != None:
            self.fName = fileName
        
            # reading the header, constructing the full mesh dtype
            # and re-reading the whole file into a mesh dtype
            dtHdr = np.dtype( self.headerFormat() )
            hdr   = np.fromfile(self.fName, dtype = dtHdr, count = 1 )
            self.hdr = hdr[0]

            # re-reading the whole file, header and the data
            dtMesh    = self.constructMeshDtype( self.hdr['nSpecs'], self.hdr['nSteps'])
            data      = np.fromfile( fileName, dtype = dtMesh, count = 1)
            self.data = data[0]
        else:
            self.data = None
            self.hdr  = None
            
        self.fig = None
        self.axs = None

    def constructMeshDtype(self, nSpecs, nSteps):
        meshFormat = self.meshFormat( nSpecs, nSteps )
        meshDtype  = np.dtype( meshFormat )
        #print meshDtype
        #adasdads
        return  meshDtype

    def meshFormat(self, nSpecs, nSteps):
         
        return [ 
                 ('hdr'    , np.dtype( self.headerFormat ()               ), 1),
                 ('state'  , np.dtype( self.stateFormat  ( nSpecs, nSteps)), 1),
                 ('heating', np.dtype( self.heatingFormat( nSteps)        ), 1),
                 ('cooling', np.dtype( self.coolingFormat( nSteps)        ), 1),
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
        n = int(nSteps)
        m = int(nSpecs)
        fmt = [
                  ('gasT'  , np.float64, (n, 1) ),
                  ('dustT' , np.float64, (n, 1) ),
                  ('Av'    , np.float64, (n, 1) ),
                  ('abun'  , np.float64, (m, n) ),
               ]
        return fmt

    def heatingFormat(self, nSteps):
        n = int(nSteps)
        return [
                  ('photo'    , np.float64, (n, 1) ),
                  ('cIon'     , np.float64, (n, 1) ),
                  ('molHydro' , np.float64, (n, 1) ),
                  ('H2pump'   , np.float64, (n, 1) ),
                  ('ggColl'   , np.float64, (n, 1) ),
                  ('visc'     , np.float64, (n, 1) ),
                  ('cr'       , np.float64, (n, 1) ),
               ]

    def coolingFormat(self, nSteps):
        n = int(nSteps)
        return [
                  ('metaStable'    , np.float64, (n, 1) ),
                  ('fineStructure' , np.float64, (n, 1) ),
                  ('roVib'         , np.float64, (n, 1) ),
                  ('recom'         , np.float64, (n, 1) ),
                  ('lymanAlpha'    , np.float64, (n, 1) ),
               ]
    
    def getData(self):
        return self.data
    
    def saveData(self, fName):
        self.data.tofile( fName )

    def toBuffer(self, fObj):
        self.data.tofile( fObj )

    def setData(self, data):
        self.data = data
    
    def setFigure(self, figObj, axObj, axRef):
        self.fig = figObj
        self.axs = axObj
        self.axsRef = axRef
        
    def plot(self, chemNet):
        
        data = self.data
        spcs = chemNet.species
        
        if self.fig == None:
            self.fig, self.axs = pyl.subplots(2, 2, sharex=True, sharey=False) 

        # subplot 0,0
        pyl.subplot(self.axsRef[0,0])
        pyl.hold(False)
        pyl.semilogy(data['state']['Av'],  data['state']['gasT'] , 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['dustT'], 'b' )
        pyl.axis([0, 20, 1, 100000])
        pyl.ylabel('$T(K)$')
        pyl.text(0.4, 1e3, '$T_{gas}$' , color='r')
        pyl.text(0.4, 1e4, '$T_{dust}$', color='b')
        pyl.text(0.8, 5e5,'$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (np.log10(data['hdr']['G0']), np.log10(data['hdr']['nGas']), np.log10(data['hdr']['gammaMech']) ) )
        pyl.setp(pyl.gca(), xticks=[])

        # subplot 0,1        
        pyl.subplot(self.axsRef[0,1])        
        pyl.hold(False)
        
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H+'].num], 'r' )
        pyl.hold()
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H'].num],   'g' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H2'].num],  'b' )
        pyl.semilogy(data['state']['Av'],  data['state']['abun'][spcs['e-'].num], 'c' )
        pyl.axis([0, 20.0, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$H^+$' , color='r')
        pyl.text(0.4, 1e-9 , '$H$'   , color='g')
        pyl.text(0.4, 1e-8 , '$H_2$' , color='b')
        pyl.text(0.4, 1e-7 , '$e^-$' , color='c')
        ax1 = pyl.gca()        
        ax2 = pyl.gca().twinx()
        pyl.setp(ax1, xticks=[])
        pyl.setp(ax2, xticks=[])
        pyl.setp(ax2, 'ylim',(1e-12, 2) )
        # redundant, but we do it just to get the right ticks on the y axis
        ax2.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H'].num],   'g' )
        ax2.semilogy(data['state']['Av'],  data['state']['abun'][spcs['H2'].num],  'b' )
        ax2.semilogy(data['state']['Av'],  data['state']['abun'][spcs['e-'].num], 'c' )
        # deleting all the ticks on the first axis
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False
        # enabling y ticks on the second axis
        for tick in ax2.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
        ax2.set_ylabel('abun')
        
        # subplot 1,0        
        pyl.subplot(self.axsRef[1,0])
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

        # subplot 1,1        
        pyl.subplot(self.axsRef[1,1])
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
        pyl.setp(pyl.gca(), yticks=[])
