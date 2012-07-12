import numpy as np
import pylab as pyl
import inspect

from ismUtils import getSlabThicknessFromAv

class mesh( ):
    """
    This class read a PDR mesh and returns an object containing all the data of that mesh
    
    the class attributes are the same as the xxxxFormats defined below.
        
        m = mesh(path)
        m.hdr     contains the info in the header         ( dtype )
        
        m.data    contains everything abt the mesh        ( dtype ) 
          data['hdr']     = m.hdr
          data['state']   = temperatures, abundances.... 
          data['thermo']  = total heating and cooling  
          data['cooling'] = detailed info about cooling
          data['heating'] = detailed info about heating
          
    templates to extract 1D arrays along the slab
         Av         = self.data['state']['Av'][0][0,:]
         abunSpecie = self.data['state']['abun'][0]
    """

    def __init__(self, fileName=None, chemNet = None, metallicity = None):
        
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
            
        if chemNet != None:
            self.set_chemNet(chemNet)
        else:
            self.chemNet = None

        if metallicity != None:
            self.set_metallicity(metallicity)
        else:
            self.metallicity = None
            
        self.figInit = 0
        self.fig     = None
        self.axs     = None
        self.axsRef  = None

        
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
                 ('therm'  , np.dtype( self.thermoFormat( nSteps)         ), 1),
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
                  ('gasT'  , np.float64, ( n) ),
                  ('dustT' , np.float64, ( n) ),
                  ('Av'    , np.float64, ( n) ),
                  ('abun'  , np.float64, (m, n) ),
               ]
        return fmt

    def thermoFormat(self, nSteps):
        n = int(nSteps)
        return [
                  ('heating', np.float64, (n) ),
                  ('cooling', np.float64, (n) ),
               ]

    def heatingFormat(self, nSteps):
        n = int(nSteps)
        return [
                  ('photo'    , np.float64, (n) ),
                  ('cIon'     , np.float64, (n) ),
                  ('molHydro' , np.float64, (n) ),
                  ('H2pump'   , np.float64, (n) ),
                  ('ggColl'   , np.float64, (n) ),
                  ('visc'     , np.float64, (n) ),
                  ('cr'       , np.float64, (n) ),
               ]

    def coolingFormat(self, nSteps):
        n = int(nSteps)
        return [
                  ('metaStable'    , np.float64, ( n) ),
                  ('fineStructure' , np.float64, ( n) ),
                  ('roVib'         , np.float64, ( n) ),
                  ('recom'         , np.float64, ( n) ),
                  ('lymanAlpha'    , np.float64, ( n) ),
               ]
    
    def getData(self):
        return self.data
    
    def saveData(self, fName):
        self.data.tofile( fName )

    def toBuffer(self, fObj):
        self.data.tofile( fObj )

    def setData(self, data):
        self.data = data
    
    def setFigureObjects(self, figObj=None, axObj=None, axRef=None):
        
        if figObj != None and axObj != None and axRef != None:
            pass
        else:
            figObj, axObj = pyl.subplots(2, 2, sharex=False, sharey=False, figsize=(12,12))
            axRef = np.array( [[221,222], [223,224]])   

        self.fig     = figObj
        self.axs     = axObj
        self.axsRef  = axRef        
        self.figInit = 1
        
    def setupFigures(self):
        
        if self.figInit == 0:
            self.setFigureObjects()
            
        # subplot 0,0        
        pyl.subplot(self.axsRef[0,0])
        pyl.hold(False)
        self.plt00tgasPlt,  = pyl.semilogy([1],  [1], 'r' )
        pyl.hold()
        self.plt00tdustPlt, = pyl.semilogy([1],  [1], 'b' )
        pyl.axis([0, 20, 1, 100000])
        pyl.ylabel('$T(K)$')
        pyl.text(0.4, 1e3, '$T_{gas}$' , color='r')
        pyl.text(0.4, 1e4, '$T_{dust}$', color='b')
        self.plt00Ttl = pyl.text(0.8, 5e5,'$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (0, 0, 0 ) )
        # enabling y ticks on the second axis
        for tick in pyl.gca().xaxis.get_major_ticks():
            tick.label1On = False
        
        
        # subplot 0,1
        pyl.subplot(self.axsRef[0,1])        
        pyl.hold(False)
        self.plt01Spec1Plt,  = pyl.semilogy([1],  [1], 'r' )
        pyl.hold()
        self.plt01Spec2Plt, = pyl.semilogy([1],  [1], 'g' )
        self.plt01Spec3Plt, = pyl.semilogy([1],  [1], 'b' )
        self.plt01Spec4Plt, = pyl.semilogy([1],  [1], 'c' )
        self.plt01Spec5Plt, = pyl.semilogy([1],  [1], 'y' )
        pyl.axis([0, 20.0, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$H^+$' , color='r')
        pyl.text(0.4, 1e-9 , '$H$'   , color='g')
        pyl.text(0.4, 1e-8 , '$H_2$' , color='b')
        pyl.text(0.4, 1e-11 ,'$e^-$' , color='c')
        pyl.text(0.4, 1e-12 ,'$He$'  , color='y')
        self.plt01ax1 = pyl.gca()
        for tick in self.plt01ax1.xaxis.get_major_ticks():
            tick.label1On = False
        self.plt01ax2 = pyl.gca().twinx()
        pyl.setp(self.plt01ax2, 'ylim',(1e-12, 2) )
        pyl.setp(self.plt01ax2, 'xlim',(0, 20) )
        # redundant, but we do it just to get the right ticks on the y axis
        self.plt01ax2.semilogy([1],  [1], 'g' )
        self.plt01ax2.semilogy([1],  [1], 'b' )
        self.plt01ax2.semilogy([1],  [1], 'c' )
        # deleting all the ticks on the first axis
        for tick in self.plt01ax1.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False
        # enabling y ticks on the second axis
        for tick in self.plt01ax2.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = True
        self.plt01ax2.set_ylabel('abun')

        # subplot 1,0
        pyl.subplot(self.axsRef[1,0])
        pyl.hold(False)
        self.plt10Spec1Plt,  = pyl.semilogy([1],  [1], 'r' )
        pyl.hold()
        self.plt10Spec2Plt, = pyl.semilogy([1],  [1], 'g' )
        self.plt10Spec3Plt, = pyl.semilogy([1],  [1], 'b' )
        self.plt10Spec4Plt, = pyl.semilogy([1],  [1], 'c' )
        pyl.axis([0, 20, 1e-12, 2])
        pyl.text(0.4, 1e-11 , '$O$'   , color='c')
        pyl.text(0.4, 1e-10, '$C^+$' , color='r')
        pyl.text(0.4, 1e-9 , '$C$'   , color='g')
        pyl.text(0.4, 1e-8 , '$CO$'  , color='b')
        pyl.xlabel('$A_V$')
        pyl.ylabel('abun')
        
        #subplot 1,1
        pyl.subplot(self.axsRef[1,1])
        pyl.hold(False)
        self.plt11Spec1Plt,  = pyl.semilogy([1],  [1], 'r' )
        pyl.hold()
        self.plt11Spec2Plt, = pyl.semilogy([1],  [1], 'g' )
        self.plt11Spec3Plt, = pyl.semilogy([1],  [1], 'b' )
        pyl.axis([0, 20, 1e-12, 2])
        pyl.text(0.4, 1e-10, '$HCN$'  , color='r')
        pyl.text(0.4, 1e-9 , '$HNC$'  , color='g')
        pyl.text(0.4, 1e-8 , '$HCO^+$', color='b')
        pyl.xlabel('$A_V$')
        for tick in pyl.gca().yaxis.get_major_ticks():
            tick.label1On = False
        #pyl.setp(pyl.gca(), yticks=[])
        
    def plot(self):
        
        if self.chemNet == None:
            strng = "Error : chemical networ object not set." 
            raise NameError(strng)
        else:
            chemNet = self.chemNet
        
        if self.figInit == 0:
            self.setupFigures()
            
        data = self.data
        spcs = chemNet.species
        
        if self.fig == None:
            self.fig, self.axs = pyl.subplots(2, 2, sharex=True, sharey=False) 

        # subplot 0,0
        self.plt00tgasPlt.set_xdata( data['state']['Av'] )
        self.plt00tgasPlt.set_ydata( data['state']['gasT'] )
        self.plt00tdustPlt.set_xdata( data['state']['Av'] )
        self.plt00tdustPlt.set_ydata( data['state']['dustT'] )
        self.plt00Ttl.set_text('$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (np.log10(data['hdr']['G0']), np.log10(data['hdr']['nGas']), np.log10(data['hdr']['gammaMech']) ) )
        
        # subplot 0,1
        self.plt01Spec1Plt.set_xdata( data['state']['Av'] )
        self.plt01Spec1Plt.set_ydata( data['state']['abun'][spcs['H+'].num] )
        self.plt01Spec2Plt.set_xdata( data['state']['Av'] )
        self.plt01Spec2Plt.set_ydata( data['state']['abun'][spcs['H'].num] )
        self.plt01Spec3Plt.set_xdata( data['state']['Av'] )
        self.plt01Spec3Plt.set_ydata( data['state']['abun'][spcs['H2'].num] )
        self.plt01Spec4Plt.set_xdata( data['state']['Av'] )
        self.plt01Spec4Plt.set_ydata( data['state']['abun'][spcs['e-'].num] )
        self.plt01Spec5Plt.set_xdata( data['state']['Av'] )
        self.plt01Spec5Plt.set_ydata( data['state']['abun'][spcs['He'].num] )
        
        # subplot 1,0
        self.plt10Spec1Plt.set_xdata( data['state']['Av'] )
        self.plt10Spec1Plt.set_ydata( data['state']['abun'][spcs['C+'].num] )
        self.plt10Spec2Plt.set_xdata( data['state']['Av'] )
        self.plt10Spec2Plt.set_ydata( data['state']['abun'][spcs['C'].num] )
        self.plt10Spec3Plt.set_xdata( data['state']['Av'] )
        self.plt10Spec3Plt.set_ydata( data['state']['abun'][spcs['CO'].num] )
        self.plt10Spec4Plt.set_xdata( data['state']['Av'] )
        self.plt10Spec4Plt.set_ydata( data['state']['abun'][spcs['O'].num] )

        # subplot 1,1
        self.plt11Spec1Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec1Plt.set_ydata( data['state']['abun'][spcs['HCN'].num] )
        self.plt11Spec2Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec2Plt.set_ydata( data['state']['abun'][spcs['HNC'].num] )
        self.plt11Spec3Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec3Plt.set_ydata( data['state']['abun'][spcs['HCO+'].num] )
        
    ## computes the average temperature, weighted by the column density of the XX
    #  specie. which collides with YY_i other species. The average density of the 
    #  YY_i species is also weighed by the column density of the specie XX in each
    #  slab. Only slabs which have abundances greater than XX_threshold are taken
    #  into account in computing N(XX) and the average collider densities and the 
    #  average temepratures.
    def getRadexParameters(self, colliderStr=None, speciesStr=None, threshold=None):

        coll = colliderStr
        spec = speciesStr
        xMin = threshold
        
        # checking the metallciity
        if self.metallicity == None:
            strng = "Error : metallicity of mesh not set"
            raise NameError(strng)
        else:
            Z = self.metallicity
            
        m   = self
        net = self.chemNet 
        
        nGas = m.data['hdr']['nGas']
        Av   = m.data['state']['Av']
        gasT = m.data['state']['gasT']

        xSpec   = m.data['state']['abun'][ net.species[spec].num ]
        xColle  = m.data['state']['abun'][ net.species['e-'].num ]
        xCollHP = m.data['state']['abun'][ net.species['H+'].num ]
        xCollH  = m.data['state']['abun'][ net.species['H'].num ]
        xCollHe = m.data['state']['abun'][ net.species['He'].num ]
        xCollH2 = m.data['state']['abun'][ net.species['H2'].num ]

        inds = np.nonzero( xCollH2  > xMin  ) # ;;; remove this later
        if len(inds[0]) == 0:                 # ;;; remove this later
            return (None, None, None)
        
        # assigning the thickness of the last slab to the one before the last one
        dx    = getSlabThicknessFromAv(Av, nGas, Z)
        dxNew = np.ndarray( len(dx)+1, dtype = np.float64 )
        dxNew[0:-1] = dx
        dxNew[-1]   = dx[-1]
        dx = dxNew

        # selecting the slabs which are usefull
        Av    = Av[inds]
        dx    = dx[inds]
        xSpec = xSpec[inds]
        gasT  = gasT[inds]
        xColle  = xColle[inds]
        xCollHP = xCollHP[inds]
        xCollH  = xCollH[inds]
        xCollHe = xCollHe[inds]
        xCollH2 = xCollH2[inds]

        #print xCollH2
        # calculating the means
        nSpec    = xSpec * nGas
        nColle   = xColle * nGas
        nCollHP  = xCollHP * nGas
        nCollH   = xCollH * nGas
        nCollHe  = xCollHe * nGas
        nCollH2  = xCollH2 * nGas
        
        NSpec     = nSpec * dx 
        N_specLVG = np.sum(NSpec)

        TMean = np.sum( NSpec*gasT    ) / N_specLVG
        nColleMean  = np.sum( NSpec*nColle  ) / N_specLVG
        nCollHPMean = np.sum( NSpec*nCollHP ) / N_specLVG
        nCollHMean  = np.sum( NSpec*nCollH  ) / N_specLVG
        nCollHeMean = np.sum( NSpec*nCollHe ) / N_specLVG
        nCollH2Mean = np.sum( NSpec*nCollH2 ) / N_specLVG
        
        nDenseColl = {'e-': nColleMean ,
                      'H+': nCollHPMean,
                      'H' : nCollHMean ,
                      'He': nCollHeMean,
                      'H2': nCollH2Mean}
            
        #print 'mesh.py:gas T at the end of the slab = ', gasT[-1] 
        return (
                TMean, 
                nDenseColl,
                N_specLVG
                )
    
    ## 
    def getRadexParameters(self, colliderStr=None, speciesStr=None, threshold=None):


    def set_chemNet(self, chemNet):
        self.chemNet = chemNet
    def set_metallicity(self, metallicity):
        self.metallicity = metallicity
    