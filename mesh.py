import numpy as np
import pylab as pyl
import inspect

from ismUtils import getSlabThicknessFromAv

class mesh( ):
    """This class read a PDR mesh and returns an object containing all the data of that mesh.
        
        :param string fileName: path of the binary file holding the data of the mesh. this is assigned to :data:`self.fileName <fName>`
        :param chemicalNetwork chemNet: the chemical network object instance. this is assigned to :data:`self.chemNet <chemNet>`
        :param numpy.float64 metallicity: the metallicity to be used in the computations for this mesh.  this is assigned to :data:`self.metallicity <metallicity>`
        :returns: mesh object

        .. code-block:: python
        
            m = mesh(path)
            m.data['hdr']  # dtype defined in self.headerFormat()
            m.data['hdr']['version']
            m.data['hdr']['G0']
            m.data['hdr']['remaining tags in self.headerFormat()']
            
            m.data['state']  # dtype defined in self.stateFormat()
            m.data['state']['gasT']
            m.data['state']['dustT']
            m.data['state']['remaining tags in self.stateFormat()']
            
            m.data['therm'] # dtype defined in self.thermoFormat()
            m.data['therm']['heating'] 
            m.data['therm']['cooling']
            
            m.data['cooling'] # dtype defined in self.coolingFormat()
            m.data['cooling']['metaStable']
            m.data['cooling']['fineStructure']
            m.data['cooling']['roVib']
            m.data['cooling']['remainig tags in self.coolingFormat()']
            
            m.data['heating'] # dtype defined in self.heatingFormat()
            m.data['heating']['photo']
            m.data['heating']['cIon']
            m.data['heating']['molHydro']
            m.data['heating']['remainig tags in self.heatingFormat()']
            
            # Only for binary files with version 2 have the follwing extra data
            m.data['metaStableCoolingCmponents'] # dtype defined in self.coolingFormaMetaStable()
            m.data['metaStableCoolingCmponents']['C']
            m.data['metaStableCoolingCmponents']['C+']
            m.data['metaStableCoolingCmponents']['Fe']
            m.data['metaStableCoolingCmponents']['remainig species in self.coolingFormatMetaStable()']
            
            m.data['fineStructureCoolingComponents'] # dtype defined in self.coolingFormatFineStructure()
            m.data['fineStructureCoolingComponents']['Si']['popDens']['0']
            m.data['fineStructureCoolingComponents']['Si']['popDens']['1']
            m.data['fineStructureCoolingComponents']['Si']['rate']['1-0']
            m.data['fineStructureCoolingComponents']['remainig species in self.coolingFormatFineStructure']
            
            m.data['selfSheilding']  # dtype defined in self.selfSheildingFormat()
            m.data['selfSheilding']['H2']
            m.data['selfSheilding']['CO']
            m.data['selfSheilding']['13CO']
    """
    def __init__(self, fileName = None, chemNet = None, metallicity = None):

        if fileName != None:
            self.fName = fileName #: path (py:string) of the binary file holding the data of the mesh.
        
            # reading the header, constructing the full mesh dtype
            dtHdr = np.dtype( self.headerFormat() ) 
            hdr   = np.fromfile(self.fName, dtype = dtHdr, count = 1 ) 
            self.hdr = hdr[0]  #: data of the header 

            # re-reading the whole file and constructing the full dtype
            # with the header info read above 
            dtMesh    = self.constructMeshDtype( self.hdr['nSpecs'], self.hdr['nSteps'], self.hdr['version'])
            data      = np.fromfile( fileName, dtype = dtMesh, count = 1)
            #print fileName
            self.data = data[0] 
        else:
            self.data = None  #: (numpy.dtype) holds all the data of the PDR
            self.hdr  = None  #: (numpy.dtype) holds the header info the PDR. This is the same as :data:`self.data['hdr'] <data>`
            
        if chemNet != None:
            self.set_chemNet(chemNet)
        else:
            self.chemNet = None #:

        if metallicity != None:
            self.set_metallicity(metallicity)
        else:
            self.metallicity = None #:
            
        self.figInit = 0     #:
        self.fig     = None  #:
        self.axs     = None  #:
        self.axsRef  = None  #:

        
    def constructMeshDtype(self, nSpecs, nSteps, version ):
        """ constructs the numpy dtype of the mesh which will be used to read the contents of the binary data fiel
        """
        meshFormat = self.meshFormat( nSpecs, nSteps, version )
        meshDtype  = np.dtype( meshFormat )
        return  meshDtype

    def meshFormat(self, nSpecs, nSteps, version):
        """ define the format of the mesh from which the dtype is contrsucted.
        
             :param np.int32 nSpecs: number of species in the chemical network in the file
             :param np.int32 nSteps: number of slabs in the file
             :param np.int32: the version of the data file. For now can handle up to version = 2.
             :returns: list
        """
        dt =  [  ('hdr'    , np.dtype( self.headerFormat ()               ), 1),
                 ('state'  , np.dtype( self.stateFormat  ( nSpecs, nSteps)), 1),
                 ('therm'  , np.dtype( self.thermoFormat( nSteps)         ), 1),
                 ('heating', np.dtype( self.heatingFormat( nSteps)        ), 1),
                 ('cooling', np.dtype( self.coolingFormat( nSteps)        ), 1),
              ]

        # if the data file version is the second version, append the rest of the 
        # file format to be read, detailed cooling and self sheilding stuff
        if version == 2:
            dt.append( ('metaStableCoolingCmponents'    , np.dtype( self.coolingFormaMetaStable(nSteps))    , 1), )
            dt.append( ('fineStructureCoolingComponents', np.dtype( self.coolingFormatFineStructure(nSteps)), 1), )
            dt.append( ('selfSheilding'                 , np.dtype( self.selfSheildingFormat( nSteps))      , 1), )
            
        return dt
    
    def headerFormat(self):
        """ define and return the format of the mesh header."""

        return [ 
                  ('version'  , np.int32   ),
                  ('G0'       , np.float64 ),
                  ('nGas'     , np.float64 ),
                  ('gammaMech', np.float64 ),
                  ('nSteps'   , np.int32   ),
                  ('nSpecs'   , np.int32   ),
               ]
        
    def stateFormat(self, nSpecs, nSteps):
        """ define the format of the state of the mesh."""
        
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
        """ define the format of the total heating and cooling."""

        n = int(nSteps)
        return [
                  ('heating', np.float64, (n) ),
                  ('cooling', np.float64, (n) ),
               ]

    def heatingFormat(self, nSteps):
        """ define the format of the heating components."""
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
        """ define the format of the cooling components."""
        n = int(nSteps)
        return [
                  ('metaStable'    , np.float64, ( n) ),
                  ('fineStructure' , np.float64, ( n) ),
                  ('roVib'         , np.float64, ( n) ),
                  ('recom'         , np.float64, ( n) ),
                  ('lymanAlpha'    , np.float64, ( n) ),
               ]

    def coolingFormaMetaStable(self, nSteps):
        """ define the format of the meta stable cooling components."""
        n = np.int32(nSteps)
        return [
                  ('C'  , np.float64, (n) ),
                  ('C+' , np.float64, (n) ),
                  ('Fe' , np.float64, (n) ),
                  ('Fe+', np.float64, (n) ),
                  ('O'  , np.float64, (n) ),
                  ('O+' , np.float64, (n) ),
                  ('S'  , np.float64, (n) ),
                  ('S+' , np.float64, (n) ),
                  ('Si' , np.float64, (n) ),
                  ('Si+', np.float64, (n) ),
               ]

    # tansitions = ['upper1-lower1', 'upper2-lower2'....etc ]
    #   for example 
    # tansitions = ['1-0', '2-1'....etc ]#
    def coolingFmtFsSpecie(self, levels, transitions, n):
        """ define the format of fine structure cooling transitions of a certain specie.
        
             :param list levels: a list of strings for the names of the levels which can be used to access the population densities. For example :

             .. code-block:: python
             
                 levels = ['0', '1', '2', '3', '3-p']
                 
             would allow for the abundances of each level to be accessed through :
             
             .. code-block:: python
             
                 ['SPECIE']['popDens']['0'] # np.float64 array of the population density of level '0' 
                 ['SPECIE']['popDens']['1'] # np.float64 array of the population density of level '1'

             :param list transitions: a list containing the transitions. For each transitions an array is created as a dtype which can be used to access the cooling rate. For example.

             .. code-block:: python
             
                 transitions = ['1-0','2-1','2-0']

             would allow for the cooling rates for the transitions to be accessed through :
             
             .. code-block:: python
             
                 ['SPECIE']['rate']['1-0'] # np.float64 array of the cooling rate for transition '1-0' 
                 ['SPECIE']['rate']['2-1'] # np.float64 array of the cooling rate for transition '2-1' 
            
             :param np.int32 n: the number of steps which will be the size of all the np.float64 arrays mentioned above             
        """
        fmtCool    = []
        fmtPopDens = []
        
        for trans in transitions:
            fmtCool.append(     (trans, np.float64, (n)),   )
        for level in levels:
            fmtPopDens.append(  (level, np.float64, (n)),   )
            
        return [           
                 ('rate'   , np.dtype(fmtCool), 1),
                 ('popDens', np.dtype(fmtPopDens), 1),
               ]
        
    def coolingFormatFineStructure(self, nSteps):
        """ defines the format of the contents of each specie's fine structure cooling info
            
            :param np.int32 nSteps: the number of steps in the mesh to be set as the size of the arrays holding the fine structure info for the cooling rate and population densitites for each transition and level respectively.
        
        """
        n = np.int32(nSteps)
        return [ ('C+' , np.dtype(self.coolingFmtFsSpecie(['0','1']    , ['1-0']            , n)), 1),
                 ('Si+', np.dtype(self.coolingFmtFsSpecie(['0','1']    , ['1-0']            , n)), 1),
                 ('C'  , np.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('O'  , np.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('S'  , np.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('Fe+', np.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('Si' , np.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                ]

    def selfSheildingFormat(self, nSteps):
        """ defines the format of the arrays which will hold the self sheilding data for the slabs in the mesh.
         
        :param np.int32 nSteps: the number of steps in the mesh to be set as the size of the arrays holding the fine structure info for the cooling rate and population densitites for each transition and level respectively.
        """
        return [ ('H2'  ,  np.float64, (nSteps)),
                 ('CO'  ,  np.float64, (nSteps)),
                 ('13CO',  np.float64, (nSteps)),
               ] 
    def getColumnDensity(self, specsStrs = None, maxAv = None):
        """returns a list of the same lenght as specsStrs containing the column densities
           of the species in the list specsStrs. It is assumed that self.chemNet, 
           self.metallicity and self.data are set.
           
           :param list specsStrs: a list which holds the string names of the species
               whose column densities are to be computed
           :param chemicalNetwork.chemicalNetwork() chemNet: an instance of the object chemicalNetwork
            that holds the info about all the chemistry of the mesh. If this is not passed, self.chemNet 
            is used (assuming it is set).
        """
        
        colDensities = []
        
        # getting mesh parameters
        m = self.data
        nDense_m = m['hdr']['nGas'] 
        abun_m   = m['state']['abun']
        Av_m = m['state']['Av']

        # setting the thickness of the last slab to the one before it
        dxSlabs          =  getSlabThicknessFromAv(Av_m, nDense_m, self.metallicity)
        dxSlabsNew       =  np.ndarray( len(dxSlabs)+1, dtype = np.float64 )
        dxSlabsNew[0:-1] =  dxSlabs
        dxSlabsNew[-1]   =  dxSlabs[-1]
        dxSlabs          =  dxSlabsNew
        
        if maxAv == None:
            slabsIdx = np.arange(m['hdr']['nSteps'])
        else:
            slabsIdx = np.where( Av_m < maxAv)
             
        # computing the column densities for the species in the list
        for specStr in specsStrs:
            specIdx = self.chemNet.species[specStr].num
            nDensSpec = nDense_m * abun_m[ specIdx ][ : ]
            colDens = np.sum( dxSlabs[slabsIdx] * nDensSpec[slabsIdx] )
            
            colDensities.append( colDens )

        return colDensities

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

    def plotMeshGeometry(self):
        
        fig = pyl.figure()
        axs = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        #Av = self.data['state']['Av']
        AvMax  = 30.0
        dAvMin = 0.001
        
        Av = 0.0
        
        print Av
        
        dAv = dAvMin
        
        while Av <= AvMax + 1e-10:
            
            Av += dAv
             
            if Av < 1.0:
                dAv = 10.0**(int(np.log10(Av)))
                
            print Av
                

        """
        while dAv < 1.0:
            
            
        Av = [dAv]
        
        Av = np.arange(31)
        for i in np.arange(30):
            if i > 0 and i <= 10: 
                Av[i] = np.float64(i)/100
            if i > 10 and i <= 20:
                Av[i] = np.float64(i)/10
            if i <= 10: 
                Av[i] = np.float64(i)/100
            
            
        y  = np.ones(Av.shape)
        #asdasd
        axs.semilogx(Av, y, 'r')
        axs.semilogx(Av, y, 'ro')

        axs.set_ylim( ymin=0, ymax=2)
        axs.set_xlim( xmin=0.01, xmax=10)
        #axs.xaxis.set_scale('log')
        pyl.show()
        """
    
    def set_chemNet(self, chemNet):
        self.chemNet = chemNet
    def set_metallicity(self, metallicity):
        self.metallicity = metallicity
    
