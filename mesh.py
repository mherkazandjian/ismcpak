import numpy
import pylab

from ismUtils import getSlabThicknessFromAv, AvToLength
from mylib.utils.misc import fetchNestedDtypeValue

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
            m.data['hdr']['nGas']
            m.data['hdr']['gammaMech']  
            m.data['hdr']['nSteps']
            m.data['hdr']['nSpecs']
            
            m.data['state']  # dtype defined in self.stateFormat()
            m.data['state']['gasT']    # shape = (nSteps)
            m.data['state']['dustT']   # shape = (nSteps)
            m.data['state']['Av']      # shape = (nSteps)
            m.data['state']['abun']    # shape = (nSpecs, nSteps)
            
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
            dtHdr = numpy.dtype( self.headerFormat() ) 
            hdr   = numpy.fromfile(self.fName, dtype = dtHdr, count = 1 ) 
            self.hdr = hdr[0]  #: data of the header 

            # re-reading the whole file and constructing the full dtype
            # with the header info read above 
            dtMesh    = self.constructMeshDtype( self.hdr['nSpecs'], self.hdr['nSteps'], self.hdr['version'])
            data      = numpy.fromfile( fileName, dtype = dtMesh, count = 1)
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
        meshDtype  = numpy.dtype( meshFormat )
        return  meshDtype

    def meshFormat(self, nSpecs, nSteps, version):
        """ define the format of the mesh from which the dtype is contrsucted.
        
             :param numpy.int32 nSpecs: number of species in the chemical network in the file
             :param numpy.int32 nSteps: number of slabs in the file
             :param numpy.int32: the version of the data file. For now can handle up to version = 2.
             :returns: list
        """
        dt =  [  ('hdr'    , numpy.dtype( self.headerFormat ()               ), 1),
                 ('state'  , numpy.dtype( self.stateFormat  ( nSpecs, nSteps)), 1),
                 ('therm'  , numpy.dtype( self.thermoFormat( nSteps)         ), 1),
                 ('heating', numpy.dtype( self.heatingFormat( nSteps)        ), 1),
                 ('cooling', numpy.dtype( self.coolingFormat( nSteps)        ), 1),
              ]

        # if the data file version is the second version, append the rest of the 
        # file format to be read, detailed cooling and self sheilding stuff
        if version == 2:
            dt.append( ('metaStableCoolingCmponents'    , numpy.dtype( self.coolingFormaMetaStable(nSteps))    , 1), )
            dt.append( ('fineStructureCoolingComponents', numpy.dtype( self.coolingFormatFineStructure(nSteps)), 1), )
            dt.append( ('selfSheilding'                 , numpy.dtype( self.selfSheildingFormat( nSteps))      , 1), )
            
        return dt
    
    def headerFormat(self):
        """ define and return the format of the mesh header."""

        return [ 
                  ('version'  , numpy.int32   ),
                  ('G0'       , numpy.float64 ),
                  ('nGas'     , numpy.float64 ),
                  ('gammaMech', numpy.float64 ),
                  ('nSteps'   , numpy.int32   ),
                  ('nSpecs'   , numpy.int32   ),
               ]
        
    def stateFormat(self, nSpecs, nSteps):
        """ define the format of the state of the mesh."""
        
        n = int(nSteps)
        m = int(nSpecs)
        fmt = [
                  ('gasT'  , numpy.float64, ( n) ),
                  ('dustT' , numpy.float64, ( n) ),
                  ('Av'    , numpy.float64, ( n) ),
                  ('abun'  , numpy.float64, (m, n) ),
               ]
        return fmt
    
    def thermoFormat(self, nSteps):
        """ define the format of the total heating and cooling."""

        n = int(nSteps)
        return [
                  ('heating', numpy.float64, (n) ),
                  ('cooling', numpy.float64, (n) ),
               ]

    def heatingFormat(self, nSteps):
        """ define the format of the heating components."""
        n = int(nSteps)
        return [
                  ('photo'    , numpy.float64, (n) ),
                  ('cIon'     , numpy.float64, (n) ),
                  ('molHydro' , numpy.float64, (n) ),
                  ('H2pump'   , numpy.float64, (n) ),
                  ('ggColl'   , numpy.float64, (n) ),
                  ('visc'     , numpy.float64, (n) ),
                  ('cr'       , numpy.float64, (n) ),
               ]

    def coolingFormat(self, nSteps):
        """ define the format of the cooling components."""
        n = int(nSteps)
        return [
                  ('metaStable'    , numpy.float64, ( n) ),
                  ('fineStructure' , numpy.float64, ( n) ),
                  ('roVib'         , numpy.float64, ( n) ),
                  ('recom'         , numpy.float64, ( n) ),
                  ('lymanAlpha'    , numpy.float64, ( n) ),
               ]

    def coolingFormaMetaStable(self, nSteps):
        """ define the format of the meta stable cooling components."""
        n = numpy.int32(nSteps)
        return [
                  ('C'  , numpy.float64, (n) ),
                  ('C+' , numpy.float64, (n) ),
                  ('Fe' , numpy.float64, (n) ),
                  ('Fe+', numpy.float64, (n) ),
                  ('O'  , numpy.float64, (n) ),
                  ('O+' , numpy.float64, (n) ),
                  ('S'  , numpy.float64, (n) ),
                  ('S+' , numpy.float64, (n) ),
                  ('Si' , numpy.float64, (n) ),
                  ('Si+', numpy.float64, (n) ),
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
             
                 ['SPECIE']['popDens']['0'] # numpy.float64 array of the population density of level '0' 
                 ['SPECIE']['popDens']['1'] # numpy.float64 array of the population density of level '1'

             :param list transitions: a list containing the transitions. For each transitions an array is created as a dtype which can be used to access the cooling rate. For example.

             .. code-block:: python
             
                 transitions = ['1-0','2-1','2-0']

             would allow for the cooling rates for the transitions to be accessed through :
             
             .. code-block:: python
             
                 ['SPECIE']['rate']['1-0'] # numpy.float64 array of the cooling rate for transition '1-0' 
                 ['SPECIE']['rate']['2-1'] # numpy.float64 array of the cooling rate for transition '2-1' 
            
             :param numpy.int32 n: the number of steps which will be the size of all the numpy.float64 arrays mentioned above             
        """
        fmtCool    = []
        fmtPopDens = []
        
        for trans in transitions:
            fmtCool.append(     (trans, numpy.float64, (n)),   )
        for level in levels:
            fmtPopDens.append(  (level, numpy.float64, (n)),   )
            
        return [           
                 ('rate'   , numpy.dtype(fmtCool), 1),
                 ('popDens', numpy.dtype(fmtPopDens), 1),
               ]
        
    def coolingFormatFineStructure(self, nSteps):
        """ defines the format of the contents of each specie's fine structure cooling info
            
            :param numpy.int32 nSteps: the number of steps in the mesh to be set as the size of the arrays holding the fine structure info for the cooling rate and population densitites for each transition and level respectively.
        
        """
        n = numpy.int32(nSteps)
        return [ ('C+' , numpy.dtype(self.coolingFmtFsSpecie(['0','1']    , ['1-0']            , n)), 1),
                 ('Si+', numpy.dtype(self.coolingFmtFsSpecie(['0','1']    , ['1-0']            , n)), 1),
                 ('C'  , numpy.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('O'  , numpy.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('S'  , numpy.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('Fe+', numpy.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                 ('Si' , numpy.dtype(self.coolingFmtFsSpecie(['0','1','2'], ['1-0','2-1','2-0'], n)), 1),
                ]

    def selfSheildingFormat(self, nSteps):
        """ defines the format of the arrays which will hold the self sheilding data for the slabs in the mesh.
         
        :param numpy.int32 nSteps: the number of steps in the mesh to be set as the size of the arrays holding the fine structure info for the cooling rate and population densitites for each transition and level respectively.
        """
        return [ ('H2'  ,  numpy.float64, (nSteps)),
                 ('CO'  ,  numpy.float64, (nSteps)),
                 ('13CO',  numpy.float64, (nSteps)),
               ] 
    def getColumnDensity(self, specsStrs = None, maxAv = None):
        """returns a list of the same lenght as specsStrs containing the column densities
           of the species in the list specsStrs. It is assumed that self.chemNet, 
           self.metallicity and self.data are set.
           
           :param list specsStrs: a list which holds the string names of the species
               whose column densities are to be computed. This MUST be a list, even if it
               is one species.
           :param chemicalNetwork.chemicalNetwork() chemNet: an instance of the object chemicalNetwork
            that holds the info about all the chemistry of the mesh. If this is not passed, self.chemNet 
            is used (assuming it is set).
            
           :todo: make use of self.compute_integrated_quantity in this method.
        """
        
        if specsStrs.__class__ is not [].__class__:
            raise TypeError('the parameter specsStrs is not a list.')
        
        colDensities = []
        
        # getting mesh parameters
        data = self.data
        nDense_m = data['hdr']['nGas'] 
        abun_m   = data['state']['abun']
        Av_m = data['state']['Av']

        # setting the thickness of the last slab to the one before it
        dxSlabs = self.compute_dx() 
        
        if maxAv == None:
            slabsIdx = numpy.arange(data['hdr']['nSteps'])
        else:
            slabsIdx = numpy.where( Av_m < maxAv)

        # computing the column densities for the species in the list
        for specStr in specsStrs:
            specIdx = self.chemNet.species[specStr].num
            nDensSpec = nDense_m * abun_m[ specIdx ][ : ]
            NSpec = dxSlabs[slabsIdx] * nDensSpec[slabsIdx]
            colDens = numpy.sum( NSpec )
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
    
    def setFigureObjects(self, figObj=None, axObj=None):
        
        if figObj != None and axObj != None:
            pass
        else:
            figObj, axObj = pylab.subplots(2, 2, sharex=False, sharey=False, figsize=(12,12))

        self.fig     = figObj
        self.axs     = axObj
        self.figInit = 1
        
    def setupFigures(self, plot_Av_range = None):
        """Setup the ranges, labels...etc..for the plots.
        
           :param list plot_Av_range: a list holding the (min, max) of the x axis.  
        """
        
        if plot_Av_range == None:
            plot_x_range = [0, 20]
        else:
            plot_x_range = plot_Av_range
            
        if self.figInit == 0:
            self.setFigureObjects()
        
        def add_NH_labels(axes):
            #setting the labels in N(H) on the twin x-axis on top
            twiny = axes.twiny()
            twiny.set_xlim( axes.get_xlim() )
            #removing the x ticks of both axes (original and the twin)
            for tick in twiny.xaxis.get_major_ticks():
                tick.label1On = False
            for tick in axes.xaxis.get_major_ticks():
                tick.label1On = False
    
            xticksloc = twiny.xaxis.get_majorticklocs()
            xticksNH = AvToLength(xticksloc, 1.0, self.metallicity)
            
            expnnt = 10.0**numpy.floor(numpy.log10(xticksNH[-1]))
            twiny.set_xlabel('N(H) [x %.1e cm^-2]' % expnnt)
            xticksNH_strs = ['%.2f' % (tickv/expnnt) for tickv in xticksNH]
            twiny.set_xticklabels(xticksNH_strs, size = 'small', rotation = 45)

        # subplot 0,0
        #----------------------------------------------------------------        
        self.plt00tgasPlt,  = self.axs[0,0].semilogy([1],  [1], 'r' )
        self.plt00tdustPlt, = self.axs[0,0].semilogy([1],  [1], 'b' )
        self.plt00_v_line,  = self.axs[0,0].semilogy([1],  [1], 'k--' )
        self.axs[0,0].set_xlim(plot_x_range[0], plot_x_range[1])
        self.axs[0,0].set_ylim(1, 100000)
        self.axs[0,0].set_ylabel('$T(K)$')
        self.axs[0,0].text(0.4, 1e3, '$T_{gas}$' , color='r')
        self.axs[0,0].text(0.4, 1e4, '$T_{dust}$', color='b')
        self.plt00Ttl = self.axs[0,0].text(-2.0, 4e6,'$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (0, 0, 0 ) )
        # enabling y ticks on the second axis
        for tick in self.axs[0,0].xaxis.get_major_ticks():
            tick.label1On = False

        #setting the labels in N(H) on the twin x-axis on top
        add_NH_labels(self.axs[0,0])
        
        # subplot 0,1
        #----------------------------------------------------------------
        self.plt01Spec1Plt, = self.axs[0,1].semilogy([1],  [1], 'r' )
        self.plt01Spec2Plt, = self.axs[0,1].semilogy([1],  [1], 'g' )
        self.plt01Spec3Plt, = self.axs[0,1].semilogy([1],  [1], 'b' )
        self.plt01Spec4Plt, = self.axs[0,1].semilogy([1],  [1], 'c' )
        self.plt01Spec5Plt, = self.axs[0,1].semilogy([1],  [1], 'y' )
        self.plt01_v_line,  = self.axs[0,1].semilogy([1],  [1], 'k--' )

        self.axs[0,1].text(0.4, 1e-10, '$H^+$' , color='r')
        self.axs[0,1].text(0.4, 1e-9 , '$H$'   , color='g')
        self.axs[0,1].text(0.4, 1e-8 , '$H_2$' , color='b')
        self.axs[0,1].text(0.4, 1e-11 ,'$e^-$' , color='c')
        self.axs[0,1].text(0.4, 1e-12 ,'$He$'  , color='y')

        #setting ranges of the original axis 
        self.plt01ax1 = self.axs[0,1]
        self.plt01ax1.set_xlim(plot_x_range[0], plot_x_range[1])
        self.plt01ax1.set_ylim(1e-12, 2)
        
        #setting ranges of the twin axis 
        self.plt01ax2 = self.axs[0,1].twinx()
        self.plt01ax2.set_xlim(plot_x_range[0], plot_x_range[1])
        self.plt01ax2.set_ylim(1e-12, 2)
        
        self.plt01ax2.semilogy([1],  [1], 'g' )
        self.plt01ax2.semilogy([1],  [1], 'b' )
        self.plt01ax2.semilogy([1],  [1], 'c' )
        self.plt01ax2.set_ylabel('abun')

        #removing the x ticks of both axes (original and the twin)
        for tick in self.plt01ax1.xaxis.get_major_ticks():
            tick.label1On = False
        for tick in self.plt01ax2.xaxis.get_major_ticks():
            tick.label1On = False
        #removing the y labels of the original axis
        for tick in self.plt01ax1.yaxis.get_major_ticks():
            tick.label1On = False

        #setting the labels in N(H) on the twin x-axis on top
        add_NH_labels(self.axs[0,1])

        # subplot 1,0
        #-------------------------------------------------------------------
        self.plt10Spec1Plt, = self.axs[1,0].semilogy([1],  [1], 'r' )
        self.plt10Spec2Plt, = self.axs[1,0].semilogy([1],  [1], 'g' )
        self.plt10Spec3Plt, = self.axs[1,0].semilogy([1],  [1], 'b' )
        self.plt10Spec4Plt, = self.axs[1,0].semilogy([1],  [1], 'c' )
        self.plt10_v_line,  = self.axs[1,0].semilogy([1],  [1], 'k--' )
        self.axs[1,0].set_xlim(plot_x_range[0], plot_x_range[1])
        self.axs[1,0].set_ylim(1e-12, 2)
        self.axs[1,0].text(0.4, 1e-11, '$O$'   , color='c')
        self.axs[1,0].text(0.4, 1e-10, '$C^+$' , color='r')
        self.axs[1,0].text(0.4, 1e-9 , '$C$'   , color='g')
        self.axs[1,0].text(0.4, 1e-8 , '$CO$'  , color='b')
        self.axs[1,0].set_xlabel('$A_V$')
        self.axs[1,0].set_ylabel('abun')

        #subplot 1,1
        #------------------------------------------------------------------
        self.plt11Spec1Plt, = self.axs[1,1].semilogy([1],  [1], 'r' )
        self.plt11Spec2Plt, = self.axs[1,1].semilogy([1],  [1], 'g' )
        self.plt11Spec3Plt, = self.axs[1,1].semilogy([1],  [1], 'b' )
        self.plt11_v_line,  = self.axs[1,1].semilogy([1],  [1], 'k--' )
        self.axs[1,1].set_xlim(plot_x_range[0], plot_x_range[1])
        self.axs[1,1].set_ylim(1e-12, 2)   
        self.axs[1,1].text(0.4, 1e-10, '$CS$'  , color='r')
        self.axs[1,1].text(0.4, 1e-9 , '$HNC$'  , color='g')
        self.axs[1,1].text(0.4, 1e-8 , '$HCO^+$', color='b')
        self.axs[1,1].set_xlabel('$A_V$')
        for tick in self.axs[1,1].yaxis.get_major_ticks():
            tick.label1On = False
        
    def plot_v_lines_used_in_chemnet(self):
        """plotting the vertical lines on the gui indicating the positions
           in the slab used for the chemistry.
        """

        self.plt00_v_line.set_xdata( [self.chemNet.Av, self.chemNet.Av] )
        self.plt00_v_line.set_ydata( [1, 100000])
        
        self.plt01_v_line.set_xdata( [self.chemNet.Av, self.chemNet.Av] )
        self.plt01_v_line.set_ydata( [1e-12, 1])
        
        self.plt10_v_line.set_xdata( [self.chemNet.Av, self.chemNet.Av] )
        self.plt10_v_line.set_ydata( [1e-12, 1])
        
        self.plt11_v_line.set_xdata( [self.chemNet.Av, self.chemNet.Av] )
        self.plt11_v_line.set_ydata( [1e-12, 1])

        
    def plot(self):
        
        if self.chemNet == None:
            strng = "Error : chemical network object not set." 
            raise NameError(strng)
        else:
            chemNet = self.chemNet
        
        if self.figInit == 0:
            self.setupFigures()
            
        data = self.data
        spcs = chemNet.species
        
        if self.fig == None:
            self.fig, self.axs = pylab.subplots(2, 2, sharex=True, sharey=False) 

        # subplot 0,0
        self.plt00tgasPlt.set_xdata( data['state']['Av'] )
        self.plt00tgasPlt.set_ydata( data['state']['gasT'] )
        self.plt00tdustPlt.set_xdata( data['state']['Av'] )
        self.plt00tdustPlt.set_ydata( data['state']['dustT'] )
        self.plt00Ttl.set_text('$\log_{10} G_0 = $ %4.2f $\log_{10} n_{gas} = $ %4.2f  $\log_{10} \Gamma_{mech} = $  %5.2f\n ' % (numpy.log10(data['hdr']['G0']), numpy.log10(data['hdr']['nGas']), numpy.log10(data['hdr']['gammaMech']) ) )

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
        self.plt11Spec1Plt.set_ydata( data['state']['abun'][spcs['CS'].num] )
        self.plt11Spec2Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec2Plt.set_ydata( data['state']['abun'][spcs['HNC'].num] )
        self.plt11Spec3Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec3Plt.set_ydata( data['state']['abun'][spcs['HCO+'].num] )
        
        """;;;remove those later;;;"""
        """
        m = data
        nDense_m = m['hdr']['nGas']
        Av_m = m['state']['Av']
        y = data['fineStructureCoolingComponents']['C+']['rate']['1-0']
        
        # setting the thickness of the last slab to the one before it
        dxSlabs = self.compute_dx()
        
        q = y
        v = numpy.sum( q*dxSlabs ) / (2.0 * numpy.pi)
        print ';;;', y
        print ';;;', dxSlabs
        print ';;;', q*dxSlabs
        
        self.axs[1,1].set_ylim(1e-10, 1e-0)   
        self.plt11Spec1Plt.set_xdata( data['state']['Av'] )
        self.plt11Spec1Plt.set_ydata( q*dxSlabs )
        """
        """;;;remove those later;;;"""

        
    ## computes the average temperature, weighted by the column density of the XX
    #  specie. which collides with YY_i other species. The average density of the 
    #  YY_i species is also weighed by the column density of the specie XX in each
    #  slab. Only slabs which have abundances greater than XX_threshold are taken
    #  into account in computing N(XX) and the average collider densities and the 
    #  average temepratures.
    def getRadexParameters(self, speciesStr = None, threshold = None, Av_range = None):
        """Returns a list (TMean, nDenseColl, N_specLVG which are 
        the weighted averaged temperature, number density of the collider species and the
        column density of speciesStr. 

        :param string speciesStr: the sting of the species whose column desnity is to be returned.
          (note use the function self.getColumnDensity()).
        :param float threshold: only slabs with an abundance bigger than threshold are considered
          in computing the column density anf the Tmean. Set this to a negative number to make sure
          all the slabs are considered.
        :return: (TMean, nDenseColl, N_specLVG, Av_range)\n  
          TMean (flaot)\n
          nDenseColl (dict) : the number density of the colliders e-, H+, H, He, H2 which are the keys
            of the dict holding the corresponding number densities.\n
          N_spcLVG (float)
        """
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

        # assigning/computing the thickness of the last slab to the one before the last one
        dx    = getSlabThicknessFromAv(Av, nGas, Z)
        dxNew = numpy.ndarray( len(dx)+1, dtype = numpy.float64 )
        dxNew[0:-1] = dx
        dxNew[-1]   = dx[-1]
        dx = dxNew
                
        if Av_range == None:
            Av_range = [0, Av.max()]

        inds_in_Av_range = numpy.where( (Av >= Av_range[0])*(Av <= Av_range[1]) )  
        
        gasT = m.data['state']['gasT'][inds_in_Av_range]
        dx = dx[inds_in_Av_range]
        #abundances of species to be used later in getting the input needed by radex
        xSpec   = m.data['state']['abun'][ net.species[spec].num ][inds_in_Av_range]
        xColle  = m.data['state']['abun'][ net.species['e-'].num ][inds_in_Av_range]
        xCollHP = m.data['state']['abun'][ net.species['H+'].num ][inds_in_Av_range]
        xCollH  = m.data['state']['abun'][ net.species['H'].num ][inds_in_Av_range]
        xCollHe = m.data['state']['abun'][ net.species['He'].num ][inds_in_Av_range]
        xCollH2 = m.data['state']['abun'][ net.species['H2'].num ][inds_in_Av_range]

        #getting the indicies of the slab which have xH2 greater than xMin
        inds = numpy.nonzero( xCollH2  > xMin  )
        if len(inds[0]) == 0:                
            return (None, None, None)
        
        
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

        # calculating the means
        nSpec    = xSpec * nGas
        nColle   = xColle * nGas
        nCollHP  = xCollHP * nGas
        nCollH   = xCollH * nGas
        nCollHe  = xCollHe * nGas
        nCollH2  = xCollH2 * nGas
        
        NSpec     = nSpec * dx
        N_specLVG = numpy.sum(NSpec)

        TMean = numpy.sum( NSpec*gasT    ) / N_specLVG
        nColleMean  = numpy.sum( NSpec*nColle  ) / N_specLVG
        nCollHPMean = numpy.sum( NSpec*nCollHP ) / N_specLVG
        nCollHMean  = numpy.sum( NSpec*nCollH  ) / N_specLVG
        nCollHeMean = numpy.sum( NSpec*nCollHe ) / N_specLVG
        nCollH2Mean = numpy.sum( NSpec*nCollH2 ) / N_specLVG
        
        nDenseColl = {'e-': nColleMean ,
                      'H+': nCollHPMean,
                      'H' : nCollHMean ,
                      'He': nCollHeMean,
                      'H2': nCollH2Mean,
                     }
        
        return (
                TMean, 
                nDenseColl,
                N_specLVG,
                Av_range,
                )

    def plotMeshGeometry(self):
        
        fig = pylab.figure()
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
                dAv = 10.0**(int(numpy.log10(Av)))
                
            print Av
                

        """
        while dAv < 1.0:
            
            
        Av = [dAv]
        
        Av = numpy.arange(31)
        for i in numpy.arange(30):
            if i > 0 and i <= 10: 
                Av[i] = numpy.float64(i)/100
            if i > 10 and i <= 20:
                Av[i] = numpy.float64(i)/10
            if i <= 10: 
                Av[i] = numpy.float64(i)/100
            
            
        y  = numpy.ones(Av.shape)
        #asdasd
        axs.semilogx(Av, y, 'r')
        axs.semilogx(Av, y, 'ro')

        axs.set_ylim( ymin=0, ymax=2)
        axs.set_xlim( xmin=0.01, xmax=10)
        #axs.xaxis.set_scale('log')
        pylab.show()
        """
    
    def set_chemNet(self, chemNet):
        self.chemNet = chemNet
    def set_metallicity(self, metallicity):
        self.metallicity = metallicity
    
    def compute_dx(self):
        """computes and returns the thickness of the slabs for this mesh using self.data.
        
           :warning: do some extensive tests for this.
           :warning: the thickness of the last slab is set to be equal to the one before it. 
               This is a good apporximation. 
        """
        
        nDense_m = self.data['hdr']['nGas'] 
        Av_m     = self.data['state']['Av']

        # setting the thickness of the last slab to the one before it
        dx          =  getSlabThicknessFromAv(Av_m, nDense_m, self.metallicity)
        dxNew       =  numpy.ndarray( len(dx)+1, dtype = numpy.float64 )
        dxNew[0:-1] =  dx
        dxNew[-1]   =  dx[-1]
        dx          =  dxNew

        return dx
        
    def compute_integrated_quantity(self, quantity, Av_range = None):
        """does the integral \sum_0^{N-1} f_i dx_i, where f_i is the quantity
        to be integrated upon and dx is the thickness of each slab.
        
        If Av_range = [Av_1, Av_2] is provided, the integration is done over that range, 
        otherwise, a range [0, Av_max] is used, where Av_max is the maximum
        Av of the PDR model.
        """

        q = fetchNestedDtypeValue(self.data, quantity)
        dxSlabs = self.compute_dx()

        if Av_range == None:
            integrated_quantity = numpy.sum(q*dxSlabs)
        else:
            Av = self.data['state']['Av']
            inds =  numpy.where((Av >= Av_range[0])*(Av < Av_range[1]))
            integrated_quantity = numpy.sum(q[inds]*dxSlabs[inds])
            
        return integrated_quantity