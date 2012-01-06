import numpy as np
import struct

import pylab as pyl
import matplotlib.pyplot as plt

import pickle  #### remove later

class mesh(  ):
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
        
        self.fileName = fileName
        self.fObj     = open(self.fileName, mode='rb')
        self.hdr      = self.readBlockFromFileObj( self.headerFormat() )
        self.state    = self.readBlockFromFileObj( self.stateFormat()  )
        self.heating  = self.readBlockFromFileObj( self.heatingFormat() )
        self.cooling  = self.readBlockFromFileObj( self.coolingFormat() )
        # hold all the info of the mesh in a dtype variable with the 
        # formats specified below
        self.data     = self.getData()
        

    def meshFormat(self):
        return [ 
                   ('hdr'    , np.dtype( self.headerFormat()   ), 1),
                   ('state'  , np.dtype( self.stateFormat()    ), 1),
                   ('heating', np.dtype( self.heatingFormat()  ), 1),
                   ('cooling', np.dtype( self.coolingFormat()  ), 1),
                ]
    
    def headerFormat(self):
        return [ 
                  ('version'  , np.int32  , 1),
                  ('G0'       , np.float64, 1),
                  ('nGas'     , np.float64, 1),
                  ('gammaMech', np.float64, 1),
                  ('nSteps'   , np.int32  , 1),
                  ('nSpecs'   , np.int32  , 1),
               ]
        
    def stateFormat(self):
        return [
                  ('gasT'  , np.float64, (self.hdr['nSteps'], 1) ),
                  ('dustT' , np.float64, (self.hdr['nSteps'], 1) ),
                  ('Av'    , np.float64, (self.hdr['nSteps'], 1) ),
                  ('abun'  , np.float64, (self.hdr['nSpecs'], self.hdr['nSteps']) ),
               ]

    def heatingFormat(self):
        return [
                  ('photo'    , np.float64, (self.hdr['nSteps'], 1) ),
                  ('cIon'     , np.float64, (self.hdr['nSteps'], 1) ),
                  ('molHydro' , np.float64, (self.hdr['nSteps'], 1) ),
                  ('H2pump'   , np.float64, (self.hdr['nSteps'], 1) ),
                  ('ggColl'   , np.float64, (self.hdr['nSteps'], 1) ),
                  ('visc'     , np.float64, (self.hdr['nSteps'], 1) ),
                  ('cr'       , np.float64, (self.hdr['nSteps'], 1) ),
               ]

    def coolingFormat(self):
        return [
                  ('metaStable'    , np.float64, (self.hdr['nSteps'], 1) ),
                  ('fineStructure' , np.float64, (self.hdr['nSteps'], 1) ),
                  ('roVib'         , np.float64, (self.hdr['nSteps'], 1) ),
                  ('recom'         , np.float64, (self.hdr['nSteps'], 1) ),
                  ('lymanAlpha'    , np.float64, (self.hdr['nSteps'], 1) ),
               ]

    def getData(self):       
        meshFmt = self.meshFormat() 
        data = np.empty((),  dtype = meshFmt)
        
        data['hdr']     = self.hdr
        data['state']   = self.state      
        data['heating'] = self.heating
        data['cooling'] = self.cooling

        return data 
        

    def saveData(self, fName):
        f = open(fName, 'wb')
        pickle.dump(self.data, f)
        f.close() 

    def readBlockFromFileObj(self, dtype):
        dtype    = np.dtype( dtype )
        dataRead = self.readDtypeDataFromFileObject( dtype, self.fObj)
        dtypeObj = self.setValuesToDtype(dataRead, dtype)
        return dtypeObj
        
    def readDataBlock(self):
        
        dataDtype = np.dtype( self.dataFormat() )
        dtype    = dataDtype
        dataRead = self.readDtypeDataFromFileObject( dtype, self.fObj)
        data     = self.setValuesToDtype(dataRead, dtype)
        return data

    def numpyToStructFormatChar(self, str ):
        str = str.split('<')
        str = str[1]
        strRet = ''
    
        if str == 'i4':
            strRet = 'I'
        if str == 'f8':
            strRet = 'd'
        
        if strRet == '':
            return -1
        else:
            return strRet
         
    def numItems(self, feild ):
        if len(feild) == 2:        
            times = 1
        else:
            shape = feild[2]
            times = 1
            for d in shape:
                times *= d
        
        return times
    
    def readDtypeDataFromFileObject(self, dtype, fObj):
        nBytes = 0
        collateFmt = []
        collateFmt.append('<')

        # generating the format string used to read and format the bytes
        for fld in dtype.descr:
            # name of the field
            name = fld[0]
            # type of the field use to unpack into a struct
            typeStr = self.numpyToStructFormatChar(fld[1])
            # size of each item in bytes
            size =  np.dtype( fld[1] ).itemsize
            # number of items in this feild
            times = self.numItems( fld )

            for t in range(times):
                collateFmt.append(typeStr)
                nBytes += size
                    
        #        print name, typeStr, size, times 
        #        print len(fld)
        #        print '---------------------'
    
        fmt   = ''.join(collateFmt)
        #print fmt
        dataRead = fObj.read(nBytes);
        dataRead = (struct.unpack(fmt, dataRead))
        #    print fmt
        #    print dataRead

        return dataRead


    def setValuesToDtype(self, dataRead, dtype):

        var = np.empty((),  dtype = dtype)

        #setting the vlaues to the attributes
        offset=0
        for fld in dtype.descr:
            name  = fld[0]
            times = self.numItems( fld )
            nextOffset = offset + times
            fldData = dataRead[offset:nextOffset]
            # casting the read data into the correct type (if it is not just in case)     
            fldData = np.array( fldData, dtype = var[name].dtype)
            offset = nextOffset
            var[name] = fldData.reshape( var[name].shape )

        return var    

    def plot(self, eSpcs):
        
        data = self.data
        
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=False)

        fig.text(0.3, 0.04, '$A_V$')
        fig.text(0.7, 0.04, '$A_V$')
        fig.text(0.04, 0.25, 'Abun', rotation = 90)
        fig.text(0.04, 0.75, 'T(K)', rotation = 90)
        fig.text(1.0 - 0.08, 0.25, 'Abun', rotation = 90)
        fig.text(1.0 - 0.08, 0.75, 'Abun', rotation = 90)


        axs[0,0].axis([0, 20, 1, 10000])
        axs[0,0].semilogy(data['state']['Av'],  data['state']['gasT'] , 'r' )
        axs[0,0].semilogy(data['state']['Av'],  data['state']['dustT'], 'b' )
        fig.text(0.4, 0.85, '$T_{gas}$' , color='r')
        fig.text(0.4, 0.82, '$T_{dust}$', color='b')

        axs[0,1].axis([0, 20, 1e-12, 2])
        axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH_P], 'r' )
        axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH],   'g' )
        axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xH2],  'b' )
        axs[0,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xe_N], 'c' )
        fig.text(0.85, 0.84, '$H^+$' , color='r')
        fig.text(0.85, 0.81, '$H$'   , color='g')
        fig.text(0.85, 0.78, '$H_2$' , color='b')
        fig.text(0.85, 0.75, '$e^-$' , color='c')


        axs[0,1].axis([0, 20, 1e-12, 2])
        axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xC_P], 'r' )
        axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xC],   'g' )
        axs[1,0].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xCO],  'b' )
        fig.text(0.4, 0.44, '$C^+$' , color='r')
        fig.text(0.4, 0.41, '$C$'   , color='g')
        fig.text(0.4, 0.38, '$CO$'  , color='b')

        axs[0,1].axis([0, 20, 1e-12, 2])
        axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHCN], 'r' )
        axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHNC],   'g' )
        axs[1,1].semilogy(data['state']['Av'],  data['state']['abun'][eSpcs.xHCO_P],  'b' )
        fig.text(0.8, 0.44, '$HCN$'   , color='r')
        fig.text(0.8, 0.41, '$HNC$'   , color='g')
        fig.text(0.8, 0.38, '$HCO^+$' , color='b')        