import numpy as np

class guess_db():
    """a module which provides functionalities to prepare/read/write database files
    which can be used by the PDR code to guess values for the surface slab.
    """
    
    def __init__(self, path = None):
        
        if path != None:
            self.set_path(path)
        else:
            self.path = None #: path to the database file
        
        self.db = None #: dictionary that will hold the data of the database    
    
    def read(self):
        """read a guess db binary file (see documentation of the typedef 
        database in vars.h in the pdr code. A dict is retured with  the
        following keys : 
        
            {
             'n'       ,  # np.int32, number of points in the db 
             'nSpecs'  ,  # np.int32, number of species per point
             'logDens' ,  # np.ndarray(n, dtype = np.float64) densities of the points
             'logG0    ,  # np.ndarray(n, dtype = np.float64) G0 of the points
             'logGammaM,  # np.ndarray(n, dtype = np.float64) gmech of the points
             'Teq'     ,  # np.ndarray(n, dtype = np.float64) eq temp of the surface slab of the points
             'abun'    ,  # np.ndarray((nSpec, n), dtype = np.float64) eq abun of the surface slabs corresponding to each point
            }  
         raises a value error if the database integrity check fails.    
        """

        db = {}
        
        #-----------------------------------------------------------------
        def check_integrity(fObj):
            offset = np.fromfile(fObj, dtype = np.int64, count = 1)

            if len(offset) == 0:
                raise ValueError('End of file reached')
            else:
                offset = offset[0]
                
            if (fObj.tell() - np.int64(0).nbytes) != offset:
                raise ValueError('Database integrity failed.') 
        #-----------------------------------------------------------------
        
        print 'reading the file %s' % self.path
        fObj = open(self.path, 'rb')

        #reading the header        
        check_integrity(fObj)
        db['n']      = np.fromfile(fObj, dtype = np.int32, count = 1)[0]
        db['nSpecs'] = np.fromfile(fObj, dtype = np.int32, count = 1)[0]
        print 'will read %d nPts each with %d species abundances' % (db['n'], db['nSpecs'])
        
        #reading the coordinates of the points n,G0,Gmech
        check_integrity(fObj)
        db['logDens'] = np.fromfile(fObj, dtype = np.float64, count = db['n'])
        check_integrity(fObj)
        db['logG0'] = np.fromfile(fObj, dtype = np.float64, count = db['n'])
        check_integrity(fObj)
        db['logGammaM'] = np.fromfile(fObj, dtype = np.float64, count = db['n'])
        
        #reading the temperatures
        check_integrity(fObj)
        db['Teq'] = np.fromfile(fObj, dtype = np.float64, count = db['n'])
        
        #reading the abundances
        db['abun'] = np.ndarray((db['nSpecs'], db['n']), dtype = np.float64)
        for i in np.arange(db['n']):
            check_integrity(fObj)
            db['abun'][:,i] = np.fromfile(fObj, dtype = np.float64, count = db['nSpecs'])
            
        #final integrity check
        check_integrity(fObj)

        #setting the db, closing the file and returning the db
        self.set_db(db)
        fObj.close()
        return self.db 
    
    def write(self, outFile = None):
        """Writes the self.db dict to a binary file. The binary file has the follwing format:
                  offset|n,nSpec|offset|logDens|offset|logG0|offset|logGammaM|offset
                  |Teq|offset|abun0|offset|abun2|offset...|abun(n-1)|offset|
           If outFile is not provided, the db is written to self.path.
        """

        #setting the file path        
        if outFile == None:
            path = self.get_path()
        else:
            path = outFile
        
        #opening the file and getting the data of the db
        fObj = open(path, 'wb')
        db = self.db
        
        #writing the header
        np.int64(fObj.tell()).tofile(fObj);  np.array((db['n'], db['nSpecs'])).tofile(fObj) # header
        
        #writing the point coordinates
        np.int64(fObj.tell()).tofile(fObj);  db['logDens'].tofile(fObj)   #nGas
        np.int64(fObj.tell()).tofile(fObj);  db['logG0'].tofile(fObj)     #G0
        np.int64(fObj.tell()).tofile(fObj);  db['logGammaM'].tofile(fObj) #gmech
        
        #writing the temperatures
        np.int64(fObj.tell()).tofile(fObj);  db['Teq'].tofile(fObj)   #T_eq

        #writing the abundances
        for i in np.arange(db['n']):
            np.int64(fObj.tell()).tofile(fObj);  db['abun'][:,i].tofile(fObj)
            
        #writing final check num
        np.int64(fObj.tell()).tofile(fObj)

        fObj.close()
        print 'wrote the guess database file sucessfully : %s' % path
        
        return True 
        
    def set_db_data(self, logDens, logG0, logGammaM, Teq, abun):
        """sets the enteris in the database. When this is called,
        self.db is deleted and the entries are replaced with the passed values.
        All of the parameters should be numpy float64 arrays/int32 dtypes.
        """

        if abun.shape[1] != logDens.size:
            raise ValueError("the shape of db['abun'] should be (nSpecs, nPts).")
        
        del self.db
        db = {}
        
        db['n']         = np.int32(abun.shape[1])
        db['nSpecs']    = np.int32(abun.shape[0])
        db['logDens']   = np.float64(logDens)
        db['logG0']     = np.float64(logG0)
        db['logGammaM'] = np.float64(logGammaM)
        db['Teq']       = np.float64(Teq)
        db['abun']      = np.float64(abun)

        self.set_db(db)
                
    
    #setters and getters
    def set_path(self, path):
        self.path = path
    def get_path(self):
        return self.path
    def set_db(self, db):
        self.db = db
    def get_db(self):
        return db