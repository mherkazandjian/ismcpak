from amuse.community.pdr.pyUtils.guess_db import guess_db
import numpy as np

#dbPath = '/home/mher/ism/database/database2.dat'
#dbPath = '/home/mher/ism/database/1.dat'
dbPath = '/home/mher/tmp/1.dat'
pdrGuessDb = guess_db( path = dbPath )
data = pdrGuessDb.read()

dbPath2 = '/home/mher/tmp/foo.dat'
pdrGuessDb.write(dbPath2)
pdrGuessDb.set_path(dbPath2)

data2 = pdrGuessDb.read()

for key in data:
    x = np.sum(data[key] - data2[key])
    
print 'total different = ', x
