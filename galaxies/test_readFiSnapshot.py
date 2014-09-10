from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor

home = '/home/mher'

filename = home + '/data1/mher/ism/runs/galaxies/coset2run4/coset-2-std/firun/fiout.000000'

data = read_set_from_file(filename, format = FiFileFormatProcessor)

gas = data[0]

p1 = gas[0]

#print p[0]
#print len(p[0]),len(p[1]),len(p[2])

print 'done'