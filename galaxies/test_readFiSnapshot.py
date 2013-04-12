from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor

home = '/home/mher'

filename = home + "/amuse/sandbox/mher/test/test.000000"

p = read_set_from_file(filename, format = FiFileFormatProcessor)
  
print p[0]
#print len(p[0]),len(p[1]),len(p[2])

print 'done'