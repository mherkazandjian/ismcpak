from amuse.community.fi.interface import FiMap

from amuse.io import read_set_from_file

from amuse.io.fi_io import FiFileFormatProcessor

from amuse.units import units,nbody_system

from matplotlib import pyplot

import numpy

import time

fname = "/data1/mher/ism/runs/galaxies/coset2run4/coset-6-std/firun/fiout.000020"
#parts=read_set_from_file("/home/mher/ism/runs/galaxies/test1/test/test.000000",format = FiFileFormatProcessor)
#parts=read_set_from_file("/data1/mher/ism/runs/galaxies/coset2run4/coset-2-std/firun/fiout.000000",format = FiFileFormatProcessor)
parts=read_set_from_file(fname, format = FiFileFormatProcessor)
#asdadsa
gas=parts[0]
convert=nbody_system.nbody_to_si(1. | units.kpc,1.e9 | units.MSun)

gas = gas[::1]

t0 = time.time()

#mapper=FiMap(convert,mode="openmp") #uses all the cores
mapper=FiMap(convert) #uses one core

mapper.parameters.image_width=20| units.kpc
mapper.parameters.projection_mode="perspective"  #as if it is being viewed from a camera perspective

mapper.parameters.projection_direction=[0,0,-1]
mapper.parameters.upvector=[0,1,0]
mapper.parameters.viewpoint= [0., 0., 10.] | units.kpc
mapper.parameters.image_size=[1024,1024]

print mapper.parameters

gas.weight=gas.mass.number

mapper.particles.add_particles(gas)

print mapper.image

im=mapper.image.pixel_value

vmax=im.max()

print 'time generating map = ', time.time() - t0

pyplot.imshow(numpy.log10(im))

pyplot.show()
