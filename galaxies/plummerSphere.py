import struct
import numpy

from collections import namedtuple

from amuse.io import base
from amuse.units import units
from amuse.units import nbody_system
from amuse.support.core import late

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor

from amuse import datamodel

from amuse.support import data
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.plummer import new_plummer_sphere
import numpy as np
import Gnuplot
import time

convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
plummer = new_plummer_sphere(1000, convert_nbody)
stars = plummer.copy()
plotter = Gnuplot.Gnuplot()
plotter.splot(plummer.position.value_in(units.parsec))

time.sleep(10)
