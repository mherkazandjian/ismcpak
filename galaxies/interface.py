import numpy
from meshUtils

class grid_interpolator():
    """returns an interpolated quantity from the pdr meshes"""
    
           
    px, py, pz = (gas[:].x.value_in(units.parsec), gas[:].y.value_in(units.parsec), gas[:].z.value_in(units.parsec))

    