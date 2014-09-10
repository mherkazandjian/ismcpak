""" plot a certain curve from multiple meshes on top of each other in the same plot"""
from mesh import *
from meshUtils import *
import numpy as np
from time import *

runDirPath =  '/home/mher/ism/runs/oneSided/runSingleMesh/'

# constructing the archive
t0 = time()
arxv = meshArxv(  )
arxv.construct( runDirPath )
arxv.plotCurvesFromMeshes(qx = ['state','Av'], qy = ['state', 'gasT'], qz = ['hdr', 'gammaMech'])
print 'done'

