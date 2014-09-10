# plots the emissions grid
#------------------------------------------------------------------------------------
import os
import matplotlib
import numpy

if 'particle3' in os.uname():
    matplotlib.use('Qt4Agg')
    
import meshUtils
#import plot_utils

home         = '/home/mher'
dirPath      = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
imageSaveDir = '/home/mher/tmp'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

arxvPDR.line_ratio_app()
