# -*- coding: utf-8 -*-
'''
<keywords>
ism, line ratios, fitting, minimization, ismcpak, one model
</keywords>

<description>
given a bunch of line ratios, finds the best fitting model from a database of PDR meshes
</description>
'''

import numpy
import lineDict
from scipy import interpolate
import scipy
import time

import os
import matplotlib
import sys
matplotlib.use('Qt4Agg')
import line_ratio_utils
import pylab    
import meshUtils

#--------------------------------------------------------------------------------------------------------------
home           = '/home/mher'
dirPath        = os.path.join(home, 'ism/runs/oneSided/dynamicMesh-z-1.0/')
image_save_dir = '/home/mher/tmp'
data_sheet     = '/home/mher/ism/marissa/NGC_253_data/NGC_253_Fluxes_Master_newsheet.xlsx' #path to the excel sheet containing the data
#--------------------------------------------------------------------------------------------------------------

ratios  = () 
#ratios += ('CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 'CO7-6/CO1-0', 'CO8-7/CO1-0', 'CO9-8/CO1-0', 'CO10-9/CO1-0', 'CO11-10/CO1-0', 'CO12-11/CO1-0', 'CO13-12/CO1-0')
ratios += ('CO4-3/CO3-2', 'CO5-4/CO3-2', 'CO6-5/CO3-2', 'CO7-6/CO3-2', 'CO8-7/CO3-2', 'CO9-8/CO3-2', 'CO10-9/CO3-2', 'CO11-10/CO3-2', 'CO12-11/CO3-2', 'CO13-12/CO3-2')
#ratios += ('CO6-5/CO5-4', 'CO7-6/CO5-4', 'CO8-7/CO5-4', 'CO9-8/CO5-4', 'CO10-9/CO5-4', 'CO11-10/CO5-4', 'CO12-11/CO5-4', 'CO13-12/CO5-4')
#ratios += ('13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',)
#ratios += ('13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO5-4/CO1-0', '13CO6-5/CO1-0',)
ratios += ('HCO+4-3/HCO+1-0', 'HCO+4-3/HCO+1-0',)
#ratios += ('CS3-2/CS2-1', 'CS4-3/CS2-1', 'CS5-4/CS2-1',)
ratios += ('HCO+4-3/CO1-0', 'HCO+4-3/CO4-3',)
ratios += ('HCN3-2/HCN1-0', 'HCN4-3/HCN1-0',)
ratios += ('HCN1-0/HNC1-0',)
ratios += ('HCO+1-0/HNC1-0',)
ratios += ('HCN4-3/HCO+4-3',)
ratios += ('HCO+4-3/HCO+1-0',)
ratios += ('HCO+7-6/HCO+1-0',)

ratios += ('13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',)
ratios += ('13CO5-4/13CO3-2', '13CO6-5/13CO3-2',)

Av_rng = [0.0, 100.0]

###########################################################################################################

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = dirPath, readDb=True)

#reading the lines from the excel sheet
obs = line_ratio_utils.read_observations(data_sheet)

#making the observed line ratio
obs_ratios = line_ratio_utils.ratios()
obs_ratios.make_ratios(obs, ratios)
    
#loading the lines info from the PDR database
specStrs, codes = obs_ratios.species_and_codes()
arxvPDR.readDbsRadex(species=specStrs, in_Av_rng=Av_rng)

#getting the Xi2 for all the models in the PDR database
data = line_ratio_utils.Xi2_line_ratios(obs_ratios, arxvPDR)

#finding some minima about the Xi2 from all the models
ind_min = data['v'].argmin()
x_min, y_min, z_min, t_min = data['x'][ind_min], data['y'][ind_min], data['z'][ind_min], data['t'][ind_min]
Xi_min = data['v'][ind_min]
print 'the global minimum value of the Xi2 = %e\n         (x,y,z,t) = (%f, %f, %f, %f)' % (Xi_min, x_min, y_min, z_min, t_min)

line_ratio_utils.plot_single_model_ratios(arxvPDR, x_min, y_min, z_min, t_min, obs_ratios)
pylab.title(pylab.gca().get_title() + 'model with the global minimum Xi2')

#getting the pdr model with minimum Xi2 withouth mechanical heating
zero_gmech = -10.0
inds_zero_gmech = numpy.where( numpy.fabs(1.0 - data['z']/zero_gmech) < 1e-6 )[0]
data_zero_gmech = data[inds_zero_gmech]
ind_min_zero_gmech = data_zero_gmech['v'].argmin()
x_min_zero_gmech, y_min_zero_gmech, z_min_zero_gmech, t_min_zero_gmech = data_zero_gmech['x'][ind_min_zero_gmech], data_zero_gmech['y'][ind_min_zero_gmech], data_zero_gmech['z'][ind_min_zero_gmech], data_zero_gmech['t'][ind_min_zero_gmech]
Xi_min_zero_gmech = data_zero_gmech['v'][ind_min_zero_gmech]
print 'the minimum (without mechanical heating) of the Xi2 = %e\n         (x,y,z,t) = (%f, %f, %f, %f)' % (Xi_min_zero_gmech, x_min_zero_gmech, y_min_zero_gmech, z_min_zero_gmech, t_min_zero_gmech)

line_ratio_utils.plot_single_model_ratios(arxvPDR, x_min_zero_gmech, y_min_zero_gmech, z_min_zero_gmech, t_min_zero_gmech, obs_ratios)
pylab.title(pylab.gca().get_title() + 'model (without gmech) with the global minimum Xi2')


# getting the interpolation function
ti = time.time()

#f_xi2_interp = interpolate.LinearNDInterpolator(
#                                                 numpy.vstack([data['x'], data['y'], data['z'], data['t']]).T, 
#                                                 data['v']
#                                                )#getting the pdr model with the minimum Xi2

f_xi2_interp = interpolate.NearestNDInterpolator(
                                                 numpy.vstack([data['x'], data['y'], data['z'], data['t']]).T, 
                                                 data['v']
                                                )

tf = time.time()
print 'constructed he interpolation function from %d points in %f seconds' % (data.size, tf-ti)

# defining the points in the uniform 2D grid                                                                                                                                                                                                 
grid_x, grid_y = numpy.mgrid[0.0:6.0:50j, 0.0:6.0:50j] 

nPts = numpy.product(grid_x.shape)
xNew = grid_x.reshape(nPts)
yNew = grid_y.reshape(nPts)
zNew = xNew.copy()
tNew = xNew.copy()
        
zNew[:] =  z_min
tNew[:] =  t_min
dataNew = numpy.array( [xNew, yNew, zNew, tNew] ).T            

ti = time.time()
tNew = f_xi2_interp(dataNew)
tf = time.time()

print 'interpolated %d points in %f seconds at a rate of %e pts/sec' % (nPts, tf-ti, nPts / (tf-ti))
tNew = numpy.reshape(tNew, grid_x.shape)

grd = tNew

#getting the location of the location of the minimum 
ind = grd.flatten().argmin()
x_min  = grid_x.flatten()[ind]
y_min  = grid_y.flatten()[ind]
Xi_min = grd.flatten()[ind]


# plotting it
pylab.figure()
pylab.subplot(111)
    
im = pylab.imshow(numpy.log10(grd.T), extent=(0.0, 6.0, 0.0, 6.0), origin='lower')

pylab.colorbar(im, shrink = 0.8, extend = 'both')

print 'minimum value of the Xi2 = %e in this 2D grid is at (x,y) = (%f, %f)' % (Xi_min, x_min, y_min)
pylab.plot(x_min + 6.0 / 50.0, y_min + 6.0 / 50.0, 'ks', markersize=25)


pylab.show()

print 'done'
