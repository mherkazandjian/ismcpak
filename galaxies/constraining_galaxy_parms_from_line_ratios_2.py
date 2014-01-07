'''
In this routine, we compute the Xi2 statistic for line ratios by using available emission data from teh PDR
grids. (no interpolationis done).

  - load pdr database
  - load computed line ratios as a function of beam size
  - load the emission for the specified line ratios to be included in the Xi2 statistic
  - compute the Xi2
  
'''
#########################################################################################################
import os, sys, time, pickle
import scipy
import numpy
import pylab

import meshUtils
from mylib.utils.interpolation import sectioned_4D_interpolator
import line_ratio_utils
import constraining
#======================================================parameters=================================================

home = '/home/mher'

params = {#'rundir': home + '/ism/runs/galaxies/coset2run4/coset-2-std', # the path of the dir containing the simulation
          'rundir': home + '/ism/runs/galaxies/coset2run4/coset-9-sol',  # the path of the dir containing the simulation
          'imres' : 100,                                                 # resolution of the maps to be produced imres x imres
          'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-1.0-low-res/',   # the path to the dir containing the PDR database
          #'pdrDb' : home + '/ism/runs/oneSided/sph-db-z-0.2-low-res/',  # the path to the dir containing the PDR database
          #'pdrDb' :  home + '/ism/runs/oneSided/dynamicMesh-z-1.0/',
                    
          'snap_index'  : numpy.arange(4, 4 + 1, 1),           

          
#          'line_ratios' : [
#                           'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 
#
#                           '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0',
#
#                           '13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO4-3/CO1-0', '13CO5-4/CO1-0', '13CO6-5/CO1-0',
#                          ],
          
          'lines'       : {
                           'include'     : [
                                            'CO1-0', 'CO2-1', 'CO3-2', 'CO4-3', 'CO5-4', 'CO6-5',                           
                                            '13CO1-0', '13CO2-1', '13CO3-2', '13CO4-3', '13CO5-4', '13CO6-5',                           
                                           ],
                           'combinations' :[
                                            'CO/CO', '13CO/13CO', '13CO/CO'
                                           ],
                         },
                        
          'beam_size'   : 2.0,
          'em_unit'     : 'fluxKkms',## plotting the line ratios and the modelled ones
          'error_bars'  : 0.2, 
          'save_info'   : False,
          'save_secies' : ['CO', '13CO'],
          #'interpolator' : scipy.interpolate.NearestNDInterpolator, 
          'interpolator' : scipy.interpolate.LinearNDInterpolator, 
          }
#############################################################################################################

## reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = params['pdrDb'], readDb=True)

if 'line_ratios' in params:
    line_ratios = params['line_ratios']
elif 'lines' in params:
    line_ratios = line_ratio_utils.line_ratio_combinations(params['lines']['include'], params['lines']['combinations'])

## loading the pickle file which holds the luminisity info
fname = os.path.join(params['rundir'],'analysis','line_luminosity_CO_13CO_snap_%d.pkl' % params['snap_index'][0])

fObj = open(fname, 'r')
    
luminosity = pickle.load(fObj )

fObj.close()

print 'read luminisities file :\n\t\t %s' % fname

## make line ratios from the mock luminosities 
obs_mock_ratios = line_ratio_utils.ratios()
beam_sizes = luminosity['beam_r']*2.0
ind = numpy.argmin(numpy.fabs(beam_sizes - params['beam_size']))
print 'available beam size used : %e' % beam_sizes[ind]
for line_ratio in line_ratios:

    line1, line2 = line_ratio_utils.lines_involved(line_ratio)
    
    v1, v2 = luminosity['lum_r'][line1][ind], luminosity['lum_r'][line2][ind]
    
    obs_mock_ratios.make_ratios(
                                {
                                  line1:{'fluxKkms': v1, 'err': params['error_bars']*v1}, 
                                  line2:{'fluxKkms': v2, 'err': params['error_bars']*v2}
                                },
                                ratios = [line_ratio],
                                em_unit = 'fluxKkms'
                               )

    print line_ratio, obs_mock_ratios[line_ratio]
    
obs_mock_ratios.species_and_codes()

## loading the emission info from all the models for all Avs (also check for the consistenscy of the 
## number of models...i.e same number of models for all the lines)
model_em = {}
for i, line in enumerate(obs_mock_ratios.codes):
    v, grid_coords = arxvPDR.get_emission_from_all_radex_dbs_for_Av_range(
                                                                          line = line, 
                                                                          Avs = 'all', 
                                                                          quantity = params['em_unit'],
                                                                          keep_nans = True,
                                                                         )
    model_em[line] = 10.0**v
    print line, v.size, grid_coords.shape
    print '----------------------'

    ## some checks of the sizes
    if v.size != grid_coords.shape[0]:
        raise ValueError('number of elements in the emission values is different from the number of modesl.')
    
    if i == 0:
        nModels = v.size
    else:
        if nModels != v.size:
            raise ValueError('the number of elements for this line differes at least from that of one of the other lines')

f = constraining.Xi2_line_ratios_single_component(obs_data = obs_mock_ratios, 
                                                  model_data = model_em, 
                                                  model_parms = grid_coords, 
                                                  line_ratios = line_ratios
                                                  )

f.compute_model_line_ratios()

f.compute_Xi2()

f.print_minima()

f.plot_results()
f.plot_results(True)

## making the interpolator for the Xi2 statistic as a function of the model coordinates
#fInterp = sectioned_4D_interpolator(grid_coords, Xi2, params['interpolator'])

#################
'''
Av =  18.0
gm = -23.0
+
n_i, g0_i = numpy.meshgrid(
                            numpy.linspace(-3.0, 4.0, 20),
                            numpy.linspace(-3.0, 4.0, 20)
                          )

nx, ny = n_i.shape
sz = n_i.size

n_i  = n_i.flatten()
g0_i = g0_i.flatten()
Av_i = numpy.ones(sz, 'f8')*Av
gm_i = numpy.ones(sz, 'f8')*gm

coords_interp = numpy.array([n_i, g0_i, gm_i, Av_i]).T
 
t0 = time.time()
Xi2_grd = fInterp(coords_interp)
dt = time.time() - t0

Xi2_grd = Xi2_grd.reshape((nx,ny))

#taking the log of the grid
Xi2_grd = numpy.log10(Xi2_grd)
mn, mx = numpy.nanmin(Xi2_grd), numpy.nanmax(Xi2_grd)
print 'min, max = ', mn, mx

im = pylab.imshow(
                  Xi2_grd,
                  extent=[-3.0, 4.0, -3.0, 4.0],
                  vmin=mn,
                  vmax=mx, 
                  interpolation='bessel', #intepolation used for imshow
                  origin='lower'
                 )

pylab.colorbar()
pylab.show()
'''
