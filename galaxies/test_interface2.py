# in this test file, we load a database of PDR meshes and extract information
# at conditions determined by the SPH simulations
import time
import sys
import os
import subprocess

import matplotlib
matplotlib.use('Qt4Agg')

from scipy import interpolate
import numpy
from numpy import log10
import pylab

from amuse.io import read_set_from_file
from amuse.io.fi_io import FiFileFormatProcessor
from amuse.units import units, constants
import meshUtils
from mylib.utils.misc  import xselect
from mylib.utils.histogram import hist_nd
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------parameters--------------------------------------

home = '/home/mher'
#snapshots run dir
dirname = home + "/ism/runs/galaxies/coset2run4/coset-9-sol"

dtime = 0.0625

#PDR database file
#pdrDatabaseDirPath =  home + '/ism/runs/oneSided/uniformSweep2-z-1.0/'
#pdrDatabaseDirPath =  home + '/ism/runs/oneSided/dynamicMeshTest1/'
pdrDatabaseDirPath =  home + '/ism/runs/oneSided/sph-db-z-1.0-tmp/'

#reading and setting up the pdr database
arxvPDR = meshUtils.meshArxv(dirPath = pdrDatabaseDirPath, readDb=True)

snaps = numpy.arange(0,5,1)

for snap in snaps:
    
    #snapshot index
    snapIndex = snap
    
    suffix = '%06d' % snapIndex
    #path to fi snapshots
    snapName = 'fiout.%s' % suffix 
    filename = dirname + '/firun/' + snapName 
    
    snap_time = float(snapIndex)*dtime
    
    #loading the sph simulation data
    print 'loading snapshot %s : ' % filename
    gas, dark, stars = read_set_from_file(filename, format = FiFileFormatProcessor)
    print 'done reading fi snapshot : %s' % filename
    print 'number of sph particles in snapshot = %d' %  len(gas)
    
        
    """
    #getting an interpolated quantity from the PDR meshes for the states of the SPH gas
    #particles.
    dataInterp = arxvPDR.grid_interpolator(gas.rho, gas.fuvheat, gas.dethdt,
                                           quantity = ['state', 'gasT'], slabIdx = 0,)
    """
    
    #plotting the temperature vs the gas density
    meanmwt=1.3|units.amu
    n_gas_cgs = (gas.rho/meanmwt).value_in( units.cm **-3 )
    gasT = gas.temperat.value_in(units.K)
    G0 = 6.54*gas.fuvheat.value_in( units.none )
    g_mech = gas.dethdt.as_quantity_in( units.erg / (units.g * units.s ) )
    # converting the mechanical heating rate from per unit mass (of H gas)     
    # to per unit volume (see notesISM.odt)
    g_mech = g_mech.value_in( g_mech.unit ) # gMech in erg / (g s )
    g_mech = g_mech * 1.6474e-24 * n_gas_cgs # gMech in erg / (cm^3 s) 
    m_cgs = gas.mass.value_in(units.g)
    
    #plotting the gmech vs n_gas of the sph particles
    """
    pylab.loglog(n_gas_cgs[0::100], g_mech[0::100],'.')
    pylab.xlim([1, 1e6])
    pylab.ylim([1e-26, 1e-16])
    pylab.ylabel(r'$\Gamma$ [erg cm$^{-3}$ s$^{-1}]')
    pylab.xlabel(r'n [cm$^{-3}$]')
    """
    #info = xselect(x = numpy.log10(n_gas_cgs), y = numpy.log10(gasT) )
    
    n_gas_pdr = 10.0**arxvPDR.grid_x
    gasT_pdr = arxvPDR.getQuantityFromAllMeshes(['state', 'gasT'], slabIdx = 0)
    
    fig1 = pylab.figure()
    pylab.loglog(n_gas_cgs, gasT, 'r.')
    pylab.loglog(n_gas_pdr, gasT_pdr, 'b.')
    #pylab.show()
    
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------------
    #plotting the galaxy spatial (particles above a certain density)
    #---------------------------------------------------------------
    #lnMin, lG0Min, lgMechMin = 0.0, 0.0, -26.0   #log 10 of n, G0, gMech of SPH particles to keep
    lnMin, lG0Min, lgMechMin = -3.0, -3.0, -50.0  #database minimum boundaries
    bsMin, bsMax = -4.0, 4.0
    showEm       = True
    tansition    = 10 #the index of hte transition, for example: 0 => 1-0
     
    #plotting parms
    nBins          = 100    #mesh resolution use in binning
    plt_every_part = 100
    plt_ln_rng     = [lnMin,4] 
    plt_lg0_rng    = [lG0Min,3] 
    plt_lgm_rng    = [-35, -22]
    plt_em_rng     = [-10.0, -2.0] 
    #------------------------------------------
    
    #=============================================
    #getting the interpolation functio for CO(1-0)
    #loading the molecular specie radex database
    arxvPDR.readDbsRadex(species=['CO','HCN','HNC','HCO+'],Av=10.0)
    
    def funcRadex(mesh, **kwargs):
        
        if mesh == None:
            return numpy.nan
        else:
            return mesh[kwargs['transitionIdx']]['fluxcgs']  #a radex quantity
    
    def funcPDR(meshObj, **kwargs):
        return (1.0/(2.0*numpy.pi))*meshObj.compute_integrated_quantity(kwargs['quantity'], Av_range = kwargs['Av_range'])
    
    """
    v = arxvPDR.apply_function_to_all_meshes(funcPDR, 
                                             func_kw={'quantity': ['fineStructureCoolingComponents','C+','rate','1-0'],
                                                      'Av_range': [0.0, 10.0]
                                                     }
                                            )
    """
    
    def get_interpolation_function_radex(specStr, transIdx):
        
        arxvPDR.use_radexDb(specStr=specStr, Av=10.0) #.. todo:: do this for different Avs too    
    
        v = arxvPDR.apply_function_to_all_radex_meshes(funcRadex, func_kw={'transitionIdx':transIdx})
    
        v = numpy.array(v)
        inds_valid = numpy.isfinite(v)
        v = v[inds_valid] 
        xGrd, yGrd, zGrd = arxvPDR.grid_x[inds_valid], arxvPDR.grid_y[inds_valid], arxvPDR.grid_z[inds_valid]
        data = numpy.array([xGrd, yGrd, zGrd], dtype = numpy.float64).T
        fInterp = interpolate.NearestNDInterpolator(data, numpy.log10(v))
        return fInterp
    
    def make_interpolation_functions():
            
        funcs = {}
        
        funcs['CO-1-0']   = get_interpolation_function_radex('CO', 0)
        funcs['CO-3-2']   = get_interpolation_function_radex('CO', 2)
        funcs['CO-7-6']   = get_interpolation_function_radex('CO', 6)
        funcs['CO-16-15'] = get_interpolation_function_radex('CO', 15)
    
        funcs['HCN-1-0']  = get_interpolation_function_radex('HCN', 0)
        funcs['HCN-3-2']  = get_interpolation_function_radex('HCN', 2)
        
        funcs['HNC-1-0']  = get_interpolation_function_radex('HNC', 0)
        funcs['HNC-3-2']  = get_interpolation_function_radex('HNC', 2)
        
        funcs['HCO+-1-0'] = get_interpolation_function_radex('HCO+', 0)
        funcs['HCO+-3-2'] = get_interpolation_function_radex('HCO+', 2)
    
        return funcs
    
    print 'getting interpolation functions of the emissions'
    funcs = make_interpolation_functions()
    print 'done getting interpolation functions of the emissions'
    
    #=============================================
    
    interp = 'bessel' #intepolation used for imshow
    
    xkpc = gas.x.value_in(units.parsec)/1e3
    ykpc = gas.y.value_in(units.parsec)/1e3
    zkpc = gas.z.value_in(units.parsec)/1e3
    
    #setting the lower bound of the mechanical heating of the sph particles of the 
    #minimum of the pdr database which is almost negligable (since some of the 
    #sph particles have a zero gmech, their log would be -inf, so we clip it to the
    #mininmum of the pdr database
    g_mech = g_mech.clip(10.0**lgMechMin)
    
    #selecting particles within a spatial range and within the 
    #density and g0 range of the database
    inds = numpy.where(
                       (numpy.log10(n_gas_cgs) >= lnMin)*
                       (numpy.log10(G0) >= lG0Min)*
                       (numpy.log10(g_mech) >= lgMechMin)*                   
                       (xkpc > bsMin)*
                       (xkpc < bsMax)*
                       (ykpc > bsMin)*
                       (ykpc < bsMax)
                      )[0]
    
    x1,y1, z1 = xkpc[inds], ykpc[inds], zkpc[inds]
    n1, g01, g_mech1 = n_gas_cgs[inds], G0[inds], g_mech[inds]
    m1 = m_cgs[inds]
    
    #computing alpha for each sph particle (amount of gm wrt surface heating for the PDR model)
    f_log_gamma_surf=arxvPDR.construct3DInterpolationFunction(quantity=['therm','heating'], slabIdx=0, log10=True, interpolator='nearest') 
    f_log_Tkin_surf=arxvPDR.construct3DInterpolationFunction(quantity=['state','gasT'], slabIdx=0, log10=True, interpolator='nearest') 
    dataNew = numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.ones(n1.shape)*lgMechMin]).T
    gammaSurf_sph_from_pdr = 10.0**f_log_gamma_surf(dataNew)
    Tkin_sph_from_pdf = 10.0**f_log_gamma_surf(dataNew)
    alpha_sph = g_mech1/gammaSurf_sph_from_pdr
    
    indsNan = numpy.where(numpy.isnan(Tkin_sph_from_pdf))[0]
    
    #selecting particles which have alpha < 1.0
    
    if showEm:
        print 'getting the emissions of each sph particle'
        #getting emissions with mechanical heating    
        emissions = {}
        data =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.log10(g_mech1)]).T
        
        emissions['CO-1-0']   = 10.0**funcs['CO-1-0'](data)
        emissions['CO-3-2']   = 10.0**funcs['CO-3-2'](data)
        emissions['CO-16-15'] = 10.0**funcs['CO-16-15'](data)
        
        emissions['HCN-1-0'] = 10.0**funcs['HCN-1-0'](data)
        emissions['HCN-3-2'] = 10.0**funcs['HCN-3-2'](data)
        
        emissions['HNC-1-0'] = 10.0**funcs['HNC-1-0'](data)
        emissions['HNC-3-2'] = 10.0**funcs['HNC-3-2'](data)
        
        emissions['HCO+-1-0'] = 10.0**funcs['HCO+-1-0'](data)
        emissions['HCO+-3-2'] = 10.0**funcs['HCO+-3-2'](data)
    
        #getting emissions without mechanical heating    
        emissions_no_gm = {}    
        data =  numpy.array([numpy.log10(n1), numpy.log10(g01), numpy.ones(n1.size)*lgMechMin]).T
    
        emissions_no_gm['CO-1-0']   = 10.0**funcs['CO-1-0'](data)
        emissions_no_gm['CO-3-2']   = 10.0**funcs['CO-3-2'](data)
        emissions_no_gm['CO-16-15'] = 10.0**funcs['CO-16-15'](data)
        
        emissions_no_gm['HCN-1-0'] = 10.0**funcs['HCN-1-0'](data)
        emissions_no_gm['HCN-3-2'] = 10.0**funcs['HCN-3-2'](data)
    
        emissions_no_gm['HNC-1-0'] = 10.0**funcs['HNC-1-0'](data)
        emissions_no_gm['HNC-3-2'] = 10.0**funcs['HNC-3-2'](data)
        
        emissions_no_gm['HCO+-1-0'] = 10.0**funcs['HCO+-1-0'](data)
        emissions_no_gm['HCO+-3-2'] = 10.0**funcs['HCO+-3-2'](data)
        print 'done getting the emissions of each sph particle'
    
    #else:
    #    CO = numpy.zeros(n1.size)
    #    CO_no_gm = numpy.zeros(n1.size)
    
    #--------------------------------------------------------------------------------
    #getting the distribution of mass, particle count as a function of number density
    #--------------------------------------------------------------------------------
    print 'getting the distributions as a function of gas density'
    nDens_hist = hist_nd(numpy.log10(n1.reshape((1,n1.shape[0]))),  loc=True,
                         nbins = nBins, mn = plt_ln_rng[0], mx = plt_ln_rng[1], reverse_indicies = True)
    fMass_vs_n = numpy.zeros((nBins), dtype=numpy.float64) 
    fLumCO_vs_n = numpy.zeros((nBins), dtype=numpy.float64) 
    fLumCO_no_gm_vs_n = numpy.zeros((nBins), dtype=numpy.float64) 
    
    for i in numpy.arange(nBins):
        inds_in_bin = nDens_hist.get_indicies([i])
    
        if inds_in_bin.size > 0:
            fMass_vs_n[i] = m1[inds_in_bin].sum()
            fLumCO_vs_n[i]  = emissions['CO-1-0'][inds_in_bin].sum()
            fLumCO_no_gm_vs_n[i] = emissions_no_gm['CO-1-0'][inds_in_bin].sum()
            
    #normalizing the computed distributions to unity
    fn = nDens_hist.f/nDens_hist.f.sum()
    fMass_vs_n /= m1.sum()
    fLumCO_vs_n /= fLumCO_vs_n.sum()
    fLumCO_no_gm_vs_n /= fLumCO_no_gm_vs_n.sum()
    print 'done getting the distributions as a function of gas density'
    
    fig2, axs2 = pylab.subplots(1, 3, sharex=True, sharey=False, figsize=(15,5),
                                subplot_kw = {'xlim':plt_ln_rng,
                                              'ylim':[-4.0,0.0],
                                              'xlabel':r'$\log_{10} n [cm^{-3}]$'})
    pylab.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.5, hspace=0.0)
    
    
    axs2[0].set_ylabel(r'$\int_0^n f(n) dn$', fontsize='xx-large')
    axs2[0].plot(nDens_hist.f.cntrd, numpy.log10(fn.cumsum()))
    
    axs2[1].set_ylabel(r'$\int_0^n f_{M}(n) dn$', fontsize='xx-large')
    axs2[1].plot(nDens_hist.f.cntrd, numpy.log10(fMass_vs_n.cumsum()))
    
    axs2[2].set_ylabel(r'$\int_0^n f_{L_{CO(%d,%d)}}(n) dn$' % (tansition+1, tansition), fontsize='xx-large')
    axs2[2].plot(nDens_hist.f.cntrd, numpy.log10(fLumCO_vs_n.cumsum()),'b')
    axs2[2].plot(nDens_hist.f.cntrd, numpy.log10(fLumCO_no_gm_vs_n.cumsum()),'r')
    #---------------------------------------------------------------------------------------
    # Done getting the distribution of mass, particle count as a function of number density
    #---------------------------------------------------------------------------------------
    print 'getting the spatial distrubutions'
    
    f_sph_part = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    f_m        = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    f_mean_n   = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    f_mean_g0  = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    f_mean_gm  = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    f_mean_em  = {}
    f_mean_em_no_gm  = {}
    for key in emissions:
        f_mean_em[key] = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    for key in emissions_no_gm:
        f_mean_em_no_gm[key] = numpy.zeros((nBins, nBins), dtype=numpy.float64)
        
    #f_mean_CO = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    #f_mean_CO_no_gm = numpy.zeros((nBins, nBins), dtype=numpy.float64)
    
    hist = hist_nd(numpy.vstack((x1,y1)), mn = bsMin, mx=bsMax, nbins=nBins, reverse_indicies=True) 
    for i in numpy.arange(nBins):
        for j in numpy.arange(nBins):
            inds_in_bin = hist.get_indicies([i,j])
    
            f_sph_part[i,j] = inds_in_bin.size
            if inds_in_bin.size > 0:
                f_m[i,j] = numpy.log10(m1[inds_in_bin].sum())
                f_mean_n[i,j]  = numpy.log10(numpy.mean(n1[inds_in_bin]))
                f_mean_g0[i,j] = numpy.log10(numpy.mean(g01[inds_in_bin]))
                f_mean_gm[i,j] = numpy.log10(numpy.mean(g_mech1[inds_in_bin]))
                for key in emissions:
                    f_mean_em[key][i,j] = max(numpy.log10(numpy.mean(emissions[key][inds_in_bin])),plt_em_rng[0])
                for key in emissions_no_gm:
                    f_mean_em_no_gm[key][i,j] = max(numpy.log10(numpy.mean(emissions_no_gm[key][inds_in_bin])),plt_em_rng[0])
                #f_mean_CO[i,j] = max(numpy.log10(numpy.mean(CO[inds_in_bin])),plt_em_rng[0])
                #f_mean_CO_no_gm[i,j] = max(numpy.log10(numpy.mean(CO_no_gm[inds_in_bin])), plt_em_rng[0])
            else:
                f_sph_part[i,j] = 1
                f_m[i,j] = 0.0
                f_mean_n[i,j]  = lnMin
                f_mean_g0[i,j] = lG0Min
                f_mean_gm[i,j] = lgMechMin
                for key in emissions:
                    f_mean_em[key][i,j] = plt_em_rng[0]
                for key in emissions_no_gm:
                    f_mean_em_no_gm[key][i,j] = plt_em_rng[0]
    
                #f_mean_CO[i,j] = plt_em_rng[0]
                #f_mean_CO_no_gm[i,j] = plt_em_rng[0]
    
    #print 'no gmech   (log10(gm) = %f): log10[CO] (min,max) = (%f,%f)' % (lgMechMin, f_mean_CO.min(),f_mean_CO.max()) 
    #print 'with gmech (from sph parts): log10[CO] (min,max) = (%f,%f)' % (f_mean_CO_no_gm.min(), f_mean_CO_no_gm.max()) 
    
    print 'done getting the spatial distributuions'
    
    fig3, axs3 = pylab.subplots(4, 4, sharex=True, sharey=True, figsize=(16,16), 
                               subplot_kw = {'xlim':[bsMin, bsMax],
                                             'ylim':[bsMin,bsMax],
                                             'aspect':'equal',
                                             'adjustable':'datalim',
                                             }
                               )
    for ax in axs3[:,0]: ax.set_ylabel('y(kpc)')
    for ax in axs3[3,:]: ax.set_xlabel('x(kpc)')
    
    pylab.subplots_adjust(left=0.05, bottom=0.05, right=0.9, top=0.95, wspace=0.3, hspace=0.3)
    
    plt00, = axs3[0,0].plot(x1[::plt_every_part], y1[::plt_every_part], '.', markersize=1)
    
    im01   = axs3[0,1].imshow(numpy.log10(f_sph_part), extent=[bsMin, bsMax, bsMin, bsMax], vmin=0, vmax=6, interpolation=interp)
    ttl01  = axs3[0,1].set_title(r'$f(N)$', size='large')
    cbar01 = pylab.colorbar(im01, ax=axs3[0,1], orientation='vertical')
    
    im02   = axs3[0,2].imshow(f_m, extent=[bsMin, bsMax, bsMin, bsMax], vmin=30, vmax=40, interpolation=interp)
    ttl02  = axs3[0,2].set_title(r'$f(m)$', size='large')
    cbar02 = pylab.colorbar(im02, ax=axs3[0,2], orientation='vertical')
    
    im03   = axs3[0,3].imshow(f_mean_n       , extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_ln_rng[0] , vmax=plt_ln_rng[1] , interpolation=interp)
    ttl03  = axs3[0,3].set_title(r'$f(\bar{n})$', size='large')
    cbar03 = pylab.colorbar(im03, ax=axs3[0,3], orientation='vertical')
    
    im10   = axs3[1,0].imshow(f_mean_g0      , extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_lg0_rng[0], vmax=plt_lg0_rng[1], interpolation=interp)
    ttl10  = axs3[1,0].set_title(r'$f(\bar{g_0})$', size='large')
    cbar10 = pylab.colorbar(im10, ax=axs3[1,0], orientation='vertical')
    
    im11   = axs3[1,1].imshow(f_mean_gm      , extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_lgm_rng[0], vmax=plt_lgm_rng[1], interpolation=interp)
    ttl11  = axs3[1,1].set_title(r'$f(\bar{\Gamma_m})$', size='large')
    cbar11 = pylab.colorbar(im11, ax=axs3[1,1], orientation='vertical')
    
    #CO 1-0
    im12   = axs3[1,2].imshow(f_mean_em['CO-1-0'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    tt12   = axs3[1,2].set_title(r'$f(L_{CO(1-0})$', size='large')
    cbar12 = pylab.colorbar(im12, ax=axs3[1,2], orientation='vertical')
    
    im13   = axs3[1,3].imshow(f_mean_em_no_gm['CO-1-0'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    ttl13  = axs3[1,3].set_title(r'$f(L_{CO(1-0)})$ $\Gamma_m = 0$', size='large')
    cbar13 = pylab.colorbar(im13, ax=axs3[1,3], orientation='vertical')
    
    #CO-3-2
    im20   = axs3[2,0].imshow(f_mean_em['CO-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    tt20   = axs3[2,0].set_title(r'$f(L_{CO(3-2})$', size='large')
    cbar20 = pylab.colorbar(im20, ax=axs3[2,0], orientation='vertical')
    
    im21   = axs3[2,1].imshow(f_mean_em_no_gm['CO-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    ttl21  = axs3[2,1].set_title(r'$f(L_{CO(3-2)})$ $\Gamma_m = 0$', size='large')
    cbar21 = pylab.colorbar(im21, ax=axs3[2,1], orientation='vertical')
    
    #HCN-3-2
    im22   = axs3[2,2].imshow(f_mean_em['HCN-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    tt22   = axs3[2,2].set_title(r'$f(L_{HCN(3-2})$', size='large')
    cbar22 = pylab.colorbar(im22, ax=axs3[2,2], orientation='vertical')
    
    im23   = axs3[2,3].imshow(f_mean_em_no_gm['HCN-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    ttl23  = axs3[2,3].set_title(r'$f(L_{HCN(3-2)})$ $\Gamma_m = 0$', size='large')
    cbar23 = pylab.colorbar(im23, ax=axs3[2,3], orientation='vertical')
    
    #HNC-3-2
    im30   = axs3[3,0].imshow(f_mean_em['HNC-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    tt30   = axs3[3,0].set_title(r'$f(L_{HNC(3-2})$', size='large')
    cbar30 = pylab.colorbar(im30, ax=axs3[3,0], orientation='vertical')
    
    im31   = axs3[3,1].imshow(f_mean_em_no_gm['HNC-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    ttl31  = axs3[3,1].set_title(r'$f(L_{HNC(3-2)})$ $\Gamma_m = 0$', size='large')
    cbar31 = pylab.colorbar(im31, ax=axs3[3,1], orientation='vertical')
    
    #HCO+-3-2
    im32   = axs3[3,2].imshow(f_mean_em['HCO+-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    tt32   = axs3[3,2].set_title(r'$f(L_{HCO+(3-2})$', size='large')
    cbar32 = pylab.colorbar(im32, ax=axs3[3,2], orientation='vertical')
    
    im33   = axs3[3,3].imshow(f_mean_em_no_gm['HCO+-3-2'], extent=[bsMin, bsMax, bsMin, bsMax], vmin=plt_em_rng[0] , vmax=plt_em_rng[1] , interpolation=interp)
    ttl33  = axs3[3,3].set_title(r'$f(L_{HCO+(3-2)})$ $\Gamma_m = 0$', size='large')
    cbar33 = pylab.colorbar(im33, ax=axs3[3,3], orientation='vertical')
    
    pylab.figtext(0.05, 0.97, '%.3f Gyr' % snap_time, size='xx-large')
    image_fname =  dirname + '/analysis/maps.%s.png' % suffix
    fig3.savefig(image_fname)
    print 'wrote image file : %s' % image_fname
    
    #subprocess.Popen(['xv', image_fname])
    #pylab.show()
