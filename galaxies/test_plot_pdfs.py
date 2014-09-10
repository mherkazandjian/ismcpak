from mylib.utils import templates
from numpy import *
from pylab import *
from galaxies import fi_utils, coset9_sol_info
from mylib import units, constants

xmin, xmax = -2.0, 8.0 

fig = templates.subplots_grid(1, 1, hspace=0.0, wspace=0.0,
                              fig = {'kwargs':{
                                               'figsize' : {3.55, 2.9}
                                               }
                                     },
                              axs = {
                                     'left' :0.2, 'bottom': 0.2, 'w':0.75, 'h':0.75
                                    },
                              )

# density bins to be plotted
log10x, dlog10x  = linspace(xmin, xmax, 100, retstep=True, endpoint=False)
n_bins = 10.0**(log10x + dlog10x*0.5)

ax = fig.sub[0,0]


#--------------------------------------------------------------------------------------
## the probability density at that density (not in log scale) (wada2001)

wada_n_cgs = (array([10.0**-0.5, 1e2 , 1e4 , 1e6])*constants.M_SUN_CGS / constants.M_PROTON_CGS)/units.PC2CM**3.0
wada_pdf   = array([5e-2       , 8e-3, 2e-4, 1e-6])

ax.loglog(wada_n_cgs, wada_pdf, 'o' )

nPDF_obj = fi_utils.density_distribution('lognormal', mu=wada_n_cgs[0], sigma=22.4)
nPDF = nPDF_obj.pdf
ax.loglog(n_bins, nPDF(n_bins)/log(10), '-.', label=r'$\log_{10}(\sigma$) = 3.7')

#--------------------------------------------------------------------------------------

ax.legend(loc=0, prop={'size':6})

'''
fig.set_xticks([1e-1, 1e1, 1e3, 1e5, 1e7], 
               labels=[r'10$^{-1}$', r'10$^1$', r'10$^3$', r'10$^5$', r'10$^7$'], 
               size='x-small')
fig.set_yticks([1e-6, 1e-4, 1e-2, 1e-0], 
               labels=[r'10$^{-6}$', r'10$^{-4}$', r'10$^{-2}$', r'10$^{0}$'], 
               size='x-small')
'''     

fig.sub_setp('xlim', [10.0**xmin, 10.0**xmax])
fig.sub_setp('ylim', [1e-8, 1.0])

show()