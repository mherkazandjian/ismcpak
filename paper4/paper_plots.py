from numpy import *
from pylab import *
from galaxies import fi_utils
from mylib.utils import templates
import lineDict

def plot_methods_fig(gas, save_figs, fig_paths):
    
    def plot_original(gas, fig_path):
        '''plotting the original distributuion'''

        fig = templates.subplots_grid(1,1, hspace=0.15, wspace=0.15,
                                      fig = {'kwargs':{
                                                       'figsize' : {3.55, 2.95}
                                                       }
                                             },
                                      axs = {
                                             'left' :0.20, 'bottom': 0.22, 'w':0.75, 'h':0.74
                                            },
                                      )

        #figure()
        ax = fig.sub[0,0]

        ######### plotting the original PDFs #############
        
        ## selecting the original gas particles
        gas_orig_amuse = gas[where(isnan(gas.parent) == True)]
        gas_orig = fi_utils.gas_set(len(gas_orig_amuse)) 
        gas_orig.copy_attr_from_amuse_set(gas_orig_amuse, ['n'])
        
        ## using the ln of the densities                
        x = log(gas_orig.n)

        ## parameters for the plots
        xmin, xmax = log(1e-8), log(1e6)
        nbins = 100 
        fit_func_rng  = [1e-2, 1e+4]
                
        ## computing the PDF of the original particles
        f, bins = histogram(x, range=[xmin, xmax], bins=nbins, normed=True)
        
        ## plotting the PDF
        plt_orig_pdf, = ax.step(bins[1::], log10(f), 'b', lw=1)
        
        ## fitting the distribution to a log-normal function
        rv, N = gas_orig.fit_interval(*fit_func_rng)        
        x = linspace(log(1e-3), log(1e+6), 100)
        plt_fit_pdf, = ax.plot(x, log10(rv.pdf(x)), 'r--', lw=1)

        ## setting the appropriate labels on the x-axis
        ticks = [1e-6, 1e-4, 1e-2, 1e-0, 1e2, 1e3, 1e4, 1e5, 1e6]
        ax.set_xticks( log(ticks) )
        ax.set_xticklabels([r'$10^{%d}$' % i for i in log10(ticks)], size='small')
        
        ax.set_ylim(log10(1e-10), log10(1))
        ticks = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0]
        ax.set_yticks( log10(ticks) )
        ax.set_yticklabels([r'$10^{%d}$' % i for i in log10(ticks)], size='small')
        
        ######### plotting the sampled PDF #############
        
        ## selecting the original gas particles
        gas_sampled_amuse = gas[where(isnan(gas.parent) == False)]
        gas_sampled = fi_utils.gas_set(len(gas_sampled_amuse)) 
        gas_sampled.copy_attr_from_amuse_set(gas_sampled_amuse, ['n'])
        
        ## using the ln of the densities                
        x = log(gas_sampled.n)

        ## parameters for the plots
        xmin, xmax = log(1e2), log(1e6)
        nbins = 30
                
        ## computing the PDF of the original particles
        f, bins = histogram(x, range=[xmin, xmax], bins=nbins, normed=True)
        
        ## plotting the PDF from which the new densities which are sampled
        f_sampled, = ax.step(bins[1::], log10(f/150), 'g--', lw=1)

        ## computing the PDF of the combined distribution of the original and sampled particles
        
        # using the ln of the densities 
        weights_matched, intervals, weights_func = gas.match_weights(20, 
                                                                     fit_func_rng, 
                                                                     1e2, 
                                                                     [1e2, gas.n.max()], 
                                                                     30, 
                                                                     0.001)
        
        gas.weights = weights_matched
        
        f_combined, bins = histogram(log(gas.n),
                                     range=[log(1e-6), log(1e6)], bins=nbins, normed=True,
                                     weights=gas.weights,
                                     )
        plt_combined, = ax.plot( (bins[0:-1] + bins[1:])/2.0, log10(f_combined), 'k+', markersize=4)


        #plt_combined, = ax.plot( (bins[0:-1] + bins[1:])/2.0, log10(f/150), 'go', markersize=3)                                     

        ax.legend([plt_orig_pdf, plt_fit_pdf, f_sampled, plt_combined], 
                  ['SPH', 'fit', 'sampled', 'combined'], loc=0, fontsize='x-small')
        ax.set_xlabel(r'$n / {\rm cm}^{-3}$', size='small')
        ax.set_ylabel('PDF', size='small')

        if save_figs == True:
            fig.fig.savefig(fig_path)
            print 'figure saved to:\n\t%s' % fig_path
        #fig.preview(env='aa')
        
    ###
    
    def plot_T_n_PDFs(gas, fig_path):

        fig = templates.subplots_grid(2,1, hspace=0.15, wspace=0.15,
                                      fig = {'kwargs':{
                                                       'figsize' : {3.55, 4.95}
                                                       }
                                             },
                                      axs = {
                                             'left' :0.17, 'bottom': 0.12, 'w':0.8, 'h':0.84
                                            },
                                      )

        ######### plotting the original PDFs of the temperature#############

        ## selecting the original gas particles
        gas_orig_amuse = gas[where(isnan(gas.parent) == True)]
        gas_orig = fi_utils.gas_set(len(gas_orig_amuse)) 
        gas_orig.copy_attr_from_amuse_set(gas_orig_amuse, ['n', 'T'])
        
        ## using the ln of the densities
        x1 = log10(gas_orig.T)

        ax = fig.sub[0,0]
        
        ax.hist(x1, 100)

        ticks = [1e1, 1e2, 1e3, 1e4, 1e5]
        ax.set_xticks( log10(ticks) )
        ax.set_xticklabels(['10', '100', '1000', r'$10^4$', r'$10^5$'], size='small')
                
        ticks = linspace(2.0e4, 1e5, 5)
        ax.set_yticks( ticks )
        ax.set_yticklabels(['2', '4', '6', '8', '10'], size='small')
        ax.set_ylim(0, 1.3e5)
        ax.set_xlim(1, 5)
        ax.set_xlabel(r'$T ({\rm K})$', size='small')
        ax.set_ylabel(r'$N$ ' + r'($\times 10^4$)', size='small')

        ax.text(2, 1.05e5, 'bound', size='x-small')
        ax.text(4, 0.95e5, 'unbound', size='x-small')
        
        ######### plotting the original PDFs of the temperature#############
        ax = fig.sub[1,0]
        
        x2 = log10(gas_orig.n)

        ax.hist(x2, 80)

        ax.hist(x2[x1 < log10(1000)], 50, label='T < 1000 K')

        ax.hist(x2[x1 < log10(500)], 50, label='T < 500 K')

        ax.hist(x2[x1 < log10(100)], 50, label='T < 100 K')

        ticks = [1e-6, 1e-4, 1e-2, 1e-0, 1e2, 1e4]
        ax.set_xticks( log10(ticks) )
        ax.set_xticklabels([r'$10^{%d}$' % i for i in log10(ticks)], size='small')
        ax.set_xlim(-6, 5)
        ax.set_xlabel(r'$n ( {\rm cm}^{-3} )$', size='small')

        ticks = linspace(2.0e4, 1e5, 5)
        ax.set_yticks( ticks )
        ax.set_yticklabels(['2', '4', '6', '8', '10'], size='small')

        ax.set_ylim(0, 1.3e5)
        ax.set_ylabel(r'$N$ ' + r'($\times 10^4$)', size='small')
        
        ax.legend(loc=0, fontsize='small')
        
        if save_figs == True:
            fig.fig.savefig(fig_path)
            print 'figure saved to:\n\t%s' % fig_path
        #fig.preview(env='aa')
    ###
    
    plot_T_n_PDFs(gas, fig_paths['fig1'])
    plot_original(gas, fig_paths['fig2'])
    
    show()


def total_luminosity(gas, params, fig_save_path):
    '''plots the total luminosity formatted figure'''
    
    fig = templates.subplots_grid(1,1, hspace=0.15, wspace=0.15,
                                  fig = {'kwargs':{
                                                   'figsize' : {3.55, 2.95}
                                                   }
                                         },
                                  axs = {
                                         'left' :0.18, 'bottom': 0.15, 'w':0.8, 'h':0.82
                                        },
                                  )
    
    ax = fig.sub[0,0]
    
    specsStrs =  ['CO', '13CO', 'HCN', 'HNC', 'HCO+', 'CS', 'SiO']
    colors    =  ['k' , 'r'   , 'g'  , 'b'  , 'c'   , 'y' , 'm']
    
    ## plotting the total luminosity weighted by number of sampled points
    weights_filename = params['rundir'] + '/firun/' + 'weights_func.%06d.npz' % params['snap_index']
    gas.use_weights(weighting='by-number', weights_filename = weights_filename)
    
    sym       = '-'
    
    for i, specStr in enumerate(specsStrs):
        
        print specStr
        
        x, y = gas.get_total_luminosity_ladder(specStr)
        
        ax.semilogy(x+1, y*1e6, colors[i] + sym, label=specStr)
    
    if True:
        ## plotting the total luminosity of the original points
        gas.use_weights(weighting='original-only')
            
        sym       = '--'
        
        
        for i, specStr in enumerate(specsStrs):
            
            print specStr
            
            x, y = gas.get_total_luminosity_ladder(specStr)
            
            ax.semilogy(x+1, y*1e6, colors[i] + sym)
        
    xlabel(r'$J_{\rm up}$', size='small')
    ylabel(r'L (K km s$^{-1}$ pc$^2$)', size='small')
    xlim([1,15])
    xticks(size='small')
    
    ylim([1e5, 1e10])
    yticks([1e5, 1e6, 1e7, 1e8, 1e9], size='small')
    legend(loc=0, prop={'size':'x-small'})
    
    show()
    
    if fig_save_path != None:
        fig.fig.savefig(fig_save_path)
        print 'saved figure to\n\t%s' % fig_save_path
        
    return fig 
    
def plot_flux_maps(gas, hist, params, fig_save_path):
    
    import matplotlib.cm as cm
    
    fig = templates.subplots_grid(7, 7, hspace=0.0, wspace=0.0,
                                  fig = {'kwargs':{
                                                   'figsize' : {7.2, 6.2}
                                                   }
                                         },
                                  axs = {
                                         'left' :0.05, 'bottom': 0.05, 'w':0.82, 'h':0.94
                                        },
                                  cbar = {'space':0.01, 'scale':0.7, 'sz':0.02,
                                          'orientation':'vertical', 'cmap':cm.OrRd},
                                  )
    
    ax = fig.sub[0,0]
    fig.set_xlabel('x (kpc)', space=-0.03, size='small')
    fig.set_ylabel('y (kpc)', space=-0.03, size='small')
    

    fig.sub_setp('xlim', [-8, 8])
    fig.sub_setp('ylim', [-8, 8])
 
    pos_labels = [-6, -3, 0, 3, 6]
    fig.set_xticks(pos_labels, labels=pos_labels, size='x-small')
    fig.set_yticks(pos_labels, labels=pos_labels, size='x-small')

    fig.clean_subplot_labels()
    fig.set_title('W (K km s$^{-1}$)', space=[+0.47, -0.4], size=10, rotation=90)

    fig.set_cbar_range([-6, 4])
    fig.set_cbar_ticks([-6, -4, -2, 0, 2, 4], 
                       [r'10$^{-6}$', r'10$^{-4}$', 0.01, 1, 100, r'10$^{4}$'])

    for r, spec_lines in enumerate(params['lines']):
        for c, line in enumerate(spec_lines):
            print 'em_fluxKkms_%s' % line
            map_data = fi_utils.make_map(gas, hist, attr='em_fluxKkms_%s' % line, as_log10=True, 
                                         func=fi_utils.mean_flux, show=True, in_ax=fig.sub[r,c],    
                                         extent=[-8, 8, -8, 8], cmap=cm.OrRd, vmin=-6, vmax=4, aspect='auto')
        print '------------------'
    
    for r, spec_lines in enumerate(params['lines']):
        for c, line in enumerate(spec_lines):
            latex_str = lineDict.lines[line]['latex']
            print latex_str
            fig.sub[r,c].text(-6, 5, latex_str, size=8)
        print '----------------'
        
    if fig_save_path != None:
        fig.fig.savefig(fig_save_path)
        print 'figure saved to:\n\t%s' % fig_save_path
    #fig.preview(env='aa')
    
    show()
    
