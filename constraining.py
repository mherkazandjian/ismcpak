import line_ratio_utils
import numpy
import pylab

class Xi2_line_ratios_single_component(object):
    
    def __init__(self, obs_data=None, model_data=None, model_parms=None, line_ratios=None):
         
        self.model_em = model_data  #: the model emission of the input PDR models (dict, line_ratio_string:values pairs)
        self.model_parms = model_parms  #: the model parameters of the input PDR models (n_models ndarray x 4 ndarray)
        self.obs_ratios = obs_data   #: the observed (measured) line ratios corresponding to self.line_ratios_use
        self.line_ratios_use = line_ratios   #: the line ratio strings (list of strings)
        
        self.model_line_ratios = None
        
        self.nModels = self.model_em.values()[0].size 
        
        self.Xi2 = None
        self.Xi2_valid = None
        self.inds_valid = None
        self.model_parms_valid = None      #: the model emission of the input PDR models which hold valid data self.model_parms
        self.model_line_ratios_valid = {}  #: the model line ratios of the input PDR models which hold valid data 
        
        self.ind_global_min = None
        
    def compute_model_line_ratios(self):
        ## computing the model line ratios
        model_ratios = {}
        for line_ratio in self.line_ratios_use:
        
            line1, line2 = line_ratio_utils.lines_involved(line_ratio)
        
            model_ratios[line_ratio] = self.model_em[line1]/self.model_em[line2]
            
        self.model_line_ratios = model_ratios
    
    def compute_Xi2(self):
        
        ## computing the Xi2 statistic
        Xi2 = numpy.zeros(self.nModels, 'f8')
        
        for line_ratio in self.obs_ratios:
        
            obs = self.obs_ratios[line_ratio]['v']
            err = self.obs_ratios[line_ratio]['e']
        
            mod = self.model_line_ratios[line_ratio]
            
            Xi2 += ((obs - mod)**2.0)/err**2 # goodness of fit
            #Xi2 += ((obs - mod)**2)/mod  # pearsons chi squared test
        
        self.Xi2 = Xi2
        
        ## keeping the finite values of the Xi2 statistic
        inds_valid = numpy.isfinite(Xi2)
        Xi2_valid = Xi2[inds_valid]
        self.Xi2_valid = Xi2_valid

        self.inds_orig = numpy.arange(Xi2.size)[inds_valid]
        
        self.model_parms_valid = self.model_parms[inds_valid, :]
        for line_ratio in self.obs_ratios:
            self.model_line_ratios_valid[line_ratio] = self.model_line_ratios[line_ratio][inds_valid]

    def model_inds_for_min_Xi2(self):
        
        ## the index of the model which minimum Xi2 for all n,G0,Av,gmech
        ind_global_min = numpy.argmin(self.Xi2_valid)
        ind_model_orig = self.inds_orig[ind_global_min]
        
        self.ind_global_min = ind_global_min
        self.ind_global_min_orig = ind_model_orig

        ## the index of the model with minimum Xi2 for all n,G0,Av,gmech=0
        inds_zero_gm = numpy.where(self.model_parms_valid[:,2] == numpy.min(self.model_parms_valid[:,2]) )[0]
        inds_orig_zero_gm = self.inds_orig[numpy.arange(self.Xi2_valid.size)[inds_zero_gm]]

        Xi2_valid_zero_gm = self.Xi2_valid[inds_zero_gm]
        ind_min_zero_gm = numpy.argmin(Xi2_valid_zero_gm)
        
        ind_model_orig = inds_orig_zero_gm[ind_min_zero_gm]
        
        self.ind_global_min_no_gmech = ind_min_zero_gm
        self.ind_global_min_no_gmech_orig = ind_model_orig

        return self.ind_global_min, self.ind_global_min_orig,\
               self.ind_global_min_no_gmech, self.ind_global_min_no_gmech_orig

    def get_model_parms_for_min_Xi2(self):

        self.model_inds_for_min_Xi2()
        
        ## parms of the model with the gloab Xi2 min
        model_parms_for_global_Xi2_min = self.model_parms_valid[self.ind_global_min,:]
        
        ## getting the params of model with the min Xi2 for gmech = 0 
        inds_zero_gm = numpy.where(self.model_parms_valid[:,2] == numpy.min(self.model_parms_valid[:,2]) )[0]
        grid_coords_valid_zero_gm = self.model_parms_valid[inds_zero_gm,:]
        Xi2_valid_zero_gm = self.Xi2_valid[inds_zero_gm]
        ind_min_zero_gm = numpy.argmin(Xi2_valid_zero_gm)
        model_parms_for_zero_gmech_Xi2_min = grid_coords_valid_zero_gm[ind_min_zero_gm,:]
        
        
        return model_parms_for_global_Xi2_min, numpy.min(self.Xi2_valid),\
               model_parms_for_zero_gmech_Xi2_min, numpy.min(Xi2_valid_zero_gm) 
    
    def get_model_em_for_min_Xi2(self):

        print 'emission intensities for the best fit model'
        print '\t     line     intensity model     intensity observed '
        for line in self.obs_ratios.species_and_codes()[1]:
            print '\t%10s :   %e         %e' % (line, 
                                                self.model_em[line][self.inds_orig[self.ind_global_min]],
                                                1
                                               )
            
        print '\tline ratio                model        observed'
        for line_ratio in self.line_ratios_use:
            print '\t%-20s : %e  %e' % (line_ratio, self.model_line_ratios[line_ratio][self.ind_global_min_orig],\
                                        self.obs_ratios[line_ratio]['v'])
        
       
    def print_minima(self):

        self.get_model_parms_for_min_Xi2()
        
        ## minimum infor for 4D parameter space
        model_parms_for_global_Xi2_min, Xi2_min_global, model_parms_for_zero_gmech_Xi2_min, Xi2_min_zero_gm = self.get_model_parms_for_min_Xi2()
        
        print 'global minimum : Xi2 = ', Xi2_min_global
        print ' \tn, g0, gm, Av = ', model_parms_for_global_Xi2_min

        """
        for line in self.obs_ratios.species_and_codes()[1]:
            print '\t%10s : %e' % (line, self.model_em[line][self.inds_orig[self.ind_global_min]])
            
        print '\tline ratio                model        observed'
        for line_ratio in self.line_ratios_use:
            print '\t%-20s : %e  %e' % (line_ratio, self.model_line_ratios[line_ratio][self.ind_global_min_orig],\
                                        self.obs_ratios[line_ratio]['v'])
        """
        
        ## minimum infor for 4D parameter space (no gmech)                
        inds_models_zero_gm = numpy.where(self.model_parms_valid[:,2] == numpy.min(self.model_parms_valid[:,2]) )[0]
        inds_models_orig_zero_gm = self.inds_orig[numpy.arange(self.Xi2_valid.size)[inds_models_zero_gm]]
        #print '--------------------------------------------------------------------------------------------------'
         
        print 'global minimum (no gmech): Xi2 = ', Xi2_min_zero_gm
        print '     n, g0, gm, Av = ', model_parms_for_zero_gmech_Xi2_min 
        
        """
        for line in self.obs_ratios.species_and_codes()[1]:
            print '\t%10s : %e' % (line, self.model_em[line][inds_models_orig_zero_gm[self.ind_global_min_no_gmech]])
        print '---------------------------------------------------------'
        print '\tline ratio                model        observed'
        for line_ratio in self.line_ratios_use:
            print '\t%-20s : %e  %e' % (line_ratio, self.model_line_ratios[line_ratio][self.ind_global_min_no_gmech_orig],\
                                        self.obs_ratios[line_ratio]['v'])                        
        """

    def plot_results(self, fig=None, ax=None, no_gmech=False):
        ## plotting the line ratios and the modelled ones

        if fig == None:
            fig = pylab.figure(figsize=(4,8))
        if ax == None:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        if no_gmech == False:
            ind_model = self.ind_global_min_orig
        else:
            ind_model = self.ind_global_min_no_gmech_orig
        
        obs = self.obs_ratios
        mod = self.model_line_ratios
        
        def plot_ratio_ladder(spec_num=None, spec_denom=None, in_denom='1-0', **kwargs):
            ## plotting the observed ratios for CO/CO line ratios
            line_ratios = obs.get_ratios_by_species(spec_num, spec_denom, **kwargs)
            print 'plotting line ratios %s/%s    :' % (spec_num, spec_denom), line_ratios
                    
            vo, eo = obs.get_values_and_error(line_ratios)
            Jup = numpy.arange(len(vo))+1
    
            ax.errorbar(Jup, vo,  yerr=eo, fmt='o', markersize=3, **kwargs)
            #ax.text(len(vo)-1, vo[-1], '%s/%s' %  (spec_num, spec_denom), size='x-small')
            
            ## plotting the model line ratios
            vm = []
            for i, line_ratio in enumerate(line_ratios):
                vm.append(mod[line_ratio][ind_model])
            ax.plot(Jup, vm, linestyle='--', label='%s/%s' %  (spec_num, spec_denom), **kwargs)
        #
        
        plot_ratio_ladder(spec_num='CO', spec_denom='CO', in_denom='1-0', color='k')
        plot_ratio_ladder(spec_num='13CO', spec_denom='13CO', in_denom='1-0', color='r')
        plot_ratio_ladder(spec_num='13CO', spec_denom='CO', in_denom='1-0', color='g')
        plot_ratio_ladder(spec_num='HCN', spec_denom='CO', in_denom='1-0', color='c')
        plot_ratio_ladder(spec_num='HNC', spec_denom='CO', in_denom='1-0', color='m')
        plot_ratio_ladder(spec_num='HCO+', spec_denom='CO', in_denom='1-0', color='y')
        

        ax.set_xlim(-1, 17)
        ax.set_ylim(1e-7, 2.0)
        
        ax.set_xticks(numpy.arange(0, 16, 2))
        ax.set_xticklabels(numpy.arange(0, 16, 2))
        ax.set_ylabel('line ratio')
        ax.set_xlabel('J = i + 1/ J = 1-0')
        ax.set_yscale('log')
        
        ax.legend(loc=0, prop={'size':'x-small'})
        #pylab.show()
        