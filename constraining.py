import line_ratio_utils
import numpy
import pylab

class Xi2_line_ratios_single_component(object):
    
    def __init__(self, obs_data, model_data, model_parms, line_ratios):
         
        self.model_em = model_data
        self.model_parms = model_parms

        self.obs_ratios = obs_data
        
        self.line_ratios_use = line_ratios
        
        self.model_line_ratios = None
        
        self.nModels = self.model_em.values()[0].size 
        
        self.Xi2 = None
        self.Xi2_valid = None
        self.inds_valid = None
        self.model_parms_valid = None
        self.model_line_ratios_valid = {}
        
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

        
    def print_minima(self):
           
        ind_global_min = numpy.argmin(self.Xi2_valid)
        ind_model_orig = self.inds_orig[ind_global_min]
        print 'global minimum : Xi2 = ', numpy.min(self.Xi2_valid)
        print ' \tn, g0, gm, Av = ', self.model_parms_valid[ind_global_min,:]

        for line in self.obs_ratios.species_and_codes()[1]:
            print '\t%10s : %e' % (line, self.model_em[line][self.inds_orig[ind_global_min]])
        print '\tline ratio                model        observed'
        for line_ratio in self.line_ratios_use:
            print '\t%-20s : %e  %e' % (line_ratio, self.model_line_ratios[line_ratio][ind_model_orig], self.obs_ratios[line_ratio]['v'])

        self.ind_global_min = ind_global_min
        self.ind_global_min_orig = ind_model_orig

        #######################
                
        inds_zero_gm = numpy.where(self.model_parms_valid[:,2] == numpy.min(self.model_parms_valid[:,2]) )[0]
        inds_orig_zero_gm = self.inds_orig[numpy.arange(self.Xi2_valid.size)[inds_zero_gm]]
        grid_coords_valid_zero_gm = self.model_parms_valid[inds_zero_gm,:]
        Xi2_valid_zero_gm = self.Xi2_valid[inds_zero_gm]
        ind_min_zero_gm = numpy.argmin(Xi2_valid_zero_gm)
        print '--------------------------------------------------------------------------------------------------'
        
        print 'global minimum (no gmech): Xi2 = ', numpy.min(Xi2_valid_zero_gm)
        print '     n, g0, gm, Av = ', grid_coords_valid_zero_gm[ind_min_zero_gm,:]
        print '     n, g0, gm, Av = ', self.model_parms[inds_orig_zero_gm[ind_min_zero_gm],:]
        
        for line in self.obs_ratios.species_and_codes()[1]:
            print '\t%10s : %e' % (line, self.model_em[line][inds_orig_zero_gm[ind_min_zero_gm]])
        print '---------------------------------------------------------'
        print '\tline ratio                model        observed'
        ind_model_orig = inds_orig_zero_gm[ind_min_zero_gm]
        for line_ratio in self.line_ratios_use:
            print '\t%-20s : %e  %e' % (line_ratio, self.model_line_ratios[line_ratio][ind_model_orig], self.obs_ratios[line_ratio]['v'])
                        
        
        self.ind_global_min_no_gmech = ind_min_zero_gm
        self.ind_global_min_no_gmech_orig = ind_model_orig

    def plot_results(self, no_gmech=False):
        ## plotting the line ratios and the modelled ones

        fig = pylab.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        if no_gmech == False:
            ind_model = self.ind_global_min_orig
        else:
            ind_model = self.ind_global_min_no_gmech_orig
        
        obs = self.obs_ratios
        mod = self.model_line_ratios
        
        ## plotting the observed ratios for CO/CO line ratios
        line_ratios = obs.get_ratios_by_species('CO','CO', in_denom='1-0', sort_num=True)
        print 'plotting line ratios:', line_ratios
                
        vo, eo = obs.get_values_and_error(line_ratios)
        Jup = numpy.arange(len(vo))+1
        ax.errorbar(Jup, vo, fmt='ko', yerr=eo)
        ax.text(len(vo)-1, vo[-1], 'CO/CO', size='small')
        ## plotting the model line ratios
        vm = []
        for i, line_ratio in enumerate(line_ratios):
            vm.append(mod[line_ratio][ind_model])
        ax.plot(Jup, vm,'k--')
        #-----------------------------------------------------
        ## plotting the observed ratios for 13CO/13CO line ratios
        line_ratios = obs.get_ratios_by_species('13CO','13CO', in_denom='1-0', sort_num=True)        
        print 'plotting line ratios:', line_ratios

        vo, eo = obs.get_values_and_error(line_ratios)
        Jup = numpy.arange(len(vo))+1
        ax.errorbar(Jup, vo, fmt='ro', yerr=eo)
        ax.text(len(vo)-1, vo[-1], '13CO/13CO', size='small')
        ## plotting the model line ratios
        vm = []
        for i, line_ratio in enumerate(line_ratios):
            vm.append(mod[line_ratio][ind_model])
        ax.plot(Jup, vm,'r--')
        #-----------------------------------------------------
        ## plotting the observed ratios for 13CO/CO line ratios
        line_ratios = obs.get_ratios_by_species('13CO','CO', in_denom='1-0', sort_num=True)
        print 'plotting line ratios:', line_ratios
        
        vo, eo = obs.get_values_and_error(line_ratios)
        Jup = numpy.arange(len(vo))
        
        ax.errorbar(Jup, vo, fmt='go', yerr=eo)
        ax.text(len(vo)-1, vo[-1], '13CO/CO', size='small')
        ## plotting the model line ratios
        vm = []
        for i, line_ratio in enumerate(line_ratios):
            vm.append(mod[line_ratio][ind_model])
        ax.plot(Jup, vm,'g--')
        #-----------------------------------------------------
        
        ax.set_xlim(-1, 6)
        ax.set_ylim(0.001, 2.0)
        
        ax.set_xticklabels(['', 'i=0', 'i=1', 'i=2', 'i=3', 'i=4', 'i=5'])
        ax.set_ylabel('line ratio')
        ax.set_xlabel('J = i + 1/ J = 1-0')
        pylab.yscale('log')
        pylab.show()
