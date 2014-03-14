# -*- coding: utf-8 -*-

import numpy
import pylab
import itertools

from xlrd import open_workbook
import collections

import lineDict

def codes_from_ratio(ratio_str):
    '''gets the line codes from a ratio string '''
    code1, code2 = ratio_str.split('/')
            
    return code1.strip(), code2.strip()

def reverse_ratio_str(ratio_str):
    '''returns the inverted ratio of the input ratio string x/y becomes y/x'''
    code1, code2 = codes_from_ratio(ratio_str)
    return '%s/%s' % (code2, code1)

def line_ratio_latex(ratio_str):
    
    code1, code2 = codes_from_ratio(ratio_str) 

    latex_ratio_str = r'%s/%s' % (lineDict.lines[code1]['latex'], lineDict.lines[code2]['latex']) 
    
    return latex_ratio_str  

class ratios(collections.OrderedDict):
    '''A class (based on a dict object) which provides utilities and storage for 
    line ratio observatios
    
        r = ratios()
        r['CO1-0'] = {'fluxcgs':1e-16, 'err':1e-17}
    '''
    
    def __init__(self, **kwargs):
        
        super(ratios, self).__init__()
        
        self.codes = None #: all the unique codes of the lines in the dict 
        self.lines = {'luminosity' : {}, # the luminosity of the lines
                      'intensity'  : {}, # the intensity  of the lines
                      } #: a dict holding the emission of the lines corresponding to self.codes
        
        self.specStrs = None #: all the unique species string in the dict

    def add_ratio(self, rStr, v=None, e=None, mn=None, mx=None):
        '''add a line ratio item given the range of the ratio (v - err/2, v+err/2)'''

        r = {}
        
        if v != None and e != None:
            if mn != None or mx != None:
                raise ValueError('mn or mx keywords can not be passed when (v,e) are passed')
            else:
                r['v'] = v
                r['e'] = e

        if mn != None and mx != None:
            if v != None or e != None:
                raise ValueError('e or v keywords can not be passed when (mn,mx) are passed')
            else:
                r['v'] = (mx + mn)/2.0
                r['e'] = (mx - mn)/2.0
        
        self[rStr] = r 
    
    def make_ratios(self, lines, ratios, em_unit=None, lum=None):
        '''takes a dict of lines ('CO1-0', 'HCN4-3'...ect) and ratio string ('HCN4-3/CO1-0',...)
        and adds the ratios as dict keys and computes the errors using propagation of error
        for ratios (assuming the errors in line[SPEC_STR]['err'] exist)
        '''
        
        if em_unit == None:
            unit = 'fluxcgs'
        else:
            unit = 'fluxKkms'
            
        for ratio in ratios:
            
            code1, code2 = codes_from_ratio(ratio)
            
            v1, e1 = lines[code1][unit], lines[code1]['err']
            v2, e2 = lines[code2][unit], lines[code2]['err'] 
             
            #the value of the line ratio
            v = v1 / v2
            
            #the error on the line ratio
            e = v * numpy.sqrt( (e1/v1)**2.0 + (e2/v2)**2.0 )
            
            self[ratio] = {}
            
            self[ratio]['v'] = v  
            self[ratio]['e'] = e 
        
            self.lines['intensity'][code1]  = lines[code1]
            self.lines['intensity'][code2]  = lines[code2]
            if lum != None:
                self.lines['luminosity'][code1] = lum[code1]
                self.lines['luminosity'][code2] = lum[code2]
        
    def species_and_codes(self):
        '''returns the species and the line codes from a list of line ration passed as strings
        
        For example :
        
        .. code-block:: python
        
            specStrs, codes = get_info_from_line_ratio(['CO4-3/CO1-0', 'HCO+1-0/HCN1-0'])
        
        the values of specStrs and codes is :
        
        .. code-block:: python
        
            specStr = ['CO', 'HCO+', 'HCN']
            codes   = ['CO4-3', 'CO1-0', 'HCO+1-0', 'HCN1-0']
        '''
        
        specStrs = {}
        codes = {}
        
        for ratio in self:
            
            #getting the line codes invloved in the ratio
            code1, code2 = codes_from_ratio(ratio)            
            codes[code1] = True
            codes[code2] = True
            
            #getting the species string for each component of the ratio
            specStr1 = lineDict.lines[code1]['specStr']
            specStr2 = lineDict.lines[code2]['specStr']            
            specStrs[specStr1] = True        
            specStrs[specStr2] = True
        
        self.specStrs = specStrs.keys()
        self.codes = codes.keys()
        
        return specStrs.keys(), codes.keys()
    
    def get_all_values(self):
        '''returns the values of all the ratios'''
        return [self[ratio_str]['v'] for ratio_str in self.keys()]
    
    def get_all_errors(self):
        '''returns the values of all the ratios'''
        return [self[ratio_str]['e'] for ratio_str in self.keys()]

    def get_all_values_and_error(self, line_ratios):
        '''returns all the line ratio values and errors for the specified ratios in the list'''
        
        return self.get_all_values(line_ratios), self.get_all_errors(line_ratios)
    
    def get_values(self, line_ratios):
        '''returns a list of line ratio values for the specified ratios in the list'''
        
        v_ratios = []
        for line_ratio in line_ratios:
            v_ratios.append(self[line_ratio]['v'])
            
        return numpy.array(v_ratios, 'f8')

    def get_errors(self, line_ratios):
        '''returns a list of line ratio errors for the specified ratios in the list'''
        
        e_ratios = []
        for line_ratio in line_ratios:
            e_ratios.append(self[line_ratio]['e'])
            
        return numpy.array(e_ratios, 'f8')
    
    def get_values_by_species(self, spec1, spec2):
        '''returns the line ratios of all the lines where the ratios involve spec1 and spec2 in the num and denom respectively'''
        
        line_ratios = self.get_ratios_by_species(spec1, spec2)
        
        v = self.get_values(line_ratios)
        
        return v

    def get_errors_by_species(self, spec1, spec2):
        '''returns the line ratios of all the lines where the ratios involve spec1 and spec2 in the num and denom respectively'''
        
        line_ratios = self.get_ratios_by_species(spec1, spec2)
        
        v = self.get_errors(line_ratios)
        
        return v

    def get_values_and_errors_by_species(self, spec1, spec2):
        '''returns the line ratios values and errors of all the lines where the ratios involve spec1 and spec2 in the num and denom respectively'''
        
        return self.get_values_by_species(spec1, spec2), self.get_errors_by_species(spec1, spec2) 
        
    def get_values_and_error(self, line_ratios):
        '''returns the line ratio values and errors for the specified ratios in the list'''
        
        return self.get_values(line_ratios), self.get_errors(line_ratios)

    def get_ratios_by_species(self, spec1, spec2, in_num=None, in_denom=None, sort_num=None, **kwargs):
        '''returns the line ratios string of the line ratios whose numerator involves spec1 and denominator spec2
        For example 
        
            get_values_by_species('13CO', 'CO')
        
        returns a list of line ratio strings such as 13CO2-1/CO1-0, 13CO2-1/CO2-1 ... 13COX_Y/COM-N    
        '''
        
        ret = []
        
        for line_ratio in self:
            
            line1, line2 = lines_involved(line_ratio)
            
            if in_num != None and in_num not in line1:
                continue 

            if in_denom != None and in_denom not in line2:
                continue 
             
            if lineDict.lines[line1]['specStr'] == spec1 and lineDict.lines[line2]['specStr'] == spec2:
                ret.append(line_ratio)
        
        ret = numpy.array(ret)
        
        if sort_num != None and sort_num == True:
            '''sort the line ratio strings with increasing transition index of the line in the numerator'''
            idx_num = []
            for line_ratio in ret:
                
                line1, line2 = lines_involved(line_ratio)
                
                idx1 = lineDict.lines[line1]
                
                idx_num.append(idx1['radexIdx'])

            inds = numpy.argsort(numpy.array(idx_num))
            ret = ret[inds]
            
        return list(ret)
    
    #def get_
    def show(self):
        
        print 'ratio                |   value   |    err   '
        print '---------------------+-----------+----------'

        for item in self:
            print '%-20s |%10.2e |%10.2e' % (item, self[item]['v'], self[item]['e'])
        
def read_observations(fname):
    '''returns info about the observed lines as a dict object from an excel sheet (Marissa's format)'''

    lines = {}
    
    wb = open_workbook(fname)
    
    for s in wb.sheets():
    
        if s.name == 'identified':
    
            print 'Sheet name : *%s*' % s.name
            
            for row in range(s.nrows):
                
                if row < 15-1:
                    continue
                
                specStr    = s.cell(row, 0).value.strip()
                transition = s.cell(row, 1).value
                fluxcgs    = s.cell(row, 41).value
                rel_err    = 0.3 #setting the relative error manually (for now)
                
                #an empty line
                if len(specStr) == 0:
                    continue
                
                #getting the 
                specStr = lineDict.get_default_specStr(specStr)
                
                #cleaning the transition string
                transition = transition.replace(' ','')
                transition = transition.replace(u'–', u'-')
                
                if specStr == None:
                    continue
                    
                if type(fluxcgs) == int:
                    continue
                
                line = specStr + transition
                line = line.encode('ascii')
                print line, fluxcgs
     
                lines[line] = {'fluxcgs': numpy.float64(fluxcgs), 'err': numpy.float64(fluxcgs)*rel_err}

    return lines
#

def Xi2_line_ratios(obs_ratios, arxvPDR):
    '''Computes the Xi2 statistic given the observed lines and a PDR arxv.'''
    
    
    allData = numpy.recarray([],[('x', 'f8'),('y', 'f8'),('z', 'f8'),('t', 'f8'),('v', 'f8'),])
    
    models = {} 
    
    specStrs, codes = obs_ratios.species_and_codes()

    #collecting all the line intensities of the ratios involved in the observations (obs_ratio)
    #from the model database. Proccessing one Av at a time...
    for i, AvStr in enumerate(arxvPDR.radexDbs):
                
        Av = numpy.float64(AvStr)

        #array which will hold the grid points and the values for this Av
        data = numpy.recarray((arxvPDR.nMeshes), allData.dtype.descr)

        #getting the emissions for each line from the PDR database for all the models for the current Av
        for code in codes:
            models[code] = 10.0**arxvPDR.get_emissions_from_databases(line={'type':'radex-lvg', 'code':code}, Av_use=Av)

        #defining the array which will hold the Xi2 for all the models for this Av
        Xi2 = numpy.zeros(arxvPDR.nMeshes, 'f8')
        
        #compute the Xi2
        for obs_ratio in obs_ratios:
            
            #the line codes invloved in this ratio
            code1, code2 = codes_from_ratio(obs_ratio)
            
            #the ratios for all the models at this Av for this particular line ratio
            model_ratio = models[code1] / models[code2] 
            
            #computing the Xi2
            f_o = obs_ratios[obs_ratio]['v']
            f_e = obs_ratios[obs_ratio]['e']
            f_m = model_ratio
            
            Xi2 += ((f_m - f_o)/f_e)**2.0
        #
        
        data.x = arxvPDR.grid_x
        data.y = arxvPDR.grid_y
        data.z = arxvPDR.grid_z
        data.t = Av
        data.v = Xi2

        allData = numpy.hstack((allData, data) )

    #removing the first entry (redundant ;;; .. todo:: fix this [low priority])
    allData = allData[1::]
    
    #filtering out the points which have Nans 
    inds_not_nan = numpy.where( numpy.isfinite(allData['v']) )

    return allData[inds_not_nan]
        
    return allData[1::]

def plot_single_model_ratios(arxvPDR, x, y, z, t, obs_ratios):

    specStrs, codes = obs_ratios.species_and_codes()

    model_data = {}
    
    #getting the emissions for each line from the PDR model
    for code in codes:
        
        specStr = lineDict.lines[code]['specStr']

        arxvPDR.use_radexDb(specStr=specStr, Av=t)
        ind_model = arxvPDR.get_mesh_index(x=x, y=y, z=z)

        radex_mesh = arxvPDR.meshesRadex[ind_model]
        
        if radex_mesh != None:  
            line_flux = radex_mesh['fluxcgs'][lineDict.lines[code]['radexIdx']]
        else:
            line_flux = 0.0
        
        model_data[code] = line_flux

    #plotting the model ratios and the observed ones    
    i_all = [] 
    fo_all = []
    fe_all = []
    fm_all = []
    lbl_all = []
    
    #getting the ratios
    for i, obs_ratio in enumerate(obs_ratios):

        #the line codes invloved in this ratio
        code1, code2 = codes_from_ratio(obs_ratio)
        
        #the ratios for the model
        try:
            model_ratio = model_data[code1] / model_data[code2]
            print 'ratio %s/%s ok' % (code1, code2)

            i_all.append(i)
            fo_all.append(obs_ratios[obs_ratio]['v'])
            fe_all.append(obs_ratios[obs_ratio]['e'])
            fm_all.append(model_ratio)
            lbl_all.append(obs_ratio)
    
        except ZeroDivisionError:
            print 'ratio %s/%s cannot be computed.' % (code1, code2)
            print '    %s = %e' % (code1, model_data[code1])
            print '    %s = %e' % (code2, model_data[code2])
             
    i_all = numpy.array(i_all)
    fo_all = numpy.array(fo_all)
    fe_all = numpy.array(fe_all)
    fm_all = numpy.array(fm_all)
    
    #computing the Xi2 
    Xi2 = ((fo_all - fm_all)/fe_all)**2.0
    
    pylab.figure()
    pylab.xticks(i_all)
    pylab.errorbar(i_all, fo_all, yerr=[fe_all, fe_all], xerr = 0, fmt='s')
    pylab.plot(i_all, fm_all, 'ro')
    
    #pylab.gca().set_xticklabels(lbl_all)
    for i, label in enumerate(lbl_all):
        pylab.text(i, fo_all[i], label) 
                   
    pylab.xlim([i_all.min()-1, i_all.max()+1])
    pylab.yscale('log')

    pylab.title('Xi2 = %f' % Xi2.sum() )
    
def same_species_ratio(line_ratio):
    '''returns true if the lines involved are for the same species. 
        13CO2-1/13CO1-0 ruturns True,
        13CO2-1/CO1-0 ruturns False
    '''
    line1, line2 = lines_involved(line_ratio)
    
    if lineDict.lines[line1]['specStr'] == lineDict.lines[line2]['specStr']:
        return True
    else:
        return False
    

def lines_involved(line_ratio):
    '''returns the lines involved in the line ratio. 
        13CO2-1/13CO1-0 ruturns 13CO2-1, 13CO1-0
    '''
    line1, line2 = line_ratio.split('/')
    
    return line1, line2

def all_lines_involved(line_ratios):
    
    lines = []

    for line_ratio in line_ratios:
        
        line1, line2 = lines_involved(line_ratio)
        
        lines.append(line1)
        lines.append(line2)

    lines = numpy.unique(lines)
    
    return lines

def species_involved(lines):
    
    specs = []
    
    for line in lines:
        
        specs.append( lineDict.lines[line]['specStr'] )

    return numpy.unique(specs)

def line_ratio_combinations(lines, combinations=None):

    x = list(itertools.permutations(lines, 2))
    
    line_ratios = {}
    
    for line1, line2 in x:
        
        spec1, spec2 = lineDict.lines[line1]['specStr'], lineDict.lines[line2]['specStr']
        
        idx1, idx2 = lineDict.lines[line1]['radexIdx'], lineDict.lines[line2]['radexIdx']
        
        if spec1 == spec2 and idx1 < idx2:
            continue
        
        if combinations != None:

            dont_consider = True
                        
            for combination in combinations:
                
                spec_num, spec_denom = combination.split('/')
                
                if spec1 == spec_num and spec2 == spec_denom:
                    dont_consider = False

            if dont_consider == True:
                continue
            
        line_ratio = '%s/%s' % (line1, line2)
        line_ratio_inv = '%s/%s' % (line2, line1)
        
        if line_ratio_inv not in line_ratios:
            line_ratios[line_ratio] = True
        
    return line_ratios.keys()

#########################################

ratio_sets = {
              ### CO, 13CO ratios in the denominator (all ladders)
              'COCO' : (
                        'CO2-1/CO1-0', 'CO3-2/CO1-0', 'CO4-3/CO1-0', 'CO5-4/CO1-0', 'CO6-5/CO1-0', 'CO7-6/CO1-0',   
                        'CO8-7/CO1-0', 'CO9-8/CO1-0', 'CO10-9/CO1-0', 'CO11-10/CO1-0', 'CO12-11/CO1-0', 'CO13-12/CO1-0', 
                        'CO14-13/CO1-0', 'CO15-14/CO1-0', 
                        ), 
              '13CO13CO' : (
                            '13CO2-1/13CO1-0', '13CO3-2/13CO1-0', '13CO4-3/13CO1-0', '13CO5-4/13CO1-0', '13CO6-5/13CO1-0', '13CO7-6/13CO1-0', 
                            '13CO8-7/13CO1-0', '13CO9-8/13CO1-0', '13CO10-9/13CO1-0', '13CO11-10/13CO1-0', '13CO12-11/13CO1-0', '13CO13-12/13CO1-0', 
                            '13CO14-13/13CO1-0', '13CO15-14/13CO1-0',
                           ),
              '13COCO' :  (
                           '13CO2-1/CO1-0', '13CO3-2/CO1-0', '13CO4-3/CO1-0', '13CO5-4/CO1-0', '13CO6-5/CO1-0', '13CO7-6/CO1-0', 
                           '13CO8-7/CO1-0', '13CO9-8/CO1-0', '13CO10-9/CO1-0', '13CO11-10/CO1-0', '13CO12-11/CO1-0', '13CO13-12/CO1-0', 
                           '13CO14-13/CO1-0', '13CO15-14/CO1-0',
                          ), 
              'HCNCO' : (
                         'HCN1-0/CO1-0', 'HCN2-1/CO1-0', 'HCN3-2/CO1-0', 'HCN4-3/CO1-0', 'HCN5-4/CO1-0',  'HCN6-5/CO1-0', 'HCN7-6/CO1-0',
                        ),

              'HNCCO' : (
                         'HNC1-0/CO1-0', 'HNC2-1/CO1-0', 'HNC3-2/CO1-0', 'HNC4-3/CO1-0', 'HNC5-4/CO1-0',  'HNC6-5/CO1-0', 'HNC7-6/CO1-0',
                        ),
              'HCO+CO' : (
                          'HCO+1-0/CO1-0', 'HCO+2-1/CO1-0', 'HCO+3-2/CO1-0', 'HCO+4-3/CO1-0', 'HCO+5-4/CO1-0',  'HCO+6-5/CO1-0', 'HCO+7-6/CO1-0',
                         ),
              'SiOCO' : (
                         'SiO1-0/CO1-0', 'SiO2-1/CO1-0', 'SiO3-2/CO1-0', 'SiO4-3/CO1-0', 'SiO5-4/CO1-0',  'SiO6-5/CO1-0', 'SiO7-6/CO1-0',
                        ),
              'CSCO' : (
                        'CS1-0/CO1-0', 'CS2-1/CO1-0', 'CS3-2/CO1-0', 'CS4-3/CO1-0', 'CS5-4/CO1-0',  'CS6-5/CO1-0', 'CS7-6/CO1-0',
                       ), 
              'CNCO' : (
                        'CN1-0/CO1-0', 'CN2-1/CO1-0', 'CN3-2/CO1-0', 'CN4-3/CO1-0', 'CN5-4/CO1-0',  'CN6-5/CO1-0', 'CN7-6/CO1-0',
                       ),
              'HCO+13CO' : (
                         'HCO+1-0/13CO1-0', 'HCO+2-1/13CO1-0', 'HCO+3-2/13CO1-0', 'HCO+4-3/13CO1-0', 'HCO+5-4/13CO1-0',  'HCO+6-5/13CO1-0', 'HCO+7-6/13CO1-0',
                         ),
              ### CO ratios in the denominator (all ladders)
              'COCO-low'     : ('CO2-1/CO1-0',), 
              '13CO13CO-low' : ('13CO2-1/13CO1-0',),
              '13COCO-low'   : ('13CO1-0/CO1-0',), 
              'HCNCO-low'    : ('HCN1-0/CO1-0',),
              'HNCCO-low'    : ('HNC1-0/CO1-0',),
              'HCO+CO-low'   : ('HCO+1-0/CO1-0',),
              'SiOCO-low'    : ('SiO1-0/CO1-0',),
              'CSCO-low'     : ('CS1-0/CO1-0',), 
              'CNCO-low'     : ('CN1-0/CO1-0',),
              ##########################################
              # ratios among high density tracers
              ##########################################
              'HCNHNC' : (
                         'HCN1-0/HNC1-0', 'HCN2-1/HNC1-0', 'HCN3-2/HNC1-0', 'HCN4-3/HNC1-0', 'HCN5-4/HNC1-0',  'HCN6-5/HNC1-0', 'HCN7-6/HNC1-0',
                         ),
              'HCNHCO+' : (
                         'HCN1-0/HCO+1-0', 'HCN2-1/HCO+1-0', 'HCN3-2/HCO+1-0', 'HCN4-3/HCO+1-0', 'HCN5-4/HCO+1-0',  'HCN6-5/HCO+1-0', 'HCN7-6/HCO+1-0',
                         ),
              'HNCHCO+' : (
                         'HNC1-0/HCO+1-0', 'HNC2-1/HCO+1-0', 'HNC3-2/HCO+1-0', 'HNC4-3/HCO+1-0', 'HNC5-4/HCO+1-0',  'HNC6-5/HCO+1-0', 'HNC7-6/HCO+1-0',
                         ),
              # J > 3-2 in the numerator
              'HNCHCNhigh' : (
                              'HNC4-3/HCN1-0', 'HNC5-4/HCN1-0',  'HNC6-5/HCN1-0', 'HNC7-6/HCN1-0',
                             ),
              'HCNHCO+high' : (
                              'HCN4-3/HCO+1-0', 'HCN5-4/HCO+1-0',  'HCN6-5/HCO+1-0', 'HCN7-6/HCO+1-0',
                              ),
              'HNCHCO+high' : (
                              'HNC4-3/HCO+1-0', 'HNC5-4/HCO+1-0',  'HNC6-5/HCO+1-0', 'HNC7-6/HCO+1-0',
                              ),
              }

ratios_set1 = ratio_sets['COCO'] + ratio_sets['13CO13CO'] + ratio_sets['13COCO'] +\
              ratio_sets['HCNCO'] + ratio_sets['HNCCO'] + ratio_sets['HCO+CO'] +\
              ratio_sets['SiOCO'] + ratio_sets['CSCO'] + ratio_sets['CNCO']
              
ratios_set2 = ratio_sets['COCO'] + ratio_sets['13CO13CO'] + ratio_sets['13COCO']

ratios_set3 = ratio_sets['COCO'] + ratio_sets['13CO13CO'] + ratio_sets['13COCO'] +\
              ratio_sets['HCNCO'] + ratio_sets['HNCCO'] + ratio_sets['HCO+CO'] 

ratios_set4 = ratio_sets['HCNCO'] + ratio_sets['HNCCO'] + ratio_sets['HCO+CO'] 

ratios_set5 = ratio_sets['HCNHNC'] + ratio_sets['HCNHCO+'] + ratio_sets['HNCHCO+'] 

ratios_set6 = ratio_sets['HNCHCNhigh'] + ratio_sets['HCNHCO+high'] + ratio_sets['HNCHCO+high'] 
             
low_ratios = ratio_sets['COCO-low'] + ratio_sets['13CO13CO-low'] + ratio_sets['13COCO-low'] +\
             ratio_sets['HCNCO-low'] + ratio_sets['HNCCO-low'] + ratio_sets['HCO+CO-low'] +\
             ratio_sets['SiOCO-low'] + ratio_sets['CSCO-low'] + ratio_sets['CNCO-low']
