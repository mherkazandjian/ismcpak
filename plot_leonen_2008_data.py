import numpy
import pylab
import line_ratio_utils

data_fname = '/home/mher/ism/observatios/loenen2008/loenen2008-table-3.2.npz' 
data = numpy.load(data_fname)['arr_0']
fig, axs = pylab.subplots(2,2, sharex=False, sharey=False)

pylab.setp(axs[0,0], 'xlim', [-0.8, 0.8], 'ylim', [-0.8, 0.8])
pylab.setp(axs[0,1], 'xlim', [-1.0, 0.4], 'ylim', [-0.8, 0.8])
pylab.setp(axs[1,0], 'xlim', [-0.8, 0.8], 'ylim', [-1.0, 0.4])
pylab.setp(axs[1,1], 'xlim', [-2.0, 1.0], 'ylim', [-2.0, 2.0])


def get_line_ratio_from_object(obj, r_str):
    '''return a line ratio for a certain object '''
    
    def get_value_and_error(obj, strng):
        r = 0.5*(obj[strng]['mx'] + obj[strng]['mx'])
        e = 0.5*(obj[strng]['mx'] - obj[strng]['mx'])
        return r, e
    #
    
    if r_str in obj.dtype.names:
        r, e = get_value_and_error(obj, r_str)
    elif line_ratio_utils.reverse_ratio_str(r_str) in obj.dtype.names:
        r, e = get_value_and_error(obj, line_ratio_utils.reverse_ratio_str(r_str))
        
        #taking the inverse of the data
        r = 1.0/r
        #computing the error (not so accurate but for now it is good)
        percent_err = numpy.fabs(1.0 - e/r)
        e = percent_err * r 
        
    else:
        raise 'line ratio %s not available in the object'

    return r, e
#
def get_useful_data(r1_str, r2_str, data):
    ''' '''
    
    x, xerr = [], []
    y, yerr = [], []
    
    for obj in data:
        
        x_obj, xerr_obj = get_line_ratio_from_object(obj, r1_str)
        y_obj, yerr_obj = get_line_ratio_from_object(obj, r2_str)
        
        if numpy.isnan(x_obj) == False and numpy.isnan(y_obj) == False:  
            x.append(x_obj)
            y.append(y_obj)
            xerr.append(xerr_obj)
            yerr.append(yerr_obj)
    
    return numpy.array(x, 'f8'), numpy.array(xerr, 'f8'), numpy.array(y, 'f8'), numpy.array(yerr, 'f8') 
#

panels = [
          ['HNC1-0/HCO+1-0' , 'HCO+1-0/HCN1-0'],
          ['HNC1-0/HCN1-0'  , 'HCO+1-0/HCN1-0'],
          ['HNC1-0/HCO+1-0' , 'HNC1-0/HCN1-0' ],
          #['HNC1-0/CO1-0'   , 'HNC1-0/HCN1-0']          
          ['HCO+1-0/CO1-0'   , 'HNC1-0/HCN1-0']          
         ]

axs = axs.flatten()

for i, panel in enumerate(panels):
        
    rstr1, rstr2 = panel 
    
    x, xerr, y, yerr = get_useful_data(rstr1, rstr2, data)
    
    axs[i].plot(numpy.log10(x), numpy.log10(y), 'ks')
    axs[i].set_xlabel(rstr1)
    axs[i].set_ylabel(rstr2)

axs[0].plot([-1,1], [1, -1])
axs[1].plot([0,0] , [-1, 1])
axs[2].plot([-1,1], [0.0, 0.0])

pylab.figtext(0.1, 0.1, r'HCO$^+$(4-3)/CO(1-0)')
pylab.figtext(0.3, 0.1, r'HCN(1-0)/HNC(1-0)')
pylab.figtext(0.5, 0.1, r'HCN(1-0)/HCO$^+$(1-0)')
pylab.figtext(0.5, 0.4, r'HNC(1-0)/HCO$^+$(1-0)')

pylab.show()
