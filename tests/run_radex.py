"""
driver script to run radex through python and read back the output
A single LVG run
"""
from radex import Radex
from time import time


# path of the radex excutable
# radex_path = '~/ism/code/radex/Radex/bin-gcc/radex'
radex_path = '/ism/ismcpak/radex/Radex/bin-gcc/radex'
mol_data_dirPath = '/ism/ismcpak/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'

# parameters that will be passed to radex single partner
infile = {
    'specStr':                  'HCO+',
    'freqRange':                [0, 0],
    'tKin':                     20.0,
    'collisionPartners':        ['H2'],
    'nDensCollisionPartners': [1e4],
    'tBack':                    2.73,
    'molnDens':                1e13,
    'lineWidth':                1.0,
    'runAnother':               0  # run another is not supported
}

"""
infile = {
    'spec_str':                  'CO',
    'freq_range':                [0, 0],
    't_kin':                     30.0,
    'collision_partners':        ['H2','H' , 'e-', 'H+', 'He'],
    'n_dens_collision_partners': [1e1 , 1e2, 1e2, 1e1 , 1e3 ],
    't_back':                    2.73,
    'mol_n_dens':                1e16,
    'line_width':                1.0,
    'run_another':               0        # run another is not supported
}
"""

# creating the radex process instance
radex_obj = Radex(radex_path, mol_data_dirPath)

t0 = time()
# setting put the parameters, running and parsing the output
radex_obj.set_infile(infile)
radex_obj.run(check_input=True, verbose=True)

if radex_obj.get_status() & radex_obj.FLAGS['SUCCESS']:

        # radex_obj.parse_output()
        t1 = time()
        print('time running radex = %f ' % (t1 - t0))
        print('number of iterations = %d' % radex_obj.get_n_iter())

        # fetch the info of the transiotion from 1->0
        hcop10_info = radex_obj.get_transition(1)
        print(hcop10_info)
        
        # printing all transitions and fluxes
        print('header')
        print('------')
        print(radex_obj.output_hdr)
        
        for transition in radex_obj.transitions:
            print(transition['upper'], transition['lower'], transition['fluxcgs'])

        radex_obj.plot_model()
        
else:
    if radex_obj.get_status() &  radex_obj.FLAGS['ITERWARN']:
        print('did not converge')
        
        print('warnings')
        print('--------')
        print(radex_obj.warnings)

pop_up = radex_obj.transitions['pop_up']
pop_down = radex_obj.transitions['pop_down']
