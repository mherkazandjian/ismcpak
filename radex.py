import os
import sys
import numpy
import pylab
import subprocess
import logging
import traceback

from misc import default_logger

radex_bin_path = os.environ['HOME'] + '/ism/code/ismcpak/radex/Radex/bin-gcc/radex'
radex_moldata_path = os.environ['HOME'] + '/ism/ismcpak/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'
print('default radex path       {}'.format(radex_bin_path))
print('defualt radex data path  {}'.format(radex_moldata_path))


class Radex(object):
    """
    Wrapper class which runs radex parses its output into objects.

    - test_radex.py is a sample script for using this wrapper
    - radex_view.py is a UI that makes use of this wrapper

    Upon a successull executation, all the transitions are stored in
    self.transitions.
     
    .. todo:: docment the order in which things must be called.
    
    .. warning:: make sure that the lower and upper string names for the
     transitions are not longer than 10 character, or change the length
     accordingly in self.generate_transition_dtype().
    """
    def __init__(self, exec_path, mol_datadir, logger=None):
        """
        Constructor

        :param exec_path:
        :param mol_datadir:
        :param logger:
        """
        self.exec_path = None if exec_path is None else os.path.expanduser(exec_path)
        """str: The path to the radex executable"""

        self.mol_datadir = None if mol_datadir is None else os.path.expanduser(mol_datadir)
        """
        str: Path to the directory containing the transition data for
        the species
        """

        self.molDataFiles = {
            'CO':   'co.dat',
            '13CO': '13co.dat',
            'HCO+': 'hco+@xpol.dat',
            'HCN':  'hcn.dat',  #hcn@xpol.dat, hcn.dat
            'HNC':  'hnc.dat',
            'CS':   'cs@xpol.dat',
            'C':    'catom.dat',
            'C+':   'c+.dat',
            'O':    'oatom.dat',
            'SiO':  'sio.dat',
            'CN':   'cn.dat'
        }
        """
        dict : dictonary for the species files, the keys of this dict are as
         follows (see radex.py for an example):
          
          SPECIE_STRING : FILENAME

        .. literalinclude:: radex.py
           :lines: 24-34
           :linenos:

        .. todo:: implement a method to generate this dict automatically from
         all the files in the directory :data:`molDataFiles`.

        .. todo:: also implement a dict for the collision partners like
         'e-': 'e', but it might not be necessary since RADEX takes both
         'e' and 'e-'
        """

        self.infile = None
        """
        dict : dictionary that holds all the input parameters. It is used to
         construct the input parameter file that is piped to the radex
         executable. It should be of the form.

          .. code-block:: python
             :linenos:

             inFile = {
              'specStr':                'CO',  
              'freqRange':              [0, 50000],
              'tKin':                   10.0,
              'collisionPartners':      ['H2'],
              'nDensCollisionPartners': [1e3],
              'tBack':                  2.73,
              'molnDens':               1e14, # column density of 'specStr'
              'lineWidth':              1.0,
              'runAnother':             1
             }

        All of the items in this dict correspond to the same parameters in
        radex except: 
           'collisionPartners'
           'nDensCollisionPartners'
           'molDataDir'
           'specStr' 

        'collisionPartners' is a list of strs of species such as 'H2' or 'H'

        .. code-block:: python

           'collisionPartners' : ['H2'] 

        or

        .. code-block:: python

           'collisionPartners' : ['H2','H+','e-']

        'nDensCollisionPartners' should be a list of the same length holding the 
        number density of each collision partner respectively. For example :
        
        .. code-block:: python
          
           'nDensCollisionPartners' : [1e3]

        or
           
        .. code-block:: python
        
           'nDensCollisionPartners' : [1e3, 1e1, 1e-2]
           
        the parameter molData in the input parameter file is constructed by appending
        the value of 'specStr' to self.molDataDir 
        """

        self.n_coll_partner = None
        """
        integer : number of collision partners. This is the length of the
        list self.inFile['collisionPartners']
        """

        self.proccess = None
        """
        subproccess.popen object: the radex process of the radex executable"""

        self.raw_output = None
        """
        str: the output of the run which would be dumped by RADEX when ran
        standalone
        """

        self.warnings = None
        """
        list of strings: A list of strings containing the dumped by radex, if
        any. If there are any warnings, the appropriate flag is set in
        self.#FLAGS
        """

        self.n_iter = None
        """integer: number of iterations used when the run finishes."""

        self.output_hdr    = None
        """string : the header of the raw output. This is useful for inspecting
        wheather the input parameters were constructed properly.
        """

        self.n_transitions = None
        """the number of transitions"""

        self.transitions = None
        """
        an ndarray of dtype self.transitionFormat holding all the info of the
        transitions. The dtype has the following keys:
        
        .. code-block:: python
            
            'upper'    : string         # the upper level of the transition     
            'lower'    : string         # the lower level of the transition 
            'E_up'     : numpy.float64  # Enery of the upper level 
            'Tex'      : numpy.float64  # computed excitation temperatur
            'tau'      : numpy.float64  # computed optical depth
            'T_R'      : numpy.float64  # computed rotational temperature
            'pop_up'   : numpy.float64  # computed pop density in the upper lvl
            'pop_down' : numpy.float64  # computed pop density in the lower lvl
            'fluxKkms' : numpy.float64  # computed flux in $K.km^{-1}.s^{-1}$
            'fluxcgs'  : numpy.float64  # computed flux in cgs
        """

        self.transition_dtype = self.generate_transition_dtype()
        """holds the numpy.dtype of the transitions."""
        
        self.FLAGS = {
            'DEFAULT': 0x000,
            'ERROR': 0x001,
            'PARMSOK': 0x002,
            'RUNOK': 0x004,
            'SUCCESS': 0x008,
            'WARNING': 0x010,
            'ITERWARN': 0x020,

            #out of range parameter flags
            'TKIN_OUT_OF_RANGE': 0x100,
            'N_DENS_COLLIOSION_PARTERS_OUT_OF_RANGE': 0x110,
            'MOL_N_DENS_OUT_OF_RANGE': 0x120,
        }
        """
        dict: Flags set when running radex which can be used to examing the
        output. The flags are :

           .. code-block:: python
           
             'DEFAULT':  Default status, should be called before each run 
             'ERROR':    failed 
             'PARMSOK':  parameters are ok 
             'RUNOK':    radex did all the iterations (but not NECCESARILY converged)
             'SUCCESS':  succeeded with no major warning (output should make sense)
             'WARNING':  success with warnings, warning are assigned to self.warnings.  
             'ITERWARN': number of iterations warning

        See the source for the value of each flag
        """ 

        self.status = self.FLAGS['DEFAULT']
        """int: the default stauts set from self.FLAGS. The flags are set in a
        bitwise fashion. To check if a flag is set, it can be done in the
        following way:
        
        .. code-block:: python
        
           stt = self.get_status() & self.FLAGS['FLAG']
           
        if the flag is set stt would be a number greater than 0, zero otherwise
        (for more details on the values of the flags, see radex.py)
        """

        # attr used for generating the plots
        self.fig = None
        """matplotlib figure: when set, the output is plotted in this figure"""

        self.axs = None
        """
        nump.ndarry matplotlib.axes: run output are plotted in these axes
        objects for a single model self.axs has the shape (4,). When nx is
        larger than 1, the shape of self.#axs is  (4,nx). In other words,
        self.#axs is the object returned by self.fig,
        self.axs = pylab.subplots(4, nx )
        """

        self.nx  = None
        """integer : number of models to run and plot"""

        self.ny  = 4
        """integer : number of horizontal division in the figure (default)"""

        self.logger = None
        
        if logger is None:
            self.logger = default_logger()
        else:
            self.logger = logger 

    def set_infile_parm(self, parm, value):
        """set values to the parameters to be passed to radex.

        :param str parm: A key to be modified or added in the self.inFile dict.
        :param arbitrary value: The value to be set to the key. [the type
         depends on the key to be modified, see :data:`inFile` documentation]
        """
        self.infile[parm] = value

    def get_infile_param(self, parm):
        """
        return a parameter from self.infile given the key in the dict as the
        argument 'parm' which is a string used to extract the vlaue from the
        dict.
        """
        return self.infile[parm]

    def set_infile(self, infile):
        self.infile = infile

    def getInFile(self):
        return self.infile

    def get_status(self):
        return self.status

    def set_flag(self, flag):
        self.status |= self.FLAGS[flag]

    def get_raw_output(self):
        return self.raw_output

    def set_status(self, status):
        self.status = status

    def set_logger(self, logger):
        self.logger = logger

    def set_default_status(self):
        """
        set the default status and self.transitions, self.raw_output and
        self.n_transitions to None
        """
        self.status = self.FLAGS['DEFAULT']
        self.transitions = None
        self.n_transitions = None
        self.raw_output = None

    def get_transition(self, idx):
        """
        Extract individual transition information from the self.#transitions
        list. For example, for CO, get_transition(0) would return the transition
        info for the 1-0 transition

        :param int32 idx: The index of the transition in the transition list.
         must be between 0 and len(self.#transitions)
        :return: dict :data:`transitions` item
        """
        return self.transitions[idx]

    def get_warnings(self):
        return self.warnings

    def get_n_iter(self):
        return self.n_iter

    def print_set_flags(self):
        """
        print in a pretty way all the flags (whether they are set or not) in
        self.status
        """
        print('             FLAG_NAME                     SET')
        print('---------------------------------------------------')
        for i, key in enumerate(self.FLAGS):
            print(
                '%-40s  %r' % (
                    key, (self.get_status() & self.FLAGS[key]) >= 1)
            )

    def set_axes(self, fig, axs):
        """set the axes and figure objects"""
        self.axs = axs
        self.fig = fig

    def filter_colliders(self):
        """
        remove colliders from the dictionary of the self.#inFile if their
        density is less than the minimum accepted value by RADEX (for now, this
        ie 1e-3 cm^{-3})
        """
        # removing the colliders which have abundances less than the one range
        # that radex accepts
        for i, nDense in enumerate(self.infile['nDensCollisionPartners']):
            if nDense <= 1e-3:
                self.logger.debug(
                    'poped %e %s' % (
                        self.infile['nDensCollisionPartners'][i],
                        self.infile['collisionPartners'][i])
                )
                self.infile['nDensCollisionPartners'].pop(i)
                self.infile['collisionPartners'].pop(i)

    def gen_input_file_content_as_str(self):
        """
        Generate the parameter file contents from self.inFile

        This string can be passed to the radex executable, as a string.

        :return: (str)
        """
        self.n_coll_partner = len(self.infile['collisionPartners'])
        
        strng = ''

        # set value for "Molecular data file ?"
        strng  += '%s/%s\n' % (
            self.mol_datadir, self.molDataFiles[ self.infile['specStr']]
        )

        # set dummy value for "Name of output file ?" (
        strng  += 'foo\n'

        # set value for "inimum and maximum output frequency [GHz] ?"
        strng  += '%d %d\n' % (
            self.infile['freqRange'][0], self.infile['freqRange'][1]
        )

        # set value for "Kinetic temperature [K] ?"
        strng  += '%f\n' % (self.infile['tKin'])

        # set value for "Number of collision partners ?"
        strng  += '%d\n' % self.n_coll_partner

        # set value for "Number of collision partners ?  1"
        # and "Type of partner 1 ?"
        # and "Density of collision partner  1 [cm^-3] ?"
        for i in numpy.arange(self.n_coll_partner):
            strng += '%s\n' % self.infile['collisionPartners'][i]
            strng += '%e\n' % self.infile['nDensCollisionPartners'][i]

        # set value for "Background temperature [K] ?"
        strng += '%f\n' % self.infile['tBack']

        # set value for "Molecular column density [cm^-2] ?"
        strng += '%e\n' % self.infile['molnDens']

        # set value for "Line width [km/s] ?"
        strng += '%f\n' % self.infile['lineWidth']

        # set value for "Another calculation [0/1] ?"
        strng += '%d'   % self.infile['runAnother']

        # DEBUG
        #print '------------\n%s\n-----------------\n' % strng
        return strng

    def check_parameters(self):
        """
        chech whether the contents of self.infile are within the ranges where
        RADEX can work.

        :raise:  exception NameException
        :attention: This method sets the value of self.status. The flag
         'PARMSOK' is set if the paramteres are ok, 'ERROR' is set otherwise.
        """
        infile = self.infile
        
        for item in infile:
            if infile[item] is None:
                strng = 'Error : missing input parameter %s ' % item
                raise NameError(strng)
            
        # check for correct range for the kinetic temperature
        if infile['tKin'] <= 0.1 or infile['tKin'] >= 1e4 :
            self.set_flag('TKIN_OUT_OF_RANGE')
            self.set_flag('ERROR')

        # check for correct range for the densities of the collion partners
        for i, collPartner in enumerate(infile['collisionPartners']):
            n_dens = infile['nDensCollisionPartners'][i]
            if n_dens < 1e-3 or n_dens > 1e12 :
                self.set_flag('N_DENS_COLLIOSION_PARTERS_OUT_OF_RANGE')
                self.set_flag('ERROR')

        # check for correct range for the species column density
        if infile['molnDens'] < 1e5 or infile['molnDens'] > 1e25:
            self.set_flag('MOL_N_DENS_OUT_OF_RANGE')
            self.set_flag('ERROR')

        if self.flag_is_set('ERROR'):
            return self.FLAGS['ERROR']
        else:
            self.set_flag('PARMSOK')
            return self.FLAGS['PARMSOK']

    def run(self, check_input=None, verbose=None):
        """
        Run the radex executable.
        
        :param bool check_input: By default this is False. In this case, the
         input parameters are not checked and the flag 'PARMSOK' (see #FLAGS) is
         set to :data:`status`. Otherwise, set this to True to force a paremter
         input check.
        :param bool verbose: By default this is False. Set it to True to to
         write the raw output to stdout.
        :return: (int) :data:`status`. Upon a successful run, 'RUNOK' and
         'SUCCESS' flags are set. If the number of iterations excceeds 10,000
         'RUNOK', 'WARNING' and 'ITERWARN' flags are set. if the 'PARMSOK' is
         true :data:`rawOutput` and :data:`transitions` are set, otherwise they
         remaine None.

        .. todo:: extract other warnings also, not just the one due to the max
        iterations.
              
        .. warning:: when running the same radex instance multiple times, make
         sure to set the status to the default before calling :data:`run` using
         :data:`set_default_status()`.
        """
        self.warnings = []

        if check_input is True:
            self.check_parameters()
        else:
            self.set_flag('PARMSOK')

        if self.flag_is_set('PARMSOK'):
            radex_input = self.gen_input_file_content_as_str()
            
            if verbose is True:
                print('-'*100)
                print('{}{}{}'.format(
                    '----------------',
                    'radex input parameters passed to the executable',
                    '------------------------'
                ))
                print('-'*100)
                print(radex_input)
                print('-'*100)

            self.proccess = subprocess.Popen(
                self.exec_path,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            radex_output, radex_err = self.proccess.communicate(input=radex_input)

            if verbose is True:
                print('{}{}{}'.format(
                    '-------------',
                    'raw output of Radex stdout',
                    '------------------'
                ))
                print(radex_output)
                print('-'*100)
                print('radex errors:')
                print(radex_err)
                print('-'*100)

            self.raw_output = radex_output
            self.parse_output()
            self.set_flag('RUNOK')
        
            if self.n_iter == 10000:  # .. todo:: get this from radex.inc
                self.set_flag('WARNING')
                self.set_flag('ITERWARN')
            else:
                self.set_flag('SUCCESS')
        else:
            if verbose is True:
                print(
                    'radex.py: parameters NOT ok. radex.status = ', self.status
                )
            
        return self.status
        
    def run_mutiple_trials(self, expected_sum_pop_dens = None,
                                 rel_pop_dens_tol = None, 
                                 change_frac_trial = None, 
                                 max_trials = None,
                                 strict = None,
                                 verbose = None):
        """Runs the same parameters set a few times by doing minor random steps around the
        input original parameters (since although radex sometimes fails for a certain
        parameter set, it does converge for the same set with tiny variations around them).
        
        :param float64 expected_sum_pop_dens: ideally sum pop dens lower should be
          1, for values less than tol the parameters are changed by 'changeFac' and 
          another attempt is done until the sum of pop dens is close to 1 up to tol.
        :param float64 pop_dens_tol: the absolute tolerence for which the solution is accepted to be
          close enough to expected_sum_pop_dens (otherwise we iterate until this is satisfied).
        :param float64 change_frac_trial: the relative amount by which the input 
          parameters to radex are changed.
        :param int max_trials: number of time random changes are done before giving up
        :param bool strict: if this is True, the shift attempts are also done to attempt to 
        get a solution which makes sense even when warnings are issued.
        """
        
        radexOutputMakesSense = False
        
        inFileOrig = self.infile.copy() # original parameters

        nTried = 0
        while not radexOutputMakesSense:
            
            self.set_default_status()
            self.run(check_input= True, verbose = verbose)

            if self.flag_is_set('PARMSOK') == False:
                self.logger.debug('parameters not ok')
                self.print_set_flags()
                break

            if nTried >= max_trials: #no meaningful output after nTried trials
                print(self.infile)
                break # here radexOutputMakesSense woulbe be False
                
            #does shifts when radex issues a warning
            if strict == True and (self.flag_is_set('WARNING') or self.flag_is_set('ITERWARN')):

                self.logger.debug('radex has a warning')
                self.logger.warn('====> radex issued a warning during trial %d' % nTried)
                self.logger.warn('====> varying slightly the input parameters')
                
                if verbose == True:
                    self.logger.debug('---------------------------------------------------------')
                    self.logger.debug('------------------------radex raw output-----------------')
                    self.logger.debug('---------------------------------------------------------')
                    print(self.raw_output)
                    self.logger.debug('---------------------------------------------------------')

                self.rand_shift_inFile_params(inFileOrig, change_frac_trial)                
                nTried += 1
                continue
            
            # checking the integrity of the solution in case it converged
            #if self.flag_is_set('SUCCESS') or self.flag_is_set('WARNING') or self.flag_is_set('ITERWARN'):
            if self.flag_is_set('SUCCESS'):

                #the checks to be done                
                do_check_for_pop_density_sum          = False
                do_check_for_pop_density_continuity   = True
                do_check_for_pop_density_positivity   = True
                do_check_for_pop_density_all_non_zero = True
                  
                #variable which holds the boolean indicating the test results
                pop_dense_check = True

                #checking for the sum of the population densities if it is close to the expected one
                if do_check_for_pop_density_sum:
                    pop_dense_sum_check = self.check_for_pop_density_sum(expected_sum_pop_dens, rel_pop_dens_tol)
                    pop_dense_check &= pop_dense_sum_check

                #checking if the transitions are smooth (i.e ladder has no abrupt order of magnitude jumops)
                if do_check_for_pop_density_continuity:                
                    continuity_check = self.check_for_pop_density_continuity()  
                    pop_dense_check &= continuity_check
                    
                #checking for pop dense positivity (since they can not be negative)
                if do_check_for_pop_density_positivity:
                    pop_density_positivity = self.check_for_pop_density_positivity()
                    pop_dense_check &= pop_density_positivity

                #checking for pop dense are all negligible
                if do_check_for_pop_density_all_non_zero:
                    pop_density_all_non_zero = self.check_for_pop_density_all_non_zero()
                    pop_dense_check &= pop_density_all_non_zero
                
                if  pop_dense_check == False:
                    
                    self.logger.warn('====> pop dens does NOT make sense, trial %d' % nTried)
                    self.logger.warn('====> total pop dense lower = %18.15f' % self.transitions['pop_down'].sum())
                    self.logger.warn('====> number of iterations = %d' % self.n_iter)
                    if do_check_for_pop_density_sum:
                        if pop_dense_sum_check == False:
                            self.logger.warn('====> pop_dense_sum_check = False')
                   
                    if do_check_for_pop_density_continuity:         
                        if continuity_check == False:
                            self.logger.warn('====> continuity_check = False')

                    if do_check_for_pop_density_positivity:
                        if pop_density_positivity == False:
                            self.logger.warn('====> pop_density_positivity = False')

                    if do_check_for_pop_density_all_non_zero:
                        if pop_density_all_non_zero == False:
                            self.logger.warn('====> pop_density_all_nonzero = False')

                    self.rand_shift_inFile_params(inFileOrig, change_frac_trial)                
                    
                    nTried += 1
                    continue
                else:
                    radexOutputMakesSense = True
                    print('====>', 'pop dense make sense, trial %d with %d iterations' % (nTried, self.n_iter))
                    break
            
            else: # radex did not succeed
                #gets here only in case of an error (or a warning when strict = False) 
                self.logger.debug('radex failed')
                self.logger.debug('---------------------------------------------------------')
                self.logger.debug('------------------------radex raw output-----------------')
                self.logger.debug('---------------------------------------------------------')
                print(self.raw_output)
                self.logger.debug('---------------------------------------------------------')
                break  # here radexOutputMakesSense woulbe be False
        #
                        
        return self.status, radexOutputMakesSense
        
    def rand_shift_inFile_params(self, inFileOrig, factor):
        """change the in input parameters 'inFileParms' within an amount randomly determined by 'factor'. 
        This might be useful since there is a numerical bug in radex which fails sometimes and 
        the output varies hugely even for parameters which are a relative distance of 1e-6
        apart. The shifted parameters are set to self.inFile
        """
        
        self.infile['tKin']     = inFileOrig['tKin'] * (1.0 + numpy.random.rand() * factor)
        self.infile['molnDens'] = inFileOrig['molnDens'] * (1.0 + numpy.random.rand() * factor)
        for i, dens in enumerate(self.infile['nDensCollisionPartners']):
            self.infile['nDensCollisionPartners'][i] = inFileOrig['nDensCollisionPartners'][i] * (1.0 + numpy.random.rand() * factor)

    def print_warnings(self):
        """prints the warning strings dumped by Radex"""
        self.logger.debug('--------------warnings---------------')
        for warning in self.warnings:
            print(warning)
            self.logger.debug('--------------warnings---------------')

    def generate_transition_dtype(self):
        """generates the trasition dtype which will be used in assigning the transition info in self.transitionsNdArry.
        """
        fmt = [
               ('upper'   , numpy.str_, 10),   
               ('lower'   , numpy.str_, 10),   
               ('E_up'    , numpy.float64),
               ('Tex'     , numpy.float64), 
               ('tau'     , numpy.float64), 
               ('T_R'     , numpy.float64), 
               ('pop_up'  , numpy.float64), 
               ('pop_down', numpy.float64), 
               ('fluxKkms', numpy.float64), 
               ('fluxcgs' , numpy.float64), 
              ] # 2 * 10b + 8 * 8b = 84b per transition 
       
        return numpy.dtype(fmt)
    
    def flag_is_set(self, flag):
        """return True 'flag' is set, false otherwise"""
        return self.status & self.FLAGS[flag]

    def _parse_output(self):
        output = self.raw_output
        lines = output.splitlines()
        n_lines = len(lines)

        # looking for the part of the output after the warning
        i = 4  # 5th line (skipping into)
        while True:
            if lines[i][0] == '*' and lines[i + 1][0] == '*':
                break
            i += 1

        self.warnings = lines[4:i]  # extracting the warning
        # print self.warnings

        # collecting the header (lines starting with a '*')
        line_num_hdr_start = i
        while True:
            if lines[i][0] == '*' and lines[i + 1][0] != '*':
                break
            i += 1
        line_num_hdr_end = i + 1

        self.output_hdr = lines[line_num_hdr_start:line_num_hdr_end]
        # print  self.outputHdr

        line_num = line_num_hdr_end
        line_splt = lines[line_num].split()
        self.n_iter = numpy.int32(line_splt[3])

        transitions = []

        # -----------------------------------------------------------------
        # parses a line containing the info of a transition line and returns
        # a dict of the into
        def parse_line_data(line):
            info = {}

            # print self.rawOutput
            line_splt = line.split()
            # print lineSplt
            info['upper'] = line_splt[0].strip()
            info['lower'] = line_splt[1].strip()
            info['E_up'] = numpy.float64(line_splt[2])
            info['Tex'] = numpy.float64(line_splt[5])
            info['tau'] = numpy.float64(line_splt[6])
            info['T_R'] = numpy.float64(line_splt[7])
            info['pop_up'] = numpy.float64(line_splt[8])
            info['pop_down'] = numpy.float64(line_splt[9])
            info['fluxKkms'] = numpy.float64(line_splt[10])
            info['fluxcgs'] = numpy.float64(line_splt[11])

            return info

        # ------------------------------------------------------------------

        # look for the first line which starts with a '1' indicating the
        # first transition line
        while True:
            line_num += 1
            if 'transition data' in lines[line_num + 1]:
                line_num += 2
                break

        # parse the output lines info
        for line in lines[line_num:n_lines - 1]:
            transition = parse_line_data(line)
            transitions.append(transition)

        # copy the content of self.transitions into self.transitionsNdarrya
        # :note: do this at one go without storing things first in
        # self.transitiosn
        # ******************************************************
        self.n_transitions = len(transitions)
        transitions_ndarray = numpy.ndarray(
            self.n_transitions,
            dtype=self.transition_dtype
        )
        for i, trans in enumerate(transitions):
            for key in trans.keys():
                transitions_ndarray[i][key] = trans[key]
        self.transitions = transitions_ndarray
        # ******************************************************

    def parse_output(self):
        """
        Once radex exectues and dumps transition information, this method is
        used to extract the line data.

        :return: None. The instance variable self.transitions is set.
           
           .. note:: in parsing the output, the transitions info is between the lines 
             containing the units:
          
               "(K)    (GHz) ..."
             
             and
         
               "Another calculation"
              
            thisway, we can get the number of transitions (the number of new lines, 
            the work on parsing the data without the need to append anythign to a 
            list..just preallocate the dtype and fill in the values.
        """
        try:
            self._parse_output()
        except (ValueError, IndexError) as exc:
            print(self.raw_output)
            err_str = 'parsing the output failed'
            print('-'*100)
            print(repr(exc))
            traceback.print_exc(file=sys.stdout)
            print('-'*100)
            raise NameError(err_str)

    def check_for_pop_density_sum(self, expected_sum, tolerence):
        """
        return True if all the population densities sum is close to one up to
        a tolerence and False otherwise
        """
        total_pop_dens_lower = numpy.sum(self.transitions[:]['pop_down'])
        diff = numpy.fabs(expected_sum - total_pop_dens_lower)

        # check if the sum of the population desity makes sense
        if diff > tolerence:
            return False
        else:
            return True

    def check_for_pop_density_positivity(self):
        """
        return True if all the population densities are positive and False
        otherwise
        """
        inds_neg = numpy.where(self.transitions[:]['pop_down'] < 0.0)[0]
        if inds_neg.size >= 1:
            return False
        else:
            return True        

    def check_for_pop_density_all_non_zero(self):
        """
        return True if all the population densities are non negligibel and False
        otherwise
        """
        inds_tiny = numpy.where(self.transitions[:]['pop_down'] < 1e-10)[0]
        if inds_tiny.size == self.transitions[:]['pop_down'].size:
            return False
        else:
            return True        
         
    def check_for_pop_density_continuity(self):
        """
        return true if the difference between the population densities of two
        consecutive transitions is not greate than 'log_max_diff'
        """
        pop_down = self.transitions['pop_down'] 

        # value of the log of the diff in consecutive pop dense above which it
        # is considered bogus usually a value below 6 sometimes returns an error
        # for good data, but above 6 (7,8,9) nail down the anomaleous ladders
        log_max_diff = 7.0
        dlogx = numpy.log10(pop_down[0:-2]) - numpy.log10(pop_down[1:-1])
        ind = numpy.where(numpy.fabs(dlogx) > log_max_diff)[0]
        if ind.size > 0:
            return False
        else:
            return True

    def setup_plot(self, nx=None, fig=None, axs=None):
        """
        Set up the object to plot the radex output

        :param integer nx: The number of vertical columns to setup (one model
         per column)
        :param matplotlib.figure fig: Set this object to use a figure which
        is pre-defined, otherwise, a new figure object is created. If this is
         not set the axs keyword is igonored and a self.#axs object is created.
        :param numpy.ndarray of matplotlib.axes object : The axes where the
         models will be plotted. Should have dimensions nx x 4
        :return  self.#fig, self.#axs: in this method, also self.nx is set to nx
        """
        self.nx = nx
        
        if fig is None and axs is None:
            self.make_axes(nx)
        else:
            self.set_axes(fig, axs)
            
        return self.fig, self.axs
            
    def make_axes(self, nx=None):
        """create and assings the self.#fig and self.#axs attributes
        
           :param nx: see self.setup_plot()
        """
        if nx is None:
            nx = 1

        # axs[0,0] is the one in the top left corner
        fig, axs = pylab.subplots(
            4, nx,
            sharex=False, sharey=False,
            figsize=(8,8)
        )
        pylab.subplots_adjust(
            left=0.1,
            bottom=0.15,
            right=0.95,
            top=0.95,
            wspace=0.0,
            hspace=0.0
        )
        self.set_axes(fig, axs)
        
    def plotModelInFigureColumn(self, allTrans = None, inAxes = None, title = None, em_unit='fluxcgs'):
        """Plots the output of a certain model in a certain column in a predefined figure
        
          :param allTrans (integer list or numpy.int32): The indicies of the transitions in self.#transitions whose
             data should be plotted, ex : numpy. 
          :param inAxes ( nump.ndarray matplorlib.axes ): The axes colomn in which the model info will
             be plotted. Fluxes are plotted in inAxes[0], T_ex, T_rot are plotted in inAxes[1], 
             optical depth in inAxes[2] and population densiities in inAxes[3]
          :param title (string): The title to be written at the top of the axes column. 
        """
        
        if allTrans is None:
            allTrans = numpy.arange(len(self.transitions))

        nTrans = len(allTrans)
        #----------------flux-------------------------
        axes = inAxes[0]
        axes.lines = []
        xticksStrs = ()
        
        xPlot = allTrans
        yPlot = numpy.ndarray((nTrans,), dtype=numpy.float64)
        for i in numpy.arange(nTrans):
    
            thisTrans = allTrans[i]
            xPlot[i] = thisTrans

            transition = self.get_transition(thisTrans)
            yThis = transition[em_unit]
             
            yPlot[i] = yThis
            
            # construncting the latex string of the transitions
            upperStr = transition['upper'].split('_')
            if len(upperStr) == 2:
                upperStr = '%s_{%s}'% (upperStr[0], upperStr[1])
            else:
                upperStr = upperStr[0]
            lowerStr = transition['lower'].split('_')
            if len(lowerStr) == 2:
                lowerStr = '%s_{%s}'% (lowerStr[0], lowerStr[1])
            else:
                lowerStr = lowerStr[0]
            
            thisTickStr = '$%s-%s$' %( upperStr, lowerStr )
            xticksStrs = xticksStrs + (thisTickStr ,)

        #plotting the intensities    
        axes.semilogy(xPlot, yPlot, 'b')
        if em_unit == 'fluxcgs':
            axes.axis([numpy.min(allTrans), numpy.max(allTrans), 1e-15, 1e-1])
        if em_unit == 'fluxKkms':
            axes.axis([numpy.min(allTrans), numpy.max(allTrans), 1e-4, 1e4])            
        axes.set_xticks( allTrans, minor = False )
        yticks = 10.0**numpy.linspace(numpy.log10(axes.get_ylim()[0]), numpy.log10(axes.get_ylim()[1]), 5)
        print(yticks)
        axes.set_yticks(yticks , minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        ##numpy.savetxt('/home/mher/ism/tmp/intensities.out', numpy.array([xPlot, yPlot]).T, '%e')
        ##print '-----------> saved the file /home/mher/intensities.out'
        
        #---------------Tex and T_R--------------------------
        axes = inAxes[1]
        axes.lines = []        

        xPlot = allTrans
        yPlot1 = numpy.ndarray(nTrans, dtype=numpy.float64)
        yPlot2 = numpy.ndarray(nTrans, dtype=numpy.float64)
        for i in numpy.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.get_transition(thisTrans)
            yThis1 = transition['Tex']
            yThis2 = transition['T_R'] 
    
            yPlot1[i] = yThis1
            yPlot2[i] = yThis2
            
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.semilogy(xPlot, yPlot2, 'r')
        if self.infile != None: #plotting the horizontal line corresponding to the input kinetic temp
            axes.semilogy([0, 1000], [self.infile['tKin'], self.infile['tKin']], 'k--')
        
        axes.axis([numpy.min(allTrans), numpy.max(allTrans), 1, 10000])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        
        #------------------optical depth----------------------------------
        axes = inAxes[2]
        axes.lines = []        

        xPlot = allTrans
        yPlot = numpy.ndarray(nTrans, dtype=numpy.float64)
        for i in numpy.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.get_transition(thisTrans)
            yThis = transition['tau'] 
    
            yPlot[i] = yThis
        axes.plot(xPlot, yPlot, 'b')
        axes.axis([numpy.min(allTrans), numpy.max(allTrans), -1, numpy.max(yPlot)])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        ##numpy.savetxt('/home/mher/ism/tmp/tau.out', numpy.array([xPlot, yPlot]).T, '%f')
        ##print '-----------> saved the file /home/mher/tau.out'
        
        #------------------population densities-----------------------------
        axes = inAxes[3]
        axes.lines = []        

        xPlot = allTrans
        yPlot1 = numpy.ndarray(nTrans, dtype=numpy.float64)
        yPlot2 = numpy.ndarray(nTrans, dtype=numpy.float64)
        for i in numpy.arange(nTrans):
    
            thisTrans = allTrans[i]
    
            xPlot[i] = thisTrans
    
            transition = self.get_transition(thisTrans)
            yThis1     = transition['pop_up']
            yPlot1[i]  = yThis1
            
        axes.semilogy(xPlot, yPlot1, 'b')
        axes.axis([numpy.min(allTrans), numpy.max(allTrans), 1e-16, 1])
        axes.set_xticks( allTrans, minor = False )
        axes.set_xticklabels( xticksStrs, rotation = -45 )
        
        self.setLabels(em_unit=em_unit)
        
    def clearCurves(self):
        """Clears the curves in an axes column. This is usually used when radex fails
        and the data in a certain column need to be removed. Here it is assumes there
        is one column."""

        for ax in self.axs:            
            ax.lines = []
                
    def setLabels(self, em_unit=''):
        """Set the appropriate labels of all the axes"""

        axs = self.axs
        def removeAll_xLabels(ax):
            for tick in ax.axes.xaxis.get_major_ticks():
                tick.label1On = False
        def removeAll_yLabels(ax):
            for tick in ax.axes.yaxis.get_major_ticks():
                tick.label1On = False
                        
        if self.nx == 1:
            axsLeft = axs.tolist()
            axsBotm = [axs[3]]
        else:
            axsLeft = axs[:,0].tolist()
            axsBotm = axs[3,:].tolist()

            # removing all x and y labels from axes to the left
            # of th zeroth column and above the bottom row
            for axRow in axs[0:-1]:
                for ax in axRow[1:]:
                    removeAll_xLabels(ax)
                    removeAll_yLabels(ax)
        
        #removing the x labeles of from axes on the zeroth column
        #excluding the one in the bottom one (lower-left corner)            
        for ax in axsLeft[0:-1]:
            removeAll_xLabels(ax)
        
        #removing the y labeles of from axes on the bottom row, to the right of the 
        #one on the lower-left corner
        for ax in axsBotm[1:]:
            removeAll_yLabels(ax)
                    
        # plotting the ylabels of the axes
        axsLeft[0].set_ylabel('Flux[%s]' % em_unit)
        axsLeft[1].set_ylabel('Tex, Trot')
        axsLeft[2].set_ylabel('tau')
        axsLeft[3].set_ylabel('pop dens')
        # plotting the xlabels
        for ax in axsBotm:
            ax.set_xlabel('Trans') #;;; use the up-down string (rotated 90 deg) to plot the trans labele

        # removing the firt and last labels from the bottom and left axes
        if self.nx > 1:
            for ax in axsBotm + axsLeft:
                xticks = ax.axes.xaxis.get_major_ticks()
                yticks = ax.axes.yaxis.get_major_ticks()
                for tick in [ xticks[0], xticks[-1], yticks[0], yticks[-1] ]:
                    tick.label1On = False
        
        #removing all labels from axes with no lines (only for multiple models)
        if self.nx > 1:
            for ax in axsBotm:
                if len(ax.lines) == 0:
                    removeAll_xLabels(ax)
                    ax.set_xlabel('')
    
    def plot_model(self):
        """plotting a single model
           todo allow for a choice of the number of transitions to be plotted
        """

        self.setup_plot(nx = 1)
        self.plotModelInFigureColumn(allTrans=None, inAxes = self.axs, title='')
        self.setLabels()   
        pylab.show()
        
    def clear(self):
        self.n_coll_partner = None
        self.raw_output = None
        self.warnings = []
        self.n_iter = None
        self.output_hdr = None
        self.n_transitions = None
        self.transitions = None
        self.status = None

    def setupLogger(self):
        """sets up the logger which will prepend info about the printed stuff. Assignes a value to self.logger."""
        
        # setting up the logger                                                                                                                                                                                                              
        # create logger                                                                                                                                                                                                                      
        self.logger = logging.getLogger('simple_example')
        self.logger.setLevel(logging.DEBUG)

        # create console handler and set level to debug                                                                                                                                                                                      
        ch = logging.StreamHandler( sys.stdout )  # setting the stream to stdout                                                                                                                                                             
        ch.setLevel(logging.DEBUG)

        # create formatter                                                                                                                                                                                                                   
        formatter = logging.Formatter('[%(asctime)s %(funcName)s() %(filename)s:%(lineno)s] %(message)s') # this was the original in the example                                                                              

        # add formatter to ch                                                                                                                                                                                                                
        ch.setFormatter(formatter)

        # add ch to logger                                                                                                                                                                                                                   
        self.logger.addHandler(ch)

        return self.logger

    def copy(self):
        radexNew = radex(self.exec_path, self.mol_datadir)
        radexNew.set_infile(self.infile.copy())
        radexNew.logger = self.logger
        
        return radexNew


try:
    from PyQt4 import QtGui, QtCore
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
    
    class radex_gui(QtGui.QWidget):
    
        def __init__(self, parent=None):
            super(radex_gui, self).__init__(parent)
    
            self.setup_radex_and_set_default_vaules()
            
            self.initUI()
    
        def initUI(self):
            
            #-----------------------------------------------------------------------
            # a figure instance to plot on
            self.figure = self.radexObj.fig
    
            # this is the Canvas Widget that displays the `figure`
            # it takes the `figure` instance as a parameter to __init__
            self.canvas = FigureCanvas(self.figure)
    
            # this is the Navigation widget
            # it takes the Canvas widget and a parent
            self.toolbar = NavigationToolbar(self.canvas, self)
            #-----------------------------------------------------------------------
    
            #---------laying out the line edit feilds-------------------------------        
            self.lbl_specStr = QtGui.QLabel('specie')
            default_text = self.radexObj.inFile['specStr']
            self.qle_specStr = QtGui.QLineEdit(default_text, self)
            self.qle_specStr.textChanged[str].connect(self.onChanged_specStr)
    
            self.lbl_freqRange = QtGui.QLabel('freqRange')
            default_text = '%d, %d' % (self.radexObj.inFile['freqRange'][0],self.radexObj.inFile['freqRange'][1])
            self.qle_freqRange = QtGui.QLineEdit(default_text, self)
            self.qle_freqRange.textChanged[str].connect(self.onChanged_freqRange)
    
            self.lbl_T = QtGui.QLabel('T')
            default_text = '%.2e' % (self.radexObj.inFile['tKin'])
            self.qle_T = QtGui.QLineEdit(default_text, self)        
            self.qle_T.textChanged[str].connect(self.onChanged_T)
    
            self.lbl_n_coll = QtGui.QLabel('n(H2)')
            default_text = '%.2e' % (self.radexObj.inFile['nDensCollisionPartners'][0])        
            self.qle_n_coll = QtGui.QLineEdit(default_text, self)
            self.qle_n_coll.textChanged[str].connect(self.onChanged_n_coll)
    
            self.lbl_t_back = QtGui.QLabel('t_back')
            default_text = '%.2e' % (self.radexObj.inFile['tBack'])        
            self.qle_t_back = QtGui.QLineEdit(default_text, self)
            self.qle_t_back.textChanged[str].connect(self.onChanged_t_back)
            
            self.lbl_molnDens = QtGui.QLabel('N(mol)')
            default_text = '%.2e' % (self.radexObj.inFile['molnDens'])
            self.qle_molnDens = QtGui.QLineEdit(default_text, self)
            self.qle_molnDens.textChanged[str].connect(self.onChanged_molnDens)
    
            self.lbl_lineWidth = QtGui.QLabel('lineWidth')
            default_text = '%.2e' % (self.radexObj.inFile['lineWidth'])        
            self.qle_lineWidth = QtGui.QLineEdit(default_text, self)
            self.qle_lineWidth.textChanged[str].connect(self.onChanged_lineWidth)
            #-----------------------------------------------------------------
            
            # Just some button connected to `plot` method
            self.button = QtGui.QPushButton('Plot')
            self.button.clicked.connect(self.plot)
            
            grid = QtGui.QGridLayout()
            grid.setSpacing(10)
            
            grid.addWidget(self.toolbar, 1, 0, 1, 3) #location 1,0 on the grid spanning 1 row and 3 columns
    
            #-----------------------------------------
            grid.addWidget(self.lbl_specStr, 1, 4)
            grid.addWidget(self.qle_specStr, 1, 5, 1, 2)
    
            grid.addWidget(self.lbl_freqRange, 2, 4)
            grid.addWidget(self.qle_freqRange, 2, 5, 1, 2)
    
            grid.addWidget(self.lbl_T, 3, 4)
            grid.addWidget(self.qle_T, 3, 5, 1, 2)
    
            grid.addWidget(self.lbl_n_coll, 4, 4)
            grid.addWidget(self.qle_n_coll, 4, 5, 1, 2)
    
            grid.addWidget(self.lbl_t_back, 5, 4)
            grid.addWidget(self.qle_t_back, 5, 5, 1, 2)
            
            grid.addWidget(self.lbl_molnDens, 6, 4)
            grid.addWidget(self.qle_molnDens, 6, 5, 1, 2)
    
            grid.addWidget(self.lbl_lineWidth, 7, 4)
            grid.addWidget(self.qle_lineWidth, 7, 5, 1, 2)        
            #----------------------------------------
            
            grid.addWidget(self.button , 15, 5)
            
            grid.addWidget(self.canvas , 2, 0, 15, 3)
    
            
            self.setLayout(grid) 
            
            self.setGeometry(300, 300, 1000, 700)
            self.setWindowTitle('Radex')    
            self.show()
    
        def onChanged_specStr(self, text):
            self.radexObj.inFile['specStr'] = str(text)
        def onChanged_freqRange(self, text):
            #range = parse text into an object [rangemin, rangemax] 
            self.radexObj.inFile['freqRange'] = [0,0]
        def onChanged_T(self, text):
            self.radexObj.inFile['tKin'] = numpy.float(text)
        def onChanged_n_coll(self, text):
            self.radexObj.inFile['nDensCollisionPartners'][0] = numpy.float(text)
        def onChanged_t_back(self, text):
            self.radexObj.inFile['tBack'] = numpy.float(text)
        def onChanged_molnDens(self, text):
            self.radexObj.inFile['molnDens'] = numpy.float(text)
        def onChanged_lineWidth(self, text):
            self.radexObj.inFile['lineWidth'] = numpy.float(text)
    
        def plot(self):
    
            radexObj = self.radexObj 
                   
            radexObj.run(check_input= True, verbose = True)
    
            if radexObj.get_status() &  radexObj.FLAGS['SUCCESS'] :
                
                    radexObj.parse_output()
                    
                    # printing all transitions and fluxes
                    print('header')
                    print('------')
                    print(radexObj.outputHdr)
                    
                    for transition in radexObj.transitions:
                        print(
                            transition['upper'],
                            transition['lower'],
                            transition['fluxcgs']
                        )
                    radexObj.plotModelInFigureColumn(inAxes = radexObj.axs, title='')
    
                    self.canvas.draw()                
            else:
                
                if radexObj.get_status() &  radexObj.FLAGS['ITERWARN']:
                    print('did not converge')
                    
                    print('warnings')
                    print('--------')
                    print(radexObj.warnings)
    
        def set_label_text(self):
            '''Sets the'''
            pass 
                
        def setup_radex_and_set_default_vaules(self):
            '''sets the default values of the radex parameters'''
            
            # path of the radex excutable
            radexPath      = radex_bin_path
            molDataDirPath = radex_moldata_path

            # parameters that will be passed to radex single partner
            inFile = { 'specStr'                : 'CO'   ,
                       'freqRange'              : [0, 0] ,
                       'tKin'                   : 100.0  ,
                       'collisionPartners'      : ['H2'] ,
                       'nDensCollisionPartners' : [1e3]  ,
                       'tBack'                  : 2.73   ,
                       'molnDens'               : 1e14   ,
                       'lineWidth'              : 1.0    ,
                       'runAnother'             : 0      }
            
            # creating the radex process instance
            self.radexObj = Radex(radexPath, molDataDirPath)
            self.radexObj.set_infile(inFile)
            self.radexObj.setup_plot(1)
                    
    def run_radex_gui():
            
        app = QtGui.QApplication(sys.argv)
        gui = radex_gui()
        sys.exit(app.exec_())

except:
    print('can not use radex in gui mode..failed to import PyQt4, or matplotlib.backends.backend_qt4agg')
    