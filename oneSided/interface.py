from amuse.community import *

class pdrInterface(CodeInterface):
    
    include_headers = ['worker.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="worker", **keyword_arguments)

    @legacy_function
    def initialize():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def add_mesh():
        """ the input values only for rho and G0 are in log10, but that for T and Lmech are not."""
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('id'    , dtype='int32'  , direction=function.OUT, description='mesh ID')
        function.addParameter('rho'   , dtype='float64', direction=function.IN,  description='gas density')
        function.addParameter('T'     , dtype='float64', direction=function.IN,  description='gas temperature')
        function.addParameter('G0'    , dtype='float64', direction=function.IN,  description='radiation field')
        function.addParameter('Lmech' , dtype='float64', direction=function.IN,  description='mechanical heating')
        function.result_type = 'int32'
        return function

    @legacy_function
    def calc_equilibrium():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_errorFlags():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('ids'     , dtype='int32', direction=function.IN     , description  = "ids of the meshes")
        function.addParameter('nMeshes' , dtype='int32', direction=function.LENGTH , description  = "number of meshes")
        function.addParameter('errFlags', dtype='int32', direction=function.OUT    , description  = "errFlags ids or flags returned")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_temperature():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('ids'     , dtype='int32'  , direction=function.IN     , description  = "ids of the meshes")
        function.addParameter('nMeshes' , dtype='int32'  , direction=function.LENGTH , description  = "number of meshes")
        function.addParameter('T'       , dtype='float64', direction=function.OUT    , description  = "temperature statisitcs")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_outputDir():
        function = LegacyFunctionSpecification()
        function.addParameter('outputDirPath', dtype='string', direction=function.IN, description  = "path of the ouptut directory.")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_species_fName():
        function = LegacyFunctionSpecification()
        function.addParameter('species_fName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_underUbundant_fName():
        function = LegacyFunctionSpecification()
        function.addParameter('underUbundant_fName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rate99_fName():
        function = LegacyFunctionSpecification()
        function.addParameter('rate99_fName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_selfSheilding_CO_fName():
        function = LegacyFunctionSpecification()
        function.addParameter('selfSheilding_CO_fName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_rotationalCooling_baseName():
        function = LegacyFunctionSpecification()
        function.addParameter('rotationalCooling_baseName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_vibrationalCooling_baseName():
        function = LegacyFunctionSpecification()
        function.addParameter('vibrationalCooling_baseName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_zeta():
        function = LegacyFunctionSpecification()
        function.addParameter('zeta', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_S_depletion():
        function = LegacyFunctionSpecification()
        function.addParameter('S_depletion', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_TTol():
        function = LegacyFunctionSpecification()
        function.addParameter('TTol', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_CTol():
        function = LegacyFunctionSpecification()
        function.addParameter('CTol', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_metalicity():
        function = LegacyFunctionSpecification()
        function.addParameter('metalicity', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_AvMax():
        function = LegacyFunctionSpecification()
        function.addParameter('AvMax', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_slabSizeCrit():
        function = LegacyFunctionSpecification()
        function.addParameter('slabSizeCrit', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_max_deltaAv():
        function = LegacyFunctionSpecification()
        function.addParameter('max_deltaAv', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_min_deltaAv():
        function = LegacyFunctionSpecification()
        function.addParameter('min_deltaAv', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_min_deltaAv():
        function = LegacyFunctionSpecification()
        function.addParameter('min_deltaAv', dtype='float64', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_maxSlabs():
        function = LegacyFunctionSpecification()
        function.addParameter('maxSlabs', dtype='int32', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_database_fName():
        function = LegacyFunctionSpecification()
        function.addParameter('database_fName', dtype='string', direction=function.IN, description  = "")
        function.result_type = 'int32'
        return function
