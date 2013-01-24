from specie import specie

def baseSpecies():
    """ Defines the elements and basic species from which all the other
        species are made.  Returns a list of :data:`specie` objects which
        define each specie. This is just a sample file, similar files
        can be defined and imported into scripts, or the specie definition
        of this file can be modified.   
        
        .. literalinclude:: baseSpeciesDefault.py
           :lines: 19-38  
           :linenos:
           
        these can be included in a script as :
        
        .. code-block:: python
        
            from baseSpeciesDefault import baseSpecies
            baseSpecs = baseSpecies.baseSpecies()
    """
    specList = [  specie('CRPHOT', specType = -1, charge=0 , init=1),
                  specie('PHOTON', specType = -1, charge=0 , init=1),
                  specie('CRP'   , specType = -1, charge=0 , init=1),
                  specie('PAH'   , specType = 0 , charge=0 , init=1),
                  specie('H2V'   , specType = 0 , charge=0 , init=1, comp = [ ['H',2] ]),
                  specie('13C'   , specType = 0 , charge=0 , init=1),
                  specie('Na'    , specType = 0 , charge=0 , init=1),
                  specie('Mg'    , specType = 0 , charge=0 , init=1),
                  specie('Si'    , specType = 0 , charge=0 , init=1),
                  specie('Cl'    , specType = 0 , charge=0 , init=1),
                  specie('Fe'    , specType = 0 , charge=0 , init=1),
                  specie('He'    , specType = 0 , charge=0 , init=1),
                  specie('H'     , specType = 0 , charge=0 , init=1),
                  specie('M'     , specType = -1, charge=0 , init=1),
                  specie('C'     , specType = 0 , charge=0 , init=1),
                  specie('N'     , specType = 0 , charge=0 , init=1),
                  specie('O'     , specType = 0 , charge=0 , init=1),
                  specie('P'     , specType = 0 , charge=0 , init=1),
                  specie('S'     , specType = 0 , charge=0 , init=1),
                  specie('e-'    , specType = 0 , charge=-1, init=1) ]
    return specList