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
    specList = [  
               specie('CRPHOT', specType = -1, charge=0 , init=1, mass = 0.0                                      ,),
               specie('PHOTON', specType = -1, charge=0 , init=1, mass = 0.0                                      ,),
               specie('CRP'   , specType = -1, charge=0 , init=1, mass = 0.0                                      ,),
               specie('PAH'   , specType = 0 , charge=0 , init=1, mass = 720.0                                    ,),
               specie('H2V'   , specType = 0 , charge=0 , init=1, mass = 2.0*1.00794, comp = [ ['H', 2, 'self'] ] ,),
               specie('13C'   , specType = 0 , charge=0 , init=1, mass = 13.0033                                  ,),
               specie('Na'    , specType = 0 , charge=0 , init=1, mass = 22.9897                                  ,),
               specie('Mg'    , specType = 0 , charge=0 , init=1, mass = 24.3050                                  ,),
               specie('Si'    , specType = 0 , charge=0 , init=1, mass = 28.0855                                  ,),
               specie('Cl'    , specType = 0 , charge=0 , init=1, mass = 35.4532                                  ,),
               specie('Fe'    , specType = 0 , charge=0 , init=1, mass = 55.8452                                  ,),
               specie('He'    , specType = 0 , charge=0 , init=1, mass = 4.00260                                  ,),
               specie('H'     , specType = 0 , charge=0 , init=1, mass = 1.00794                                  ,),
               specie('M'     , specType = -1, charge=0 , init=1, mass = 0.0                                      ,),
               specie('C'     , specType = 0 , charge=0 , init=1, mass = 12.0107                                  ,),
               specie('N'     , specType = 0 , charge=0 , init=1, mass = 14.0067                                  ,),
               specie('O'     , specType = 0 , charge=0 , init=1, mass = 15.9994                                  ,),
               specie('P'     , specType = 0 , charge=0 , init=1, mass = 30.9737                                  ,),
               specie('S'     , specType = 0 , charge=0 , init=1, mass = 32.0655                                  ,),
               specie('e-'    , specType = 0 , charge=-1, init=1, mass = 5.48576e-5                               ,), 
               ]
    
    return specList