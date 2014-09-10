import chemistry.core
import numpy
from numpy import nan

###declare just one species
s1 = chemistry.core.species()
print s1

#declare and intialize two species
s2 = chemistry.core.species(20,
                            name    = ['CRPHOT', 'PHOTON', 'CRP', 'PAH', 'H2V', '13C', 'Na', 'Mg', 'Si', 'Cl', 'Fe', 'He', 'H', 'M', 'C', 'N', 'O', 'P', 'S', 'e-'],
                            derived = [    0   ,    0    ,   0  ,  0   ,  1   ,   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0  ],  
                            charge  = [    0   ,    0    ,   0  ,  0   ,  0   ,   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0  ],
                            init    = [    1   ,    1    ,   1  ,  1   ,  1   ,   1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1  ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1  ],
                            passive = [    1   ,    1    ,   1  ,  0   ,  0   ,   0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0  ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0  ],
                            comp    = { 'H2': [[2,'H']] }
                           ) 


s2.pprint()
