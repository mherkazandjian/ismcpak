from amuse.units import units

lines = {
         'C609'      : {'ismcpak' : ['fineStructureCoolingComponents','C','rate','1-0'], 
                        'latex'   : r'[CI]609\mu m', 
                        'nu'      : 492e9 | units.Hz,
                       },
         'C369'      : {'ismcpak' : ['fineStructureCoolingComponents','C','rate','2-1'], 
                        'latex'   : r'[CI]369\mu m', 
                        'nu'      : 492e9 | units.Hz,
                       },
         'CO1-0'     : {'radexIdx'  : 0,  'latex' : 'CO(1-0)'                              ,},
         'CO2-1'     : {'radexIdx'  : 1,  'latex' : 'CO(2-1)'                              ,},
         'CO3-2'     : {'radexIdx'  : 2,  'latex' : 'CO(3-2)'                              ,},
         'CO4-3'     : {'radexIdx'  : 3,  'latex' : 'CO(4-3)'                              ,},
         'CO5-4'     : {'radexIdx'  : 4,  'latex' : 'CO(5-4)'                              ,},
         'CO6-5'     : {'radexIdx'  : 5,  'latex' : 'CO(6-5)'                              ,},
         'CO7-6'     : {'radexIdx'  : 6,  'latex' : 'CO(7-6)'                              ,},
         'CO8-7'     : {'radexIdx'  : 7,  'latex' : 'CO(8-7'                               ,},
         'CO9-8'     : {'radexIdx'  : 8,  'latex' : 'CO(9-8)'                              ,},
         'CO10-9'    : {'radexIdx'  : 9,  'latex' : 'CO(10-9)'                             ,},
         'CO11-10'   : {'radexIdx'  : 10, 'latex' : 'CO(11-10)'                            ,},
         'CO12-11'   : {'radexIdx'  : 11, 'latex' : 'CO(12-11)'                            ,},
         'CO13-12'   : {'radexIdx'  : 12, 'latex' : 'CO(13-12)'                            ,},
         'CO14-13'   : {'radexIdx'  : 13, 'latex' : 'CO(14-13)'                            ,},
         'CO15-14'   : {'radexIdx'  : 14, 'latex' : 'CO(15-14)'                            ,},
         #------------------------------------------------------------------------------------------
         '13CO1-0'     : {'radexIdx'  : 0,  'latex' : '13CO(1-0)'                              ,},
         '13CO2-1'     : {'radexIdx'  : 1,  'latex' : '13CO(2-1)'                              ,},
         '13CO3-2'     : {'radexIdx'  : 2,  'latex' : '13CO(3-2)'                              ,},
         '13CO4-3'     : {'radexIdx'  : 3,  'latex' : '13CO(4-3)'                              ,},
         '13CO5-4'     : {'radexIdx'  : 4,  'latex' : '13CO(5-4)'                              ,},
         '13CO6-5'     : {'radexIdx'  : 5,  'latex' : '13CO(6-5)'                              ,},
         '13CO7-6'     : {'radexIdx'  : 6,  'latex' : '13CO(7-6)'                              ,},
         #---------------------------------------------------------------------------------------------
         'HCO+1-0'   : {'radexIdx'  : 0,  'latex' : 'HCO+(1-0)'                            ,},
         'HCO+4-3'   : {'radexIdx'  : 3,  'latex' : 'HCO+(4-3)'                            ,},
         'HCO+7-6'   : {'radexIdx'  : 6,  'latex' : 'HCO+(7-6)'                            ,},
         'HCN1-0'    : {'radexIdx'  : 0,  'latex' : 'HCN(1-0)'                            ,},
         'HCN3-2'    : {'radexIdx'  : 2,  'latex' : 'HCN(3-2)'                            ,},
         'HCN4-3'    : {'radexIdx'  : 3,  'latex' : 'HCN(4-3)'                            ,},
        }

'''
    quantity = ['fineStructureCoolingComponents','C','rate','1-0'] # CI 609um
    flux['C609'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
    quantity = ['fineStructureCoolingComponents','C','rate','2-1'] # CI 369um
    flux['C369'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
    quantity = ['fineStructureCoolingComponents','C+','rate','1-0'] # CII 158um
    flux['C+158'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
'''