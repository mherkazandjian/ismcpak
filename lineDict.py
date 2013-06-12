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
         'CO1-0'     : {'radexIdx'  : 0, 'latex' : 'CO(1-0)'                              ,},
         'CO10-9'    : {'radexIdx'  : 9, 'latex' : 'CO(10-9)'                             ,},
         'HCO+1-0'   : {'radexIdx'  : 0, 'latex' : 'HCO+(1-0)'                            ,},
         'HCO+4-3'   : {'radexIdx'  : 3, 'latex' : 'HCO+(4-3)'                            ,},
         'HCO+7-6'   : {'radexIdx'  : 6, 'latex' : 'HCO+(7-6)'                            ,},
        }

'''
    quantity = ['fineStructureCoolingComponents','C','rate','1-0'] # CI 609um
    flux['C609'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
    quantity = ['fineStructureCoolingComponents','C','rate','2-1'] # CI 369um
    flux['C369'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
    quantity = ['fineStructureCoolingComponents','C+','rate','1-0'] # CII 158um
    flux['C+158'] = (1.0/(2.0*numpy.pi))*pdrMeshObj.compute_integrated_quantity(quantity, Av_range = [0.0, Av_max])
'''