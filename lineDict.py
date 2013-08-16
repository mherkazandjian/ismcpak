from amuse.units import units
from leidenLambda import molData

lines = {
         'C+158'     : {'ismcpak' : ['fineStructureCoolingComponents','C+','rate','1-0'],
                        'latex'  : r'[CII] 158$\mu m$',  
                        'latex2'  : r'C$^+$ 158$\mum$',
                        'nu'      : 000e0 | units.Hz,
                        'type'    : 'pdr',
                        'specStr' : 'C+',
                        'radexIdx': 0,
                        'u'       : {'P':'2', 'J':'3/2'}, #;;; .. todo:: check if the QN are correct 
                        'l'       : {'P':'2', 'J':'1/2'}, #;;; .. todo:: check if the QN are correct 
                       },
         'C609'      : {'ismcpak' : ['fineStructureCoolingComponents','C','rate','1-0'], 
                        'latex'   : r'[CI] 609$\mu m$', 
                        'nu'      : 000e0 | units.Hz,
                        'type'    : 'pdr',
                        'specStr' : 'C',
                        'radexIdx': 0,
                        'u'       : {'J':'1'}, #;;; .. todo:: check if the QN are correct 
                        'l'       : {'J':'0'}, #;;; .. todo:: check if the QN are correct 
                       },
         'C369'      : {'ismcpak' : ['fineStructureCoolingComponents','C','rate','2-1'], 
                        'latex'   : r'[CI] 369$\mu m$', 
                        'nu'      : 000e0 | units.Hz,
                        'type'    : 'pdr',
                        'specStr' : 'C',
                        'radexIdx': 2,                        
                        'u'       : {'J':'2'}, #;;; .. todo:: check if the QN are correct 
                        'l'       : {'J':'1'}, #;;; .. todo:: check if the QN are correct     
                       },
         'O63'       : {'ismcpak' : ['fineStructureCoolingComponents','O','rate','1-0'], 
                        'latex'   : r'[OI] 63$\mu m$', 
                        'nu'      : 000e0 | units.Hz,
                        'type'    : 'pdr',
                        'specStr' : 'O',
                        'radexIdx': 0,                        
                        'u'       : {'P':'2', 'J':'3'}, #;;; .. todo:: check if the QN are correct 
                        'l'       : {'P':'1', 'J':'3'}, #;;; .. todo:: check if the QN are correct 
                       },
         #----------------------------------------------------------------------------------------------------------
         'CO1-0'     : {'radexIdx'  : 0,  'latex' : 'CO(1-0)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO2-1'     : {'radexIdx'  : 1,  'latex' : 'CO(2-1)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO3-2'     : {'radexIdx'  : 2,  'latex' : 'CO(3-2)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO4-3'     : {'radexIdx'  : 3,  'latex' : 'CO(4-3)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO5-4'     : {'radexIdx'  : 4,  'latex' : 'CO(5-4)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO6-5'     : {'radexIdx'  : 5,  'latex' : 'CO(6-5)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO7-6'     : {'radexIdx'  : 6,  'latex' : 'CO(7-6)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO8-7'     : {'radexIdx'  : 7,  'latex' : 'CO(8-7'   , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO9-8'     : {'radexIdx'  : 8,  'latex' : 'CO(9-8)'  , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO10-9'    : {'radexIdx'  : 9,  'latex' : 'CO(10-9)' , 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO11-10'   : {'radexIdx'  : 10, 'latex' : 'CO(11-10)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO12-11'   : {'radexIdx'  : 11, 'latex' : 'CO(12-11)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO13-12'   : {'radexIdx'  : 12, 'latex' : 'CO(13-12)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO14-13'   : {'radexIdx'  : 13, 'latex' : 'CO(14-13)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO15-14'   : {'radexIdx'  : 14, 'latex' : 'CO(15-14)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         'CO16-15'   : {'radexIdx'  : 15, 'latex' : 'CO(16-15)', 'type'    : 'radex-lvg' , 'specStr': 'CO'  , 'specStr-other' : ['12CO'], },
         #------------------------------------------------------------------------------------------------------------
         '13CO1-0'     : {'radexIdx'  : 0,  'latex' : r'$^{13}$CO(1-0)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO2-1'     : {'radexIdx'  : 1,  'latex' : r'$^{13}$CO(2-1)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO3-2'     : {'radexIdx'  : 2,  'latex' : r'$^{13}$CO(3-2)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO4-3'     : {'radexIdx'  : 3,  'latex' : r'$^{13}$CO(4-3)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO5-4'     : {'radexIdx'  : 4,  'latex' : r'$^{13}$CO(5-4)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO6-5'     : {'radexIdx'  : 5,  'latex' : r'$^{13}$CO(6-5)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO7-6'     : {'radexIdx'  : 6,  'latex' : r'$^{13}$CO(7-6)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO8-7'     : {'radexIdx'  : 7,  'latex' : r'$^{13}$CO(8-7)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO9-8'     : {'radexIdx'  : 8,  'latex' : r'$^{13}$CO(9-8)'  , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO10-9'    : {'radexIdx'  : 9,  'latex' : r'$^{13}$CO(10-9)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO11-10'   : {'radexIdx'  : 10, 'latex' : r'$^{13}$CO(10-10)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO12-11'   : {'radexIdx'  : 11, 'latex' : r'$^{13}$CO(10-11)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO13-12'   : {'radexIdx'  : 12, 'latex' : r'$^{13}$CO(10-12)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO14-13'   : {'radexIdx'  : 13, 'latex' : r'$^{13}$CO(10-13)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO15-14'   : {'radexIdx'  : 14, 'latex' : r'$^{13}$CO(10-14)' , 'type' : 'radex-lvg' , 'specStr': '13CO' ,},
         '13CO16-15'   : {'radexIdx'  : 15, 'latex' : r'$^{13}$CO(16-15)', 'type' : 'radex-lvg' , 'specStr': '13CO' ,},                  
         #------------------------------------------------------------------------------------------------------------
         'HCO+1-0'   : {'radexIdx'  : 0,  'latex' : r'HCO$^+$(1-0)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+2-1'   : {'radexIdx'  : 1,  'latex' : r'HCO$^+$(2-1)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+3-2'   : {'radexIdx'  : 2,  'latex' : r'HCO$^+$(3-2)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+4-3'   : {'radexIdx'  : 3,  'latex' : r'HCO$^+$(4-3)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+5-4'   : {'radexIdx'  : 4,  'latex' : r'HCO$^+$(5-4)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+6-5'   : {'radexIdx'  : 5,  'latex' : r'HCO$^+$(6-5)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+7-6'   : {'radexIdx'  : 6,  'latex' : r'HCO$^+$(7-6)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+8-7'   : {'radexIdx'  : 7,  'latex' : r'HCO$^+$(8-7)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+9-8'   : {'radexIdx'  : 8,  'latex' : r'HCO$^+$(9-8)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+10-9'  : {'radexIdx'  : 9,  'latex' : r'HCO$^+$(10-9)' , 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+11-10' : {'radexIdx'  : 10, 'latex' : r'HCO$^+$(11-10)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+12-11' : {'radexIdx'  : 11, 'latex' : r'HCO$^+$(12-11)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+13-12' : {'radexIdx'  : 12, 'latex' : r'HCO$^+$(13-12)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+14-13' : {'radexIdx'  : 13, 'latex' : r'HCO$^+$(14-13)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+15-14' : {'radexIdx'  : 14, 'latex' : r'HCO$^+$(15-14)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         'HCO+16-15' : {'radexIdx'  : 15, 'latex' : r'HCO$^+$(16-15)', 'type'   : 'radex-lvg' , 'specStr': 'HCO+'  ,},
         #------------------------------------------------------------------------------------------------------------
         'HCN1-0'   : {'radexIdx'  : 0,  'latex' : r'HCN(1-0)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN2-1'   : {'radexIdx'  : 1,  'latex' : r'HCN(2-1)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN3-2'   : {'radexIdx'  : 2,  'latex' : r'HCN(3-2)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN4-3'   : {'radexIdx'  : 3,  'latex' : r'HCN(4-3)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN5-4'   : {'radexIdx'  : 4,  'latex' : r'HCN(5-4)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN6-5'   : {'radexIdx'  : 5,  'latex' : r'HCN(6-5)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN7-6'   : {'radexIdx'  : 6,  'latex' : r'HCN(7-6)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN8-7'   : {'radexIdx'  : 7,  'latex' : r'HCN(8-7)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN9-8'   : {'radexIdx'  : 8,  'latex' : r'HCN(9-8)'  , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN10-9'  : {'radexIdx'  : 9,  'latex' : r'HCN(10-9)' , 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN11-10' : {'radexIdx'  : 10, 'latex' : r'HCN(11-10)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN12-11' : {'radexIdx'  : 11, 'latex' : r'HCN(12-11)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN13-12' : {'radexIdx'  : 12, 'latex' : r'HCN(13-12)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN14-13' : {'radexIdx'  : 13, 'latex' : r'HCN(14-13)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN15-14' : {'radexIdx'  : 14, 'latex' : r'HCN(15-14)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         'HCN16-15' : {'radexIdx'  : 15, 'latex' : r'HCN(16-15)', 'type'   : 'radex-lvg' , 'specStr': 'HCN'  ,},
         #------------------------------------------------------------------------------------------------------------
         'HNC1-0'   : {'radexIdx'  : 0,  'latex' : r'HNC(1-0)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC2-1'   : {'radexIdx'  : 1,  'latex' : r'HNC(2-1)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC3-2'   : {'radexIdx'  : 2,  'latex' : r'HNC(3-2)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC4-3'   : {'radexIdx'  : 3,  'latex' : r'HNC(4-3)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC5-4'   : {'radexIdx'  : 4,  'latex' : r'HNC(5-4)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC6-5'   : {'radexIdx'  : 5,  'latex' : r'HNC(6-5)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC7-6'   : {'radexIdx'  : 6,  'latex' : r'HNC(7-6)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC8-7'   : {'radexIdx'  : 7,  'latex' : r'HNC(8-7)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC9-8'   : {'radexIdx'  : 8,  'latex' : r'HNC(9-8)'  , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC10-9'  : {'radexIdx'  : 9,  'latex' : r'HNC(10-9)' , 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC11-10' : {'radexIdx'  : 10, 'latex' : r'HNC(11-10)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC12-11' : {'radexIdx'  : 11, 'latex' : r'HNC(12-11)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC13-12' : {'radexIdx'  : 12, 'latex' : r'HNC(13-12)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC14-13' : {'radexIdx'  : 13, 'latex' : r'HNC(14-13)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC15-14' : {'radexIdx'  : 14, 'latex' : r'HNC(15-14)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         'HNC16-15' : {'radexIdx'  : 15, 'latex' : r'HNC(16-15)', 'type'   : 'radex-lvg' , 'specStr': 'HNC'  ,},
         #------------------------------------------------------------------------------------------------------------         
         'CN1_0.5-0_0.5' : {'radexIdx'  : 0,  'latex' : r'CN(1$_{1/2}$-0$_{1/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},
         'CN1_1.5-0_0.5' : {'radexIdx'  : 1,  'latex' : r'CN(1$_{3/2}$-0$_{1/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},
         'CN2_1.5-1_1.5' : {'radexIdx'  : 2,  'latex' : r'CN(2$_{3/2}$-1$_{3/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},         
         'CN2_1.5-1_0.5' : {'radexIdx'  : 3,  'latex' : r'CN(2$_{3/2}$-1$_{1/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},         
         'CN2_2.5-1_1.5' : {'radexIdx'  : 4,  'latex' : r'CN(2$_{5/2}$-1$_{3/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},
         'CN3_2.5-2_2.5' : {'radexIdx'  : 5,  'latex' : r'CN(3$_{5/2}$-2$_{5/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},
         'CN3_2.5-2_1.5' : {'radexIdx'  : 6,  'latex' : r'CN(3$_{5/2}$-2$_{3/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},
         'CN3_3.5-2_2.5' : {'radexIdx'  : 7,  'latex' : r'CN(3$_{7/2}$-2$_{5/2}$)', 'type' : 'radex-lvg' , 'specStr': 'CN' ,},         
         #------------------------------------------------------------------------------------------------------------         
         'CS1-0' : {'radexIdx'  : 0,  'latex' : r'CS(1-0)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS2-1' : {'radexIdx'  : 1,  'latex' : r'CS(2-1)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS3-2' : {'radexIdx'  : 2,  'latex' : r'CS(3-2)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS4-3' : {'radexIdx'  : 3,  'latex' : r'CS(4-3)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS5-4' : {'radexIdx'  : 4,  'latex' : r'CS(5-4)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS6-5' : {'radexIdx'  : 5,  'latex' : r'CS(6-5)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS7-6' : {'radexIdx'  : 6,  'latex' : r'CS(7-6)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS8-7' : {'radexIdx'  : 7,  'latex' : r'CS(8-7)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         'CS9-8' : {'radexIdx'  : 8,  'latex' : r'CS(9-8)', 'type' : 'radex-lvg' , 'specStr': 'CS' ,},
         #------------------------------------------------------------------------------------------------------------         

        }

def get_default_specStr(specStr):
    '''returns formal string representation (according to this dictionary) of the specie by trying to 
    match it to the pre-defined other names of the specie. For example if specStr is '12CO' then 'CO'
    is returned.
    '''

    for line in lines:
        if specStr == lines[line]['specStr']:
            return specStr
        else:
            if 'specStr-other' in lines[line]:
                if specStr in lines[line]['specStr-other']:
                    return lines[line]['specStr']
    
    return None

