import rxn_rates.collisions

r = rxn_rates.collisions.from_cross_sections(
                                             m1=1.00739*1.6602e-24, #mass of H in cgs 
                                             m2=1.00849*1.6602e-24, #mass of H in cgs
                                             T_rng=[1e3, 1e6],
                                             T_res=1000.0,
                                             E_rng=[0.0, 100.0] #energy range used for the integration in eV
                                            )
