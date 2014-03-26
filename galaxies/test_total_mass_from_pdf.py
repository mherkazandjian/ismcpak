"""
<keywords>
pdf, mass, total, number of clouds
</keywords>
<description>
test mass_from_pdf
</description>
"""

from galaxies import fi_utils, coset9_sol_info

# radius density scaling function, r(kpc) = R(n[cm-3]) (kpc) 
r_func = coset9_sol_info.r_sph_kpc


## the probability density at that density (not in log scale
nPDF = coset9_sol_info.nPDF_sph.pdf
#nPDF = coset9_sol_info.nPDF_sph_test.pdf

nmin, nmax = 1e-6, 1e12

fi_utils.mass_from_pdf(nPDF, r_func, nmin, nmax, res=100.0)
'''
fi_utils.plot_luminosity_from_pdf(nPDF, r_func, arxvPDR, F, params,  nmin, nmax)
'''