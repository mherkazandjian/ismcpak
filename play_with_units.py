"""
<keywords>
units, conversion, ism
</keywords>
<description>
test unit conversions
</description>
"""

from mylib import units, constants


SM_si = 1.98e30
SM_cgs = SM_si * 1000

mass_H_si = 1.67e-27
mass_H_cgs = mass_H_si*1e-3 
 
# converting from SM / pc^3 to H_atoms / cm^3
density_SM_pc3 = 10.0**(-0.5)
density_H_atoms_pc3 = (density_SM_pc3 * constants.M_SUN_CGS / constants.M_PROTON_CGS)/units.PC2CM**3.0

print 'density %e SM/pc^3 = %e H_atoms.cm^-3' % (density_SM_pc3, density_H_atoms_pc3)