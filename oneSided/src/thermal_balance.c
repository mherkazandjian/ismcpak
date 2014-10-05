#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

double thermal_balance(mesh *msh, const double T, const double Td, const int k, const double AvL, const double *ymol) 
{
  double flux, heating, cooling, fTot;
  
  heating=cooling=0.0;
  
  /* --------------------------------- COOLING ------------------------------------------------------------------------ */
  cooling += metaStableCooling(msh, T, k, ymol);          CHECK_ERR;                       /* COMPUTING METASTABLE COOLING    */
  cooling += rovVibCooling(msh, T, k, ymol[31], ymol[46], ymol[93], ymol[3], ymol[29]);    /* ROTATIONAL-VIBRATIONAL COOLING  */
  cooling += recombinationCooling(msh, T, k, AvL, ymol[3], ymol[29]);                      /* RECOMBINATION COOLING           */
  cooling += OIandHI_LymanAlphaCooling(msh, T, k, ymol[13], ymol[3], ymol[31], ymol[29]);  /* OI and HI laymen_\alpha cooling */
  
  /* --------------------------------- HEATING -------------------------------------------------------------- */
  heating += heating_PAH(msh, T, AvL, k, ymol[29]);                         /* BY SMALL GRAPHITIC GRAINS AND PAH's  */
  heating += heating_carbon(msh, AvL, k, ymol[8]);                          /* DUE TO CARBON IONIZATION             */
  heating += heating_molecular_hydrogen(msh, AvL, k, ymol[31]);             /* DUE TO H2 PHOTODISSOCIATION          */
  heating += heating_viscous(msh, T, AvL, k, ymol[29], ymol[7]);            /* VISCOUS HEATING                      */
  heating += heating_cosmic_rays(msh, k, ymol[31]);                         /* HEATING BY COSMIC RAYS               */
  heating += heating_H2_pumping(msh, T, k, ymol[3], ymol[31], ymol[299]);   /* DUE TO H_2 COLLISIONAL DE-EXCITATION */  
  heating += pv.gamma_mech;                                                 /* MECHANICAL HEATING                   */
    
  /* ----------------------------------HEATING or COOLING ----------------------------------------------------*/
  flux = fineStructureCooling(msh, T, k, ymol);       CHECK_ERR;            /* FINE STRUCTURE COOLING               */
  if( flux > 0 ) 
      cooling += fabs(flux);
  else
      heating += fabs(flux);
  
  if( T < Td ) 
    heating += fabs(heating_gas_grain_collisions(msh, T, Td, k));           
  else
    cooling += fabs(heating_gas_grain_collisions(msh, T, Td, k));
  
  if( verbose & VERB4_MASK) fprintf(outputFd, "%scooling rate  = ", INDENT4);
  if( verbose & VERB4_MASK) fprintf(outputFd, "%+e\n", cooling);
  if( verbose & VERB4_MASK) fprintf(outputFd, "%sheating rate  = ", INDENT4);
  if( verbose & VERB4_MASK) fprintf(outputFd, "%+e\n", heating);
  if( verbose & VERB4_MASK) fprintf(outputFd, "%s--------------------------------------\n", INDENT4);
  
  if( heating < 0.0 || cooling < 0.0) errMsg(11); CHECK_ERR;

  msh->heating.total[k]=heating;
  msh->cooling.total[k]=cooling;
  
  fTot=heating - cooling;
  return fTot;
}
 
