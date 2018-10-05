#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/*** L is determined for the vibrational cooling of H2 ******************/
double L_Vib_H2(double T, double H2dens) 
{
  double Llte, L0, L;

  Llte = 1.10e-18 * exp(-6744e0 / T);
  L0 = 1.19e-24 * sqrt(T) * exp(-18100e0/ (T + 1190e0)) * exp(-5987e0 / T);  
  L = pow(1e0 / L0 + H2dens / Llte, -1e0);

  /* fprintf(stdout,"%le %le %le %le \n", T, 1e0 / L, 1e0 / L0, H2dens / Llte);*/
  
  return L;

}
/*** Interpolates in the rotational cooling tables for a given column
     density, which determines Llte, L0, alpha and n1/2, and with
     these parameters L is determined. **********************************/
double L_Rot_H2(double T, double H2dens) 
{
  int i, coolspec, dummy;
  double dumtemp1, dumtemp2, dumT, fac;
  double L, L0, Llte, ncrit, alpha;
  
  fac=dummy=L=L0=Llte=ncrit=alpha=0.0;
  
  if (T < 100e0) { 
    coolspec = H2_low_T;
    dumT = T;
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      if (dumT >= pow(10e0,coollines[coolspec].paramH2_rot[0][i]) && dumT < pow(10e0,coollines[coolspec].paramH2_rot[0][i+1])) {
	dummy = i;
	dumtemp1 = pow(10e0,coollines[coolspec].paramH2_rot[0][i]);
	dumtemp2 = pow(10e0,coollines[coolspec].paramH2_rot[0][i+1]);
	fac=(dumT - dumtemp1) / (dumtemp2 - dumtemp1);
      } 
    }
    L0 = pow(10e0,-coollines[coolspec].paramH2_rot[1][dummy]) + 
      fac * (pow(10e0,-coollines[coolspec].paramH2_rot[1][dummy+1]) - pow(10e0,-coollines[coolspec].paramH2_rot[1][dummy]));
    L0 = L0 * exp(-509e0 / T);
    Llte = pow(10e0,-coollines[coolspec].paramH2_rot[2][dummy]) + 
      fac * (pow(10e0,-coollines[coolspec].paramH2_rot[2][dummy+1]) - pow(10e0,-coollines[coolspec].paramH2_rot[2][dummy]));
    Llte = Llte * exp(-509e0 / T);
    ncrit =  pow(10e0,coollines[coolspec].paramH2_rot[3][dummy]) + 
      fac * (pow(10e0,coollines[coolspec].paramH2_rot[3][dummy+1]) - pow(10e0,coollines[coolspec].paramH2_rot[3][dummy]));
    alpha =  coollines[coolspec].paramH2_rot[4][dummy] + 
      fac * (coollines[coolspec].paramH2_rot[4][dummy+1] - coollines[coolspec].paramH2_rot[4][dummy]);
  }
  if (T >= 100e0) { 
    coolspec = H2;
    dumT = log10(T);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      if (dumT >= coollines[coolspec].paramH2_rot[0][i] && dumT < coollines[coolspec].paramH2_rot[0][i+1]) {
	dummy = i;
	dumtemp1 = coollines[coolspec].paramH2_rot[0][i];
	dumtemp2 = coollines[coolspec].paramH2_rot[0][i+1];
	fac=(dumT - dumtemp1) / (dumtemp2 - dumtemp1);
      } 
    }
    L0 = coollines[coolspec].paramH2_rot[1][dummy] + 
      fac * (coollines[coolspec].paramH2_rot[1][dummy+1] - coollines[coolspec].paramH2_rot[1][dummy]);
    L0 = pow(10e0, -L0);
    Llte = coollines[coolspec].paramH2_rot[2][dummy] + 
      fac * (coollines[coolspec].paramH2_rot[2][dummy+1] - coollines[coolspec].paramH2_rot[2][dummy]);
    Llte = pow(10e0, -Llte);
    ncrit =  coollines[coolspec].paramH2_rot[3][dummy] + 
      fac * (coollines[coolspec].paramH2_rot[3][dummy+1] - coollines[coolspec].paramH2_rot[3][dummy]);
    ncrit = pow(10e0, ncrit);
    alpha =  coollines[coolspec].paramH2_rot[4][dummy] + 
      fac * (coollines[coolspec].paramH2_rot[4][dummy+1] - coollines[coolspec].paramH2_rot[4][dummy]); 
  }

  if (T < pow(10e0,1.9)) {
    L = pow(1e0 / L0 + H2dens / Llte, -1e0);
  } else {
    L = pow(1e0 / L0 + H2dens / Llte + 1e0 / L0 * pow(H2dens / ncrit, alpha) * (1e0 - ncrit * L0 / Llte), -1e0);
  }
  
  /*  fprintf(stdout,"%le %le %le %le %le \n",T, 1 / L, 1e0 / L0, H2dens / Llte, 
      1e0 / L0 * pow(H2dens / ncrit, alpha) * (1e0 - ncrit * L0 / Llte) );*/

  return L;


}
/*** Interpolates in the vibrational cooling tables for a given column
     density and temperature and then L is calculated ************************/
double L_Vib_H2O_and_CO(double T, double Ntilde, double H2dens, int coolspec) 
{
  int i,dummy1, dummy2;
  double dumtemp1, dumtemp2, dumdens1, dumdens2;
  double fac1, fac2;
  double Llte1, Llte2, Llte, L, L0;

  L0=0.0;
  dummy1 = coollines[coolspec].rows_vib - 2;
  dummy2 = coollines[coolspec].columns_vib - 2;
  fac1 = 1.0;
  fac2 = 1.0;
  
  for (i=0; i < coollines[coolspec].rows_vib; i++) {
    if (Ntilde >= coollines[coolspec].ColumnDens_vib[i] && Ntilde < coollines[coolspec].ColumnDens_vib[i+1]) {
      dummy1 = i;
      dumdens1 = coollines[coolspec].ColumnDens_vib[i];
      dumdens2 = coollines[coolspec].ColumnDens_vib[i+1];
      fac1=(Ntilde - dumdens1) / (dumdens2 - dumdens1);
    } 
  }
  for (i=0; i < coollines[coolspec].columns_vib; i++) {
    if (T >= coollines[coolspec].Temperature_vib[i] && T < coollines[coolspec].Temperature_vib[i+1]) {
      dummy2 =i;
      dumtemp1 = coollines[coolspec].Temperature_vib[i];
      dumtemp2 = coollines[coolspec].Temperature_vib[i+1];
      fac2=(T - dumtemp1) / (dumtemp2 - dumtemp1);
    } 
  }

  Llte1 = coollines[coolspec].Llte_vib[dummy2][dummy1] + 
    (coollines[coolspec].Llte_vib[dummy2][dummy1 + 1] - 
     coollines[coolspec].Llte_vib[dummy2][dummy1]) * fac1;
  Llte2 = coollines[coolspec].Llte_vib[dummy2 + 1][dummy1] + 
    (coollines[coolspec].Llte_vib[dummy2 + 1][dummy1 + 1] - 
     coollines[coolspec].Llte_vib[dummy2 + 1][dummy1]) * fac1;
  Llte = Llte1 + (Llte2 - Llte1) * fac2;

  if (coolspec == H2O) {
    Llte = pow(10e0,-Llte) * exp(-2325e0 / T);
  }
  if (coolspec == CO) {
    Llte = pow(10e0,-Llte) * exp(-3080e0 / T);
  }

  if (coolspec == H2O) {
    L0 = 1.03e-26 * T * exp(-47.5 / pow(T, 1e0/3e0)) * exp(-2325e0 / T);
  }
  if (coolspec == CO) {
    L0 = 1.83e-26 * T * exp(-68e0 / pow(T, 1e0/3e0)) * exp(-3080e0 / T);
  }
  
  L = pow(1e0 / L0 + H2dens / Llte, -1e0);  
  /*
  fprintf(stdout,"%le %le %le %le \n",T ,1e0 / L, 1e0 / L0, H2dens / Llte);
  */

  return L;

}
/*** Interpolates in the tables for a given column density and 
     temperature *******************************************************/
double L_Rot_H2O_and_CO(double T, double Ntilde, double H2dens, int coolspec) 
{
  int i, dummy1, dummy2;
  double dumtemp1, dumtemp2, dumdens1, dumdens2; 
  double fac1, fac2;
  double L, L0, Llte, Llte1, Llte2, ncrit, ncrit1, ncrit2, alpha, alpha1, alpha2;

  dummy1 = coollines[coolspec].rows_rot - 2;
  dummy2 = coollines[coolspec].columns_rot - 2;
  fac1 = 1.0;
  fac2 = 1.0;

  for (i=0; i < coollines[coolspec].rows_rot; i++) {
    if (Ntilde >= coollines[coolspec].ColumnDens_rot[i] && Ntilde < coollines[coolspec].ColumnDens_rot[i+1]) {
      dummy1 = i;
      dumdens1 = coollines[coolspec].ColumnDens_rot[i];
      dumdens2 = coollines[coolspec].ColumnDens_rot[i+1];
      fac1=(Ntilde - dumdens1) / (dumdens2 - dumdens1);
    } 
  }
  for (i=0; i < coollines[coolspec].columns_rot; i++) {
    if (T >= coollines[coolspec].Temperature_rot[i] && T < coollines[coolspec].Temperature_rot[i+1]) {
      dummy2 =i;
      dumtemp1 = coollines[coolspec].Temperature_rot[i];
      dumtemp2 = coollines[coolspec].Temperature_rot[i+1];
      fac2=(T - dumtemp1) / (dumtemp2 - dumtemp1);
    } 
  }
  L0 = coollines[coolspec].L0_rot[dummy2] + 
    (coollines[coolspec].L0_rot[dummy2 + 1] - coollines[coolspec].L0_rot[dummy2]) * fac2;
  L0 = pow(10e0,-L0);

  Llte1 = coollines[coolspec].Llte_rot[dummy2][dummy1] + 
    (coollines[coolspec].Llte_rot[dummy2][dummy1 + 1] - 
     coollines[coolspec].Llte_rot[dummy2][dummy1]) * fac1;
  Llte2 = coollines[coolspec].Llte_rot[dummy2 + 1][dummy1] + 
    (coollines[coolspec].Llte_rot[dummy2 + 1][dummy1 + 1] - 
     coollines[coolspec].Llte_rot[dummy2 + 1][dummy1]) * fac1;
  Llte = Llte1 + (Llte2 - Llte1) * fac2;
  Llte = pow(10e0, -Llte);

  ncrit1 = coollines[coolspec].nhalf_rot[dummy2][dummy1] + 
    (coollines[coolspec].nhalf_rot[dummy2][dummy1 + 1] - 
     coollines[coolspec].nhalf_rot[dummy2][dummy1]) * fac1;
  ncrit2 = coollines[coolspec].nhalf_rot[dummy2 + 1][dummy1] + 
    (coollines[coolspec].nhalf_rot[dummy2 + 1][dummy1 + 1] - 
     coollines[coolspec].nhalf_rot[dummy2 + 1][dummy1]) * fac1;
  ncrit = ncrit1 + (ncrit2 - ncrit1) * fac2;
  ncrit = pow(10e0, ncrit); 

  alpha1 = coollines[coolspec].alpha_rot[dummy2][dummy1] + 
    (coollines[coolspec].alpha_rot[dummy2][dummy1 + 1] - 
     coollines[coolspec].alpha_rot[dummy2][dummy1]) * fac1;
  alpha2 = coollines[coolspec].alpha_rot[dummy2 + 1][dummy1] + 
    (coollines[coolspec].alpha_rot[dummy2 + 1][dummy1 + 1] - 
     coollines[coolspec].alpha_rot[dummy2 + 1][dummy1]) * fac1;
  alpha = alpha1 + (alpha2 - alpha1) * fac2;

  L = pow(1e0 / L0 + H2dens / Llte + 1e0 / L0 * pow(H2dens / ncrit, alpha) * (1e0 - ncrit * L0 / Llte),-1e0);

  return L;

}
