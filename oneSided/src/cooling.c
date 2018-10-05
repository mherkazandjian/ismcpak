#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/*----------------------------------------------------------------------------------------------------*/
/* depends on paramters : T, k                                                                        */
/* depends on globals   : ymol[3,7,8,12,13,18,19,23,24,27,28,29]                                      */
/* modifies globals     : coolingMS.CI[], coolingMS.CII[], coolingMS.FeI[],                           */
/*                        coolingMS.FeII[], coolingMS.OI[], coolingMS.OII[],                          */
/*                        coolingMS.SI[], coolingMS.SII[], coolingMS.SiI[]                            */
/*                        coolingMS.SiII[]                                                            */
/*                        coolingProcs.metaStable[]                                                   */
/*----------------------------------------------------------------------------------------------------*/
double metaStableCooling(mesh *msh, const double T, const int k, const double *ymol)
{
  coolingMetaStable *coolingMS=&(msh->cooling.ms);
  double cooling=0.0;
  
  if (T > 1e2) {
    coolingMS->CI[k]   = metastable_cooling_CI  (ymol[8],  ymol[3], ymol[29], T);  cooling += coolingMS->CI[k];
    coolingMS->CII[k]  = metastable_cooling_CII (ymol[7],  ymol[3], ymol[29], T);  cooling += coolingMS->CII[k];
    coolingMS->FeI[k]  = metastable_cooling_FeI (ymol[28], ymol[3], ymol[29], T);  cooling += coolingMS->FeI[k];
    coolingMS->FeII[k] = metastable_cooling_FeII(ymol[27], ymol[3], ymol[29], T);  cooling += coolingMS->FeII[k];
    coolingMS->OI[k]   = metastable_cooling_OI  (ymol[13], ymol[3], ymol[29], T);  cooling += coolingMS->OI[k];
    coolingMS->OII[k]  = metastable_cooling_OII (ymol[12], ymol[3], ymol[29], T);  cooling += coolingMS->OII[k];
    coolingMS->SI[k]   = metastable_cooling_SI  (ymol[24], ymol[3], ymol[29], T);  cooling += coolingMS->SI[k];
    coolingMS->SII[k]  = metastable_cooling_SII (ymol[23], ymol[3], ymol[29], T);  cooling += coolingMS->SII[k];
    coolingMS->SiI[k]  = metastable_cooling_SiI (ymol[19], ymol[3], ymol[29], T);  cooling += coolingMS->SiI[k];
    coolingMS->SiII[k] = metastable_cooling_SiII(ymol[18], ymol[3], ymol[29], T);  cooling += coolingMS->SiII[k];
  } else {
    coolingMS->CI[k]   = 0.0;
    coolingMS->CII[k]  = 0.0;
    coolingMS->FeI[k]  = 0.0;
    coolingMS->FeII[k] = 0.0;
    coolingMS->OI[k]   = 0.0;
    coolingMS->OII[k]  = 0.0;
    coolingMS->SI[k]   = 0.0;
    coolingMS->SII[k]  = 0.0;
    coolingMS->SiI[k]  = 0.0;
    coolingMS->SiII[k] = 0.0;
    cooling = 0.0;
  }
  
  coolingMS->total[k]=cooling;
  return cooling;
}
/*----------------------------------------------------------------------------------------------------*/
/* depends on paramters : T, k                                                                        */
/* depends on globals   : ymol[2 3 7 8 13 18 19 24 27 29 31]                                          */
/* modifies globals     : coolingProcs.fineStruct[]                                                   */
/*----------------------------------------------------------------------------------------------------*/
double fineStructureCooling(mesh *msh, const double T, const int k, const double *ymol)
{
  double cooling=0.0;
  coolingFineStruct *fs=&(msh->cooling.fs);
  populationDensitiesFineStruct *popDens=&(msh->cooling.fs.popDens);

  cooling += finestruct_cooling_CP (T, k, ymol[7 ], ymol[29], ymol[3], ymol[31],          msh->dx, popDens->CP,  fs->CP);    
  cooling += finestruct_cooling_SiP(T, k, ymol[18], ymol[29], ymol[3], ymol[31],          msh->dx, popDens->SiP, fs->SiP);
  cooling += finestruct_cooling_CA (T, k, ymol[8] , ymol[29], ymol[3], ymol[31], ymol[2], msh->dx, popDens->CA,  fs->CA);
  cooling += finestruct_cooling_OA (T, k, ymol[13], ymol[29], ymol[3], ymol[31], ymol[2], msh->dx, popDens->OA,  fs->OA);
  cooling += finestruct_cooling_SA (T, k, ymol[24], ymol[3],  ymol[2], ymol[31],          msh->dx, popDens->SA,  fs->SA);
  cooling += finestruct_cooling_FeP(T, k, ymol[27], ymol[29], ymol[3], ymol[31],          msh->dx, popDens->FeP, fs->FeP);
  cooling += finestruct_cooling_SiA(T, k, ymol[19], ymol[2] , ymol[3], ymol[31],          msh->dx, popDens->SiA, fs->SiA);

  msh->cooling.fs.total[k]=cooling;
  return cooling;
}
/*----------------------------------------------------------------------------------------------------*/
/* depends on paramters : T, k                                                                        */
/* depends on           : pv.dens0, DELTAV, H2O_ortho, H2O_para, CO_low_T, H2O, CO                       */
/* modified global      : cooling_rovib[][], coolingProcs.roVib[]                                     */
/*----------------------------------------------------------------------------------------------------*/
double rovVibCooling(mesh *msh, const double T, const int k, 
		     const double H2Dens, const double CODens, const double H2ODens, const double HDens, const double elecDens)
{
  int i;
  double NH2_L, NCO_L, NH2O_L, NtildeH2, NtildeCO, NtildeH2O;
  double ncol_H2, sigH, sigH2, vel, ncol_CO1, k_e, k_H2, ncol_H2O1, ii, bb, dd, beta, deltaE, CC, ncol_CO2, L_e, L_0, ncol_H2O2;
  double cooling, **Lambda_rovibcool, coolingRovibEach[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double d=1.0;

  cooling = 0.0;
  
  NH2_L  = msh->NH2_L[k];          /* setting colon densities */
  NCO_L  = msh->NCO_L[k];
  NH2O_L = msh->NH2O_L[k];

  NtildeH2  = NH2_L  / DELTAV;    /* LET OP EENHEDEN */
  NtildeCO  = NCO_L  / DELTAV;
  NtildeH2O = NH2O_L / DELTAV;
  
  ncol_H2 = H2Dens + 7e0 * HDens + 16e0 * elecDens;  /* eq-5.41 rowin thesis pp 59 */
  
  /*---------*/  
  sigH = 2.3e-15;
  sigH2 = 3.3e-16 * pow(T/1e3, -0.25);
  vel = 1.03e4 * pow(T,0.5);
  ncol_CO1 = H2Dens + 1.414 * HDens * (sigH / sigH2) + 1.3e-8 * elecDens / (sigH2 * vel);   /* eq-5.42 rowin thesis pp 59 */
  L_e = 1.03e-10 * pow(T/3e2, 0.938); // * exp(-3080e0 / T);                                /* eq-5.47 rowin thesis pp 59 */
  L_0 = 1.14e-14 * T * exp(-68.0 / pow(T, 1e0 / 3e0)); // * exp(-3080e0 / T);               /* eq-5.47 rowin thesis pp 59 */
  ncol_CO2 = H2Dens + 50e0 * HDens + elecDens * L_e / L_0;                                  /* eq-5.46 rowin thesis pp 59 */
  /*---------*/  
  
  ii = 1e0;
  bb = 20e0;
  dd = 1.9;
  beta = 11600e0 / T;
  deltaE = 2.48e-4 * bb * (ii + 1e0);
  if (bb <= 1.53) CC = 9.08e3 / (bb * (ii + 1e0));
  if (bb > 1.53) CC = 1.93e4 / (dd * bb * (ii + 1e0)) * exp(-1.18 / pow(dd, 3e0));
  k_e = (3.56e-6 * SQ(d)) / (pow(T, 0.5) * (2e0 - 1e0 / (ii + 1e0))) *          /* eq-5.44 rowin thesis pp 59 */
    exp(-beta * deltaE) * log(CC * deltaE + CC / beta * exp(-0.577 / (1e0 + 2e0 * beta * deltaE)));
  k_H2 = 7.4e-12 * pow(T, 0.5);  
  ncol_H2O1 = H2Dens + 10e0 * HDens + elecDens * k_e / k_H2;                                /* eq-5.43 rowin thesis pp 59 */
  L_e = 2.6e-6 * pow(T, -0.5); // * exp(-2325e0 / T);                           /* eq-5.49 rowin thesis pp 59 */
  L_0 = 0.64e-14 * T * exp(-47.5 / pow(T, 1e0 / 3e0)); // * exp(-2325e0 / T);   /* eq-5.49 rowin thesis pp 59 */
  ncol_H2O2 = H2Dens + 10e0 * HDens + elecDens * L_e / L_0;                                 /* eq-5.48 rowin thesis pp 59 */
  /*---------*/
  
  /* computing the rovib cooling rates */
  MAKE_2D_ARRAY(Lambda_rovibcool, 3, 2, double);
  if (T > 10e0 && T < 100e0) {
    Lambda_rovibcool[0][0] = 0e0; 
    coolingRovibEach[0] = Lambda_rovibcool[0][0];
    if (T > pow(10e0,1.6)) {
      Lambda_rovibcool[0][1] = 0.5 * ncol_H2 * H2Dens * L_Rot_H2(T, ncol_H2);
      coolingRovibEach[1] = Lambda_rovibcool[0][1];
    } else {
      Lambda_rovibcool[0][1] = 0e0;
      coolingRovibEach[1] = Lambda_rovibcool[0][1];
    }
    Lambda_rovibcool[1][0] = 0e0;
    coolingRovibEach[2] = Lambda_rovibcool[1][0];
    if (log10(NtildeH2O) >= 10e0) {
      Lambda_rovibcool[1][1] = 0.5 * ncol_H2O1 * H2ODens * 
	(0.75 * L_Rot_H2O_and_CO(T, log10(NtildeH2O), ncol_H2O1, H2O_ortho) +
	 0.25 * L_Rot_H2O_and_CO(T, log10(NtildeH2O), ncol_H2O1, H2O_para));
      coolingRovibEach[3] = Lambda_rovibcool[1][1];
    } else {
      Lambda_rovibcool[1][1] = 0e0;
      coolingRovibEach[3] = Lambda_rovibcool[1][1];
    }
    Lambda_rovibcool[2][0] = 0e0;
    coolingRovibEach[4] = Lambda_rovibcool[2][0];
    if (log10(NtildeCO) >= 14.5) {
      Lambda_rovibcool[2][1] = 0.5 * ncol_CO1 * CODens * L_Rot_H2O_and_CO(T, log10(NtildeCO), ncol_CO1, CO_low_T); 
      coolingRovibEach[5] = Lambda_rovibcool[2][1];
    } else {
      Lambda_rovibcool[2][1] = 0e0;
      coolingRovibEach[5] = Lambda_rovibcool[2][1];
    }
  }
  
  if (T >= 100) {
    Lambda_rovibcool[0][0] = 0.5 * ncol_H2 * H2Dens * L_Vib_H2(T, ncol_H2);
    coolingRovibEach[0] = Lambda_rovibcool[0][0];
    Lambda_rovibcool[0][1] = 0.5 * ncol_H2 * H2Dens * L_Rot_H2(T, ncol_H2);
    coolingRovibEach[1] = Lambda_rovibcool[0][1];
    if (log10(NtildeH2O) >= 13e0) {
      Lambda_rovibcool[1][0] = 0.5 * ncol_H2O2 * H2ODens * L_Vib_H2O_and_CO(T, log10(NtildeH2O), ncol_H2O2, H2O);
      coolingRovibEach[2] = Lambda_rovibcool[1][0];
    } else {
      Lambda_rovibcool[1][0] = 0e0;
      coolingRovibEach[2] = Lambda_rovibcool[1][0];
    }
    if (log10(NtildeH2O) >= 10e0) {
      Lambda_rovibcool[1][1] = 0.5 * ncol_H2O1 * H2ODens * L_Rot_H2O_and_CO(T, log10(NtildeH2O), ncol_H2O1, H2O);
      coolingRovibEach[3] = Lambda_rovibcool[1][1];
    } else {    
      Lambda_rovibcool[1][1] = 0e0;
      coolingRovibEach[3] = Lambda_rovibcool[1][1];
    }
    if (log10(NtildeCO) >= 13e0) {
      Lambda_rovibcool[2][0] = 0.5 * ncol_CO2 * CODens * L_Vib_H2O_and_CO(T, log10(NtildeCO), ncol_CO2, CO);
      coolingRovibEach[4] = Lambda_rovibcool[2][0];
    } else {
      Lambda_rovibcool[2][0] = 0e0;
      coolingRovibEach[4] = Lambda_rovibcool[2][0];
    }
    if (log10(NtildeCO) >= 14e0) {
      Lambda_rovibcool[2][1] = 0.5 * ncol_CO1 * CODens * L_Rot_H2O_and_CO(T, log10(NtildeCO), ncol_CO1, CO);
      coolingRovibEach[5] = Lambda_rovibcool[2][1];
    } else {
      Lambda_rovibcool[2][1] = 0e0;
      coolingRovibEach[5] = Lambda_rovibcool[2][1];
    }
  }
  
  for (i=0; i < 6; i++) cooling +=  coolingRovibEach[i];
  
  msh->cooling.roVib[k]=cooling;
  return cooling;
}
/*----------------------------------------------------------------------------------------------------*/
/* depends on paramters : T, k, AvL                                                                    */
/* depends on globals   : pv.G0, K_UV, pv.metalicity                                                        */
/* modifies globals     : cooling_recom[]                                                             */
/*----------------------------------------------------------------------------------------------------*/
double recombinationCooling(mesh *msh, const double T, const int k, const double AvL, const double HDens, const double elecDens)
{
  double BET, cooling;
  
  BET = 0.735 / pow(T, 0.068); 
  cooling = 3.49e-30 * pow(T, 0.944) *    /* eq-5.38 rowin thesis pp 58   */
           pow((pv.G0 * exp(-K_UV * AvL) * sqrt(T)/elecDens),BET) *
           elecDens * HDens * pv.metalicity;

  msh->cooling.recom[k]=cooling;
  return cooling;
}
/*----------------------------------------------------------------------------------------------------*/
/* depends on paramters : T, k                                                                        */
/* modifies globals     : coolingProcs.lymanAlpha[]                                                   */
/*----------------------------------------------------------------------------------------------------*/
double OIandHI_LymanAlphaCooling(mesh *msh, const double T, const int k, 
				 const double ODens, const double HDens, const double H2Dens, const double elecDens)
{
  double Ftot=0.0, cooling;

  cooling = 1.8e-24 * ODens * (HDens + H2Dens) * exp(-22800e0 / T);
  /* Ftot += cooling;*/      

  cooling = 7.3e-19 * elecDens * HDens * exp(-118000e0 / T);
  Ftot +=  cooling;
  
  msh->cooling.lymanAlpha[k]=Ftot;
  return Ftot;
}
