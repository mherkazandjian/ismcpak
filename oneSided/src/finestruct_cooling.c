#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

void getCollisionalTransitions(const int levels, double *F, double **DFDX, const int transitions, 
			       const double *weights, const int *lower, const int *upper,
			       const double *energy, const double *CulArr, const double T        )
{
  double lw, uw; /* lower, upper wieghts        */
  int    li, ui; /* lower, upper wight indicies */
  int i;
  double Cul, Clu;

  for (i=1; i <= transitions; i++) {    /* Collisional transitions */
    li = lower[i];  lw = weights[li];
    ui = upper[i];  uw = weights[ui];
    Cul = CulArr[i];
    Clu = (uw/lw)*Cul*exp(-(energy[ui]-energy[li])/ T); /* use dE=energy[ui]-energy[li] */
    DFDX[li][li] -= Clu;
    DFDX[li][ui] += Cul;
    DFDX[ui][li] += Clu;
    DFDX[ui][ui] -= Cul;
  }
}
void radiativeTransitions(const int transitions, const double *weights, const int *lower, const int *upper, 
			  const double *EinsA  , const double *energy , double *nu, double *betatau       ,
			  const int slabInd    , const double *dx     , double **popDensThis , double *Pnu, double **DFDX  )
{
  int i, k;
  double lw, uw;   /* lower, upper wieghts        */
  int    li, ui;   /* lower, upper wight indicies */
  double Aul, lambda, background, Bul, Blu, tau;

  for (i=1; i <= transitions; i++) {    /* Radiative transitions */
    li = lower[i];  lw = weights[li];
    ui = upper[i];  uw = weights[ui];
    Aul = EinsA[i];
    nu[i] = KB * (energy[ui] - energy[li]) / HP;
    lambda = CLIGHT / nu[i];
    background = calc_background_radiation(lambda);
    Bul = pow(CLIGHT, 2e0) / (2e0 * HP * pow(nu[i], 3e0)) * Aul;
    Blu = uw / lw * Bul;

    tau = 0e0;      /* Calculate the optical depth */
    for (k=0; k < slabInd; k++) { 
	tau += max(0e0, (Aul * pow(CLIGHT, 3e0)) / (8e0 * M_PI * pow(nu[i], 3e0)) * popDensThis[ui][k] * 
		   ((popDensThis[li][k] * uw) / (popDensThis[ui][k] * lw) - 1e0) * dx[k] / (DELTAV * 1e5));
    }
    betatau[i] = min(0.5, beta_tau_finestructure(tau));
 
    Pnu[i] = background;
    DFDX[li][li] -= Blu * background * betatau[i];
    DFDX[li][ui] += (Aul + Bul * background) * betatau[i]; 
    DFDX[ui][li] += Blu * background * betatau[i]; 
    DFDX[ui][ui] -= (Aul + Bul * background) * betatau[i]; 
  }
}

double computeCoolingRate(const int levels, double *F, double **DFDX, double **popDensThis, const double ymolThis, 
		    const int slabInd , const int transitions, const double *weights, const int *lower, 
		    const int *upper, double *Snu, const double *nu, const double *Pnu,
		    const double *EinsA, double **coolingFsThis, const double *betatau                      )
{
  int i;
  double cool=0.0;
  double Aul;
  double lw, uw;   /* lower, upper wieghts        */
  int    li, ui;   /* lower, upper wight indicies */
  int *indx;						

  for (i=1; i <= levels; i++) {
    if (i == levels) F[i]=1e0; else F[i]=0e0;
    DFDX[levels][i] = 1e0;
  }
  
  MAKE_1D_ARRAY(indx, chemNet.nSpec+1, int);		
  ludcmp(DFDX, levels, indx);				
  lubksb(DFDX, F, levels, indx);			
  FREE_1D_ARRAY(indx,chemNet.nSpec+1, int);

  for (i=1; i <= levels; i++) {
    if ( F[i]==0 ) {
      F[i] = 1e-15;
      printf("WARNING : set population density to machine precision\n");
    }
    popDensThis[i][slabInd] = F[i] * ymolThis; 
  }

  cool = 0e0;
  for (i=1; i <= transitions; i++) {
    li = lower[i];  lw = weights[li];
    ui = upper[i];  uw = weights[ui];
    Snu[i] = 2e0*HP*pow(nu[i],3e0)/SQ(CLIGHT)*pow((uw*ymolThis * F[li])/(lw*ymolThis*F[ui])-1e0,-1e0);
    //printf("aaa %e %e %e\n", ymolThis, F[li], F[ui]);
    Aul = EinsA[i];
    coolingFsThis[i][slabInd] = betatau[i]*ymolThis*F[ui]*Aul*HP*nu[i]*((Snu[i]-Pnu[i])/Snu[i]);
    cool += coolingFsThis[i][slabInd]; 
  }
  FREE_1D_ARRAY(F, levels+1, double);
  FREE_2D_ARRAY(DFDX, levels+1, levels+1, double);
  
  return cool;
}
/*----------------------------------------------------------------------------------------------------*/
double getCoolingRate(const int levels, const int transitions, const double *weights, const int *lower,
		      const int *upper, const double *energy, const double T, const double *EinsA     ,
		      const int slabInd , const double ymolThis, double **popDensFsThis, double **coolingFsThis,
		      const double *CulArr, const double *dx)
{
  double *F, **DFDX, cool;
  double nu[ transitions+1], Pnu[ transitions+1], Snu[ transitions+1], betatau[ transitions+1];
  
  MAKE_1D_ARRAY(F, levels+1, double);
  MAKE_2D_ARRAY(DFDX, levels+1, levels+1, double);
  
  getCollisionalTransitions(levels, F, DFDX, transitions, weights, lower, upper, energy, CulArr, T);
  radiativeTransitions(transitions, weights, lower, upper, EinsA, energy, nu, betatau, slabInd, dx, popDensFsThis, Pnu, DFDX);
  cool=computeCoolingRate(levels, F, DFDX, popDensFsThis, ymolThis, slabInd , transitions, weights, lower, upper, Snu, nu, Pnu, EinsA, coolingFsThis, betatau);
  return cool;
}
/*----------------------------------------------------------------------------------------------------*/
double finestruct_cooling_CP(const double T, const int slabInd, const double ymolThis, const double elecDens,
			     const double HDens, const double H2Dens, const double *dx, double **popDensFsThis, double **coolingFsThis) 
{
#define levels      2
#define transitions 1
  const double weights[ levels+1      ] = { 0.0, 2.0   , 4.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1            };
  const int    upper  [ transitions+1 ] = { 0  , 2            };
  const double energy [ levels+1 ]      = { 0.0, 0e0   , 92e0 };
  const double EinsA  [ transitions+1 ] = { 0.0, 2.4e-6       };
  const double col_e[transitions+1]     = { 0.0, 1.4e-6  * pow(T,-0.37)};
  const	double col_H[transitions+1]     = { 0.0, 5.8e-10 * pow(T,0.02) };
  const	double col_H2[transitions+1]    = { 0.0, 3.1e-10 * pow(T, 0.1) };
  
  int i;
  double CulArr[transitions+1], cool;
  
  for (i=1; i<=transitions; i++) CulArr[i]=elecDens*col_e[i] + HDens*col_H[i] + H2Dens*col_H2[i];

  cool = getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
			slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;
#undef levels 
#undef transitions
}
double finestruct_cooling_SiP(const double T, const int slabInd, const double ymolThis, const double elecDens, 
			      const double HDens, const double H2Dens, const double *dx, double **popDensFsThis, double **coolingFsThis)
{
#define levels      2
#define transitions 1
  const double weights[ levels+1      ] = { 0.0, 2.0   , 4.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1            };
  const int    upper  [ transitions+1 ] = { 0  , 2            };
  const double energy [ levels+1 ]      = { 0.0, 0e0   , 414e0};
  const double EinsA  [ transitions+1 ] = { 0.0, 2.1e-4       };
  const double col_e[transitions+1]     = { 0.0, 1.2e-5 * pow(T,-0.5) };
  const	double col_H[transitions+1]     = { 0.0, 6.5e-10 };
  const	double col_H2[transitions+1]    = { 0.0, 6.5e-10 };

  int i;
  double CulArr[transitions+1], cool;

  for (i=1; i<=transitions; i++) CulArr[i]=elecDens*col_e[i] + HDens*col_H[i] + H2Dens*col_H2[i];

  cool=getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
		      slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
double finestruct_cooling_CA(const double T, const int slabInd, const double ymolThis, const double elecDens, const double HDens,
			     const double H2Dens, const double HpDens, const double *dx, double **popDensFsThis, double **coolingFsThis) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 3.0, 5.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1, 2, 1           };
  const int    upper  [ transitions+1 ] = { 0  , 2, 3, 3           };
  const double energy [ levels+1 ]      = { 0.0, 0e0,  24.0, 63.0  };
  const double EinsA  [ transitions+1 ] = { 0.0, 7.9e-8, 2.7e-7, 2.0e-14};
        double col_e[transitions+1];
  const double col_Hp[transitions+1]    = { 0.0, 1.2e-9 * pow(T/3e2, 0.24), 5.4e-9 * pow(T/3e2, 0.35), 8.0e-10 * pow(T/3e2, 0.57)};
  const double col_H[transitions+1]     = { 0.0, 1.3e-10 * pow(T, 0.045), 7.8e-11 * pow(T, 0.035), 2.0e-10 * pow(T, 0.084) };
	double col_paraH2[transitions+1];
	double col_orthoH2[transitions+1];

  int i;
  double CulArr[transitions+1], cool;
  double fraction_ortho, fraction_para;

  if (T < 1000e0) {
    col_e[1] = exp(-0.925141e1 - 0.773782e0 * log(T) + 0.361184e0 * pow(log(T), 2e0) -
		   0.150892e-1 * pow(log(T), 3e0) - 0.656325e-3 * pow(log(T), 4e0));
    col_e[2] = exp(-0.743870e1 - 0.574430e0 * log(T) + 0.358264e0 * pow(log(T), 2e0) -
		   0.418166e-1 * pow(log(T), 3e0) - 0.235272e-2 * pow(log(T), 4e0)); 
    col_e[3] = exp(-0.769735e1 - 0.130743e1 * log(T) + 0.697638e0 * pow(log(T), 2e0) - 
		   0.111338e0 * pow(log(T), 3e0) + 0.705277e-2 * pow(log(T), 4e0));    
  }
  if (T >= 1000e0) {
    col_e[1] = exp(0.444600e3 - 0.227913e3 * log(T) + 0.425952e2 * pow(log(T), 2e0) -
		   0.347620e1 * pow(log(T), 3e0) + 0.105085e0 * pow(log(T), 4e0));
    col_e[2] = exp(0.386186e3 - 0.202192e3 * log(T) + 0.385049e2 * pow(log(T), 2e0) -
		   0.319268e1 * pow(log(T), 3e0) + 0.978573e-1 * pow(log(T), 4e0));
    col_e[3] = exp(0.350609e3 - 0.187474e3 * log(T) + 0.361803e2 * pow(log(T), 2e0) - 
		   0.303238e1 * pow(log(T), 3e0) + 0.938138e-1 * pow(log(T), 4e0));
  }

  if (T < 150e0) {
    col_paraH2[1] = 1.61e-10 * pow(T,-0.19); 
    col_orthoH2[1] = 1.05e-10 * pow(T,-0.08);
    col_paraH2[2] = 2.09e-10 * pow(T,-0.04);
    col_paraH2[3] = 1.18e-10 * pow(T,-0.07);
  }
  if (T >= 150e0 && T <= 1000) {
    col_paraH2[1] = 3.50e-11 * pow(T,0.13);
    col_orthoH2[1] = 3.72e-11 * pow(T,0.12);
    col_paraH2[2] = 6.71e-11 * pow(T,0.20);
    col_paraH2[3] = 4.22e-11 * pow(T,0.13);
  }
  if (T > 1000e0) {
    col_paraH2[1] = 8.10e-11;
    col_orthoH2[1] = 8.50e-11;
    col_paraH2[2] = 2.61e-10;
    col_paraH2[3] = 1.03e-10;
  }

  col_orthoH2[2] = 4.50e-11 * pow(T, 0.27);
  col_orthoH2[3] = 3.10e-11 * pow(T, 0.18);

  fraction_para = ortho_para_ratio(T);
  fraction_ortho = 1e0 - fraction_para;
  for (i=1; i<=transitions; i++) CulArr[i] = elecDens*col_e[i] + HpDens*col_Hp[i] + HDens*col_H[i] 
				           + H2Dens*(fraction_para*col_paraH2[i] + fraction_ortho*col_orthoH2[i]);

  cool = getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA, 
		       slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
double finestruct_cooling_OA(const double T, const int slabInd, const double ymolThis, const double elecDens, const double HDens,
			     const double H2Dens, const double HpDens, const double *dx, double **popDensFsThis, double **coolingFsThis)  
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ]   = { 0.0, 5.0   , 3.0, 1.0  };
  const int    lower  [ transitions+1 ]   = { 0  , 1, 2, 1           };
  const int    upper  [ transitions+1 ]   = { 0  , 2, 3, 3           };
  const double energy [ levels+1 ]        = { 0.0, 0e0, 228e0, 326e0 };
  const double EinsA  [ transitions+1 ]   = { 0.0, 9.0e-5, 1.7e-5, 1.0e-10 };
  const double col_e[transitions+1]       = { 0.0, 5.8e-12 * pow(T, 0.67), 4.1e-12 * pow(T, 0.69), 3.3e-12 * pow(T, 0.71) };
  const double col_Hp[transitions+1]      = { 0.0, 4.2e-11 * pow(T, 0.5) * exp(-693e0 / T), 7.5e-12 * pow(T, 0.5) * exp(-450e0 / T), 7.5e-12 * pow(T, 0.5) * exp(-1000e0 / T) };
  const double col_H[transitions+1]       = { 0.0, 4.2e-12 * pow(T, 0.67), 1.5e-11 * pow(T, 0.4), 1.1e-12 * pow(T, 0.44) };
  const	double col_paraH2[transitions+1]  = { 0.0, 3.4e-11 * pow(T, 0.32),  3.34e-15 * pow(T, 1.36), 5.77e-11 * pow(T, 0.30)};
  const double col_orthoH2[transitions+1] = { 0.0, 2.45e-11 * pow(T, 0.38), 2.74e-14 * pow(T, 1.06), 4.09e-11 * pow(T, 0.37)};

  int i;
  double CulArr[transitions+1], cool;
  double fraction_ortho, fraction_para;

  fraction_para = ortho_para_ratio(T);
  fraction_ortho = 1e0 - fraction_para;
  for (i=1; i<=transitions; i++) CulArr[i]= elecDens*col_e[i] + HpDens*col_Hp[i] + HDens*col_H[i] + H2Dens*(fraction_para*col_paraH2[i] + fraction_ortho*col_orthoH2[i]);

  cool=getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
		      slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
double finestruct_cooling_SA(const double T, const int slabInd, const double ymolThis, const double HDens, const double HpDens,
			     const double H2Dens, const double *dx, double **popDensFsThis, double **coolingFsThis) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 5.0   , 3.0, 1.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1, 2, 1           };
  const int    upper  [ transitions+1 ] = { 0  , 2, 3, 3           };
  const double energy [ levels+1 ]      = { 0.0, 0e0, 571e0, 826e0 };
  const double EinsA  [ transitions+1 ] = { 0.0, 1.4e-3, 3.0e-4, 7.1e-8 };
  const double col_H2[transitions+1]    = { 0.0, 7.5e-10, 4.2e-10, 7.1e-10 };
  const double col_Hp[transitions+1]    = { 0.0, 3.3e-8, 1.2e-8, 3.3e-8 };
  const double col_H[transitions+1]     = { 0.0, 7.5e-10, 4.2e-10, 7.1e-10 };

  int i;
  double CulArr[transitions+1], cool;

  for (i=1; i<=transitions; i++) CulArr[i]= HpDens*col_Hp[i] + HDens*col_H[i] + H2Dens*col_H2[i];

  cool=getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
		      slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
double finestruct_cooling_FeP(const double T, const int slabInd, const double ymolThis, const double elecDens, const double HDens,
			      const double H2Dens, const double *dx, double **popDensFsThis, double **coolingFsThis) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 10.0   , 8.0, 6.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1, 2, 1            };
  const int    upper  [ transitions+1 ] = { 0  , 2, 3, 3            };
  const double energy [ levels+1 ]      = { 0.0, 0e0, 554e0, 961e0  };
  const double EinsA  [ transitions+1 ] = { 0.0, 2.5e-3, 1.6e-3, 1.5e-9};
  const double col_e[transitions+1]     = { 0.0, 1.8e-6 * pow(T/100e0, -0.5), 8.7e-7 * pow(T/100e0, -0.5), 1.8e-6 * pow(T/100e0, -0.5) };
  const	double col_H2[transitions+1]    = { 0.0, 9.5e-10, 4.7e-10, 5.7e-10 };
  const double col_H[transitions+1]     = { 0.0, 9.5e-10, 4.7e-10, 5.7e-10 };

  int i;
  double CulArr[transitions+1], cool;

  for (i=1; i<=transitions; i++) CulArr[i]=elecDens*col_e[i] + HDens*col_H[i] + H2Dens*col_H2[i];

  cool=getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
		      slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
double finestruct_cooling_SiA(const double T, const int slabInd, const double ymolThis, const double HpDens, const double HDens, 
			      const double H2Dens, const double *dx, double **popDensFsThis, double **coolingFsThis)
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 3.0, 5.0  };
  const int    lower  [ transitions+1 ] = { 0  , 1     , 2  , 1    };
  const int    upper  [ transitions+1 ] = { 0  , 2     , 3  , 3    };
  const double energy [ levels+1 ]      = { 0.0, 0e0, 110e0, 320e0 };
  const double EinsA  [ transitions+1 ] = { 0.0, 8.4e-6, 4.2e-5, 2.4e-10};
  const double col_Hp[transitions+1]    = { 0.0, 7.2e-9, 2.2e-8, 7.2e-9 };
  const double col_H[transitions+1]     = { 0.0, 3.5e-10 * pow(T/100e0, -0.03), 5.0e-10 * pow(T/100e0, 0.17), 1.7e-10 * pow(T/100e0, 0.17) };
  const double col_H2[transitions+1]    = { 0.0, 3.5e-10 * pow(T/100e0, -0.03), 5.0e-10 * pow(T/100e0, 0.17), 1.7e-10 * pow(T/100e0, 0.17) };

  int i;
  double CulArr[transitions+1], cool;

  for (i=1; i<=transitions; i++) CulArr[i]= HpDens*col_Hp[i] + HDens*col_H[i] + H2Dens*col_H2[i];

  cool=getCoolingRate(levels, transitions, weights, lower, upper, energy, T, EinsA,
		      slabInd, ymolThis, popDensFsThis, coolingFsThis, CulArr, dx);
  return cool;

#undef levels 
#undef transitions
}
/* Calculates escape probabilities for finetructure lines, beta_tau */
double beta_tau_finestructure(double tau) 
{
  if (tau == 0.0e0)
    return 0.5e0;
  else 
    if (tau > 0e0 && tau < 7e0)
      return (1.0e0 - exp(-2.34 * tau))/(4.68 * tau);
    else
      return pow(4e0 * tau * sqrt(log(tau/sqrt(M_PI))),-1e0);
}
/* Calculates the background radiation */
double calc_background_radiation(double lambda) 
{
  double nu, taud;
  nu = CLIGHT / lambda;
  taud = tau100 * (0.01e0 / lambda);
  return blackbody(nu,2.7e0) + (0.42 - log(taud)) * taud * blackbody(nu,T0);
}
/* Calculates blackbody radiation for given temperature */
double blackbody(double nu, double T) 
{
  return 2e0 * HP * pow(nu,3e0) / (exp((HP * nu)/(KB * T)) - 1e0) / SQ(CLIGHT);
}
/* Calculates ortho para ratio of H2 */
double ortho_para_ratio(double T) 
{
  double weights[6], sum_ortho=0.0, sum_para=0.0;
  int i;
  const double energy[6] = { 0e0, 170.5, 509.9, 1015.1, 1681.7, 2504e0 };
  
  for (i=0; i <= 5; i++ ) weights[i] = 2e0 * (double)i + 1e0;
  for (i=0; i <= 4; i+=2) sum_para  += weights[i] * exp(-energy[i] / T);
  for (i=1; i <= 5; i+=2) sum_ortho += weights[i] * exp(-energy[i] / T);

  return sum_para / (3e0 * sum_ortho + sum_para);
}
