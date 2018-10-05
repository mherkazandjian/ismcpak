#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

void computeTransitions1( const int transitions, const double *weights, const int *lower, const int *upper,
			  const double *energy,  const double *elec_a, const double *elec_b,
			  const double Hdens, const double elecDens, const double T, double **DFDX )
{
  int i;
  double coeff, Cul, Clu;
  int li, ui;     /* lower, upper wight indicies */
  double lw, uw;  /* lower, upper wieghts        */

  for (i=1; i <= transitions; i++) {
    
    li = lower[i];      ui = upper[i];
    lw = weights[li];   uw = weights[ui];
    
    coeff = elec_a[i] * pow(T/1e4,elec_b[i]);
    Cul = elecDens * coeff + Hdens * 1e-12;
    Clu = (uw / lw) * Cul * exp( -energy[i] / T);

    DFDX[li][li] -= Clu;
    DFDX[li][ui] += Cul;
    DFDX[ui][li] += Clu;
    DFDX[ui][ui] -= Cul;
  }
}
void computeTransitions2(const int levels, const int transitions, const double *weights, const int *lower, const int *upper,
			 const double *energy, const double *EinsA, const double *elec_a, const double *elec_b,  
			 const double H_dens, const double elecDens, const double T,  
			 double *F, double **DFDX, double *nu, double *Pnu  )
{
  int i, *indx;
  double Aul, lambda, background, Bul, Blu;
  int li, ui;     /* lower, upper wight indicies */
  double lw, uw;  /* lower, upper wieghts        */
  
  for (i=1; i <= transitions; i++) {
    
    li = lower[i];     ui = upper[i];
    lw = weights[li];  uw = weights[ui]; 
    
    Aul = EinsA[i];
    nu[i] = KB * energy[i] / HP;
    lambda = CLIGHT / nu[i];
    background = 0e0; /*calc_background_radiation(lambda, G0);*/ 
    Bul = 0e0; /* pow(c,2e0) / (2e0 * HP * pow(nu[i], 3e0)) * Aul; */ 
    Blu = 0e0; /* uw / lw * Bul; */ 
    Pnu[i] = 0e0; /* background; */ 		
    
    DFDX[li][li] -= Blu * background * 0.5;   
    DFDX[li][ui] += (Aul + Bul * background) * 0.5;        
    DFDX[ui][li] += Blu * background * 0.5;       
    DFDX[ui][ui] -= (Aul + Bul * background) * 0.5;   
  }

  for (i=1; i <= levels; i++) {
    if (i == levels) F[i] = 1e0;  else F[i] = 0e0;  
    DFDX[levels][i] = 1e0;   
  }

  MAKE_1D_ARRAY(indx, chemNet.nSpec+1, int);
  ludcmp(DFDX, levels, indx);
  lubksb(DFDX, F, levels, indx);
  FREE_1D_ARRAY(indx,chemNet.nSpec+1, int);
}
/*----------------------------------------------------------------------------------------------------*/
double computeMetastableCoolingRate(const int levels, const int transitions, const double *weights, const int *lower, const int *upper,
				    const double *energy,  const double *EinsA,	const double *elec_a, const double *elec_b, 
				    const double popDensThis, const double Hdens, const double elecDens, const double T)
{
  double cool, Aul, nu[levels+1], Pnu[levels+1], Snu[levels+1], *F, **DFDX;
  int i;
  int li, ui;     /* lower, upper wight indicies */
  double lw, uw;  /* lower, upper wieghts        */
  
  MAKE_1D_ARRAY(F,levels+1,double);
  MAKE_2D_ARRAY(DFDX,levels+1,levels+1,double);

  computeTransitions1(transitions, weights, lower, upper, energy, elec_a, elec_b, Hdens, elecDens, T, DFDX);
  computeTransitions2(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, Hdens, elecDens, T, F, DFDX, nu, Pnu);

  cool = 0e0;
  for (i=1; i <= transitions; i++) {
    li = lower[i];     ui = upper[i];
    lw = weights[li];  uw = weights[ui];
    Snu[i] = 2e0 * HP * pow(nu[i], 3e0) / SQ(CLIGHT) *            /* eq-5.26 rowin thesis pp 56 */
      pow((uw * popDensThis * F[li]) / 
	  (lw * popDensThis * F[ui]) - 1e0, -1e0);
    nu[i] = KB * energy[i] / HP;
    Aul = EinsA[i];  
    cool += 0.5 * popDensThis * F[ui] * Aul * HP * nu[i];// * (Snu[i] - Pnu[i]) / Snu[i] * 0.5;   /* eq-5.25 rowin thesis pp 56 */
  }
  
  FREE_1D_ARRAY(F,levels+1,double);
  FREE_2D_ARRAY(DFDX,levels+1,levels+1,double);
  
  return cool;
}
/***** Cooling done by [CI] emission ***********************************/
double metastable_cooling_CI(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 5.0   , 9.0    }; /* 1S_{0} weight = 1, 1D_{2} weight=5, 3P_{0,1,2} weight=1+3+5 = 9 } */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 1.5e4 , 1.7e4 , 3.1e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 3.4e-4, 0.5   , 2.6e-3 };
  const double elec_a [ transitions+1 ] = { 0.0, 2.7e-8, 1.2e-8, 1.3e-8 };
        double elec_b [ transitions+1 ];

  if (T < 1e4)  { elec_b[1] = 0.57;  elec_b[2] = 0.57; elec_b[3] = 0.57; }
  if (T >= 1e4) { elec_b[1] = -0.13; elec_b[2] = 0e0;  elec_b[3] = 0e0;  }
  
  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [CII] emission ***********************************/
double metastable_cooling_CII(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      2
#define transitions 1
  const double weights[ levels+1      ] = { 0.0, 12.0   , 6.0  }; /* 4P_{5/2,3/2,1/2} weight=12, 2P_{3/2,1/2} weight=6 */
  const int    lower  [ transitions+1 ] = { 0  , 1      };
  const int    upper  [ transitions+1 ] = { 0  , 2      };
  const double energy [ transitions+1 ] = { 0.0, 6.2e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 3.6    };
  const double elec_a [ transitions+1 ] = { 0.0, 1.9e-8 };
  const double elec_b [ transitions+1 ] = { 0.0, -0.5   };

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [FeI] emission ***********************************/
double metastable_cooling_FeI(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 9.0   , 11.0  , 9.0    }; /* 5F_{4} weight=9, 5F_{5} weight=11, 5D_{4} weight=9 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 9.9e3 , 6.4e2 , 1.1e4  }; 
  const double EinsA  [ transitions+1 ] = { 0.0, 2.0e-3, 3.6e-3, 1.5e-3 };
  const double elec_a [ transitions+1 ] = { 0.0, 2.0e-7, 1.5e-7, 1.0e-7 };
        double elec_b [ transitions+1 ];

  if (T < 1e4)  { elec_b[1] = 0.57;  elec_b[2] = 0.0; elec_b[3] = 0.57; }
  if (T >= 1e4) { elec_b[1] = -0.13; elec_b[2] = 0e0; elec_b[3] = 0e0;  }

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/****** Cooling done by [FeII] emission ********************************/
double metastable_cooling_FeII(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 8.0   , 10.0   , 10.0    }; /* 4D_{7/2} weights=8, 6F_{9/2} weights=10, 6D_{9/2} weights=10 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 2.7e3 , 8.6e3  , 1.1e4   };
  const double EinsA  [ transitions+1 ] = { 0.0, 2.8e-4, 1.9e-3 , 5.6e-3  };
  const double elec_a [ transitions+1 ] = { 0.0, 2.2e-8, 2.5e-8 , 5.2e-8  };
  const double elec_b [ transitions+1 ] = { 0.0, -0.50 , -0.50  , -0.50   };

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [OI] emission ***********************************/
double metastable_cooling_OI(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 5.0   , 9.0    };  /* 1S_{0} weight=1, 1D_{2} weight=5, 3P_{0,1,2} weight=1+3+5=9 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 2.3e4 , 2.6e4 , 4.9e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 6.7e-3, 1.3   , 6.7e-2 };
  const double elec_a [ transitions+1 ] = { 0.0, 5.1e-9, 5.2e-9, 2.5e-9 };
        double elec_b [ transitions+1 ];

  if (T < 1e4)  { elec_b[1] = 0.57; elec_b[2] = 0.57; elec_b[3] = 0.57; }
  if (T >= 1e4) { elec_b[1] = 0.17; elec_b[2] = 0.15; elec_b[3] = 0.13; }

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [OII] emission ***********************************/
double metastable_cooling_OII(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 4.0   , 6.0   , 4.0    }; /* 2D_{3/2} weight=4, 2D_{5/2} weight=6, 4S_{3/2} weight=4 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 3.9e4 , 3e1   , 3.9e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 5.1e-5, 1.3e-7, 1.7e-4 };
  const double elec_a [ transitions+1 ] = { 0.0, 1.2e-8, 2.9e-8, 1.2e-8 };
  const double elec_b [ transitions+1 ] = { 0.0, -0.5  , -0.5  , -0.5   };

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [SI] emission ***********************************/
double metastable_cooling_SI(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 5.0   , 9.0    };   /* 1S_{0} weight=1, 1D_{2} weight=5, 3P_{0,1,2} weight=1+3+5=9 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 1.3e4 , 1.9e4 , 3.2e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 3.6e-2, 1.8   , 0.36   };
  const double elec_a [ transitions+1 ] = { 0.0, 5.1e-8, 2.3e-8, 2.5e-8 };
        double elec_b [ transitions+1 ];

  if (T < 1e4)  { elec_b[1] = 0.57;  elec_b[2] = 0.57; elec_b[3] = 0.57; }
  if (T >= 1e4) { elec_b[1] = -0.13; elec_b[2] = 0.0;  elec_b[3] = 0.0;  }

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [SII] emission ***********************************/
double metastable_cooling_SII(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 6.0   , 4.0   , 4.0    }; /* 2D_{5/2} weight=6, 2D_{3/2} weight=4, 4S_{3/2} weight=4 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 2.1e4 , 4.5e1 , 2.1e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 1.8e-3, 3.3e-7, 4.7e-4 };
  const double elec_a [ transitions+1 ] = { 0.0, 5.6e-8, 8e-8  , 5.6e-8 };
  const double elec_b [ transitions+1 ] = { 0.0, -0.5  , -0.5  , -0.5   };

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [SiI] emission ***********************************/
double metastable_cooling_SiI(const double densThis, const double Hdens, const double elecDens, const double T) 
{
#define levels      3
#define transitions 3
  const double weights[ levels+1      ] = { 0.0, 1.0   , 5.0   , 9.0    }; /* 1S_{0} weight=1, 1D_{2} weight=5, 3P_{2,1,0} weight=1+3+5=9 */
  const int    lower  [ transitions+1 ] = { 0  , 2    , 1    , 1     };
  const int    upper  [ transitions+1 ] = { 0  , 3    , 2    , 3     };
  const double energy [ transitions+1 ] = { 0.0, 9.1e3 , 1.3e4 , 2.2e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 3.6e-2, 0.8   , 3.7e-2 };
  const double elec_a [ transitions+1 ] = { 0.0, 6.2e-8, 2.8e-8, 3e-8   };
        double elec_b [ transitions+1 ];

  if (T < 1e4)  { elec_b[1] = 0.57;  elec_b[2] = 0.57; elec_b[3] = 0.57; }
  if (T >= 1e4) { elec_b[1] = -0.13; elec_b[2] = 0.0;  elec_b[3] = 0.0;  }

  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
/***** Cooling done by [SiII] emission ***********************************/
double metastable_cooling_SiII(const double densThis, const double Hdens, const double elecDens, const double T)
{
#define levels      2
#define transitions 1
  const double weights[ levels+1      ] = { 0.0, 12.0   , 6.0 }; /* 4P_{5/2,3/2,1/2} weight=12, 2P_{3/2,1/2} weight=6 */
  const int    lower  [ transitions+1 ] = { 0  , 1      };
  const int    upper  [ transitions+1 ] = { 0  , 2      };
  const double energy [ transitions+1 ] = { 0.0, 6.2e4  };
  const double EinsA  [ transitions+1 ] = { 0.0, 6.4e3  };
  const double elec_a [ transitions+1 ] = { 0.0, 3.9e-8 };
  const double elec_b [ transitions+1 ] = { 0.0, -0.50  };
  
  return computeMetastableCoolingRate(levels, transitions, weights, lower, upper, energy, EinsA, elec_a, elec_b, densThis, Hdens, elecDens, T );
#undef levels 
#undef transitions
}
