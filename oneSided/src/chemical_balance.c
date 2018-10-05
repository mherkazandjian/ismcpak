#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/*-------------------------------------- Chemical balance --------------------------------------------*/
/* depends on parameters :                                                                            */
/*   AvL, T, k                                                                                        */
/* depends on global variables :                                                                      */
/*   uses as input :                                                                                  */
/*      chemNet.ymol[], betaCO, betaH2, beta13CO                                                      */
/*   uses CONSTANT variables :                                                                        */
/*      ALPHA[], BETA[], GAMMA[], molno[][], indBaseSpec[][], amax[], abun[], specCode[][], pv.G0, pv.zeta  */
/*      ALBEDO, pv.dens0, pv.metalicity, PHI_PAH, nSpec, nReact, nBaseSpec, TINY2                           */
/*   modifies as output                                                                               */
/*      chemNet.ymol[], nNewton                                                                               */
/*----------------------------------------------------------------------------------------------------*/
void chemical_balance(const mesh *msh, const double T, const int k) 
{
  double Td, AvL=msh->AvL[k];
  if( verbose & VERB4_MASK) fprintf(outputFd, "%scomputing chemical balance for T=%.2f\n", INDENT4, T);

  if( verbose & VERB5_MASK) fprintf(outputFd, "%sInitializing chemical network reaction rates for this step.\n", INDENT5);
  Td=getDustTemp(AvL);
  initRatesNetwork( AvL, T, Td, msh->betaH2_L[k], msh->betaCO_L[k], msh->beta13CO_L[k]);

  if( verbose & VERB5_MASK) fprintf(outputFd, "%sgetting equilibirum abundances.\n", INDENT5);
  iterateChemicalNetwork( AvL ); CHECK_ERRV;
  //iterateChemicalNetworkOptimizedCorrect( AvL); CHECK_ERRV;
} 
/*----------------------------------------------------------------------------------------------------*/
/* converges to the equilibirum chemical abundances ( modifies only ymol which is used outside her )  */
/*----------------------------------------------------------------------------------------------------*/
void iterateChemicalNetwork(const double AvL )
{
  int i, gx, niterations, count;
  double yx;
  
  niterations = 0;

  if( verbose & VERB7_MASK) fprintf(outputFd, "%s chemical root refinement : ", INDENT7);

  /*Multi-Dim Newton Raphson, solving for equilibrium absolute abundances, ymol (equlibrium)   */
  gx = 1;
  while (gx == 1) { 

    gx=count=0;
    
    niterations++;
    nNewton++;

    assemble_linear_system(chemNet.F, chemNet.DFDX);                         /* compute F and DFDX          */
    solve_lu(chemNet.DFDX, chemNet.F, chemNet.X, chemNet.nSpec); CHECK_ERRV; /* compute J^-1 V (see notes)  */
    
    for (i = 1; i <= chemNet.nSpec; i++) {
      
      yx = chemNet.ymol[i] - chemNet.X[i];  /* better guess for i^th specie density */
      if (dabs(yx / chemNet.ymol[i] - 1.0e0) > pv.CTol && chemNet.ymol[i] > TINY2) {
	gx = 1;
	count++;
	//printf("%d %e %e %e\n", i, chemNet.ymol[i], chemNet.X[i], ics.maxAbun[i]);
      }
      
      if (yx < TINY2) yx = TINY2;
      if (AvL <= 0.1) chemNet.ymol[i] = min(ics.maxAbun[i], 0.5e0 * yx + 0.5e0 * chemNet.ymol[i]);
      if (AvL > 0.1) chemNet.ymol[i] = 0.5*yx + 0.5*chemNet.ymol[i];
    }
    //if( verbose & VERB7_MASK) fprintf(outputFd, "%03d,%d | ", niterations, count);
    if( verbose & VERB7_MASK) fprintf(outputFd, "%d| ", count);
    if( verbose & VERB7_MASK) if( (niterations % 30)==0 ) fprintf(outputFd, "\n%s                            ", INDENT7);

    if( niterations >= MAX_REFINE_ROOTS) { 
      fprintf(outputFd,"\n");
      errMsg(16); 
      CHECK_ERRV;
    }
  }
  
  //for (i=1; i <= chemNet.nSpec; i++) chemNet.ymol[i] /= pv.dens0;
  //writeArrayToBinaryFile("5.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
  //exit(-1);
  
  if( verbose & VERB7_MASK) fprintf(outputFd, "\n");
  if( verbose & VERB6_MASK) fprintf(outputFd, "%sNewton-Raphson Iterations until convergence : %03d \n", INDENT6, niterations);
}
void iterateChemicalNetworkOptimizedCorrect(const double AvL )
{
  int i, niterations, count;
  double yx, yxMin, yxMax;
  
  niterations = 0;

  if( verbose & VERB7_MASK) fprintf(outputFd, "%s chemical root refinement (optimized): ", INDENT7);

  /*Multi-Dim Newton Raphson, solving for equilibrium absolute abundances, ymol (equlibrium)   */
  while(1) { 
    
    niterations++;
    nNewton++;

    assemble_linear_system(chemNet.F, chemNet.DFDX);                         /* compute F and DFDX          */
    solve_lu(chemNet.DFDX, chemNet.F, chemNet.X, chemNet.nSpec); CHECK_ERRV; /* compute J^-1 V (see notes pp5)  */
    
    count=0;
    for (i = 1; i <= chemNet.nSpec; i++) {
      
      yx    = chemNet.ymol[i] - chemNet.X[i];  /* eq-1 next better guess for i^th specie density */
      yxMin = TINY2;                           /* minimum allowed value                          */
      yxMax = ics.maxAbun[i];                  /* maximum allowed value                          */
      
      if( yx < yxMin )
	yx = yxMin; 
      else
	if( yx > yxMax )
	  yx = yxMax;
	else
	  if ( dabs(yx / chemNet.ymol[i] - 1.0e0) > pv.CTol ) 
	    count++;
      
      chemNet.ymol[i] = yx;
    }
    
    //if( verbose & VERB7_MASK) fprintf(outputFd, "%03d,%d | ", niterations, count);
    if( verbose & VERB7_MASK) fprintf(outputFd, "%d| ", count);
    if( verbose & VERB7_MASK) if( (niterations % 30)==0 ) fprintf(outputFd, "\n%s                            ", INDENT7);
    
    if( count == 0 ) break;
 
    if( niterations >= MAX_REFINE_ROOTS) { 
      fprintf(outputFd,"\n");
      errMsg(16); 
      CHECK_ERRV;
    }

  }

  //for (i=1; i <= chemNet.nSpec; i++) chemNet.ymol[i] /= pv.dens0;
  //writeArrayToBinaryFile("5.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
  //exit(-1);
  
  if( verbose & VERB7_MASK) fprintf(outputFd, "\n");
  if( verbose & VERB6_MASK) fprintf(outputFd, "%sNewton-Raphson Iterations until convergence : %03d \n", INDENT6, niterations);
}
/*----------------------------------------------------------------------------------------------------*/
/* assembles the linear system which is solved to get the chemical equilibrium in other words,        */
/* assembles F[] and DFDX[][]                                                                         */
/*----------------------------------------------------------------------------------------------------*/
void assemble_linear_system(double *F, double **DFDX) 
{
  int i,j,k,l, **molno=chemNet.molno; 
  double vel, xmax; 

  /* Initialize F and DFDX to zero */
  for (i=0; i <= chemNet.nSpec; i++) { 
    F[i]=0.0; 
    for (j=0; j <= chemNet.nSpec; j++) 
      DFDX[i][j]=0.0;
  }
  
  /* get the rate for each specie for chemical reaction network */
  for (i=1; i <= chemNet.nReact; i++) { 
    
    /* Find the reaction rate from the reaction constant and the cenentrations */
    /* v = k [A][B][C]                                                         */
    vel=chemNet.rates[i];
    for(k=1; k <= 3; k++) if(molno[i][k] != 0) vel *= chemNet.ymol[molno[i][k]];
    
    /* Fill the DFDX */
    for(k=1; k <= 3; k++)
      if(molno[i][k] != 0) {
	for(l=1; l <= 3; l++) DFDX[molno[i][l]][molno[i][k]] -= (vel / chemNet.ymol[molno[i][k]]);
	for(l=4; l <= 7; l++) DFDX[molno[i][l]][molno[i][k]] += (vel / chemNet.ymol[molno[i][k]]);
      }
    
    /* collect the total destruction/production rates for the species in this reaction on the RHS */
    for (k=1; k <= 3; k++) if(molno[i][k] != 0) F[molno[i][k]] -= vel; 
    for (k=4; k <= 7; k++) if(molno[i][k] != 0) F[molno[i][k]] += vel;
  }
  
  /* Replace some equations with elemental equations, including the electron and grain */
  for (i=1; i <= chemNet.nBaseSpec; i++) {
    F[indBaseSpec[i]] = -ics.abun[i] * pv.dens0;  /*;;; what is being done here ? */
    for (k=1; k <= chemNet.nSpec; k++) {  
      F[indBaseSpec[i]] += specCode[k][i] * chemNet.ymol[k];
      DFDX[indBaseSpec[i]][k] = specCode[k][i];
    }
  }

  /* Scale the matrix */
  for (i = 1; i <= chemNet.nSpec; i++) {
    
    /* getting the max for each row in the linear system */
    xmax = dabs(F[i]);
    for (j = 1; j <= chemNet.nSpec; j++)
      if (DFDX[i][j] > xmax)
	xmax = dabs(DFDX[i][j]);
    xmax = max(xmax,1e-20);

    /* scaling */
    F[i] /= xmax;
    for (j = 1; j <= chemNet.nSpec; j++) DFDX[i][j] /= xmax;
  }
}
/*----------------------------------------------------------------------------------------------------*/
/* Computes the reaction rates for all the reactions in the chemical network                          */
/*----------------------------------------------------------------------------------------------------*/
void initRatesNetwork(const double AvL, const double T, const double Td, const double betaH2_L, 
		      const double betaCO_L, const double beta13CO_L)
{
  const double *ALPHA=chemNet.alpha, *BETA=chemNet.beta, *GAMMA=chemNet.gamma;
  double *rates=chemNet.rates;
  int i;

  /* initializing rates to zero */
  for (i=1; i <= chemNet.nReact; i++) rates[i] = 0.0;

  /* comute reaction rates : Reactions 1-3760 are two body reactions: (1) neutral-neutral reactions, (2) ion-neutral reactions,
                            (3) charge exchange reactions, (4) ion-ion reactions, (5) dissociative recombinations, 
			    (6) radiative recombinations, (7) associative detachments, (8) radiative associations */

  /* Two body reactions, eq-1 UMIST paper */
  for (i=1; i <= 3760; i++) rates[i]    = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* Reactions 3761-3916 are photoprocesses, eq-4*/
  for (i=3761; i <= 3916; i++) rates[i] = ALPHA[i] * (pv.G0/1.71) * exp(-GAMMA[i] * AvL);  ///TS
  /* Reactions 3917-3927 are cosmic-ray ionizations, eq-2 UMIST paper with CR scaling ( the pv.zeta term )*/
  for (i=3917; i <= 3927; i++) rates[i] = ALPHA[i] * pv.zeta / 1.3e-17;
  /* Reactions 3928-4059 are cosmic-ray induced photoreactions, eq-3 UMIST paper with CR scaling ( the pv.zeta term ) */
  for (i=3928; i <= 4059; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * GAMMA[i] / (1e0 - ALBEDO) * pv.zeta / 1.3e-17;
  /* Reactions 4060-4077 are collider reactions, eq-1 UMIST paper */
  for (i=4060; i <= 4077; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* Reactions 4078-4107 are thermonuclear reactions, eq-1 UMIST paper ( with extra pv.dens0 term ) */
  for (i=4078; i <= 4107; i++) rates[i] = pv.dens0 * ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* Reactions 4108-4113 are the sundries, eq-1 UMIST paper */
  for (i=4108; i <= 4113; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* H2 formation on dust */
  rates[4114] = ALPHA[4114] * pow(T/3.0e2, BETA[4114]) * H2_stickingFactor(T, Td) * pv.dens0 * H2_efficiency(T, Td) * pv.metalicity; 
  /* H2 photo-dissocociation */
  rates[4115] = ALPHA[4115] * pv.G0 * exp(-GAMMA[4115] * AvL); ///TS
  /* Reactions with PAHs */
  rates[4116] = ALPHA[4116] * pv.G0 * exp(-GAMMA[4116] * AvL); ///TS
  rates[4117] = ALPHA[4117] * pv.G0 * exp(-GAMMA[4117] * AvL); ///TS
  for (i=4118; i <= 4143; i++) rates[i] = ALPHA[i] * pow(T/1e2, BETA[i]) * PHI_PAH * GAMMA[i];
  
  /* Reactions with H2V, i.e. excited H2 (additional network, there are some more with 13C at the end) */
  /*--------------------------------------------------------------------------------------------------*/
  /* H2V neutral-neutral, ion-molecule etc. */
  for (i=4144; i <= 4250; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* H2V Cosmic-Ray reactions */
  for (i=4251; i <= 4254; i++) rates[i] = ALPHA[i] * pv.zeta / 1.3e-17;
  /* H2V Collider reactions */
  for (i=4255; i <= 4262; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* H2V misc reactions */
  rates[4263] = ALPHA[4263] * pow(T/3e2, BETA[4263]) * exp(-GAMMA[4263] / T);
  rates[4264] = ALPHA[4264] * pow(T    , BETA[4264]) * exp(-GAMMA[4264] / T);
  rates[4265] = ALPHA[4265] * pow(T    , BETA[4265]) * exp(-18100e0     / (T + 1200e0));
  rates[4266] = ALPHA[4266]; /* this would have a rate of zero always...check */
  rates[4267] = ALPHA[4267] * pv.G0 * exp(-GAMMA[4267] * AvL); ///TS
  rates[4268] = ALPHA[4268] * pv.G0 * exp(-GAMMA[4268] * AvL); ///TS
  /* done setting H2V reaction constants */
  /*-------------------------------------*/

  /* Take selfshlding into account for H2 and CO */
  rates[3811] *= betaCO_L; ///TS
  rates[4115] *= betaH2_L; ///TS
  rates[4268] *= betaH2_L; ///TS
  
  /* zeroing out reactions, i.e Find the equation with the right temperature range for reactions with multiple entries */
  /* ----------------------------------------------------------------------------------------------------------------- */
  if (T < 7e2)  rates[6] = 0e0;     /* eqn. 5,6 (checked) (done)*/
  if (T >= 7e2) rates[5] = 0e0;

  if (T < 2e3)  rates[104] = 0e0;   /* eqn. 103,104 (checked) (done)*/
  if (T >= 2e3) rates[103] = 0e0;

  if (T <  4.15e3) rates[241] = 0e0;   /* eqn. 240,241 (checked) */
  if (T >= 4.15e3) rates[240] = 0e0;

  if (T < 2e1)             rates[541] = rates[542] = 0e0;   /* eqn. 540-542 (checked) */
  if (T >= 2e1 && T < 3e2) rates[540] = rates[542] = 0e0;
  if (T >= 3e2)            rates[540] = rates[541] = 0e0;

  rates[799] = 0e0;   /* eqn. 798-799 (checked) removing this reaction since it trange overlaps with 798 (done) */

  if (T < 1e3)  rates[2610] = 0e0;   /* eqn. 2609,2610 (checked) */
  if (T >= 1e3) rates[2609] = 0e0;

  if (T < 1e4)  rates[2617] = 0e0;   /* eqn. 2616,2617 (checked) */
  if (T >= 1e4) rates[2616] = 0e0;

  if (T < 3.1e2)                rates[3189] = rates[3190] = 0e0;   /* eqn. 3188-3190 (checked) */
  if (T >= 3.1e2 && T < 1.24e3) rates[3188] = rates[3190] = 0e0;
  if (T >= 1.24e3)              rates[3188] = rates[3189] = 0e0;

  if (T < 3.1e2)                rates[3192] = rates[3193] = 0e0;     /* eqn. 3191-3193 (checked) */
  if (T >= 3.1e2 && T < 1.24e3) rates[3191] = rates[3193] = 0e0;
  if (T >= 1.24e3)              rates[3191] = rates[3192] = 0e0;

  if (T < 6.2e2)                rates[3205] = rates[3206] = 0e0;   /* eqn. 3204-3206 (checked) */
  if (T >= 6.2e2 && T < 4.46e3) rates[3204] = rates[3206] = 0e0;
  if (T >= 4.46e3)              rates[3204] = rates[3205] = 0e0;

  if (T < 6.2e2)                rates[3208] = rates[3209] = 0e0;   /* eqn. 3207-3209 (checked) */
  if (T >= 6.2e2 && T < 4.46e3) rates[3207] = rates[3209] = 0e0;
  if (T >= 4.46e3)              rates[3207] = rates[3208] = 0e0;

  if (T < 6.2e2)                rates[3211] = rates[3212] = 0e0;   /* eqn. 3210-3212 (checked) */
  if (T >= 6.2e2 && T < 4.46e3) rates[3210] = rates[3212] = 0e0;
  if (T >= 4.46e3)              rates[3210] = rates[3211] = 0e0;

  if (T < 2e3)                  rates[3608] = rates[3609] = 0e0;   /* eqn. 3607-3609 (checked) */
  if (T >= 2e3 && T < 6e3)      rates[3607] = rates[3609] = 0e0;
  if (T >= 6e3)                 rates[3607] = rates[3608] = 0e0;

  if (T < 7.95e3)                 rates[3614] = rates[3615] = 0e0;   /* eqn. 3613-3615 (checked) */
  if (T >= 7.95e3 && T < 2.114e4) rates[3613] = rates[3615] = 0e0;
  if (T >= 2.114e4)               rates[3613] = rates[3614] = 0e0;

  if (T < 1.5e4)    rates[3617] = 0e0;   /* eqn. 3616,3617 (chcked) */
  if (T >= 1.5e4)   rates[3616] = 0e0;

  if (T < 4e3)      rates[3680] = 0e0;   /* eqn. 3679-3681 (checked) */
  if (T >= 4e3)     rates[3679] = 0e0;

  if (T < 3e2)      rates[3707] = 0e0;   /* eqn. 3707,3708 (checked) (we have removed 3707 earlier, where are we keeping 3708??) */
  if (T >= 3e2)     rates[3708] = 0e0;

  if (T < 3e2)      rates[3710] = 0e0;   /* eqn. 3709,3710 (checked) */
  if (T >= 3e2)     rates[3709] = 0e0;

  if (T < 3e2)      rates[3732] = 0e0;   /* eqn. 3731,3732 (checked) */
  if (T >= 3e2)     rates[3731] = 0e0;

  if (T < 1e2)      rates[4066] = 0e0;    /* eqn. 4065,4066 (checked) */
  if (T >= 1e2)     rates[4065] = 0e0;

  if (T < 3e2)      rates[4082] = 0e0;    /* eqn. 4081,4082 (checked) (these are already removed when we remove the TR reactions*/
  if (T >= 3e2)     rates[4081] = 0e0;

  /* Watch out for negative gamma coefficients (enforcing a minimum temperature when computing the */
  /* reaction constants */
  rates[63]   = ALPHA[63]   * pow(max(T,1e1)   /3e2   , BETA[63])   * exp(-GAMMA[63]   / max(T, 1e1));
  rates[64]   = ALPHA[64]   * pow(max(T,1e1)   /3e2   , BETA[64])   * exp(-GAMMA[64]   / max(T, 1e1));
  rates[80]   = ALPHA[80]   * pow(max(T,3e2)   /3e2   , BETA[80])   * exp(-GAMMA[80]   / max(T, 3e2));
  rates[131]  = ALPHA[131]  * pow(max(T,1e2)   /3e2   , BETA[131])  * exp(-GAMMA[131]  / max(T, 1e2));
  rates[133]  = ALPHA[133]  * pow(max(T,3e2)   /3e2   , BETA[133])  * exp(-GAMMA[133]  / max(T, 3e2));
  rates[196]  = ALPHA[196]  * pow(max(T,3e2)   /3e2   , BETA[196])  * exp(-GAMMA[196]  / max(T, 3e2));
  rates[210]  = ALPHA[210]  * pow(max(T,2e2)   /3e2   , BETA[210])  * exp(-GAMMA[210]  / max(T, 2e2));
  rates[234]  = ALPHA[234]  * pow(max(T,2e2)   /3e2   , BETA[234])  * exp(-GAMMA[234]  / max(T, 2e2));
  rates[236]  = ALPHA[236]  * pow(max(T,2e2)   /3e2   , BETA[236])  * exp(-GAMMA[236]  / max(T, 2e2));
  rates[273]  = ALPHA[273]  * pow(max(T,2e2)   /3e2   , BETA[273])  * exp(-GAMMA[273]  / max(T, 2e2));
  rates[275]  = ALPHA[275]  * pow(max(T,3e2)   /3e2   , BETA[275])  * exp(-GAMMA[275]  / max(T, 3e2));
  rates[298]  = ALPHA[298]  * pow(max(T,2e2)   /3e2   , BETA[298])  * exp(-GAMMA[298]  / max(T, 2e2));
  rates[337]  = ALPHA[337]  * pow(max(T,2.5e2) /3e2   , BETA[337])  * exp(-GAMMA[337]  / max(T, 2.5e2));
  rates[343]  = ALPHA[343]  * pow(max(T,2e2)   /3e2   , BETA[343])  * exp(-GAMMA[343]  / max(T, 2e2));
  rates[348]  = ALPHA[348]  * pow(max(T,5e2)   /3e2   , BETA[348])  * exp(-GAMMA[348]  / max(T, 5e2));
  rates[354]  = ALPHA[354]  * pow(max(T,8e1)   /3e2   , BETA[354])  * exp(-GAMMA[354]  / max(T, 8e1));
  rates[360]  = ALPHA[360]  * pow(max(T,2e2)   /3e2   , BETA[360])  * exp(-GAMMA[360]  / max(T, 2e2));
  rates[363]  = ALPHA[363]  * pow(max(T,3e2)   /3e2   , BETA[363])  * exp(-GAMMA[363]  / max(T, 3e2));
  rates[366]  = ALPHA[366]  * pow(max(T,2e2)   /3e2   , BETA[366])  * exp(-GAMMA[366]  / max(T, 2e2));
  rates[374]  = ALPHA[374]  * pow(max(T,2.5e1) /3e2   , BETA[374])  * exp(-GAMMA[374]  / max(T, 2.5e1));
  rates[379]  = ALPHA[379]  * pow(max(T,1.5e2) /3e2   , BETA[379])  * exp(-GAMMA[379]  / max(T, 1.5e2));
  rates[387]  = ALPHA[387]  * pow(max(T,3e2)   /3e2   , BETA[387])  * exp(-GAMMA[387]  / max(T, 3e2));
  rates[390]  = ALPHA[390]  * pow(max(T,1.85e2)/3e2   , BETA[390])  * exp(-GAMMA[390]  / max(T, 1.85e2));
  rates[393]  = ALPHA[393]  * pow(max(T,1.3e1) /3e2   , BETA[393])  * exp(-GAMMA[393]  / max(T, 1.3e1));
  rates[400]  = ALPHA[400]  * pow(max(T,3e2)   /3e2   , BETA[400])  * exp(-GAMMA[400]  / max(T, 3e2));
  rates[414]  = ALPHA[414]  * pow(max(T,2e2)   /3e2   , BETA[414])  * exp(-GAMMA[414]  / max(T, 2e2));
  rates[425]  = ALPHA[425]  * pow(max(T,2e2)   /3e2   , BETA[425])  * exp(-GAMMA[425]  / max(T, 2e2));
  rates[432]  = ALPHA[432]  * pow(max(T,2e2)   /3e2   , BETA[432])  * exp(-GAMMA[432]  / max(T, 2e2));
  rates[3706] = ALPHA[3706] * pow(max(T,2e3)   /3e2   , BETA[3706]) * exp(-GAMMA[3706] / max(T, 2e3));
  rates[4087] = ALPHA[4087] * pow(max(T,2e2)   /3e2   , BETA[4087]) * exp(-GAMMA[4087] / max(T, 2e2));
  rates[4095] = ALPHA[4095] * pow(max(T,5e3)   /3e2   , BETA[4095]) * exp(-GAMMA[4095] / max(T, 5e3));
  rates[4097] = ALPHA[4097] * pow(max(T,2e3)   /3e2   , BETA[4097]) * exp(-GAMMA[4097] / max(T, 2e3));
  rates[4098] = ALPHA[4098] * pow(max(T,2e3)   /3e2   , BETA[4098]) * exp(-GAMMA[4098] / max(T, 2e3));
  rates[4099] = ALPHA[4099] * pow(max(T,2e3)   /3e2   , BETA[4099]) * exp(-GAMMA[4099] / max(T, 2e3));
  rates[4103] = ALPHA[4103] * pow(max(T,2e2)   /3e2   , BETA[4103]) * exp(-GAMMA[4103] / max(T, 2e2));
  rates[4106] = ALPHA[4106] * pow(max(T,2e2)   /3e2   , BETA[4106]) * exp(-GAMMA[4106] / max(T, 2e2));
  rates[4107] = ALPHA[4107] * pow(max(T,2e2)   /3e2   , BETA[4107]) * exp(-GAMMA[4107] / max(T, 2e2));

  /* disabling reactions of C + O+, C + O, CO + M respectively (done) */
  rates[3706] = rates[3707] = rates[4076] = 0e0;
  /* disabling reaction with type TR (done) */
  for (i=4078; i <= 4106; i++) rates[i] = 0e0;   

  /* Additional network for 13CO */
  /*-----------------------------*/
  /* 13C NB reactions */
  for (i=4269; i <= 4416; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  /* 13C photoreactions */
  for (i=4417; i <= 4426; i++) rates[i] = ALPHA[i] * (pv.G0/1.71) * exp(-GAMMA[i] * AvL); ///TS
  /* 13C reactions with CRP */
  for (i=4427; i <= 4428; i++) rates[i] = ALPHA[i] * pv.zeta / 1.3e-17;
  /* 13C reactions with CRPHOT */
  for (i=4429; i <= 4434; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * GAMMA[i] / (1e0 - ALBEDO) * pv.zeta / 1.3e-17;
  /* 13C NB reactions */
  for (i=4435; i <= 4438; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);    
  /* 13C NB reactions with PAHs*/
  for (i=4439; i <= 4440; i++) rates[i] = ALPHA[i] * pow(T/1e2, BETA[i]) * PHI_PAH * GAMMA[i]; // PAH reactions
  /* 13C NB reactions */
  for (i=4441; i <= 4453; i++) rates[i] = ALPHA[i] * pow(T/3e2, BETA[i]) * exp(-GAMMA[i] / T);
  
  rates[4426] *= beta13CO_L; ///TS

  rates[4278] = ALPHA[4278] * pow(max(T,3e2)/3e2, BETA[4278]) * exp(-GAMMA[4278] / max(T, 3e2)); //eq 80.

  if (T < 2e3)  rates[4280] = 0e0; /* 4280, 4281 checked */
  if (T >= 2e3) rates[4281] = 0e0;

  if (T < 7.95e3)                 rates[4404] = rates[4405] = 0e0; /* eqn. 4403-4405 (checked) */
  if (T >= 7.95e3 && T < 2.114e4) rates[4403] = rates[4405] = 0e0;
  if (T >= 2.114e4)               rates[4403] = rates[4404] = 0e0;

  rates[4412] = rates[4413] = 0e0; /* disabling the reactions of 13CO with O+ and O (T > 300K) respectivelt (done) */
}
/*----------------------------------------------------------------------------------------------------*/
/* Computes the timescales of the variation of all species                                            */
/*----------------------------------------------------------------------------------------------------*/
void chemSpecTimescale(chemicalNetwork *CN, int idx)
{
  int i, k, nShortTime;
  int nSpec=CN->nSpec, nReact=CN->nReact, **molno=CN->molno;
  double *F=CN->F, *rates=CN->rates, *ymol=CN->ymol;
  double vel, tau;
  char fName[512];
  FILE *fd;

  /* Initialize F[] to zero */
  for (i=0; i <= nSpec; i++) F[i]=0.0; 

  /* get the rate for each specie for chemical reaction network */
  for (i=1; i <= nReact; i++) { 
    
    /* Find the reaction rate from the reaction constant and the cenentrations */
    /* v = k [A][B][C]                                                         */
    vel=rates[i];
    for(k=1; k <= 3; k++) if(molno[i][k] != 0) vel *= ymol[molno[i][k]];
    
    /* collect the total destruction/production rates for the species in this reaction on the RHS */
    for (k=1; k <= 3; k++) if(molno[i][k] != 0) F[molno[i][k]] -= vel; 
    for (k=4; k <= 7; k++) if(molno[i][k] != 0) F[molno[i][k]] += vel;
  }

  /* Replace some equations with elemental equations, including the electron and grain */
  for (i=1; i <= CN->nBaseSpec; i++) {
    F[indBaseSpec[i]] = -ics.abun[i] * pv.dens0;  /*;;; what is being done here ? */
    for (k=1; k <= nSpec; k++)
      F[indBaseSpec[i]] += specCode[k][i] * ymol[k];
  }

  /*
  double *timeScale;
  int i;
  */

  //MAKE_1D_ARRAY(timeScale, chemNet.nSpec+1, double); /* local copies of the abudances to use as guesses */
  //FREE_1D_ARRAY(timeScale, chemNet.nSpec+1, double);
  //for(i=1; i<=chemNet.nSpec; i++)
  
  if(0) {
  sprintf(fName, "%schemTimeScale%05d.out",pv.outputDir, idx);
  if(  (fd=fopen(fName,"w"))==NULL ) {
      fprintf(stderr, "%s:%d : failed to open %s to dump the chemical timescale\n", __FILE__, __LINE__, fName);
      exit(-1);
  }

  nShortTime=0;
  for(i=1; i<nSpec; i++) {
    //tau = log10((ymol[i]/fabs(F[i]))/YEAR );
    tau = (ymol[i]/fabs(F[i]))/YEAR;
    fprintf(fd, "%03d %+e %+8.3f %+e\n",  i, ymol[i]/pv.dens0, log10(tau), tau);
    //    if( log10(tau) < 8 && fabs(ymol[i]/pv.dens0) > 1e-2 ) {
    if( log10(tau) < 6 && ymol[i]/pv.dens0 > 1e-3 ) {
	nShortTime++;
	//printf("xxx %03d %+e %+8.3f %+e\n",  i, ymol[i]/pv.dens0, log10(tau), tau);
    }
  }

  if(fd!=NULL)
      fclose(fd);
  }

  if( nShortTime != 0) {
    //printf("xxx n with short timescale = %d -------------------\n", nShortTime);
    //printf("xxx -----------------------------------------------\n");
  }
  
  
}
