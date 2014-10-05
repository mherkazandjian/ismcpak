#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

void allocMesh( mesh *m, int nSpecs )
{
    /* allocating mesh stuff */
    MAKE_1D_ARRAY(m->dustT, m->nx, double);
    MAKE_1D_ARRAY(m->gasT,  m->nx, double);

    MAKE_1D_ARRAY(m->xs, m->nx, double);
    MAKE_1D_ARRAY(m->xe, m->nx, double);
    MAKE_1D_ARRAY(m->xc, m->nx, double);
    MAKE_1D_ARRAY(m->dx, m->nx, double);

    MAKE_1D_ARRAY(m->AvL, m->nx, double);
    MAKE_1D_ARRAY(m->AvR, m->nx, double);

    MAKE_2D_ARRAY(m->abun, nSpecs, m->nx, double);  /* nSpec rows, nx colons */

    /* allocating bufferes to be used for storing stuff related to heating processes*/
    MAKE_1D_ARRAY(m->heating.photo   , m->nx, double);
    MAKE_1D_ARRAY(m->heating.cIon    , m->nx, double);
    MAKE_1D_ARRAY(m->heating.molHydro, m->nx, double);
    MAKE_1D_ARRAY(m->heating.H2pump  , m->nx, double);
    MAKE_1D_ARRAY(m->heating.ggColl  , m->nx, double);
    MAKE_1D_ARRAY(m->heating.visc    , m->nx, double);
    MAKE_1D_ARRAY(m->heating.CR      , m->nx, double);
    MAKE_1D_ARRAY(m->heating.total   , m->nx, double);

    /* allocating bufferes to be used for storing stuff related to cooling processes*/
    MAKE_1D_ARRAY(m->cooling.ms.total, m->nx, double); 
      MAKE_1D_ARRAY(m->cooling.ms.CI  , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.CII , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.FeI , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.FeII, m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.OI  , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.OII , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.SI  , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.SII , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.SiI , m->nx, double);
      MAKE_1D_ARRAY(m->cooling.ms.SiII, m->nx, double);
      MAKE_1D_ARRAY(m->cooling.fs.total,  m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.CP , 2, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.CP , 3, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.SiP, 2, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.SiP, 3, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.CA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.CA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.OA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.OA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.SA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.SA , 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.FeP, 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.FeP, 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.SiA, 4, m->nx, double);
      MAKE_2D_ARRAY(m->cooling.fs.popDens.SiA, 4, m->nx, double);
    MAKE_1D_ARRAY(m->cooling.roVib     , m->nx, double); 
    MAKE_1D_ARRAY(m->cooling.recom     , m->nx, double); 
    MAKE_1D_ARRAY(m->cooling.lymanAlpha, m->nx, double); 
    MAKE_1D_ARRAY(m->cooling.total     , m->nx, double); 
    
    /* allocating buffers which will store colon densities to the left and right */
    MAKE_1D_ARRAY(m->NH2All_L, m->nx, double);
    MAKE_1D_ARRAY(m->NH2_L   , m->nx, double);
    MAKE_1D_ARRAY(m->NCO_L   , m->nx, double);
    MAKE_1D_ARRAY(m->NC_L    , m->nx, double);
    MAKE_1D_ARRAY(m->NH2O_L  , m->nx, double);

    MAKE_1D_ARRAY(m->NH2All_R, m->nx, double);
    MAKE_1D_ARRAY(m->NH2_R   , m->nx, double); 
    MAKE_1D_ARRAY(m->NCO_R   , m->nx, double);
    MAKE_1D_ARRAY(m->NC_R    , m->nx, double);
    MAKE_1D_ARRAY(m->NH2O_R  , m->nx, double);

    /* allocate buffers to store self-sheilding due to colons on the left and right */
    MAKE_1D_ARRAY(m->betaH2_L,   m->nx, double); /* left */
    MAKE_1D_ARRAY(m->betaCO_L,   m->nx, double);
    MAKE_1D_ARRAY(m->beta13CO_L, m->nx, double);
    MAKE_1D_ARRAY(m->betaH2_R,   m->nx, double); /* right */
    MAKE_1D_ARRAY(m->betaCO_R,   m->nx, double);
    MAKE_1D_ARRAY(m->beta13CO_R, m->nx, double);
}
void freeMesh( mesh *m, int nSpecs )
{
    /* allocating mesh stuff */
    FREE_1D_ARRAY(m->dustT, m->nx, double);
    FREE_1D_ARRAY(m->gasT,  m->nx, double);

    FREE_1D_ARRAY(m->xs, m->nx, double);
    FREE_1D_ARRAY(m->xe, m->nx, double);
    FREE_1D_ARRAY(m->xc, m->nx, double);
    FREE_1D_ARRAY(m->dx, m->nx, double);

    FREE_1D_ARRAY(m->AvL, m->nx, double);
    FREE_1D_ARRAY(m->AvR, m->nx, double);

    FREE_2D_ARRAY(m->abun, nSpecs, m->nx, double);  

    /* allocating bufferes to be used for storing stuff related to heating processes*/
    FREE_1D_ARRAY(m->heating.photo   , m->nx, double);
    FREE_1D_ARRAY(m->heating.cIon    , m->nx, double);
    FREE_1D_ARRAY(m->heating.molHydro, m->nx, double);
    FREE_1D_ARRAY(m->heating.H2pump  , m->nx, double);
    FREE_1D_ARRAY(m->heating.ggColl  , m->nx, double);
    FREE_1D_ARRAY(m->heating.visc    , m->nx, double);
    FREE_1D_ARRAY(m->heating.CR      , m->nx, double);
    FREE_1D_ARRAY(m->heating.total   , m->nx, double);

    /* allocating bufferes to be used for storing stuff related to cooling processes*/
    FREE_1D_ARRAY(m->cooling.ms.total, m->nx, double); 
      FREE_1D_ARRAY(m->cooling.ms.CI  , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.CII , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.FeI , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.FeII, m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.OI  , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.OII , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.SI  , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.SII , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.SiI , m->nx, double);
      FREE_1D_ARRAY(m->cooling.ms.SiII, m->nx, double);
    FREE_1D_ARRAY(m->cooling.fs.total,  m->nx, double); 
      FREE_2D_ARRAY(m->cooling.fs.CP , 2, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.CP , 3, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.SiP, 2, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.SiP, 3, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.CA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.CA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.OA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.OA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.SA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.SA , 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.FeP, 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.FeP, 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.SiA, 4, m->nx, double);
      FREE_2D_ARRAY(m->cooling.fs.popDens.SiA, 4, m->nx, double);
    FREE_1D_ARRAY(m->cooling.roVib     , m->nx, double); 
    FREE_1D_ARRAY(m->cooling.recom     , m->nx, double); 
    FREE_1D_ARRAY(m->cooling.lymanAlpha, m->nx, double); 
    FREE_1D_ARRAY(m->cooling.total     , m->nx, double); 
    
    /* allocating buffers which will store colon densities to the left and right */
    FREE_1D_ARRAY(m->NH2All_L, m->nx, double);
    FREE_1D_ARRAY(m->NH2_L   , m->nx, double);
    FREE_1D_ARRAY(m->NCO_L   , m->nx, double);
    FREE_1D_ARRAY(m->NC_L    , m->nx, double);
    FREE_1D_ARRAY(m->NH2O_L  , m->nx, double);

    FREE_1D_ARRAY(m->NH2All_R, m->nx, double);
    FREE_1D_ARRAY(m->NH2_R   , m->nx, double); 
    FREE_1D_ARRAY(m->NCO_R   , m->nx, double);
    FREE_1D_ARRAY(m->NC_R    , m->nx, double);
    FREE_1D_ARRAY(m->NH2O_R  , m->nx, double);

    /* allocate buffers to store self-sheilding due to colons on the left and right */
    FREE_1D_ARRAY(m->betaH2_L,   m->nx, double); /* left */
    FREE_1D_ARRAY(m->betaCO_L,   m->nx, double);
    FREE_1D_ARRAY(m->beta13CO_L, m->nx, double);
    FREE_1D_ARRAY(m->betaH2_R,   m->nx, double); /* right */
    FREE_1D_ARRAY(m->betaCO_R,   m->nx, double);
    FREE_1D_ARRAY(m->beta13CO_R, m->nx, double);
    //fprintf(outputFd, "mesh freed\n");
}
/*----------------------------------------------------------------------------------------------------*/
double getDustTemp(double AvL) /* computes equilib dust temp due to incident FUV alone                */
{
  double factor1, factor2, tau100, T0; /* see thesis eq 5.28 and 5.29, 5.19 and hollenbach 1991 */

  T0     = 12.2 * pow(pv.G0, 0.2);                    /* see eq-5.28 rowin thesis pp 56    */
  tau100 = 2.7e2 * pv.G0 / pow(T0, 5e0) * pv.metalicity; /* see eq-5.29 rowin thesis pp 56 */
  factor1 = 8.9e-11 * NU0 * (TEFF/3e5) * pv.G0 * exp(-1.8 * (TEFF/3e5) * AvL);
  factor2 = 3.4e-2 * (0.42 - log(3.5e-2 * tau100 * T0)) * tau100 * pow(T0,6e0);
  return max(10e0, pow(( factor1 + pow(2.7e0,5e0) + factor2),0.2) );
}
/*----------------------------------------------------------------------------------------------------*/
/* Update column densities and self sheilding factors of the current slab                             */
/* The column densities are to computed by incrementing the colum densities to the left of the        */
/* previous slab by the amount of colum densities inside the previous slab                            */
/* Once the colum densities to the left fo the current slab 'k' are available, they are used to get   */
/* the self sheilding factors which are a function of colum densities to the left of the slab k       */
/*----------------------------------------------------------------------------------------------------*/
void updateColDensAndSelfSheildingDueToPreviousSlab(mesh *msh, const int k)
{
  const int km1=k-1;     /* index of the previous slab                                                 */
  double N_PreviousSlab; /* will hold the value of the column density of X specie in the presious slab */
  
  /* updating colon densities for some species ( colon density += abundance*slabThickness ) */
  N_PreviousSlab = (msh->abun[31-1][km1] + msh->abun[299-1][km1])*pv.dens0*msh->dx[km1]; /* 31 => nH2, 299 => nH2V */
  msh->NH2All_L[k] = msh->NH2All_L[km1] + N_PreviousSlab;
  
  N_PreviousSlab = msh->abun[31-1][km1]*pv.dens0*msh->dx[km1];  /* 31 => nH2 */
  msh->NH2_L[k]    = msh->NH2_L[km1] + N_PreviousSlab;
  
  N_PreviousSlab = msh->abun[46-1][km1]*pv.dens0*msh->dx[km1];  /* 46 => nCO */
  msh->NCO_L[k]    = msh->NCO_L[km1] + N_PreviousSlab;

  N_PreviousSlab = msh->abun[8 -1][km1]*pv.dens0*msh->dx[km1];  /* 8  => nC  */
  msh->NC_L [k]    = msh->NC_L [km1] + N_PreviousSlab;

  N_PreviousSlab = msh->abun[93-1][km1]*pv.dens0*msh->dx[km1];  /* 93 => nH2O */
  msh->NH2O_L[k]   = msh->NH2O_L[km1] + N_PreviousSlab;
  
  /* self sheilding factors due to all slabs before the current slab of index k1 */
  msh->betaH2_L[k]   = Self_shlding( msh->NH2All_L[k] );                                               /* self sheilding H2   */
  msh->betaCO_L[k]   = self_shlding_CO( max(1e0,msh->NCO_L[k]), max(1e0,msh->NH2All_L[k]), _12CO);      /* self sheilding CO   */
  msh->beta13CO_L[k] = self_shlding_CO( max(1e0,msh->NCO_L[k]), max(1e0,msh->NH2All_L[k]), _13CO); /* self sheilding 13CO */
}
/*----------------------------------------------------------------------------------------------------*/
/* compute the equilibrium get the equilibrium temperature due to coolong/heating and chemistry       */
/*----------------------------------------------------------------------------------------------------*/
int find_mesh_equilibrium(mesh *mesh)
{
  int k, j;
  double dAvTry, T_prev, T_this;
  const double min_dAv   = pv.min_deltaAv;
  const double max_dAv   = pv.max_deltaAv;
  const double relTolMax = pv.slabSizeCrit;
  const double AvMax     = pv.AvMax;
  
  /* equilibirum of the surface slab, k=0*/
  k = 0;
  getSystemEquilibrium(mesh, k); CHECK_ERR;  

#ifdef FIRST_SLAB_ONLY
  mesh->nxFilled = 1;
  printf("SKIPPING - MAX Av Reached\n");
  return 0;
#endif

  /* solve for all slabs, k >= 1, one slab at a time */
  for (k=1; k<mesh->nx; k++) {       

	  /* update col dens and self sheilding due to previous slab k-1 to */
	  /* be used in this slab k                                         */
	  updateColDensAndSelfSheildingDueToPreviousSlab(mesh, k);

	  T_prev = mesh->gasT[k-1];

	  /* setting the width of the next slab (non-adaptive) */
	  dAvTry = mesh->AvL[k] - mesh->AvL[k-1];

	  mesh->dx[k]  = Av2Position(dAvTry, pv.dens0, pv.metalicity);
	  mesh->xs[k]  = mesh->xs[k-1] + mesh->dx[k-1];
	  mesh->xc[k]  = mesh->xs[k]   + mesh->dx[k]/2.0;
	  mesh->xe[k]  = mesh->xs[k]   + mesh->dx[k];

	  /* equilibirum of the this slab */
	  getSystemEquilibrium(mesh, k);
	  if( errCode==17 ) {
		  ;
	  }
	  else
		  CHECK_ERR;

	  /* in case the temperature has been set manually bec no root was bracketted (to the minium allowed )*/
	  if( errCode==17 ) {
		  if( verbose & VERB4_MASK) fprintf(stdout, "%sused k=%03d Av=%f dAv=%f x=%e T=%f dT_Rel=%f %%\n",
				  INDENT4, k, mesh->AvL[k], dAvTry, mesh->dx[k], mesh->gasT[k],
				  100.0*(1.0-mesh->gasT[k]/T_prev));
		  break; /* while loop - adaptive temperature lookup */
	  }

	  T_this = mesh->gasT[k];

	  /* checking if no root has been found in the bracketting */
	  if( errCode==17 ) {
		  break; /* for loop, advancing throught the descreate slabs */
	  }

	  /* checking if the maximum allowed Av has been reached */
	  if(mesh->AvL[k] > AvMax) {
		  printf("SKIPPING - MAX Av Reached\n");
		  break; /* for loop, advancing throught the descreate slabs */
	  }

  } /* end for loop */
  
  mesh->nxFilled = k; /* i think this should be k+1 */

  /* manually setting the tempeature of the slabs */
  if( errCode == 17) {

	  for(k = mesh->nxFilled; k < mesh->nx; k++) {       /* solve for all slabs, one slab at a time                        */

		  T_prev = mesh->gasT[k-1];

		  //dAvTry = max_dAv;  #for a dynamic mesh
		  dAvTry = mesh->AvL[k] - mesh->AvL[k-1];  //for a pre-descretized mesh

		  mesh->dx[k]  = Av2Position(dAvTry, pv.dens0, pv.metalicity);
		  mesh->xs[k]  = mesh->xs[k-1] + mesh->dx[k-1];
		  mesh->xc[k]  = mesh->xs[k]   + mesh->dx[k]/2.0;
		  mesh->xe[k]  = mesh->xs[k]   + mesh->dx[k];
		  //mesh->AvL[k] = position2Av(mesh->xs[k], pv.dens0, pv.metalicity);  //use for an adaptive mesh

		  /* equilibirum of the this slab */
		  getSystemEquilibrium(mesh, k);
		  if( errCode==17 ) {
			  ;
		  }
		  else
			  CHECK_ERR;

		  T_this = mesh->gasT[k];

		  /* update col dens and self sheilding due to previous slab k-1 to */
		  /* be used in this slab k                                         */
		  updateColDensAndSelfSheildingDueToPreviousSlab(mesh, k);

		  /* checking if the maximum allowed Av has been reached */
		  if(mesh->AvL[k] > AvMax) {
			  printf("SKIPPING - MAX Av Reached\n");
			  break;
		  }
	  } /* end loop manually solving for fixed temeprature */

	  mesh->nxFilled = k; /* i think this should be k+1 */
  }
  
  return 0;
}
/*----------------------------------------------------------------------------------------------------*/
/* get the equilibrium temperature due to coolong/heating and chemistry                               */
/*----------------------------------------------------------------------------------------------------*/
int getSystemEquilibrium(mesh *msh, int k)
{
  double T, AvL;
  double Td, Teq=-1, Fa, Ta, Tb, Fb=666.66;
  int i, found=0;
  double Tmin=10.0, Tmax=1.0e5;
  int    maxIntervalLookups=20;
  double *buffer_a, *buffer_b;

  if( verbose & VERB1_MASK) fprintf(outputFd, "%s----------------------------------------\n", INDENT1);
  
  AvL = position2Av(msh->xc[k], pv.dens0, pv.metalicity);
  Td  = getDustTemp(AvL);  

  /* SETTING THE GUESSES FOR THE TEMPERATURE AND THE DENSITIES */
  if( verbose & VERB1_MASK) fprintf(outputFd, "%ssetting the guesses for the temperature and the densities\n", INDENT1);
  setEquilibriumGuesses( msh, k, &T); /* sets T and chemNet.ymol */

  /* dump the flux for the whole temperature range of interest for each slab */
  
  MAKE_1D_ARRAY(buffer_a, chemNet.nSpec+1, double); /* local copies of the abudances to use as guesses */
  MAKE_1D_ARRAY(buffer_b, chemNet.nSpec+1, double); /* for each side of the bracket                    */
  
  if( verbose & VERB1_MASK) fprintf(outputFd, "%sThermal/chemical balance for slab k=%04d Av=%f\n", INDENT1, k, AvL);

  /* BRACKETTING THE ROOT */
  found=0;

  if( verbose & VERB2_MASK) fprintf(outputFd, "%sLooking around the guess temperature of the previous slab\n", INDENT2);


  if( errCode == 17 ) { /* if the previous root bracketting failed, set an eq temp manually */
    if( verbose & VERB2_MASK) fprintf(outputFd, "%sWARNING: going to set the temperature manually.\n", INDENT2);
    goto GETMANUALSOLUTION;
  }

  for(i = 1; i < maxIntervalLookups; i++) { /* looking around the guess root for slabs other than the first slab*/
    
	  setRootBrackets(k, i, Tmin, Tmax, T, &Ta, &Tb);

	  if( verbose & VERB3_MASK) fprintf(outputFd, "%sChecking for roots in interval Ta=%f Tb=%f (lookup %d)\n", INDENT3, Ta, Tb, i);

	  if( i >= 2) memcpy( &(chemNet.ymol[0]), &(buffer_a[0]), (chemNet.nSpec+1)*sizeof(double) );  // ymol=buffer_a
	  chemical_balance(msh, Ta, k);                            CHECK_ERR;
	  Fa = thermal_balance(msh, Ta, Td, k, AvL, chemNet.ymol); CHECK_ERR;
	  memcpy( &(buffer_a[0]), &(chemNet.ymol[0]), (chemNet.nSpec+1)*sizeof(double) );              // buffer_a=ymol
	  if( i >= 2) if(Fa*Fb < 0 ) { found=1; break; }
	
	  if( i >= 2) memcpy( &(chemNet.ymol[0]), &(buffer_b[0]), (chemNet.nSpec+1)*sizeof(double) );  // ymol=buffer_b
	  chemical_balance(msh, Tb, k);                            CHECK_ERR;
	  Fb = thermal_balance(msh, Tb, Td, k, AvL, chemNet.ymol); CHECK_ERR;
	  memcpy( &(buffer_b[0]), &(chemNet.ymol[0]), (chemNet.nSpec+1)*sizeof(double) );              // buffer_b=ymol
	  if(Fa*Fb < 0 ) { found=1; break; }
  }
  
  /* doing a test for the root bracketting, if a root has been found or not */
  if( i >= maxIntervalLookups ) {
	  GETMANUALSOLUTION:
	  if( T < 50.0 ) {
		  setEquilibriumGuesses( msh, k, &T); /* sets T and chemNet.ymol */
		  Teq = msh->gasT[k-1];
		  chemical_balance(msh, Teq, k); CHECK_ERR;
		  Fa = thermal_balance(msh, Teq, Td, k, AvL, chemNet.ymol); CHECK_ERR;
		  errCode=17;
		  if( verbose & VERB2_MASK) fprintf(outputFd, "%sWARNING: set temperature manuall.\n", INDENT2);
		  goto  UPDATESLABDATA;
	  } else {
		  errMsg(17); CHECK_ERR;
	  }
  }
    
  if( found==0 ) errMsg(14); CHECK_ERR;
  
  /* CONVERGING TO THE ROOT */
  if( verbose & VERB2_MASK) fprintf(outputFd, "%sFount root in interval Ta=%f Tb=%f ( nLookups = %d) \n", INDENT2, Ta, Tb, i);

  if (T >= Tmin) {
    convergeToTemperatureRoot(msh, Ta, Tb, Fa, Fb, &Teq, AvL, chemNet.ymol, Td, k); CHECK_ERR;
  }
  else {
    Teq = 10e0;
    chemical_balance(msh, Teq, k); CHECK_ERR;
    Fa = thermal_balance(msh, Teq, Td, k, AvL, chemNet.ymol); CHECK_ERR;
  }

  if( verbose & VERB1_MASK) fprintf(outputFd, "%sTeq = %f\n", INDENT2, Teq);
  chemSpecTimescale(&chemNet, k);

  /* updating mesh cell data */    //if not root is found in the sweep, found may not be incremented and we would come here
  UPDATESLABDATA:
  msh->gasT[k]  = Teq;  
  msh->dustT[k] = Td;
  for (i = 1; i <= chemNet.nSpec; i++) msh->abun[i-1][k] = chemNet.ymol[i]/pv.dens0;

  return 0;
}
void convergeToTemperatureRoot(mesh *msh, double Ta, double Tb, double Fa, double Fb, double *Tx, double AvL, double *ymol, double Td, int k )
{
    double Fx=-666;
    int nIterations=-1, nBisectLocal=0, nSecantLocal=0;
    
    if( verbose & VERB3_MASK) fprintf(outputFd, "%sConverging to the root\n", INDENT3);
    
    if( (Tb-Ta)/(0.5*(Ta+Tb)) < pv.TTol ) {
	*Tx=0.5*(Ta+Tb);
	nIterations=0;
	if( verbose & VERB3_MASK) fprintf(outputFd, "%sBracketed root in accurate enough.\n", INDENT3);
    } else {
	
	while( 1 ) {
	 
	  if( verbose & VERB4_MASK) fprintf(outputFd, "%snSecant = %03d , nBisect = %03d , [Ta,Tb]=[%+.2e, %+.2e]\n", INDENT4,
					    nSecantLocal, nBisectLocal, Ta, Tb);
	  if( verbose & VERB4_MASK) fprintf(outputFd, "%s                                [Fa,Fb]=[%+.2e, %+.2e]\n", INDENT4,
					                                Fa, Fb);

	  /* setting the temperature for the next root lookup */
	  if( nSecantLocal > 10 ) {
	    *Tx=0.5*(Ta+Tb);               /* using bisection method */
	    nBisectLocal++;
	  }
	  else {
	    *Tx=Tb - Fb*(Tb-Ta)/(Fb-Fa);   /* using secant method    */
	    nSecantLocal++;
	  }
	  
	  chemical_balance(msh, *Tx, k);
	  Fx=thermal_balance(msh, *Tx, Td, k, AvL, chemNet.ymol); CHECK_ERRV;
	  
	  if(Fx*Fa<0.0) { Tb=*Tx; Fb=Fx; } else { Ta=*Tx; Fa=Fx; }

	  /* convergence criterie */
	     /* heating + cooling == 0 */
	     if( Fa == 0 ) { *Tx = Ta; Tb=Ta; Fb=Fa; Fx=0.0; break; }  
	     if( Fb == 0 ) { *Tx = Tb; Ta=Tb; Fa=Fb; Fx=0.0; break; }
	     /* acheived predefined tolerence for the temperatures */
	     if( (Tb-Ta)/(0.5*(Ta+Tb)) < pv.TTol  ) break;
	}
	
	nIterations=nBisectLocal+nSecantLocal;

	if( verbose & VERB3_MASK) {
	  fprintf(outputFd, "%sn calls to thermal root convergence          = %d\n", INDENT3, nIterations);
	  fprintf(outputFd, "%sTa = %f  Tb = %f    relTol = %e\n", INDENT3, Ta, Tb, (Tb-Ta)/(0.5*(Ta+Tb)) );
	  fprintf(outputFd, "%sFa = %e Fb = %e\n", INDENT3, Fa, Fb);
	  fprintf(outputFd, "%sthermal timescale = %e Gyr\n", INDENT3, (((3.0/2.0)*pv.dens0*KB*(*Tx))/fabs(Fx))/(YEAR*1e9));
	  
	}
    }

    if( nIterations==-1 ) errMsg(13); /* check values for Fx and/or nIterations for debugging */
}
/*----------------------------------------------------------------------------------------------------*/
/* dump output to file in binary format                                                               */
/*----------------------------------------------------------------------------------------------------*/
void writeDataOutputBinary( mesh *m, int nSpec, char *outFile )
{
  FILE *fd;
  int i;
  /* nw : is the number of slabs to include in the binary file */
  
  if( (fd = fopen(outFile, "w"))==NULL) {
    fprintf(stderr, "%s:%d: FATAL ERROR : failed to open file : %s.\n",
	    __FILE__, __LINE__, outFile);
    exit(-1);
  }

  fwrite(&dataVersion    , sizeof(int)   , 1           , fd);
  fwrite(&pv.G0          , sizeof(double), 1           , fd);
  fwrite(&pv.dens0       , sizeof(double), 1           , fd);
  fwrite(&pv.gamma_mech  , sizeof(double), 1           , fd);
  fwrite(&m->nxFilled    , sizeof(int)   , 1           , fd);
  fwrite(&nSpec          , sizeof(int)   , 1           , fd);
  fwrite(&(m->gasT[0])   , sizeof(double), m->nxFilled      , fd);
  fwrite(&(m->dustT[0])  , sizeof(double), m->nxFilled      , fd);
  fwrite(&(m->AvL[0])    , sizeof(double), m->nxFilled      , fd);
  
  for(i=0;i<nSpec;i++)
    fwrite(&(m->abun[i][0]), sizeof(double), m->nxFilled, fd);

  fwrite(&(m->heating.total[0]), sizeof(double), m->nxFilled      , fd);
  fwrite(&(m->cooling.total[0]), sizeof(double), m->nxFilled      , fd);

  fwrite(&(m->heating.photo[0])   , sizeof(double), m->nxFilled      , fd);  /* 1 */
  fwrite(&(m->heating.cIon[0])    , sizeof(double), m->nxFilled      , fd);  /* 2 */
  fwrite(&(m->heating.molHydro[0]), sizeof(double), m->nxFilled      , fd);  /* 3 */
  fwrite(&(m->heating.H2pump[0])  , sizeof(double), m->nxFilled      , fd);  /* 4 */
  fwrite(&(m->heating.ggColl[0])  , sizeof(double), m->nxFilled      , fd);  /* 5 */
  fwrite(&(m->heating.visc[0])    , sizeof(double), m->nxFilled      , fd);  /* 6 */
  fwrite(&(m->heating.CR[0])      , sizeof(double), m->nxFilled      , fd);  /* 7 */

  fwrite(&(m->cooling.ms.total[0])      , sizeof(double), m->nxFilled, fd); /* 1 */
  fwrite(&(m->cooling.fs.total[0])      , sizeof(double), m->nxFilled, fd); /* 2 */
  fwrite(&(m->cooling.roVib[0])         , sizeof(double), m->nxFilled, fd); /* 3 */
  fwrite(&(m->cooling.recom[0])         , sizeof(double), m->nxFilled, fd); /* 4 */
  fwrite(&(m->cooling.lymanAlpha[0])    , sizeof(double), m->nxFilled, fd); /* 5 */

  /* writing the components of the ms cooling */
  fwrite(&(m->cooling.ms.CI[0])  , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.CII[0]) , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.FeI[0]) , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.FeII[0]), sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.OI[0])  , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.OII[0]) , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.SI[0])  , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.SII[0]) , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.SiI[0]) , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->cooling.ms.SiII[0]), sizeof(double), m->nxFilled, fd);

  /* writing the details of the fs cooling, rates and pop densities */
  for(i=1; i<=1; i++) fwrite(&(m->cooling.fs.CP[i][0])         , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=2; i++) fwrite(&(m->cooling.fs.popDens.CP[i][0]) , sizeof(double), m->nxFilled, fd);

  for(i=1; i<=1; i++) fwrite(&(m->cooling.fs.SiP[i][0])         , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=2; i++) fwrite(&(m->cooling.fs.popDens.SiP[i][0]) , sizeof(double), m->nxFilled, fd);

  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.CA[i][0])        , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.popDens.CA[i][0]), sizeof(double), m->nxFilled, fd);

  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.OA[i][0])        , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.popDens.OA[i][0]), sizeof(double), m->nxFilled, fd);

  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.SA[i][0])        , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.popDens.SA[i][0]), sizeof(double), m->nxFilled, fd);

  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.FeP[i][0])        , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.popDens.FeP[i][0]), sizeof(double), m->nxFilled, fd);

  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.SiA[i][0])        , sizeof(double), m->nxFilled, fd);
  for(i=1; i<=3; i++) fwrite(&(m->cooling.fs.popDens.SiA[i][0]), sizeof(double), m->nxFilled, fd);

  /* writing the self sheilding stuff */
  fwrite(&(m->betaH2_L[0])  , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->betaCO_L[0])  , sizeof(double), m->nxFilled, fd);
  fwrite(&(m->beta13CO_L[0]), sizeof(double), m->nxFilled, fd);

  fclose(fd);
}
/*----------------------------------------------------------------------------------------------------*/
void errMsg(const short int x)
{
  int fatal=-1;

  fflush(stdout);
  
  fprintf(stderr, "error code %d : ", x);
  
  switch(x)  
    {
    case 1  : { fprintf(stderr,"in function main() : error reading/handling parameter file");      fatal=1; } break;
    case 4  : { fprintf(stderr,"in function readDataHeaderSpecies() : failed to open file");       fatal=1; } break;
    case 5  : { fprintf(stderr,"in function readDataHeaderUnderUbundant() : failed to open file"); fatal=1; } break;
    case 6  : { fprintf(stderr,"in function initReactionRateBuffers() : failed to open file");     fatal=1; } break;
    case 7  : { fprintf(stderr,"in function read_selfshlding_CO() : failed to open file");         fatal=1; } break;
    case 8  : { fprintf(stderr,"in function read_rotational_data() : failed to open file");        fatal=1; } break;
    case 9  : { fprintf(stderr,"in function read_rotational_data() : failed to open file");        fatal=1; } break;
    case 10 : { fprintf(stderr,"in function read_vibrational_data() : failed to open file");       fatal=1; } break;
    case 12 : { fprintf(stderr,"in function set_diagnostics_parameters() failed to open file");    fatal=1; } break;

    case 2  : { fprintf(stderr,"in function maxRow() : code must not get here.");                  fatal=1; } break;
    case 3  : { fprintf(stderr,"in function maxRow() : singular matrix - the whole row is zeros"); fatal=0; } break;
    case 11 : { fprintf(stderr,"in function thermal_balance() : either heating or cooling < 0");   fatal=0; } break;
    case 13 : { fprintf(stderr,"execution should not get here");                                   fatal=0; } break;
    case 14 : { fprintf(stderr,"in function getSystemEquilibrium() no root bracketted");           fatal=0; } break;
    case 15 : { fprintf(stderr,"in MACRO MAKE_xD_ARRAY : error allocating memory");                fatal=1; } break;
    case 16 : { fprintf(stderr,"in function iterateChemicalNetwork : maximum number of ");         
                fprintf(stderr,"iterations in refining checmical balance excceeded.\n");           fatal=0; } break;
    case 17 : { fprintf(stderr,"in function getSystemEquilibrium() : maximum number of ");         
                fprintf(stderr,"iterations in expanding bracketing interval excceeded.\n");        fatal=0; } break;
      
    default : { fprintf(stderr,"ERROR\n");                                                         fatal=1; } break;
    }

  nErrors++;
  errCode=x;
  
  fflush(stdout);
  fflush(stderr);

  if( (verbose & ERRHND_MASK) && fatal ) {
    fprintf(stderr, " This is a FATAL ERROR. Exitting program\n");
    exit(-1);
  } else {

  }
}
/*----------------------------------------------------------------------------------------------------*/
/* replaces a specific character in a string by another charachter                                    */
void replaceChar( char *str, char a, char b, int n )
{
    int i;
    for(i=0; i<n ;i++) 
	if(str[i]==a)
	    str[i]=b;
}
/*----------------------------------------------------------------------------------------------------*/
/* H2 formation efficiency and sticking factor ( eq-5.50 rowin thesis )                               */
/*----------------------------------------------------------------------------------------------------*/
double H2_stickingFactor(const double T, const double Td)
{
  return pow((1.0e0 + 0.4e0 * sqrt(T + Td) + 2.0e-3 * T + 8.0e-6 * SQ(T)),-1.0e0);
}
double H2_efficiency(const double T, const double Td)
{
  const double mu=5.0e-3, FF=1.0e-15,  Eh2=3.2e2, Ehp=6.0e2, Ehc=1.0e4, Es=2.0e2, nuH2=3.0e12, nuHc=1.3e13;
  double bbH2, betal, ksi;

  bbH2 = nuH2 * exp(-Eh2 / Td);
  betal = 2.5e-1 * SQ(1.0e0 + sqrt((Ehc - Es)/(Ehp - Es))) * exp(-Es / Td);
  ksi = pow((1.0e0 + nuHc / (2.0e0 * FF) * exp(-1.5e0 * Ehc / Td) *    
	     SQ(1.0e0 + sqrt((Ehc - Es) / (Ehp - Es)))),-1.0e0);	  	      

  return  pow((1.0e0 + mu * FF /(2.0e0 * bbH2) + betal),-1.0e0) * ksi;   
}
/*----------------------------------------------------------------------------------------------------*/
/* set diagnostics output file descriptor and verbosity level flag                                    */
/*----------------------------------------------------------------------------------------------------*/
void set_diagnostics_parameters( char *outFname )
{
  if( !(verbose & OUTPUT_MASK) ) /* make sure if first bit is not set, clear all except the forth */
    verbose = verbose & ERRHND_MASK;
  
  if( !(verbose & OUTFILE_MASK) ) { 

    outputFd=stdout;
    if( verbose & VERB1_MASK) fprintf(outputFd, "%sDumping to stdout\n", INDENT1);
    
  } else {
    
    if( (outputFd=fopen(outFname, "w"))==NULL ) 
      errMsg(12);
    else {
      if( verbose & VERB1_MASK) fprintf(outputFd, "%sDumping to file %s\n", INDENT1, outFname);
    }
    
  }
  
  if( verbose & STDERR_MASK ) /* redirecting stderr to stdout */
    stderr=outputFd;
}
/*----------------------------------------------------------------------------------------------------*/
/* finalizing the run                                                                                 */
/*----------------------------------------------------------------------------------------------------*/
void finalize( mesh *msh)
{
  TOC(tmr1);
  if( verbose & VERB1_MASK) fprintf(outputFd, "%stotal time runnning model        = %f \n", INDENT1, tmr1.dt);
  if( verbose & VERB1_MASK) fprintf(outputFd, "%stotal time LU decomposing        = %f \n", INDENT1, tmr2.tot);
  //if( verbose & VERB1_MASK) fprintf(outputFd, "%saverage time per LU decompostion = %f \n", INDENT1, tmr2.tot/((double)nNewton) );
  if( verbose & VERB1_MASK) fprintf(outputFd, "%sTotal number of Newton-raphson for chemical balance  = %d \n", INDENT1, nNewton);
  if( verbose & VERB1_MASK) fprintf(outputFd, "%sTotal number of Bisection for thermal balance        = %d \n", INDENT1, nBisect);
  //if( verbose & VERB1_MASK) fprintf(outputFd, "%sAverage number of Newton Raphson per Termal Balance  = %f \n", INDENT1, (double)nNewton/(double)nBisect );

  fflush(outputFd);
  //fclose(outputFd);
}
/*----------------------------------------------------------------------------------------------------*/
/* convert from Av to position inside the cloud and vise versa (assuming all of the gas is in H       */
/* H2 (the errors due to this assumption are quite small  ~ 1e-3 with the tests i did so far          */
/*----------------------------------------------------------------------------------------------------*/
double Av2Position(const double Av, const double dens, const double Z)
{
    return Av/(5.34e-22*Z*dens);
}
double position2Av(const double x, const double dens, const double Z)
{
    const double NH=x*dens; /* H column density */
    return 5.34e-22*NH*Z;
}
/*----------------------------------------------------------------------------------------------------*/
/*-----------------------MALLOC 2D ARRAY--------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
double ** malloc2dArray_d(int l, int m)
{
    double *data, **x2d;
    int i;
    
    data = (double *  )malloc(sizeof(double   )*l*m); 
    x2d  = (double ** )malloc(sizeof(double * )*l  );
    
    if( data==NULL || x2d==NULL )
	exit(-1);
    
    for(i=0; i<l; i++) x2d[i] = &(data[i*m]);
  
    for(i=0; i<l*m; i++)
      data[i]=(double)i;
    
    return x2d;
}
int free2dArray_d(double **x2d)
{
  free(&x2d[0][0]);
  free(&x2d[0]);

  return 0;
}
/*----------------------------------------------------------------------------------------------------*/
/*-----------------------MALLOC 3D ARRAY--------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------*/
double *** malloc3dArray_d(int l, int m, int n)
{
    double *data, **zColPtr, ***x3d;
    int i, j;
    
    data    = (double *  )malloc(sizeof(double   )*l*m*n); 
    zColPtr = (double ** )malloc(sizeof(double * )*l*m  );
    x3d     = (double ***)malloc(sizeof(double **)*l    );
    
    if( data==NULL || zColPtr==NULL || x3d==NULL )
	exit(-1);
    
    for(i=0; i<l; i++) x3d[i] = &(zColPtr[i*m]);
    for(i=0; i<l; i++) 
	for(j=0; j<m; j++)
	  x3d[i][j] = &(data[ (i*m + j)*n ]);
    
    return x3d;
}
int free3dArray_d(double ***x3d)
{
  free(&x3d[0][0][0]);
  free(&x3d[0][0]);
  free(&x3d[0]);

  return 0;
}
int writeArrayToBinaryFile(const char *fName, const double *vecPtr, const int n)
{
  FILE *fd;

  if( (fd=fopen(fName, "w"))!=NULL ) {
    return fwrite((void *)vecPtr, sizeof(double), n, fd);
  } else {
    fprintf(stderr, "fialed to open file %s.\n", fName);
    exit(-1);
  }
}
int readArrayFromBinaryFile(const char *fName, double *vecPtr, const int n)
{
  FILE *fd;

  if( (fd=fopen(fName, "r"))!=NULL ) {
    return fread((void *)vecPtr, sizeof(double), n, fd);
  } else {
    fprintf(stderr, "fialed to open file %s.\n", fName);
    exit(-1);
  }
}
void dumpFluxVsTSweep(mesh *msh, const int k, const double T, const double Td, const double AvL)
{
  FILE *fdFlux;
  double Ts, Fs;
  char fname[1024];
  double *originalGuessAbundances;
  
  if( k % 2 == 0 ) {
    MAKE_1D_ARRAY(originalGuessAbundances, chemNet.nSpec+1, double);
    memcpy( &(originalGuessAbundances[0]), &(chemNet.ymol[0]), (chemNet.nSpec+1)*sizeof(double) );
    
    sprintf(fname, "%sflux-%02d.out", pv.outputDir, k);
    if( (fdFlux=fopen(fname,"w"))==NULL  ) {
      fprintf(stderr, "%s:%d: failed to open file : %s.\n",
	      __FILE__, __LINE__, fname);
      exit(-1);
    }
    
    /* sweeping from guess T to T_min */
    Ts=T;
    memcpy( &(chemNet.ymol[0]), &(originalGuessAbundances[0]), (chemNet.nSpec+1)*sizeof(double) );
    while(Ts > 10.0 ) {
      chemical_balance(msh, Ts, k);
      Fs = thermal_balance(msh, Ts, Td, k, AvL, chemNet.ymol); CHECK_ERRV;
      fprintf(fdFlux, "%f %f %e\n", AvL, Ts, Fs);
      fprintf(stdout, "%d %f %e\n", k, Ts, Fs);
      if( Ts > 10000.0 && Ts < 100000.0 ) Ts -= 1000.0;
      if( Ts > 1000.0  && Ts < 10000.0  ) Ts -= 100.0;
      if( Ts > 100.0   && Ts < 1000.0   ) Ts -= 10.0;
      if( Ts > 10.0    && Ts < 100.0    ) Ts -= 1.0;
    }
     
    /* sweeping from guess T to T_max */
    Ts=T;
    memcpy( &(chemNet.ymol[0]), &(originalGuessAbundances[0]), (chemNet.nSpec+1)*sizeof(double) );
    while(Ts < 1.0e5 ) {
      chemical_balance(msh, Ts, k);
      Fs = thermal_balance(msh, Ts, Td, k, AvL, chemNet.ymol); CHECK_ERRV;
      fprintf(fdFlux, "%f %f %e\n", AvL, Ts, Fs);
      fprintf(stdout, "%d %f %e\n", k, Ts, Fs);
      if( Ts > 10000.0 && Ts < 100000.0 ) Ts += 1000.0;
      if( Ts > 1000.0  && Ts < 10000.0  ) Ts += 100.0;
      if( Ts > 100.0   && Ts < 1000.0   ) Ts += 10.0;
      if( Ts > 10.0    && Ts < 100.0    ) Ts += 1.0;
    }
    fclose(fdFlux);
    
    memcpy( &(chemNet.ymol[0]), &(originalGuessAbundances[0]), (chemNet.nSpec+1)*sizeof(double) );
  }
}
/* sets the new guesses for where the root might be bracketted */
int setRootBrackets(const int k, const int i, const double Tmin, const double Tmax, const double T, double *Ta, double *Tb)
{
  double dT=1.0;

  dT *= pow(2.0,i);
  *Ta = T - dT;
  *Tb = T + dT;
  
  *Ta = max(Tmin, T - dT); 
  *Tb = min(Tmax, T + dT); 
}
/* sets the temperature to the mesh struct for slab of index k */
void setSlabTemperatureMesh( const int k, const double T, mesh *msh)
{
  msh->gasT[k] = T;
}
/* sets the abundances of the mesh struct for slab of index k     */
/* The abundances are normalized to the gas number density which  */
/* (the gas number density is in terms of the total number of H   */
/* nuclei                                                         */ 
void setSlabAbundancesMesh( const int k, const double *abun, const double nDens, mesh *msh)
{
  int i;
  for (i=1; i <= chemNet.nSpec; i++) msh->abun[i-1][k] = abun[i]/nDens;
}
/*----------------------------------------------------------------------------------*/
/* compute the locations of the end of the descretized sub slabs of the 1D semi     */
/* infinite PDR slab in this routine, as long as Av is less than 1, the slab width  */
/* is mutliplied  by 10 evey dex, and above that the widths are set in a pre-defined*/
/* fashion for ex, if dAvMin = 0.001  (the zeroth slab is excluded since it         */
/* is zero                                                                          */
/*  Av = [ 0.001, 0.002,...0.01, 0.02, 0.03, .., 0.1, 0.2, .. 1, 1.2, 1.4..         */
/*         1.6..10, 10.5, 11, 11.5...., 30                                          */
/* it returns the number of slabs, the computed location of the slab endoints       */
/* are stored into AvTmp, it has to be large enough to hold all the slabs           */
/* see genPreDefined1DSlab.py                                                       */
/*--------------------------------------------------------------------------------- */
long descretizeSlab1( double *AvTmp )
{
	double epsDsc = 1e-10;
	int maxSlabs = 5000;
	int i;
	double dAv, Av;

	double dAv1 = 0.2;  /* delta Av when   1  < Av < 10 */
	double dAv2 = 1.0;  /* delta Av when   10 < Av < 30 */

	Av = 0.0;

	i = 0;  dAv = pv.min_deltaAv; Av += dAv;
	AvTmp[i] = Av;

	/* setting the thicknesses for Av < 1.0 */
	while( Av < 1.0 - epsDsc )
	{
		i += 1; Av += dAv; AvTmp[i] = Av;
		if ( fmod( fabs(log10(Av) ), 1.0  ) <= 1e-14 ) dAv *= 10.0;
	}

	/* setting the thicknesses for 1.0 < Av < 10.0 */
	while( Av < 10.0 - epsDsc )
	{
		i += 1; Av += dAv1; AvTmp[i] = Av;
	}

	/* setting the thicknesses for 10.0 < Av < AvMax */
	while( Av < pv.AvMax - epsDsc )
	{
		i += 1; Av += dAv2; AvTmp[i] = Av;
	}

	return i+1;
}
