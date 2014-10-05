#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

void init(int argc, char **argv)
{
  char outFname[256];
  sprintf(outFname, "%smesh.out", pv.outputDir);
  set_diagnostics_parameters( outFname );
  
  /* intialize root finding counters and timers */
  nNewton=nBisect=0;
  tmr1.dt=tmr2.dt=tmr3.dt=tmr4.dt=tmr1.tot=tmr2.tot=tmr3.tot=tmr4.tot=0.0;
  TIC(tmr1);
    
  /* assigning the parameters to independent variables */
  pv.dens0       = pow(10.0,pv.dens0);
  pv.G0          = pow(10.0,pv.G0);
  
  T0     = 12.2 * pow(pv.G0, 0.2);                    /* see eq-5.28 rowin thesis pp 56 */
  tau100 = 2.7e2 * pv.G0 / pow(T0, 5e0) * pv.metalicity; /* see eq-5.29 rowin thesis pp 56 ;;; veranderen bij metalliciteit */
  
  chemical_network_data();   /* read and alloc all info about the chemical network         */
  read_selfshlding_CO();     /* read and alloc stuff needed for sheilding                  */
  read_vibrational_data();   /* read data of vibrational modes                             */
  read_rotational_data();    /* read data of rotational modes                              */
  dumpVersionFile();         /* copying the version file                                   */

  init_database();
}
/*----------------------------------------------------------------------------------------------------*/
void initConstantsAndAbundances( )
{
  int i;
  
  /* setting up the intial guesses for the abundances of base species */
  for (i=5; i <= chemNet.nBaseSpec; i++) ics.abun[i] *= pv.metalicity;
  /* specific changes in abundance */
  if( pv.metalicity >= 2.0) {
    printf("WARNING : Metalicit = %f\n", pv.metalicity);
    printf("WARNING : Manually setting atomic carbon abundance equal to atomic Oxygen abundance.\n");
    ics.abun[4]  = ics.abun[6];                   // C = O
  }
  ics.abun[15] = ics.abun[4] / 40e0;              // 13C = 12C/40
  ics.abun[11] = ics.abun[11] / pv.S_depletion;   // deplete S
}
/*----------------------------------------------------------------------------------------------------*/
void setMaxAllowedDensities()
{
    int i;

    /* setting the maximum allowed absolute abundances */
    for (i=1; i <= chemNet.nSpec; i++) {
      if (i <= 163            ) ics.maxAbun[i] = pv.dens0; 
      if (i > 163  && i < 296 ) ics.maxAbun[i] = pv.dens0 * 1.00e-8;
      if (i >= 296 && i <= 298) ics.maxAbun[i] = pv.dens0 * 3.72e-7;
      if (i >= 300 && i <= 309) ics.maxAbun[i] = pv.dens0 / 40e0;
    }
    ics.maxAbun[299] = pv.dens0 * 1e-3;
}
/*----------------------------------------------------------------------------------------------------*/
void dumpVersionFile( void )
{
  char command[512];
  sprintf(command, "cp VERSION %sVERSION", pv.outputDir);
  system( command );
}
/*----------------------------------------------------------------------------------------------------*/
void setupAllocMesh( mesh *m, int nSlabs, int nSpecs )
{
  m->nx = nSlabs;
  allocMesh(m, nSpecs);
}
/*----------------------------------------------------------------------------------------------------*/
void initMesh( mesh *m )
{

  double *AvTmp; // temporary array holding the ending positions in Av of all the slabs
  AvTmp = (double *)malloc(sizeof(double)*5000);

  m->nx = descretizeSlab1( AvTmp );
  m->nx += 1; // taking the zeroth slab into account

  setupAllocMesh( m, m->nx, chemNet.nSpec); /* set the size and allocate the mesh   */

  m->AvL[0]      = 0.0; /* setting the position of the first slab to 0               */
  m->xs[0]       = 0.0;
  m->xc[0]       = 0.0;
  m->xe[0]       = 0.0;
  m->dx[0]       = 0.0;

  m->NH2All_L[0] = 0.0; /* setting the initial colon densities due to the slabs on the  */
  m->NH2_L[0]    = 0.0; /* left to zero                                                 */
  m->NCO_L[0]    = 0.0;      
  m->NC_L[0]     = 0.0;
  m->NH2O_L[0]   = 0.0;

  /* set the self sheilding factor of the first slabs on the right and left */
  m->betaH2_L[0]= m->betaCO_L[0] = m->beta13CO_L[0] = 1.0;
  m->betaH2_R[0]= m->betaCO_R[0] = m->beta13CO_R[0] = 1.0;

  for(int i = 0; i < chemNet.nSpec; i++)
	  for( int j = 0; j < m->nx; j++)
		  m->abun[i][j] = 0.0;

  /* setting the precomputed locations of Av from the temporary one */
  for( int i = 1; i < m->nx; i++)
	  m->AvL[i] = AvTmp[i-1];
  //for( int i = 0; i < m->nx; i++)
  //printf("i = %d m->AvL[i] = %f, AvTmp[i-1] = %f\n", i, m->AvL[i], AvTmp[i-1]);
  //printf(" total number of slabs = %d\n", m->nx);
  free(AvTmp);

  if( verbose & VERB1_MASK) fprintf(outputFd, "%s%s",INDENT1, "setup and initialized the mesh\n");
}
/*----------------------------------------------------------------------------------------------------*/
void setSpeciesDensityGuess( const mesh *m, const int k )
{
  int i;
  
  if(k==0) { /* setting the guess from values read from species.inp */
    //for (i=1; i <= chemNet.nSpec;     i++) chemNet.ymol[      i         ] = MIN_DENS;
    //for (i=1; i <= chemNet.nBaseSpec; i++) chemNet.ymol[ indBaseSpec[i] ] = ics.abun[i]*pv.dens0;
      
    //for (i=1; i <= chemNet.nBaseSpec; i++) chemNet.ymol[i] /= pv.dens0;
    //writeArrayToBinaryFile("dens-z-2-s-200.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
    
    //readArrayFromBinaryFile("dens-z-2-s-200.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
    //readArrayFromBinaryFile("/home/mher/ism/runs/twoSided/tsTest2/dens-z-2-s-200-slab-0-eq.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
    
    //readArrayFromBinaryFile("src/dens-z-2-s-200-slab-0-eq.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
    //readArrayFromBinaryFile("5.dat", &(chemNet.ymol[0]), chemNet.nSpec+1 );
    
    getNearestAbundancesGuessFromDatabase(pv.dens0, pv.G0, pv.gamma_mech); /* copies guess abundances from the database */
    if( verbose & VERB2_MASK) printf("%sSet the guess abundances from the database\n", INDENT2); /* to chemNet.ymol for */
    
    for (i=1; i <= chemNet.nSpec; i++) chemNet.ymol[i] *= pv.dens0;
    
  } else  /* setting chemNet.ymol[*] to m->abun[*][previous slab] i.e of index k-1 */
      for(i=1; i <= chemNet.nSpec; i++)
	chemNet.ymol[i]=m->abun[i-1][k-1]*pv.dens0; /* ymol has 1 indexing, m->abun has 0 indexing*/
}
/*----------------------------------------------------------------------------------------------------*/
void setGasTemperatureGuess( const mesh *msh, const int k, double *Tptr )
{
  if( k==0 ) {
    *Tptr = getNearestTemperateGuessFromDatabase(pv.dens0, pv.G0, pv.gamma_mech);
    if( verbose & VERB2_MASK) printf("%sSet the guess temperature to %.2f from the database\n", INDENT2, *Tptr);
  }
  else {
    *Tptr = msh->gasT[k-1];
  }

}
/*----------------------------------------------------------------------------------------------------*/
void setEquilibriumGuesses( const mesh *msh, const int k, double *Tptr)
{
  setGasTemperatureGuess(msh, k, Tptr);    
  setSpeciesDensityGuess(msh, k);  
 
}
