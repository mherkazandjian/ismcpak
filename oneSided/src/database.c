#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */
 
/*----------------------------------------------------------------------------------------------------*/
#define CHECK_DATABASE_INTEGRITY if( (bytesRead - sizeof(long)) != offset ) {\
    if( verbose & VERB1_MASK) fprintf(outputFd, "%s%s:%d:database integrity check failed\n", INDENT1, __FILE__, __LINE__); \
    exit(-1);}

/* reads the database from a binary file, and puts all the read data into the struct
 * 'database' which is a global variable.
 * The amount of memory required for the database is:
 *   dbSize = 4 + 4 + [(4+nSpec)*8]*n + 4 (bytes)
 *   n      = number of points in the database
 *   nSpec  = number of species
 * For example : for n = 1000 and nSpec = 309
 *               dbSize ~ 2.3MB
 *
 * In the binary data file, each block is saperated by a checkNum which is
 * the offset in bytes from the beginning of the binary file.
 *
 * depends on : pv.database_fName, VERB1_MASK, INDENT1, verbose, outputfd
 * modified   : database
 */
int init_database( )
{
	  FILE *fd;
	  long i, offset, nRead, bytesRead;

	  bytesRead = 0;

	  if( (fd=fopen(pv.database_fName, "r"))!=NULL ) {

	    if( verbose & VERB1_MASK) fprintf(outputFd, "%sloading database\n", INDENT1);

	      /* reading the header */
	      nRead = fread(&offset, sizeof(long), 1, fd);  // printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;

	      fread(&database.n    , sizeof(int), 1, fd);
	      fread(&database.nSpec, sizeof(int), 1, fd);
	      bytesRead += 2*sizeof(int);

	      /* allocating space for the database */
	      database.logDens   = (double *)malloc(sizeof(double)*database.n);
	      database.logG0     = (double *)malloc(sizeof(double)*database.n);
	      database.logGammaM = (double *)malloc(sizeof(double)*database.n);
	      database.T_eq      = (double *)malloc(sizeof(double)*database.n);
	      database.abunBlk   = (double *)malloc(sizeof(double)*database.n*database.nSpec);

	      if( verbose & VERB1_MASK) fprintf(outputFd,
						"%swill load info about %d points, number of species = %d \n",
						INDENT1, database.n, database.nSpec);

	      /* reading the densities */
	      nRead = fread(&offset, sizeof(long), 1, fd);     //  printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;
	      nRead = fread(&(database.logDens[0]), sizeof(double), database.n, fd);
	      bytesRead += nRead*sizeof(double);

	      /* reading the G0s' */
	      nRead = fread(&offset, sizeof(long), 1, fd);     //  printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;
	      nRead = fread(&(database.logG0[0]), sizeof(double), database.n, fd);
	      bytesRead += nRead*sizeof(double);

	      /* reading the mechanical heating */
	      nRead = fread(&offset, sizeof(long), 1, fd);       // printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;
	      nRead = fread(&(database.logGammaM[0]), sizeof(double), database.n, fd);
	      bytesRead += nRead*sizeof(double);

	      /* reading the equilibrium temperatures */
	      nRead = fread(&offset, sizeof(long), 1, fd);      // printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;
	      nRead = fread(&(database.T_eq[0]), sizeof(double), database.n, fd);
	      bytesRead += nRead*sizeof(double);

	      /* reading the abundances */
	      for(i=0; i<database.n; i++) {
	    	  nRead = fread(&offset, sizeof(long), 1, fd);    //   printf("offset = %ld\n", offset);
	    	  bytesRead += nRead*sizeof(long);
	    	  CHECK_DATABASE_INTEGRITY;
	    	  nRead = fread(&(database.abunBlk[i*database.nSpec]), sizeof(double), database.nSpec, fd);
	    	  bytesRead += nRead*sizeof(double);
	      }

	      nRead = fread(&offset, sizeof(long), 1, fd);      // printf("offset = %ld\n", offset);
	      bytesRead += nRead*sizeof(long);
	      CHECK_DATABASE_INTEGRITY;

	      /*
	      fread(&(database.abunBlk[0]), sizeof(double), database.n*database.nSpec, fd);
	      fread(&i, sizeof(int), 1, fd);

	      printf("check int = %d\n", i);
	      */
	      if( verbose & VERB1_MASK) fprintf(outputFd, "%sloaded database %s.\n", INDENT1, pv.database_fName);

	      fclose(fd);
	  } else {
	    printf("failed to open database file %s\n", pv.database_fName);
	    exit(-1);
	  }

	  return 0;
}
/*----------------------------------------------------------------------------------------------------*/
int getIndexNearestPointInDatabase(const double dens0, const double G0, const double Lmech)
{
  int i, iRet=-1;
  double x, y, z, l2Norm, lMin=-1;
  
  for(i=0; i<database.n; i++) {
    
    x = fabs( log10(dens0) - database.logDens[i] );
    y = fabs( log10(G0)    - database.logG0[i]   );
    z = fabs( log10(Lmech) - database.logGammaM[i]);

    l2Norm = sqrt( x*x + y*y + z*z );

    if(i==0) {
      lMin=l2Norm;
      iRet=0; }
    else {
      if ( l2Norm < lMin) {
	lMin = l2Norm;
	iRet = i;
      }
    }
    
  }
    
  if( verbose & VERB2_MASK) fprintf(outputFd, 
				    "%sclosest index in database = %d, [log(G0),log(dens)]=[%.2f,%.2f] l2Distance = %.2e\n", 
				    INDENT2, iRet, database.logG0[iRet], database.logDens[iRet], lMin );
  
  if(lMin<0 || iRet==-1) {
    fprintf(stderr, "closest point not found !! this is not supposed to happen\n");
    exit(-1);
  }
  
  return iRet;
}
/*----------------------------------------------------------------------------------------------------*/
double getNearestTemperateGuessFromDatabase(const double dens0, const double G0, const double gammaM)
{
  int idx;
  
  /* getting the index of the closest point in the database */
  idx=getIndexNearestPointInDatabase(dens0, G0, gammaM);
  return database.T_eq[idx];
}
/*----------------------------------------------------------------------------------------------------*/
void getNearestAbundancesGuessFromDatabase(const double dens0, const double G0, const double gammaM)
{
  int i, idx;
  
  idx=getIndexNearestPointInDatabase(dens0, G0, gammaM);
  
  for(i=0; i<chemNet.nSpec; i++) 
    chemNet.ymol[i+1] = database.abunBlk[ idx*chemNet.nSpec + i ];
}
