#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

#include "read_reaction_data.h"

/*----------------------------------------------------------------------------------------------------*/
/* read and initialize all the chemical network                                                       */
/*----------------------------------------------------------------------------------------------------*/
void chemical_network_data(void)
{
  if( verbose & VERB1_MASK) fprintf(outputFd, "%sReading the chemical network and the species\n", INDENT1);

  if( verbose & VERB2_MASK) fprintf(outputFd, "%sReading/initializing Species file and initial abundances.\n", INDENT2);
  readDataHeaderSpecies();          /* reading the header of species.inp                      */            
  readInitialAbundanceSpecies();    /* read initial abundances of initially abundant elements */
  initSpecies();                    /* read the specie index and its chemical formula as a string */

  if( verbose & VERB2_MASK) fprintf(outputFd, "%sReading/initializing Underabundant Species.\n", INDENT2);
  readDataHeaderUnderUbundant();    /* reading underAbundant.inp */
  readUnderAbundantData();

  if( verbose & VERB2_MASK) fprintf(outputFd, "%sReading/initializing reaction network and rates data.\n", INDENT2);
  initReactionRateBuffers();        /* reading chemical reation rates */
  readReactionRates();
}
/*----------------------------------------------------------------------------------------------------*/
/* initialize the species data struct and read key header paramters                                   */
void readDataHeaderSpecies()
{
  int i, j;
  char fName[256];  /* spec[0].name has a size nSpec + 1, so allocate it dynamically */

  spec     =(species *)malloc(sizeof(species)*310);
  underAbun=(species *)malloc(sizeof(species)*103); 

  /* initialize the element strings to null characters */
  for(i=0;i<310;i++) for(j=0;j<10;j++) spec[i].element[j]='\0';
  for(i=0;i<103;i++) for(j=0;j<10;j++) underAbun[i].element[j]='\0';

  sprintf(fName, "%s", pv.species_fName); /* filename of input species */

  if( verbose & VERB1_MASK ) fprintf(outputFd, "     Reading %s\n",fName); 
  if ((spec[0].FILE=fopen(fName,"r"))==NULL) errMsg(4);
  
  fscanf(spec[0].FILE, "%d", &chemNet.nSpec);             /* Number of species                        */
  fscanf(spec[0].FILE, "%d", &chemNet.nBaseSpec);         /* Number of species making up all the rest */
  fscanf(spec[0].FILE, "%d", &chemNet.nReact);            /* Number of reactions                      */

  if( verbose & VERB1_MASK ) {
    //fprintf(outputFd,"       N species in the network      = %d \n", chemNet.nSpec);
    //fprintf(outputFd,"       N elements initially abundant = %d \n", chemNet.nBaseSpec);
    //fprintf(outputFd,"       N reactions in network        = %d \n", chemNet.nReact);
    //fprintf(outputFd,"\n");
  }

  MAKE_1D_ARRAY(indBaseSpec, chemNet.nBaseSpec+1, int);       /* Replaced equations            */
  MAKE_1D_ARRAY(ics.abun   , chemNet.nBaseSpec+1, double);    /* Elemental relative abundances */
}
/*----------------------------------------------------------------------------------------------------*/
/* read the indicies of the basic species in species.inp ( second line ) and their corresponding      */
/* abundances relative to H ( third line )                                                            */ 
void readInitialAbundanceSpecies()
{
  int i;
  
  for (i=1; i <= chemNet.nBaseSpec; i++) fscanf(spec[0].FILE, "%d",  &indBaseSpec[i]); 
  for (i=1; i <= chemNet.nBaseSpec; i++) fscanf(spec[0].FILE, "%le", &ics.abun[i]); 

  /*
  if( verbose & VERB1_MASK ) {
    for (i=1; i<=chemNet.nBaseSpec; i++) fprintf(outputFd,"       element [%3d], abundance = %+15.8le\n", indBaseSpec[i], ics.abun[i]);
    fprintf(outputFd, "\n"); 
  }
  */

  MAKE_1D_ARRAY(indSpec,  chemNet.nSpec+1, int);
  MAKE_2D_ARRAY(specCode, chemNet.nSpec+1, chemNet.nBaseSpec+1, int); 
}
/*----------------------------------------------------------------------------------------------------*/
/* initializes the species indecies and numerical representation                                      */
void initSpecies( void ) 
{
  int i;
  
  for (i=1; i <= chemNet.nSpec; i++) {
    readSpecie(i);      /* read the specie index and its chemical formula as a string */
    setup_specCode(i);  /* fill in the number of elements of each basic specie        */
  }
}
/*----------------------------------------------------------------------------------------------------*/
/* reads one species line in species.inp, the element index and the species chemical formula          */
/* the specie formula string should be 10 characters wide in the file!!!                              */
void readSpecie(int i)
{
  fscanf(spec[0].FILE,"%d", &indSpec[i]);
  READ_SPACE(2, spec[0].FILE);
  fgets(spec[i].element, 10, spec[0].FILE);
  replaceChar(spec[i].element, ' ', '\0', 10);
}
/*----------------------------------------------------------------------------------------------------*/
/* sets up the variable specCode, which is an array that holds the number of atoms of each basic      */
/* specie in the molecular formula of the species in the input file                                   */

void setup_specCode( int i) 
{
  int j;
  char *elemStr;
  
  for (j=0; j < 9; j++) {
    
      elemStr=&(spec[i].element[0]);

      if (elemStr[j] == '-') specCode[i][1] = -1;
      if (elemStr[j] == '+') specCode[i][1] = 1;
      if (elemStr[j] == 'H' && 
	  elemStr[j+1] != 'e' && elemStr[j+1] != '2' && elemStr[j+1] != '3' && 
	  elemStr[j+1] != '4' && elemStr[j+1] != '5' && elemStr[j+1] != '6' && 
	  elemStr[j+1] != '7' && elemStr[j+1] != '8' && elemStr[j+1] != '9' && 
	  elemStr[j+1] != '1') specCode[i][2] = specCode[i][2] + 1;

      if (elemStr[j] == 'H' && elemStr[j+1] == '2') specCode[i][2] = specCode[i][2] + 2;
      if (elemStr[j] == 'H' && elemStr[j+1] == '3') specCode[i][2] = specCode[i][2] + 3;
      if (elemStr[j] == 'H' && elemStr[j+1] == '4') specCode[i][2] = specCode[i][2] + 4;
      if (elemStr[j] == 'H' && elemStr[j+1] == '5') specCode[i][2] = specCode[i][2] + 5;
      if (elemStr[j] == 'H' && elemStr[j+1] == '6') specCode[i][2] = specCode[i][2] + 6;
      if (elemStr[j] == 'H' && elemStr[j+1] == '7') specCode[i][2] = specCode[i][2] + 7;
      if (elemStr[j] == 'H' && elemStr[j+1] == '8') specCode[i][2] = specCode[i][2] + 8;
      if (elemStr[j] == 'H' && elemStr[j+1] == '9') specCode[i][2] = specCode[i][2] + 9;
      if (elemStr[j] == 'H' && elemStr[j+1] == '1' && elemStr[j+2] == '0') specCode[i][2] = specCode[i][2] + 10;
      if (elemStr[j] == 'H' && elemStr[j+1] == '1' && elemStr[j+2] == '3') specCode[i][2] = specCode[i][2] + 1;
      if (elemStr[j] == 'H' && elemStr[j+1] == 'e') specCode[i][3] = 1;
      if (elemStr[j] == 'C' && 
	  elemStr[j+1] != 'l' && elemStr[j+1] != '2' && elemStr[j+1] != '3' && 
	  elemStr[j+1] != '4' && elemStr[j+1] != '5' && elemStr[j+1] != '6' && 
	  elemStr[j+1] != '7' && elemStr[j+1] != '8' && elemStr[j+1] != '9' && 
	  elemStr[j+1] != '1') specCode[i][4] = specCode[i][4] + 1;
      if (elemStr[j] == 'C' && elemStr[j+1] == '2') specCode[i][4] = specCode[i][4] + 2;
      if (elemStr[j] == 'C' && elemStr[j+1] == '3') specCode[i][4] = specCode[i][4] + 3;
      if (elemStr[j] == 'C' && elemStr[j+1] == '4') specCode[i][4] = specCode[i][4] + 4;
      if (elemStr[j] == 'C' && elemStr[j+1] == '5') specCode[i][4] = specCode[i][4] + 5;
      if (elemStr[j] == 'C' && elemStr[j+1] == '6') specCode[i][4] = specCode[i][4] + 6;
      if (elemStr[j] == 'C' && elemStr[j+1] == '7') specCode[i][4] = specCode[i][4] + 7;
      if (elemStr[j] == 'C' && elemStr[j+1] == '8') specCode[i][4] = specCode[i][4] + 8;
      if (elemStr[j] == 'C' && elemStr[j+1] == '9') specCode[i][4] = specCode[i][4] + 9;
      if (elemStr[j] == 'C' && elemStr[j+1] == '1' && elemStr[j+2] == '0') specCode[i][4] = specCode[i][4] + 10;
      if (elemStr[j] == 'N' && elemStr[j+1] != 'a' && elemStr[j+1] != '2') specCode[i][5] = specCode[i][5] + 1; 
      if (elemStr[j] == 'N' && elemStr[j+1] == '2') specCode[i][5] = specCode[i][5] + 2;
      if (elemStr[j] == 'O' && elemStr[j+1] != '2') specCode[i][6] = specCode[i][6] + 1;
      if (elemStr[j] == 'O' && elemStr[j+1] == '2') specCode[i][6] = specCode[i][6] + 2;
      if (elemStr[j] == 'N' && elemStr[j+1] == 'a') specCode[i][7] = 1;
      if (elemStr[j] == 'M' && elemStr[j+1] == 'g') specCode[i][8] = 1;
      if (elemStr[j] == 'S' && elemStr[j+1] == 'i') specCode[i][9] = 1;
      if (elemStr[j] == 'P' && elemStr[j+3] != 'T') specCode[i][10] = 1;
      if (elemStr[j] == 'S' && elemStr[j+1] != 'i' && elemStr[j+1] != '2') specCode[i][11] = specCode[i][11] + 1; 
      if (elemStr[j] == 'S' && elemStr[j+1] == '2') specCode[i][11] = specCode[i][11] + 2;
      if (elemStr[j] == 'C' &&  elemStr[j+1] == 'l') specCode[i][12] = 1;
      if (elemStr[j] == 'F' && elemStr[j+1] == 'e') specCode[i][13] = 1;
      if (elemStr[j] == 'P' && elemStr[j+1] == 'A' && elemStr[j+2] == 'H') {
	specCode[i][14] = specCode[i][14] + 1;
	specCode[i][10] = specCode[i][10] - 1;
	specCode[i][2] = specCode[i][2] - 1;
      }
      if (elemStr[j] == '1' && elemStr[j+1] == '3' && elemStr[j+2] == 'C') {
	specCode[i][15] = specCode[i][15] + 1;
	specCode[i][4] = specCode[i][4] - 1;
      }
  }
}
/*----------------------------------------------------------------------------------------------------*/
/* initilize and read the header of underubundant.inp                                                 */
void readDataHeaderUnderUbundant( void )
{
  char fName[256];
  
  sprintf(fName, "%s", pv.underUbundant_fName);
  //if( verbose & VERB1_MASK ) fprintf(outputFd, "\n     Reading %s\n", fName);   

  if ((underAbun[0].FILE=fopen(fName,"r"))==NULL) errMsg(5);

  fscanf(underAbun[0].FILE, "%d", &nUnabundant); /* number of underAbundant species */  
  //if( verbose & VERB1_MASK ) fprintf(outputFd, "       number of under-abundunt species = %d\n", nUnabundant);

  MAKE_1D_ARRAY(indUnabundant, nUnabundant+1, int);
}
/*----------------------------------------------------------------------------------------------------*/
/* read the data of underubundant.inp                                                                 */
void readUnderAbundantData( void )
{
    int i;

    for (i = 1; i <= nUnabundant; i++) { 
      fscanf(underAbun[0].FILE,"%d", &indUnabundant[i]);
      READ_SPACE(2, underAbun[0].FILE);
      fgets(underAbun[i].element, 10, underAbun[0].FILE);
      replaceChar(underAbun[i].element, ' ', '\0', 10);
    }
}
/*----------------------------------------------------------------------------------------------------*/
/* initialize the chemical reaction data arrays and open the input rate file, rate99.inp              */
void initReactionRateBuffers( void )
{
  char fName[256];

  sprintf(fName, "%s", pv.rate99_fName);
  //if( verbose & VERB1_MASK ) fprintf(outputFd, "\n     Reading %s\n", fName);   

  if ( (chemNet.fd=fopen(fName,"r"))==NULL ) errMsg(6);

  MAKE_2D_ARRAY(chemNet.molno   , chemNet.nReact+1, 8, int);  
  MAKE_1D_ARRAY(chemNet.indReact, chemNet.nReact+1, int   );
  MAKE_1D_ARRAY(chemNet.alpha,    chemNet.nReact+1, double);
  MAKE_1D_ARRAY(chemNet.beta ,    chemNet.nReact+1, double);
  MAKE_1D_ARRAY(chemNet.gamma,    chemNet.nReact+1, double);
  MAKE_1D_ARRAY(ics.maxAbun, chemNet.nSpec+1, double);
  /* workspace buffers for chemical balance computations */
  MAKE_1D_ARRAY(chemNet.rates,    chemNet.nReact+1, double);
  MAKE_1D_ARRAY(chemNet.ymol,     chemNet.nSpec+1, double);
  MAKE_1D_ARRAY(chemNet.F, chemNet.nSpec+1,double); 
  MAKE_1D_ARRAY(chemNet.X, chemNet.nSpec+1,double); 
  MAKE_2D_ARRAY(chemNet.DFDX, chemNet.nSpec+1, chemNet.nSpec+1,double);
}
/*----------------------------------------------------------------------------------------------------*/
/* - read the reaction rates of all the reactions                                                     */
/* - setup the numerical representation of the chemical network                                       */ 
void readReactionRates( void )
{
  int i;  
  for (i = 1; i <= chemNet.nReact; i++) parseReactionRateLine(i);
}
/*----------------------------------------------------------------------------------------------------*/
/* extracts all neccessary info for a certain line of a reaction read from the input file             */
void parseReactionRateLine(int rxnLineNum)
{
  int i, k, **molno=chemNet.molno;

  fscanf(chemNet.fd,"%d", &chemNet.indReact[rxnLineNum]);  /* reaction number */
  
  READ_SPACE(2, chemNet.fd);
  /* setting up the entries of chemNet.molno for each reaction */
  readSpecieFld(rxnLineNum, 1, 10); /* reactant 1 */  
  readSpecieFld(rxnLineNum, 2, 10); /* reactant 2 */
  readSpecieFld(rxnLineNum, 3, 10); /* reactant 3 */
  readSpecieFld(rxnLineNum, 4, 10); /* product 1 */
  readSpecieFld(rxnLineNum, 5, 10); /* product 2 */
  readSpecieFld(rxnLineNum, 6, 6);  /* product 3 */
  readSpecieFld(rxnLineNum, 7, 6);  /* product 4 */
  
  /* strike out the underAbunant reactions */
  for (k=1; k <= 7; k++)
    if (molno[rxnLineNum][k] == -1)
      for (i=1; i <= 7; i++) molno[rxnLineNum][i] = 0;
  
  fscanf(chemNet.fd,"%le %le %le", &chemNet.alpha[rxnLineNum], &chemNet.beta[rxnLineNum], &chemNet.gamma[rxnLineNum]);
  READ_SPACE(17, chemNet.fd);
}
/*----------------------------------------------------------------------------------------------------*/
/* reads one specie fld in the reaction data file, i.e it reads the formula of a species whether it   */
/* is reactant or product                                                                             */
/* fldWidth : number of character to be read for this fld                                             */
void readSpecieFld(int j, int idx, int fldWidth)
{
  char fld[10];

  fgets(fld, fldWidth,  chemNet.fd);
  replaceChar(fld, ' ', '\0', fldWidth);
  
  chemNet.molno[j][idx] = getSpecieIndex(fld);
}
/*----------------------------------------------------------------------------------------------------*/
/* Give the species a numbers                                                                         */
double getSpecieIndex(char str[10])
{
  int i;

  /* if the specie belongs to species.inp then return its index in species.inp */
  for (i=1; i <= chemNet.nSpec; i++)        if( strcmp(spec[i].element, str)==0 )      return(indSpec[i]);
  /* if the specie belongs to underAbundant.inp then return -1                 */
  for (i=1; i <= nUnabundant; i++)  if( strcmp(underAbun[i].element, str)==0 ) return(-1);
  /* if it does not belong to both, return 0                                   */
  return(0); 
}
