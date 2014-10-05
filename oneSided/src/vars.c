#include "prototypes.h"
#include "vars.h"

/*--------------------------------------- CHEMISTRY---------------------------------------------------*/
chemicalNetwork chemNet;
initialConditions ics;
guessDatabase database;

int *indBaseSpec,      /* index of the base species making up all the other species                   */
    *indSpec,          /* unique index of the elements in the species file                            */
    **specCode;        /* numerical representation of the species                                     */ 

int nUnabundant;       /* number of un-abundant species                                               */  
int *indUnabundant;    /* index of the unabundant species                                             */

species *spec;
species *underAbun;
/*-----------------------------------------PARAMETERS-------------------------------------------------*/
parameters pv;

/*set by user,thesese shouldbe identical to the tag names of pv*/
char tagNames[NTAGS][MAXSTRLEN] = { "dens0"                       ,
				    "G0"                          ,
				    "gamma_mech"                  ,
				    "T_init"                      ,
				    "outputDir"                   ,  
				    "species_fName"               ,
				    "underUbundant_fName"         ,
				    "rate99_fName"                ,
				    "selfSheilding_CO_fName"      ,
				    "rotationalCooling_baseName"  ,
				    "vibrationalCooling_baseName" ,
				    "zeta"                        ,
				    "S_depletion"                 ,
				    "TTol"                        ,
				    "CTol"                        ,
				    "metalicity"                  ,
				    "AvMax"                       ,
				    "slabSizeCrit"                ,
				    "max_deltaAv"                 ,
				    "min_deltaAv"                 ,
				    "maxSlabs"                    ,
				    "database_fName"              };
/*---------------------------------------------MISC STUFF---------------------------------------------*/
double tau100;         /* emmision optical depth at 100 micro-meters                                  */
double T0;             /* the equilibirum dust temperature at the slab surface du to the              */
                       /* incident FUV field alone                                                    */
                       /* see hollenbach 1991                                                         */
int nNewton;           /* number of newton-raphson steps in chemical network equilibrium comp         */
int nBisect;           /* number of bisection steps in thermal balance                                */
int heapUse;           /* amount of bytes used in allocating all bufferes in the code                 */
FILE *outputFd;        /* file discriptor pointer used to dump the diagnostic output                  */

//unsigned int verbose = DUMP_OUTPUT | TO_FILE | REDIRECT_STDERR | SMOOTH_EXIT | FULL_VERBOSE;
//unsigned int verbose = DUMP_OUTPUT | TO_FILE | REDIRECT_STDERR | SMOOTH_EXIT | HALF_VERBOSE;

unsigned int verbose = DUMP_OUTPUT | REDIRECT_STDERR | SMOOTH_EXIT | VERB1_MASK | VERB2_MASK;
//unsigned int verbose = DUMP_OUTPUT | REDIRECT_STDERR | SMOOTH_EXIT | HALF_VERBOSE;
//unsigned int verbose = DUMP_OUTPUT | REDIRECT_STDERR | SMOOTH_EXIT | FULL_VERBOSE;
/*---------------------------------------- workspace buffers----------------------------------------- */
char __SPACE_BUFFER__[100];  /* dummy buffer used to read blank space                                 */
/*------------------------------struct holding data read from files-----------------------------------*/
self_shlding self_shld[MOL];       /* self sheilding data read from a table used for interpolation    */
cooling_lines coollines[COOLSPEC]; /* rotationcal cooling data from read file                         */
/*-----------------------------------------TIMERS-----------------------------------------------------*/
timer tmr1,            /* tmr1: wall time from start of simulation                                    */
      tmr2,            /* tmr2: chemical balance timer                                                */
      tmr3,            /* tmr3: timing for bracketing                                                 */
      tmr4;            /* tmr3: misc                                                                  */
struct timeval tmrTmp;
int tmrStatus;         /* status flag for the timer                                                   */

int nErrors=0;
int errCode=-1;
int dataVersion=2;
