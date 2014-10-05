#ifndef __VARS_H
#define __VARS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <ctype.h>
#include <fcntl.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#ifdef ENABEL_MPI
#include <mpi.h>
#endif

#include "definitions.h" /* Service macros   */
#include "constants.h"   /* Constants        */

/*-------------------------------------COOLING-HEATING------------------------------------------------*/
/* workspace buffers used in computing fine structure line cooling */
typedef struct { /* 56 B = 7*8  */
  double **CP; 
  double **SiP;
  double **CA;
  double **OA;
  double **SA;
  double **FeP;
  double **SiA;
} populationDensitiesFineStruct;

/* fine structure cooling utility buffer, cooling rates for some species */
typedef struct { /* 120 B = 8*8 + 56, all the arrays use 1 based indexing*/
  double **CP;   
  double **SiP;
  double **CA;
  double **OA;
  double **SA;
  double **FeP;
  double **SiA;
  double *total;
  populationDensitiesFineStruct popDens;
} coolingFineStruct;

/* storage arrays for metastable cooling for each slab */
typedef struct { /* 88 B = 11*8 */
  double *CI;
  double *CII;
  double *FeI;
  double *FeII;
  double *OI;
  double *OII;
  double *SI;
  double *SII;
  double *SiI;
  double *SiII;
  double *total;
} coolingMetaStable;

/* array where the heating due to process X is stored for each slab */ 
typedef struct { /* 64 B = 8*8 */
  double *photo;                   /* dust grains               */
  double *cIon;                    /* carbon ionization         */
  double *molHydro;                /* molecular hydrogen        */
  double *H2pump;                  /* pumping of H2 molecules   */
  double *ggColl;                  /* gas grain collisions      */
  double *visc;                    /* viscouse heating          */
  double *CR;                      /* cosmic ray heating        */
  double *total;                   /* total heating, sum of all */
} heatingProcesses;

/* array where the cooling due to process X is stored for each slab */ 
typedef struct { /* 240 B = 88 + 120 + 4*8 */
  coolingMetaStable ms;            /* metastable cooling        */
  coolingFineStruct fs;            /* fine structure cooling    */
  double *roVib;                   /* rotationa-vibrational     */
  double *recom;                   /* recombination cooling     */
  double *lymanAlpha;              /* OI and HI Lyman alpha     */
  double *total;                   /* total cooling, sum of all */
} coolingProcesses;

/* stores all rotationcal cooling data */
typedef struct {
  FILE *FILE;
  char name[32];
  int rows_vib, columns_vib;
  double *ColumnDens_vib, *Temperature_vib, *L0_vib;
  double **Llte_vib;
  int rows_rot, columns_rot;
  double *ColumnDens_rot, *Temperature_rot, *L0_rot;
  double **Llte_rot, **nhalf_rot, **alpha_rot;
  double **paramH2_rot;
} cooling_lines;
extern cooling_lines coollines[COOLSPEC];

/*------------------------------------------------THE MESH--------------------------------------------*/
typedef struct { /* 512 B = 4 + 2*8 + 4*8 + 2*8 + 1*8 + 240 + 64 + 5*8 + 5*8 + 6*8 */
  int nx;                     /* number of slabs in the x direction                                   */
  int nxFilled;               /* last index of the slab with stores data                              */
  double *gasT;               /* gas temperature                                                      */
  double *dustT;              /* dust temperature                                                     */

  double *xs;                 /* starting x coordinate of the slab from the origin                    */
  double *xe;                 /* ending   x coordinate of the slab from the origin                    */
  double *xc;                 /* centroid x coordinate of the slab from the origin                    */
  double *dx;                 /* the thickness of the slab                                            */
  
  double *AvL;                /* Av of the slabs due to the source on the left                        */
  double *AvR;                /* Av of the slabs due to the source on the right                       */

  double **abun;              /* relative adundances of all the species [specie index][slab index ]   */
  
  heatingProcesses heating;   /* heating rate for each slab                                           */
  coolingProcesses cooling;   /* cooling rate for each slab                                           */

  double *NH2All_L;           /* colon densities to the left of the start of each slab for some species */
  double *NH2_L;              /* when looking from the origin                                           */
  double *NCO_L;
  double *NC_L;
  double *NH2O_L;

  double *NH2All_R;           /* colon densities to the right of the end of each slab for some species */
  double *NH2_R;              /* when looking from the origin                                          */
  double *NCO_R;              /* alloc all these                                                       */
  double *NC_R;
  double *NH2O_R;

  double *betaH2_L;             /* self sheilding factors due to colons to the left                     */
  double *betaCO_L;
  double *beta13CO_L;
  double *betaH2_R;             /* self sheilding factors due to colons to the right                    */
  double *betaCO_R;
  double *beta13CO_R;
} mesh;
/*--------------------------------------- CHEMISTRY---------------------------------------------------*/
extern int *indBaseSpec,      /* index of the base species making up all the other species            */
           *indSpec,          /* unique index of the elements in the species file                     */
           **specCode;        /* numerical representation of the species                              */ 

extern int nUnabundant;       /* number of un-abundant species                                        */  
extern int *indUnabundant;    /* index of the unabundant species                                      */

typedef struct {              /* typedef for struct holding info about the species                    */
  FILE *FILE;                 /* file discreptor                                                      */
  char element[10];           /* string name of the specie                                            */
} species;
extern species *spec;
extern species *underAbun;

typedef struct {         /* structure holding all the info on the chemical network and species          */
  int     nSpec;         /* number of species                                                           */
  int     nBaseSpec;     /* number of specieis initially with non-zero abundances                       */
  int     nReact;        /* total number of reactions in the chemical network                           */
  int    *indReact;      /* index of the reaction in the reaction rate info file                        */
  double *alpha;         /* reaction rate constants in UMIST database (see Table 3 in the paper)        */
  double *beta;       
  double *gamma;      
  double *rates;         /* array holding all the reaction rates ( see eq-1 in the paper )              */
  FILE   *fd;            /* file discriptor for reading the chemical reaction ascii file                */
  int    **molno;        /* array holing the numerical representation of the reactions in the network   */
  double *ymol;          /* absolute abundance of all the species                                       */
  double *F;
  double *X;
  double **DFDX;
} chemicalNetwork;
extern chemicalNetwork chemNet;


typedef struct {
  double *abun;         /* array which holds the abundance of the initially abundant species           */
                        /* relative to the abundance of H                                              */
  double *maxAbun;      /* maximum allowed absolute abundances of all the species                      */
} initialConditions;
extern initialConditions ics;

typedef struct {    /* size = 4 + 4 + [(4+nSpec)*8]*n + 4 bytes                       */
  int     n;        /* number of points in the database                   */
  int     nSpec;    /* number of species in the abun block for each point */
  double *logDens;  /* log10 of density, G0, and mechanical heating       */
  double *logG0;
  double *logGammaM;
  double *T_eq;
  double *abunBlk;  /* abundances block for all the points                */
  int status;
} guessDatabase;
extern guessDatabase database;

/*-----------------------------------------PARAMETERS-------------------------------------------------*/
#define MAXSTRLEN 1000
#define NTAGS 22                   /* number of tags in the parameters struct ,set by user */
typedef struct{                    /* parameters read from parameter file, set by use      */
    double dens0;
    double G0;
    double gamma_mech;
    double T_init;
    char outputDir[MAXSTRLEN];
    char species_fName[MAXSTRLEN];
    char underUbundant_fName[MAXSTRLEN];
    char rate99_fName[MAXSTRLEN];
    char selfSheilding_CO_fName[MAXSTRLEN];
    char rotationalCooling_baseName[MAXSTRLEN];
    char vibrationalCooling_baseName[MAXSTRLEN];
    double zeta;
    double S_depletion;
    double TTol;
    double CTol;
    double metalicity;
    double AvMax;
    double slabSizeCrit;
    double max_deltaAv;
    double min_deltaAv;
    int    maxSlabs;
    char database_fName[MAXSTRLEN];
} parameters;
extern parameters pv;
extern char tagNames[NTAGS][MAXSTRLEN];
/*---------------------------------------------MISC STUFF---------------------------------------------*/
extern double tau100;         /* emmision optical depth at 100 micro-meters                           */
extern double T0;             /* the equilibirum dust temperature at the slab surface du to the       */
                              /* incident FUV field alone                                             */
                              /* see hollenbach 1991                                                  */

extern int nNewton;           /* number of newton-raphson steps in chemical network equilibrium comp  */
extern int  nBisect;          /* number of bisection steps in thermal balance                         */

extern int heapUse;           /* amount of bytes used in allocating all bufferes in the code using    */
                              /* MAKE_1D_ARRAY, MAKE_2D_ARRAY...etc                                   */
extern FILE *outputFd;        /* file discriptor pointer used to dump the diagnostic output           */
extern unsigned int verbose;  /* used to check what to dump to the file outputFd                      */

/* workspace buffers */
extern char __SPACE_BUFFER__[100];        /* dummy buffer used to read blank space                    */
/*----------------------------------------SELF SHEILDING----------------------------------------------*/
/* holds self sheilding data info which is read from a table and interpolated upon                    */
typedef struct {
  /* the following three tags are used only in CO_s */
  FILE *FILE;            /* file pointer to the data file                                            */
  int npointCO;          /* number of points in interpolation table of CO self sheilding             */
  int npointH2;          /* number of points in interpolation table of H2 self sheilding             */
  
  char name[32];            //not used anywhere.
  double *logNCO, *logNH2;  //grid points in NCO and NH2 in the table. (number of columns, number of rows in the table
  double **function;        //2D table holding the self sheilding values.
} self_shlding;
extern self_shlding self_shld[MOL];
/*-----------------------------------------TIMERS-----------------------------------------------------*/
typedef struct {              /* typedef for timekeeping                                              */
    double ti;                /* intial time                                                          */
    double tf;                /* final time                                                           */
    double dt;                /* time difference, tf - ti                                             */
    double tot;               /* agregated vaules of dt                                               */
} timer;
extern timer tmr1;            /* tmr1: wall time from start of simulation                             */
extern timer tmr2;            /* tmr2: chemical balance timer                                         */
extern timer tmr3;            /* tmr3: timing for bracketing                                          */
extern timer tmr4;            /* tmr3: misc                                                           */

extern struct timeval tmrTmp;
extern int tmrStatus;         /* status flag for the timer                                            */

/* macro definitions */
/* ----------------- */
#define TIC(tmr)  tmrStatus=gettimeofday(&tmrTmp,NULL);\
                  tmr.ti=((double)((long long)tmrTmp.tv_sec*1000000ll + (long long)tmrTmp.tv_usec))/1e6;

#define TOC(tmr)  tmrStatus=gettimeofday(&tmrTmp,NULL);\
                  tmr.tf=((double)((long long)tmrTmp.tv_sec*1000000ll + (long long)tmrTmp.tv_usec))/1e6;\
		  tmr.dt=tmr.tf-tmr.ti;\
		  tmr.tot+=tmr.dt;

#define TRUE  1
#define FALSE 0

#define OUTPUT_MASK  0x8000
#define OUTFILE_MASK 0x4000
#define STDERR_MASK  0x2000
#define ERRHND_MASK  0x1000

#define VERB1_MASK    0x0800
#define VERB2_MASK    0x0400
#define VERB3_MASK    0x0200
#define VERB4_MASK    0x0100
#define VERB5_MASK    0x0080
#define VERB6_MASK    0x0040
#define VERB7_MASK    0x0020
#define VERB8_MASK    0x0010
#define HALF_VERBOSE  0x0F00
#define FULL_VERBOSE  0x0FF0

#define DUMP_OUTPUT     0x8000
#define TO_FILE         0x4000
#define REDIRECT_STDERR 0x2000
#define SMOOTH_EXIT     0x1000

#define INDENT1 ""
#define INDENT2 "     "
#define INDENT3 "          "
#define INDENT4 "               "
#define INDENT5 "                    "
#define INDENT6 "                         "
#define INDENT7 "                              "
#define INDENT8 "                                   "

extern int nErrors;
extern int errCode;

extern int dataVersion;
#endif
