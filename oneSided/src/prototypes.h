#ifndef __PROTOTYPES_H
#define __PROTOTYPES_H

#include "vars.h"       

int main();

/* prototypes in chemistry.c */
void updateColDensAndSelfSheildingDueToPreviousSlab(mesh *msh, const int k);

/* prototyxopes of init.c */
void init(int argc, char **argv);
void initConstantsAndAbundances();
void setFirstSlabAbundanceGuesses();
void setMaxAllowedDensities();
void dumpVersionFile( void );
void setupAllocMesh( mesh *, const int, const int );
void allocMesh( mesh *, const int );
void freeMesh( mesh *, int );
void initMesh( mesh * );
void initSlabSizes( mesh * );
void setEquilibriumGuesses( const mesh *msh, const int k, double *Tptr);
void setSpeciesDensityGuess( const mesh *msh, const int k);
void setGasTemperatureGuess( const mesh *msh, const int k, double *Tptr );
void setSlabTemperatureMesh( const int k, const double T, mesh *msh);
void setSlabAbundancesMesh( const int k, const double *abun, const double nDens, mesh *msh);

/* prototypes of read_reaction_data.c */
void chemical_network_data(void);
void read_selfshlding_CO(void);
void read_rotational_data(void);
void read_vibrational_data(void);

/* prototypes of chemical_balance.c */
void chemical_balance(const mesh *, const double T, const int k); 
  void initRatesNetwork( const double, const double, const double, const double, const double, const double );
  void iterateChemicalNetwork(const double AvL );
  void iterateChemicalNetworkOptimizedCorrect(const double AvL );
     void assemble_linear_system(double *F, double **DFDX); 
     void chemSpecTimescale(chemicalNetwork *CN, int k);

/* prototypes of finestruct_cooling.c */
double finestruct_cooling_CP( const double, const int, const double, const double, const double, const double              , const double *, double **, double **);
double finestruct_cooling_SiP(const double, const int, const double, const double, const double, const double              , const double *, double **, double **);
double finestruct_cooling_CA( const double, const int, const double, const double, const double, const double, const double, const double *, double **, double **);
double finestruct_cooling_OA( const double, const int, const double, const double, const double, const double, const double, const double *, double **, double **);
double finestruct_cooling_SA( const double, const int, const double, const double, const double, const double              , const double *, double **, double **);
double finestruct_cooling_FeP(const double, const int, const double, const double, const double, const double              , const double *, double **, double **);
double finestruct_cooling_SiA(const double, const int, const double, const double, const double, const double              , const double *, double **, double **);
double beta_tau_finestructure(double tau);
double calc_background_radiation(double lambda);
double blackbody(double nu, double T);
double ortho_para_ratio(double T);

/* prototypes of Matrix_inversion.h */
void solve_lu(double **A, double *b, double *x, int n);
void ludcmp(double **, int , int *);
void lubksb(double **, double *, int , int *);
void mprove(double **A, double **Al, int n, int *indx, double *b, double *bl, int nImprove);

/* prototypes of metastable_cooling.c */
double metastable_cooling_CI(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_CII(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_FeI(const double densThis, const double Hdens, const double eldens, const double T); 
double metastable_cooling_FeII(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_OI(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_OII(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_SI(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_SII(const double densThis, const double Hdens, const double eldens, const double T);
double metastable_cooling_SiI(const double densThis, const double Hdens, const double eldens, const double T); 
double metastable_cooling_SiII(const double densThis, const double Hdens, const double eldens, const double T); 

/* prototypes of rovibrational_cooling.c */
double L_Vib_H2(double T, double H2dens);
double L_Rot_H2(double T, double H2dens); 
double L_Vib_H2O_and_CO(double T, double Ntilde, double H2dens, int coolspec);
double L_Rot_H2O_and_CO(double T, double Ntilde, double H2dens, int coolspec);

/* prototypes of self_shlding.c */
double Self_shlding(double tau); 
double small_tau(double tau); 
double big_tau(double tau);
double self_shlding_CO(double NCO, double NH2, int what_iso);

/* prototypes of thermal_balance.c */
double thermal_balance(mesh *, const double T, const double Td, const int steps, const double AvL, const double *ymol);

/* prototypes of heating.c */
double heatingEfficiency(const double T, const double AvL, const double eldens);
double heating_PAH(mesh *msh, const double T, const double AvL, const int k, const double elecDens);
double heating_carbon(mesh *msh, const double AvL, const int k, const double carbonDens);
double heating_molecular_hydrogen(mesh *msh, const double AvL, const int k, const double H2Dens);
double heating_gas_grain_collisions(mesh *msh, const double T, const double Td, const int k);
double heating_viscous(mesh *msh, const double T, const double AvL, const int k, const double elecDens, const double CIIDens);
double Gfunct(const double y);
double heating_cosmic_rays(mesh *msh, const int k, const double H2Dens);
double heating_H2_pumping(mesh *msh, const double T, const int k, const double HDens, const double H2Dens, const double H2exDens);

/* prototypes of cooling.c */
double metaStableCooling(mesh *, const double T, const int k, const double *ymol);
double fineStructureCooling(mesh *msh, const double T, const int k, const double *ymol);
double rovVibCooling(mesh *msh, const double , const int , const double , const double , const double , const double , const double );
double recombinationCooling(mesh *msh, const double T, const int k, const double AvL, const double, const double);
double OIandHI_LymanAlphaCooling(mesh *msh, const double , const int , const double , const double , const double , const double );

/* prototypes of util.c */
void allocWorkspaceBuffers( void );
double getDustTemp(double AvL);
void writeDataOutputBinary( mesh *m, int, char * );
void set_diagnostics_parameters( char *outFname );
void finalize( mesh * );

int find_mesh_equilibrium(mesh *); 
int getSystemEquilibrium(mesh *, int k); 
    void convergeToTemperatureRoot(mesh *msh, double Ta, double Tb, double Fa, double Fb, double *Tx, double AvL, double *ymol, double Td, int k );
void errMsg(const short int errNo);
int readPfile(char *fnameStr,parameters *parms);

void replaceChar( char *str, char a, char b, int n );

double H2_stickingFactor(const double , const double );
double H2_efficiency(const double T, const double Td);

double ** malloc2dArray_d(int l, int m);
int free2dArray_d(double **x2d);

double *** malloc3dArray_d(int l, int m, int n);
int free3dArray_d(double ***x3d);

int writeArrayToBinaryFile (const char *fName, const double *vecPtr, const int n);
int readArrayFromBinaryFile(const char *fName,       double *vecPtr, const int n);
void dumpFluxVsTSweep(mesh *msh, const int k, const double T, const double Td, const double AvL);
double Av2Position(const double Av, const double dens, const double Z);
double position2Av(const double x, const double dens, const double Z);
long descretizeSlab1( double *AvTmp );

/* prototypes of database.c */
int init_database( );
int getIndexNearestPointInDatabase(const double dens0, const double G0, const double gammaM);
double getNearestTemperateGuessFromDatabase(const double dens0, const double G0, const double gammaM);
void getNearestAbundancesGuessFromDatabase(const double dens0, const double G0, const double gammaM);
int setRootBrackets(const int k, const int i, const double Tmin, const double Tmax, const double T, double *Ta, double *Tb);

#endif
