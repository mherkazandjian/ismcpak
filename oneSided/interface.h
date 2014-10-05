/*-----allocate the pointer array to the mesh structs  ------*/
int initialize();
int add_mesh(int *id, double rho, double T0, double G0_FUV, double Lmech);
//int get_temperature(double *, int);

int setInitialGuessAbundances( );
int evolve_mesh( int i);
int remove_mesh( int id );
int calc_equilibrium();
int isMyMesh( int id );
/* seter functions which are otherwise read from the parameter file */
int set_outputDir( char *dirName );
int set_species_fName( char *fName );
int set_underUbundant_fName( char *fName );
int set_rate99_fName( char *fName );
int set_selfSheilding_CO_fName( char *fName );
int set_rotationalCooling_baseName( char *baseName );
int set_vibrationalCooling_baseName( char *baseName );

int set_zeta( double zeta );
int set_S_depletion( double S_depletion );
int set_TTol( double TTol );
int set_CTol( double CTol );
int set_metalicity( double metalicity );
int set_AvMax( double AvMax );
int set_slabSizeCrit( double slabSizeCrit );
int set_max_deltaAv( double max_deltaAv );
int set_min_deltaAv( double min_deltaAv );
int set_maxSlabs( int maxSlabs );
