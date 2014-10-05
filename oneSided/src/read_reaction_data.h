/* species stuff */
void readDataHeaderSpecies( void );
void readInitialAbundanceSpecies( void );

void initSpecies( void );
   void readSpecie( int );
   void setup_specCode( int );
         
/* under-ubdundant stuff */
void readDataHeaderUnderUbundant( void );
void readUnderAbundantData( void );
    
/* chem reaction rates stuff */
void initReactionRateBuffers( void );
void readReactionRates( void );
   void parseReactionRateLine( int j);
void readSpecieFld( int j, int idx, int fldWidth);
   double getSpecieIndex(char SK[10]);
