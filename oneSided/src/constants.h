#ifndef __CONSTANTS_H
#define __CONSTANTS_H

/* some constants */
#ifndef M_PI        
   #define M_PI       3.1415926536
#endif

#define MSOLAR       1.9892e33         /* solar mass in CGS                                                          */
#define RSOLAR       6.95997e+10       /* solar radius in CGS                                                        */
#define YEAR         3.1558149984e7    /* year in seconds in CGS                                                     */
#define KB           1.381e-16         /* boltzmann constant in cgs ( ergs / k )                                     */
#define HP           6.62618e-27       /* planck's constant in cgs                                                   */
#define CLIGHT       2.997925e10       /* speed of light in CGS                                                      */
#define K_UV         1.8e0             /* fitting constant attenuation of G0 due to dust, eq-5.2, rowin thesis pp 52 */
#define ESU          4.0832e-10        /* election charge in CGS                                                     */
#define m_c          1.99442e-23       /* mass of carbon nucleus                                                     */
#define m_e          9.1095e-28        /* electom mass in cgs                                                        */
#define A_MIN        1e1               /* minium radius of a dust grain in Angstroms                                 */
#define ERG          1.60206e-12       /* 1 eV in cgs                                                                */
#define NU_H         13.6e0
#define DELTAV       1.0e0             /* terbulent velocity dispersion (original value was 2.7)                     */
#define DUSTSIZE     1e-6
#define TEFF         3.0e5             /* effective temperature of the black body radiation field. Usually this       */
                                       /* has an effect only on the dust, so it is set such that it has no            */
                                       /* effect ( see note 5 on page 7 notebook 2)                                   */
#define NU0          3e15              /* frequency beyond which dust absorbtion efficiency is 1                      */
#define ALBEDO       5.0e-1            /* ratio of reflected to incident light ( for dust grains i guess )            */
#define PHI_PAH      0.5               /* scaling factor for the collision rates of ions and electrons                */

#define DEPSILON     1.0e-60           /* numerical thresholds                                                        */
#define SMALL        1.0e-30
#define TINY         1.0e-25           /* original value was 1e-20                                                    */
#define TINY2        1.0e-25           /* original value was 1e-25                                                    */
#define MIN_DENS     1.0e-25

#define MAXPOW       14                /* maximum power ??                                                            */
#define LINT         32                /* number of bits in an integer                                                */

#define MAX_REFINE_ROOTS 2000
#endif
