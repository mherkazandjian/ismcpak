#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/* Compute the self shlding factor for the molecules H2 and CO. */
double Self_shlding(double NH2All) 
{
  const double tau=1.2e-14*NH2All/DELTAV;     /* eq-5.9 rowin thesis pp 53 */
  
  if (tau <= 10e0) 
    return small_tau(tau);
  else 
    return big_tau(tau);
}
/* Self shlding factor for tau <= 10, eq-5.54 rowin thesis pp 63 */
double small_tau(double tau) 
{
    int n;  
    double nfacul, beta; 
    
    beta = 0e0; 
    nfacul = 1e0; 
    n=0; 
    beta = pow(-1e0,(double)n) * pow(tau,(double)n) /  
	(nfacul * sqrt((double)n+1e0) * pow(M_PI, (double)n/2e0)); 
    for (n=1; n < 100; n++) { 
	nfacul = nfacul * (double)n;   
	beta = beta + pow(-1e0,(double)n) * pow(tau,(double)n) /  
	    (nfacul * sqrt((double)n+1e0) * pow(M_PI, (double)n/2e0));    
    }
    return beta; 
}
/* Self shlding factor for tau > 10, eq-5.54 rowin thesis pp 63 */
double big_tau(double tau) {

  double bb, v1, beta; 

  bb = 9.2e-3 / DELTAV; 
  v1 = 5.0e2  / DELTAV; 

  beta = (pow(tau,-1e0) * pow(log(tau / sqrt(M_PI)),-5.0e-1) +  
	  pow((bb / tau), 5.0e-1)) * 
    erfc(sqrt(tau * bb / (M_PI * pow(v1,2e0))));  

  return beta;
}
/* Self-shlding for CO is retreived through interpolation in table 5 of Van Dishoeck & Black, 1988
 * table-5. The self sheilding is returned based on which istope of CO is provided by the what_is
 *  parameter. NCO and NH2 are the column densities (not in log10) of CO (which isotope????) and H2 respectively.
 */
double self_shlding_CO(double NCO, double NH2, int what_iso) 
{
  int i, row, column;
  double fac1, fac2, dumshld1, dumshld2, beta;

  row = column = 1;
  fac1 = fac2 = 0.0;
  
  for (i=1; i < self_shld[_12CO].npointCO; i++) {
      if (log10(NCO) >= self_shld[_12CO].logNCO[i-1] && log10(NCO) < self_shld[_12CO].logNCO[i]) {
	  column = i;
	  fac1 = (log10(NCO) - self_shld[_12CO].logNCO[i-1]) /
	      (self_shld[_12CO].logNCO[i] - self_shld[_12CO].logNCO[i-1]);
      }
  }
  
  for (i=1; i < self_shld[_12CO].npointH2; i++) {
    if (log10(NH2) >= self_shld[_12CO].logNH2[i-1] && log10(NH2) < self_shld[_12CO].logNH2[i]) {
      row = i;
      fac2 = (log10(NH2) - self_shld[_12CO].logNH2[i-1]) /
	(self_shld[_12CO].logNH2[i] - self_shld[_12CO].logNH2[i-1]);
    }
  }
  
  dumshld1 = self_shld[what_iso].function[column-1][row-1] + 
    fac1 * (self_shld[what_iso].function[column][row-1] - self_shld[what_iso].function[column-1][row-1]);
  
  dumshld2 = self_shld[what_iso].function[column-1][row] +
    fac1 * (self_shld[what_iso].function[column][row] - self_shld[what_iso].function[column-1][row]);

  beta = dumshld1 + fac2 * (dumshld2 - dumshld1);
  
  return beta;
}
