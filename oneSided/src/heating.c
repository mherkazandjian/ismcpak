#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/* --> Heating efficiency, eq-5.4 Rowin Thesis pp 53 */
double heatingEfficiency(const double T, const double AvL, const double eldens)
{
    double epsilon, factor1, factor2, ratio, G0p;
    
    G0p   = pv.G0 * exp(-K_UV * AvL);   /* local radiation feild attenuated by dust */
    ratio = G0p * sqrt(T) / eldens; /* ration of ionization and recombination   */
    
    factor1 = 4.87e-2 / (1e0 + 4e-3 * pow(ratio, 0.73));
    factor2 = 3.65e-2 * pow((T/1e4), 0.7)/(1e0 + 2e-4 * ratio );
    epsilon = factor1 + factor2;
    return epsilon;
}
/* --> Calculate Heating done by very small graphitic grains and PAH's */
double heating_PAH(mesh *msh, const double T, const double AvL, const int k, const double elecDens) 
{
    double Heating, epsilon;
    
    epsilon  = heatingEfficiency(T, AvL, elecDens);
    Heating  = 1e-24 * epsilon * pv.G0 * exp(-K_UV * AvL) * pv.dens0;   /* eq-5.2 rowin thesis pp 52 */
    Heating *= pv.metalicity;
    
    msh->heating.photo[k] = Heating; /* zero indexing */
    return Heating;
}
/* --> Calculates the heating due to carbon ionization */
double heating_carbon(mesh *msh, const double AvL, const int k, const double carbonDens) 
{
  double b, v1, gammaCarbon, tauC, tauH2;
  double factor1, factor2, factor3;

  if(k==0) {
    tauH2 = tauC = 0.0;
  } else {
    tauH2 = 1.2e-14*msh->NH2All_L[k]/DELTAV;     /* eq-5.9 rowin thesis pp 53 */
    tauC  = 1.1e-17*msh->NC_L[k];                /* eq-5.8 rowin thesis pp 53 */
  }

  v1    = 5e2/DELTAV;                       /* eq-5.11 rowin thesis, pp 53 */
  b     = 9.2e-3/DELTAV;                    /* eq-5.10 rowin thesis, pp 53 */
  
  factor1 = 2.79e-22 * carbonDens * pv.G0;
  factor2 = exp(-2.4 * AvL - tauC - tauH2 * b/ (M_PI * pow(v1, 2e0)));
  factor3 = pow(1e0 + tauH2 * b / (M_PI * pow(v1, 2e0)), -1e0);

  gammaCarbon = factor1* factor2 * factor3;  /* eq-5.6 rowin thesis, pp 53 */
  
  //printf("%e %e\n",gammaCarbon, carbonDens);
  msh->heating.cIon[k] = gammaCarbon; /* zero indexing */
  return gammaCarbon;
}
/* --> Calcalates the heating due to H2 photodissociation */
double heating_molecular_hydrogen(mesh *msh, const double AvL, const int k, const double H2Dens) 
{
    double gamma_H2, betatau;
  
    betatau = Self_shlding( msh->NH2All_L[k] );
    gamma_H2 = 2.23e-23 * H2Dens * pv.G0 * exp(-2.5 * AvL) * betatau;  /* eq-5.14 rowin thesis pp 54 */
    
    msh->heating.molHydro[k] = gamma_H2; /* zero indexing */
    return gamma_H2;
}
/* --> Calculates heating due to gas-grain collisions */
/*     Hollenbach & McKee 1989                        */
double heating_gas_grain_collisions(mesh *msh, const double T, const double Td, const int k) 
{
    double gamma_gas_grain;
    double factor1, factor2;
    
    factor1 = 1.2e-31 * SQ(pv.dens0) * sqrt(T/1e3) * sqrt(100e0 / A_MIN);
    factor2 = (1e0 - 0.8e0 * exp(-75e0/ T)) * (Td - T);
    gamma_gas_grain =  factor1 * factor2;  /* eq-5.18 rowin thesis pp 55   */
    gamma_gas_grain *= pv.metalicity;      
    
    msh->heating.ggColl[k] = gamma_gas_grain; /* zero indexing */
    return gamma_gas_grain;
}

/* --> Calculates viscous heating                    */
double heating_viscous(mesh *msh, const double T, const double AvL, const int k, const double elecDens, const double CIIDens) 
{
    double lambda, vth_e, vth_c, v_d, Zd, x, gamma, gamma_viscous, nd;
    double factor1, factor2;
    
    v_d   = 1e2;                        /* grain local drift velocity */
    vth_e = sqrt(3e0 * KB * T / m_e);   /* thermal speed of elections */
    vth_c = sqrt(3e0 * KB * T / m_c);   /* thermal speed of C+        */
    nd    = pv.dens0 * 1.9e-8;          /* dust density               */
    
    /* computing the grain charge */
    gamma = 0.1 * (pv.G0/pv.dens0) * sqrt(T) * exp(-K_UV * AvL);
    x = 1e0 - 0.28 / (gamma + 0.5);
    Zd = (x - 0.44) * NU_H * ERG * DUSTSIZE * pow(ESU, -2e0);   /* grain charge */

    /* computing the viscous heating */
    lambda = 1.5 * pow(Zd, -1e0) * pow(ESU,-3e0) * pow(KB * T, 1.5e0) * pow(M_PI * elecDens, -0.5e0); /* eq-5.21 rowin thesis pp 55 */
    factor1 = 8e0 * M_PI * pow(ESU, 4e0) * nd * SQ(Zd) * pow(KB * T, -1e0) * log(dabs(lambda)) * v_d;
    factor2 = CIIDens * Gfunct(v_d/vth_c) + elecDens * Gfunct(v_d/vth_e);
    gamma_viscous = factor1 * factor2; 
    gamma_viscous *= pv.metalicity; 
    
    msh->heating.visc[k] = gamma_viscous;
    return gamma_viscous;
}
/* --> Gfunct eq-5.22 rowin thesis pp 55 */
double Gfunct(const double y) 
{
    double Gfunct;
    Gfunct = 1e0 / (2e0 * pow(y, 2e0)) * (erf(y) - 2e0 / sqrt(M_PI) * y * exp(-pow(y, 2e0)));
    return Gfunct;
}
/* --> Heating by cosmic rays */
double heating_cosmic_rays(mesh *msh, const int k, const double H2Dens) 
{
    double cosmic_ray_heating;

    cosmic_ray_heating = 1.5e-11 * pv.zeta * H2Dens;  /* eq-5.24 rowin thesis pp 56 */
    msh->heating.CR[k] = cosmic_ray_heating;
    return cosmic_ray_heating;
}
/* --> Calculates the heating due the pumping of H2 molecules */
/*     also known as H_2 collisional de-excitation heating    */
double heating_H2_pumping(mesh *msh, const double T, const int k, const double HDens, const double H2Dens, const double H2exDens) 
{
    double Eexcited, gamma_h2_pump, gamma_H, gamma_H2;
    
    Eexcited      = 2.6e0 * ERG;        // KB * 5860e0;
    gamma_H       = 1e-12 * sqrt(T) * exp(-1000e0 / T);                   /* eq-5.16 rowin thesis pp 54 */
    gamma_H2      = 1.4e-12 * sqrt(T) * exp(-(18100e0 / (T + 1200e0)));  /* eq-5.17 rowin thesis pp 54 */
    gamma_h2_pump = (HDens * gamma_H + H2Dens * gamma_H2) * H2exDens * Eexcited;  /* eq-5.15 rowin thesis pp 54 */
    
    msh->heating.H2pump[k] = gamma_h2_pump;
    return gamma_h2_pump;
}
