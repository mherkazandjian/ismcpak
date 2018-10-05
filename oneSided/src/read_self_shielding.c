#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/******************** Read selfshlding data for CO, 13CO, C18O, 13CO18O *******************/
/*  Each file should have the following structure :
 *  line1 : n horizontal gridPoints 'nc' in N(X^C Y^O) (number of columns in the table) followed by the grid values for the columns.
 *  line2 : n vertical gridpoints 'nr' in N(H2) (number of rows in the table) followed by the grid values for the row.
 *  TABLE (nc columns x nr rows) of tabulated self-sheilding values.
 */
void read_selfshlding_CO(void) 
{
  int i,j,k;
  char fName[256]; /* this replaced self_shld[CO_s].name */
  
  sprintf(fName,"%s", pv.selfSheilding_CO_fName);
  //if( verbose & VERB4_MASK ) fprintf(outputFd, "     Reading %s\n", fName);
  
  //openning the file containing the tabulated self sheilding values
  if ((self_shld[_12CO].FILE=fopen(fName,"r"))==NULL) errMsg(7);
  
  fscanf(self_shld[_12CO].FILE, "%d", &self_shld[_12CO].npointCO );
  /*fprintf(outputFd, "%ld ", npointCO);*/

  //reading the N(CO) grid coordinates
  MAKE_1D_ARRAY(self_shld[_12CO].logNCO, self_shld[_12CO].npointCO, double);
  for (i=0; i < self_shld[_12CO].npointCO; i++) {
    fscanf(self_shld[_12CO].FILE, "%le", &self_shld[_12CO].logNCO[i]);
    /*fprintf(outputFd, "%le ", self_shld[_12CO].logNCO[i]);*/
  }
  /*fprintf(outputFd,"\n"); */

  //reading the N(H2) grid coordinates
  fscanf(self_shld[_12CO].FILE, "%d", &self_shld[_12CO].npointH2);
  /*fprintf(outputFd, "%ld ", npointH2);*/
  MAKE_1D_ARRAY(self_shld[_12CO].logNH2, self_shld[_12CO].npointH2, double);
  for (i=0; i < self_shld[_12CO].npointH2; i++) {
    fscanf(self_shld[_12CO].FILE, "%le", &self_shld[_12CO].logNH2[i]);
    /*fprintf(outputFd, "%le ", self_shld[_12CO].logNH2[i]);*/
  }
  /*fprintf(outputFd,"\n"); */


  for (k=_12CO; k < MOL; k++) {

    MAKE_2D_ARRAY(self_shld[k].function, self_shld[_12CO].npointCO, self_shld[_12CO].npointH2, double);
    for (j=0; j < self_shld[_12CO].npointH2; j++) {
      for (i=0; i < self_shld[_12CO].npointCO; i++) {
	fscanf(self_shld[_12CO].FILE, "%le", &self_shld[k].function[i][j]);
	/*fprintf(outputFd, "%le ", self_shld[_12CO].function[i][j]);*/
      }
      /*fprintf(outputFd,"\n"); */
    }
    /*fprintf(outputFd,"\n"); */
  }
}
/********************* Read rotational cooling data ********************/
void read_rotational_data(void) 
{
  char fName[256]; /* this replaced coollines[coolspec].name*/
  int coolspec,i,j;

  for (coolspec = H2; coolspec <= H2_low_T; coolspec++) {

    sprintf(fName,"%s%d%s", pv.rotationalCooling_baseName, coolspec,".inp");  
    //fprintf(outputFd, "     Reading %s\n", fName);

    if ((coollines[coolspec].FILE=fopen(fName,"r"))==NULL) errMsg(8);

    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].rows_rot);
    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].columns_rot);
    MAKE_2D_ARRAY(coollines[coolspec].paramH2_rot,coollines[coolspec].columns_rot,
		  coollines[coolspec].rows_rot, double);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      for (j=0; j < coollines[coolspec].columns_rot; j++) {
	fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].paramH2_rot[j][i]);	
      }
    }
  }
 
  for (coolspec=H2O; coolspec < COOLSPEC; coolspec++) {
    sprintf(fName,"%s%d%s", pv.rotationalCooling_baseName, coolspec,".inp");  
    /*    fprintf(outputFd, "Reading %s\n", fName); */
    if ((coollines[coolspec].FILE=fopen(fName,"r"))==NULL) errMsg(9);
    
    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].rows_rot);
    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].columns_rot);
    MAKE_1D_ARRAY(coollines[coolspec].ColumnDens_rot,coollines[coolspec].rows_rot,double);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].ColumnDens_rot[i]);
    }
    MAKE_1D_ARRAY(coollines[coolspec].Temperature_rot,coollines[coolspec].columns_rot,double);
    for (i=0; i < coollines[coolspec].columns_rot; i++) {
      fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].Temperature_rot[i]);
    }
    MAKE_1D_ARRAY(coollines[coolspec].L0_rot,coollines[coolspec].columns_rot,double);
    for (i=0; i < coollines[coolspec].columns_rot; i++) {
      fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].L0_rot[i]);
    }
    MAKE_2D_ARRAY(coollines[coolspec].Llte_rot,coollines[coolspec].columns_rot,
		  coollines[coolspec].rows_rot,double);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      for (j=0; j < coollines[coolspec].columns_rot; j++) {
	fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].Llte_rot[j][i]);	
      }
    }
    MAKE_2D_ARRAY(coollines[coolspec].nhalf_rot,coollines[coolspec].columns_rot,
		  coollines[coolspec].rows_rot,double);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      for (j=0; j < coollines[coolspec].columns_rot; j++) {
	fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].nhalf_rot[j][i]);	
      }
    }
    MAKE_2D_ARRAY(coollines[coolspec].alpha_rot,coollines[coolspec].columns_rot,
		  coollines[coolspec].rows_rot,double);
    for (i=0; i < coollines[coolspec].rows_rot; i++) {
      for (j=0; j < coollines[coolspec].columns_rot; j++) {
	fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].alpha_rot[j][i]);	
      }
    }
  }
}
/********************* Read vibrational cooling data *******************/
void read_vibrational_data(void) 
{
  char fName[256]; /* this replaced coollines[coolspec].name*/
  int coolspec,i,j;

  for (coolspec=H2O; coolspec < COOLSPEC; coolspec = coolspec + 3) {

    sprintf(fName,"%s%d%s", pv.vibrationalCooling_baseName, coolspec,".inp");  
    //fprintf(outputFd, "     Reading %s\n", fName);
    if ((coollines[coolspec].FILE=fopen(fName,"r"))==NULL) errMsg(10);

    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].rows_vib);
    fscanf(coollines[coolspec].FILE, "%d", &coollines[coolspec].columns_vib);
    MAKE_1D_ARRAY(coollines[coolspec].ColumnDens_vib,coollines[coolspec].rows_vib,double);
    for (i=0; i < coollines[coolspec].rows_vib; i++) {
      fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].ColumnDens_vib[i]);
    }
    MAKE_1D_ARRAY(coollines[coolspec].Temperature_vib,coollines[coolspec].columns_vib,double);
    for (i=0; i < coollines[coolspec].columns_vib; i++) {
      fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].Temperature_vib[i]);
    }
    MAKE_2D_ARRAY(coollines[coolspec].Llte_vib,coollines[coolspec].columns_vib,coollines[coolspec].rows_vib,double);
    for (i=0; i < coollines[coolspec].rows_vib; i++) {
      for (j=0; j < coollines[coolspec].columns_vib; j++) {
	fscanf(coollines[coolspec].FILE, "%le", &coollines[coolspec].Llte_vib[j][i]);	
      }
    }
  }
}
