#include "prototypes.h"  /* function headers */
#include "vars.h"        /* global variables */

/*-------------------- NR-2 LU SOLVER ----------------------*/
void solve_lu(double **A, double *b, double *x, int n) 
{
    int i,j, *indx;
    double *bl, **Al, *res;

    for(i=0;i<=n;i++) A[i][0]=A[0][i]=0.0;

    /* making local copies of the linear system */
    MAKE_1D_ARRAY(bl, n+1, double);
    for(i=1;i<=n;i++) bl[i]=b[i];
    MAKE_2D_ARRAY(Al, n+1, n+1, double);
    for(i=1;i<=n;i++)
	for(j=1;j<=n;j++)
	    Al[i][j]=A[i][j];
    
    /* allocating workspace */
    MAKE_1D_ARRAY(indx, n+1, int);
    MAKE_1D_ARRAY(res , n+1, double);
    
    /* solving with LU */
    TIC(tmr2);
    ludcmp(Al, n, indx); CHECK_ERRV;
    TOC(tmr2);
    lubksb(Al, bl, n, indx);
    mprove(A, Al, n, indx, b, bl, 1);
    

    for(i=1;i<=n;i++) x[i]=bl[i];

    /* freeing workspace */
    FREE_2D_ARRAY(Al, n+1, n+1, double);
    FREE_1D_ARRAY(indx, n+1, double);
    FREE_1D_ARRAY(res , n+1, double);
    FREE_1D_ARRAY(bl  , n+1, double);
}
/********************* LU decompositie *********************************/
/* --------- get the absolute max of a row --------- */
double maxRow( double **A, int row, int n)
{
  int i;
  double temp, big=0.0;
  
  for (i=1;i<=n;i++) 
    if ((temp=dabs(A[row][i])) > big) 
      big=temp;
  
  if (big == 0.0) { /* big must not be zero which mean the whole row is zeros */
    errMsg(3);
    CHECK_ERR;
  }
  else
    return big;

  errMsg(2); /* must not reach here */
  return 0; 
}
/* --------- get the pivot scales --------- */
void getPivots( double **A, double *pv, int n)
{
  int i;
  for (i=1;i<=n;i++) {
    pv[i]=1.0/maxRow(A,i,n); CHECK_ERRV;
  }
}
/* --------- reduce ------ */
double reduceElement( double **A, int i, int j, int n, int l)
{
  int k;
  double sum;
  
  sum=A[i][j];
  for (k=1;k<l;k++)
    sum -= A[i][k]*A[k][j];
  return sum;
}
double reduceElement2( double **A, int i, int j, int n)
{
  int k, l;
  double sum;
  
  sum=A[i][j];
  l=min(i,j)-1;
  for (k=1;k<=l;k++) sum -= A[i][k]*A[k][j];
  return sum;
} 
/* swapping macro, replace a with b and b with a */
#define SWAP(a,b){\
    double tmp;	  \
    tmp=a;\
    a=b;\
    b=tmp;\
}
/* get the max colon pivot ?? */
int getColonPivot(double **A, double *vv, int j, int n)
{
    double big, tmp;
    int i, imax;
    
    imax=1.0;
    big=0.0; /* 3 - elements on and below giagonal */
    for (i=j;i<=n;i++) {
	tmp=vv[i]*dabs(A[i][j]);
	//printf("%-+10.5f ", tmp);
	if ( tmp >= big) { 
	    big=tmp;
	    imax=i;
	}
    }
    return(imax);
}
/* swap the rows with pivot row?? */
void swapRow(double **A, double *vv, int j, int imax, int n)
{
    int k;
    
    if (j != imax)  {  /* 4 - put pivot max on the diagonal */
	for (k=1;k<=n;k++)
	    SWAP(A[imax][k],A[j][k]);
	vv[imax]=vv[j];
    }
}
/* scale the colon with the diagonal element */
void scaleColon(double **A, int j, int n)
{
    int i;
    double tmp;
    
    if (j != n) {
      tmp=1.0/(A[j][j]); /* diagonal element */
      for (i=j+1;i<=n;i++) A[i][j] *= tmp;
    }
}
/* -------- main LU routine --------- */
void ludcmp(double **A, int n, int *indx) 
{
  int i, j;
  double *vv;
  
  MAKE_1D_ARRAY(vv,n+1,double);
  getPivots(A, vv, n);                                     /* 1 */
  
  for (j=1;j<=n;j++) {
    for (i=1;i<=n;i++) A[i][j]=reduceElement2(A, i, j, n); /* 2 */
    indx[j] = getColonPivot(A, vv, j, n); CHECK_ERRV;      /* 3 */
    swapRow(A, vv, j, indx[j], n);                         /* 4 */
    if (A[j][j] == 0.0) A[j][j]=TINY;                      /* 5 */
    scaleColon(A, j, n);                                   /* 6 */
  }
  FREE_1D_ARRAY(vv,n+1,double); 
}
/********************* LU terugsubstitutie *********************************/
void lubksb(double **A, double *b, int n, int *indx) 
{ 
  int i,ii=0,ip,j;
  double sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= A[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= A[i][j]*b[j];
    b[i]=sum/A[i][i];
  }
}
/********************* iterative refinenement **********************************/
void mprove(double **A, double **Al, int n, int *indx, double *b, double *bl, int nImprove)
{
    int i, j, k;
    double sdp, *res;
    
    MAKE_1D_ARRAY(res , n+1, double);
    
    for(k=0; k<nImprove; k++) {
	for (i=1; i<=n; i++) {
	    sdp = -b[i];
	    for (j=1; j<=n; j++) sdp += A[i][j]*bl[j];
	    res[i]=sdp;
	}
	lubksb(Al, res, n, indx);     //my version
	
	for (i=1; i<=n; i++) bl[i] -= res[i];
    }

    FREE_1D_ARRAY(res , n+1, double);
}
