#include <R.h>
#include <math.h>

// simple linear interpolation using log(p)
// No extrapolation
// p_in, p_out MUST be in ascending order!!!
void fa_interp1(double* p_in, double *v_in, int *n_in, 
                double* p_out, double *v_out, int *n_out) {

  int i,j;

  i=0;
  for (j=0; j < *n_out; j++) {
    while (p_in[i] < p_out[j] && i < *n_in) i++;
    if (i==0) v_out[j] = (p_in[i] > p_out[j]) ? NA_REAL : v_in[i] ; 
    else if (i== *n_in && p_in[i] < p_out[j]) v_out[j] = NA_REAL;
    else  v_out[j] = v_in[i-1] + (v_in[i] - v_in[i-1])/(p_in[i] - p_in[i-1]) * (p_out[j] - p_in[i-1]);
  }
}

void fa_pressures(double *A, double *B, double *pref, int *nlev, double *psurf, double *pressure) {
  int i,j;
// output is the pressure values for all hybrid levels based on the given surf pressure at 1 location
// ref: FullPos 'scientist guide' (Ryad El Katib, 2002)
  double ph1, ph0, logph1, logph0, alpha;

// ph1 is usually exactly 0, but not when we pass a subset of levels
  if (A[0]==0 && B[0]==0) {
    ph1=0;
    logph1=0;
  }
  else {
    ph1 = A[0] * *pref + *psurf * B[0];
    logph1 = log(ph1);
  }

  for (i=0; i < *nlev; i++) {
    ph0 = ph1;
    logph0 = logph1;
    ph1 = A[i+1] * *pref + *psurf * B[i+1];
    logph1 = log(ph1);

//    alpha = (ph0==0) ? 1 : 1 - ph0/(ph1-ph0)*(logph1-logph0);
    alpha = 1 - ph0/(ph1-ph0)*(logph1-logph0);

//  we should return log(pressure) !!!
//  pp[i] = ph1 * exp(-alpha);
    pressure[i] = logph1 - alpha ;
  }
}


void fa_interp2(double *A, double *B, double *pref, 
                double *psurf, int * nlev, double * v_in, int *npoints,
                double *p_out, int *n_out, double *v_out) {

  int i,j;
  double *pp;
// p_out contains logarithm of output pressures!
// ALL PRESSURES MUST BE IN ASCENDING ORDER!!!
// pp should have length nlev : it contains the (log)pressure values at one location
  pp = (double*) R_alloc(*nlev, sizeof(double));
  for (i=0; i < *npoints; i++) {
    fa_pressures(A, B, pref, nlev, psurf + i, pp);
    fa_interp1(pp, v_in + i*(*nlev), nlev, p_out, v_out + i*(*n_out), n_out);
  }
}


