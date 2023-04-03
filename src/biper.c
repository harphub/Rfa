// a simple cubic spline extension
// from p[nx] to p[maxx+1]=p[1]

// nx, ny: boundaries of the physical domain
// maxx, maxy : boundaries of the extended domain (typically maxx = nx + 11)
//
// boundary conditions: either S''=0 in p(nx-2) and p(2)
// or calculate S'' by finite differences S''(2) = S(3) - 2S(2) + S(1)

#include <math.h>
#include <R.h>

void fit_spline(double p[6], double a[4], int np, int bc) {
  double M, le, da, db, ddp1, ddp2, ddp3, ddp4 ;

  // distance between p[2] and p[3] is np + 1:
  M = np + 1. ;
  if ( ! bc ) {
    ddp1 = ddp4 = 0. ;
  }
  else {
    // use 2 extra points to calculate second derivatives
    ddp1 = p[0] - 2*p[1] + p[2] ;
    ddp4 = p[3] - 2*p[4] + p[5] ;
  }
  da = 3./(3.*M + 2.) * ( (p[4] - p[3]) - (p[2] - p[1]) - (ddp1 + ddp4)/6. ) ;
  db = 3./(M + 2.)  *   (-(p[4] - p[3]) - (p[2] - p[1]) - (ddp1 - ddp4)/6. + 2.*(p[3] - p[2])/M ) ;

  ddp2 =  (da + db) ;
  ddp3 =  (da - db) ;

  a[0] = p[2] ;
  a[1] = (p[3] - p[2])/M - (ddp3 + 2.*ddp2)*M/6.;
  a[2] = 0.5 * ddp2 ;
  a[3] = (ddp3 - ddp2)/6./M ;
}

//  data is coming from R, so matrix ordering is column-major! 
#define data(i, j) data[ (i) + (j) * *maxx ]
void biper(double* data, int* nx, int* ny, int* maxx, int* maxy, int* bc) {
  // nx, ny : original dimensions (must be >= 3)
  // maxx, maxy : new dimensions (initialised to zero)
  // bc : boundary condition type (0 : S''=0, 1 : S'' = from finite difference)
  int i, j;
  double a[4], p[6] ;

  if (*maxy > *ny) {
  for (i=0 ; i < *nx ; i++) {
    // we pass 3 points at every "side"
    // if we use S''=0 we don't need p0 and p5
    p[0] = data(i, *ny - 3) ;
    p[1] = data(i, *ny - 2) ;
    p[2] = data(i, *ny - 1) ;
    p[3] = data(i, 0) ;
    p[4] = data(i, 1) ;
    p[5] = data(i, 2) ;
    fit_spline(p, a, *maxy - *ny, *bc);

    for (j = 1 ; j <= *maxy - *ny ; j++) {
      data(i,*ny - 1 + j) = a[0] + j * (a[1] + j * ( a[2] +  j * a[3] )) ;
    }
  }
  }

  if (*maxx > *nx) {
  for (j=0 ; j < *maxy ; j++) {
    p[0] = data(*nx - 3, j) ;
    p[1] = data(*nx - 2, j) ;
    p[2] = data(*nx - 1, j) ;
    p[3] = data(0, j) ;
    p[4] = data(1, j) ;
    p[5] = data(2, j) ;
    fit_spline(p, a, *maxx - *nx, *bc);

    for (i=1; i <= *maxx - *nx ; i++) {
      data(*nx - 1 + i, j) = a[0] + i * ( a[1] + i * ( a[2] + i * a[3] )) ;
    }
  }
  }

}
