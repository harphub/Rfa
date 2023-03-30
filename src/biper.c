// a simple cubic spline extension
// given f(0), f''(0), f(1), f''(1), cspline is defined by
// f(x) = f(0) + [f(1)-f(0)-(f''(1)+2f''(0))/6] x +
//        [f''(0)/2] x^2 + [(f''(1) - f''(0))/6 ] x^3

// here, spline from p[nx] to p[maxx+1]=p[1]
//       f'(0) ~ p[nx]-p[nx-1]
//       f''(0) ~ 

// nx, ny: boundaries of the physical domain
// maxx, maxy : boundaries of the extended domain (typically maxx = nx + 11)
//
// the simplest estimate of f''
// p'[nx] ~ p[nx] - p[nx-1]
// p''[nx] ~ (p[maxx]-p[nx])/(maxx-nx) - (p[nx] - p[nx-1])

#include <math.h>
#include <R.h>

void fit_spline(double p[4], double a[4], int np) {
  double ext, de, le, da, db, ddp1, ddp2 ;

  ext = np + 1. ;
  de = ext / (ext + 1.) ;
  le = 4. - de*de ;
    // our "biper" zone is between points p1 and p2
    // in future, we may improve p'' estimate by using more points
  da = ( (p[2] - p[1])/ext - (p[1] - p[0])) * 6./(ext + 1) ;
  db = ( (p[3] - p[2]) - (p[2] - p[1])/ext) * 6./(ext + 1) ;

  ddp1 = (2.*da - de*db) / le ;
  ddp2 = (2.*db - de*da) / le ;

  a[0] = p[1] ;
  a[1] = (p[2] - p[1])/ext - (ddp2 + 2.*ddp1)*ext/6.;
  a[2] = ddp1/2. ;
  a[3] = (ddp2 - ddp1)/6./ext ;
}

//  data is coming from R, so matrix ordering is column-major! 
#define data(i, j) data[ (i) + (j) * *maxx ]
void biper(double* data, int* nx, int* ny, int* maxx, int* maxy) {
  int i, j, k;
  double a[4], p[4] ;

  if (*maxy > *ny) {
  for (i=0 ; i < *nx ; i++) {
    p[0] = data(i, *ny - 2) ;
    p[1] = data(i, *ny - 1) ;
    p[2] = data(i, 0) ;
    p[3] = data(i, 1) ;
    fit_spline(p, a, *maxy - *ny);

    for (j = 1 ; j <= *maxy - *ny ; j++) {
      // use (j/ext) ?
      data(i,*ny - 1 + j) = a[0] + j * (a[1] + j * ( a[2] +  j * a[3] )) ;
    }
  }
  }

  if (*maxx > *nx) {
  for (j=0 ; j < *maxy ; j++) {
    p[0] = data(*nx - 2, j) ;
    p[1] = data(*nx - 1, j) ;
    p[2] = data(0, j) ;
    p[3] = data(1, j) ;
    fit_spline(p, a, *maxx - *nx);

    for (i=1; i <= *maxx - *nx ; i++) {
      data(*nx - 1 + i, j) = a[0] + i * ( a[1] + i * ( a[2] + i * a[3] )) ;
    }
  }
  }

}
