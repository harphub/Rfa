// a simple linear smoother (stencil [1 2 1 | 2 4 2 | 1 2 1] )
// nx, ny: boundaries of the physical domain
// maxx, maxy : boundaries of the extended domain (typically maxx = nx + 11)
// This code is mostly "in situ"
// Disadvantage: more difficult to change the filter.

#include <math.h>
#include <R.h>

//  data is coming from R, so matrix ordering is column-major! 
#define data(i, j) data[ (i) + (j) * *maxx ]
void smooth_extension(double* data, int* nx, int* ny, int* maxx, int* maxy) {
  int i, j, k;
  int jp1,jm1,ip1,im1;
  int imin, jmin, ex, ey, iter;
  double p0, p1, p2;

  ex = *maxx - *nx ;
  ey = *maxy - *ny ;
  iter = floor( ((ex>ey?ex:ey)+1.)/2.) ;

  for (k = 0 ; k < iter ; k++) {
    if (*nx < *maxx) {
    for (i = *nx + k ; i < *maxx - k ; i++) {
      p0 = data(i,0) ;
      p1 = data(i,*maxy-1) ;
      for (j=0 ; j < *maxy -1 ; j++) {
        p2 = data(i, j);
        data(i, j) += .5*(p1 + data(i, j+1)) ;
        p1 = p2;
      }
      data(i, *maxy-1) += .5*(p1 + p0) ;
    }
    for (j=0 ; j < *maxy ; j++) {
      jm1 = (j == 0) ? *maxy - 1 : j - 1 ;
      jp1 = (j == *maxy - 1) ? 0 : j + 1 ;
      ip1 = (k==0) ? 0 : *maxx - k ;

      p1 = data(*nx+k-1 , j) + 0.5*(data(*nx+k-1, jm1) + data(*nx+k-1, jp1)) ;
      p0 = data(ip1, j) + 0.5*(data(ip1, jm1) + data(ip1, jp1));

      for (i = *nx + k ; i < *maxx - k - 1; i++) {
        p2 = data(i,j);
        data(i, j) += .5*(p1 + data(i+1, j)) ;
        data(i, j) /= 4. ;
        p1 = p2;
      }
      data(*maxx - k -1, j) += .5*(p1 + p0);
      data(*maxx - k -1, j) /= 4. ;
    }
    }

    if (*ny < *maxy) {
    for (j = *ny + k ; j < *maxy - k ; j++) {
      p0 = data(0, j) ;
      p1 = data(*maxx-1, j) ;
      for (i=0 ; i < *maxx -1 ; i++) {
        p2 = data(i, j);
        data(i, j) += .5*(p1 + data(i+1, j)) ;
        p1 = p2;
      }
      data(*maxx-1, j) += .5*(p1 + p0) ;
    }
    for (i=0 ; i < *maxx ; i++) {
      im1 = (i == 0) ? *maxx - 1 : i - 1 ;
      ip1 = (i == *maxx - 1) ? 0 : i + 1 ;
      jp1 = (k==0) ? 0 : *maxy - k ;

      p1 = data(i, *ny+k-1 ) + 0.5*(data(im1,*ny+k-1) + data(ip1,*ny+k-1)) ;
      p0 = data(i, jp1) + 0.5*(data(im1, jp1) + data(ip1, jp1));

      for (j = *ny + k ; j < *maxy - k - 1; j++) {
        p2 = data(i,j);
        data(i, j) += .5*(p1 + data(i, j+1)) ;
        data(i, j) /= 4. ;
        p1 = p2;
      }
      data(i,*maxy - k -1) += .5*(p1 + p0);
      data(i,*maxy - k -1) /= 4. ;
    } }
  }
}
