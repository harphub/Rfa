#include "rfa.h"
// Read big-endian 64bit integers.
// unfortunately, INT8 doesn't work like INT4 because you only
// get automatic promotion to int, not to long or more.
// if you cast every byte to long long before shifting, 
// the macro becomes terribly long & obfuscated
// We need INT8 signed, because some values may be <0.
// probably very sub-optimal, but OK for now
// to work with sign bit, make sure you are exactly 64bit !!!
int64_t INT8(unsigned char* x){
  int64_t ll;
  int i;
  ll=0;
  for(i = 0;i<8;) ll = (ll<<8) + x[i++];
  return(ll);
}

// Reading big-endian doubles
// on a big-endian platform this is not the most efficient, I guess,
// but it is very portable
double DBL8(unsigned char * x){
  int64_t i;
  i = INT8(x);
  return(*(double*)&i);
}

// byteswap a whole vector of n numbers (integer or float)
// size is the number of bytes per number, n the number of values
void byteswap(void* data,int size,int n){
  unsigned char b;
  int i,j;
  unsigned char* buff;
  buff= (unsigned char*) data;
  for(j=0;j<n;j++){
    for(i=0;i<size/2;i++) { b=buff[i];buff[i]=buff[size-1-i];buff[size-1-i]=b;}
    buff+=size;
  }
}

// given an elliptic truncation nmsmax,nsmax and some value n
// calculate the maximum value of m
int fa_mmax(int n,int nsmax,int nmsmax){
  double nf;
  int result;
// the fabs() and 1e-10 are there in case of numerical rounding errors
// Not sure whether they are still needed in this version
  nf = ((double) n)/nsmax;
  result = (int) floor((nmsmax)*sqrt(fabs(1.-nf*nf))+1e-10);
  return(result);
}


// In Data Assimilation (backgrond error covariance calculation), some files come
// with spectral components in a different order (n <-> m)
// for the inverse operation (should you ever want it), just call with nmsmax <-> nsmax
// NEVER validated for nsmax != nmsmax
void fa_rawreorder(double* fieldin, double* fieldout,int* nmsmax, int* nsmax)
{
  int i,m,n,mmax;
  int* mbuf;

  mbuf = (int*) malloc(sizeof(int)*(*nmsmax+1));
  mbuf[0]=0;
  for(m=1;m<=*nmsmax;m++) {
    mbuf[m] = mbuf[m-1] + fa_mmax(m-1,*nmsmax,*nsmax) + 1;
//    printf("%d: %d\n",m,mbuf[m]);
  }
  
  for(n=0;n<=(*nsmax);n++){
    mmax = fa_mmax(n,*nsmax,*nmsmax);
    for(m=0;m<=mmax;m++ ){
      for(i=0;i<4;i++){
        *(fieldout++) = fieldin[4*(mbuf[m] + n) + i];
      }
    }
  }
  free(mbuf);
}


// count the number of values expected in a spectral GRIB message
// i.e. only the part that is actually compacted:  M+N > trunc
// call with sptrunc=-1 to get the total number of components in the data stream
// sptrunc=0 still means that you don't compactify the m=n=0 component!
void fa_countval_spec(int* nmsmax,int *nsmax,int* sptrunc,int*result){
  int mmax,msplit,buf,n;
  buf=0;
  for(n=0;n<= *nsmax;n++) {
    mmax = fa_mmax(n,*nsmax,*nmsmax);
    if(*sptrunc < 0) msplit= -1;
    else msplit = (n==0) ? mmax : MAX(0,*sptrunc - n) ;
    buf += 4*(mmax-msplit);
  }
  *result = buf;
}

// speed up integer powers
double fast_pow (double x, int n) {
  // found online, also e.g. in wgrib
  double result=1.;
  int sign=0;
  if (n<0) {
    n = -n;
    sign=1;
  }
  while (n) {
    if (n & 1) result *= x;
    n >>= 1;
    x *= x;
  }
  if (sign) result = 1./result;
  return result;
}
  

//The reference value (R) uses IBM single precision floating point format.
//    sAAAAAAA BBBBBBBB BBBBBBBB BBBBBBBB
//    s = sign bit, encoded as 0 means positive, 1 means negative
//    A...A = 7-bit binary integer representing the exponent/characteristic
//    B...B = 24-bit binary integer, the mantissa.
// The appropriate formula to recover the value of R is:
//    R = (-1)^s * 2^(-24) * B * 16^(A-64)
//
double IBMfloat(unsigned char*x){
  int sign, power;
  unsigned int abspower;
  long int mant;
  double value, exp;
// last 3 bytes are mantissa
  mant = INT3(x+1);
  if (mant == 0) {value=0; return(value);};

  sign = x[0] & 128; // first bit is sign bit
  power = (int) (x[0] & 127) - 64 - 6; // the other 7 bits are exponent (base 16), with 64 added
                                       // the -6 is the rescaling by 2^(-24)
//  abspower = power > 0 ? power : -power;

// calculate 16^abspower:
//  exp = 16.0;
//  value = 1.0;
//  while (abspower) {
//    if (abspower & 1) value *= exp;
//    exp *= exp;
//    abspower >>= 1;
//  }

//  if (power < 0) value = mant / value;
//  else value = mant * value;

  value = mant * fast_pow(16.0, power);
  if (sign) value = -value;
  return(value);
}
 

