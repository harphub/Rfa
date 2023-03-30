#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<inttypes.h>
#include<R.h>

// Some utilities

#define MAX(x,y) (x>y ? x : y)

// a simple test for big/little endian platforms
//#define LITTLE_ENDIAN ( *(uint16_t*)"a" < 255)

// macros for reading big-endian integers from a bit stream
// platform-independent (32bit int)
// up to four bytes, we can rely on automatic integer promotion
#define INT2(x) ( ((x)[0]<<8) + (x)[1])
#define INT3(x) ( ((x)[0]<<16) + ((x)[1]<<8) + (x)[2])
#define INT4(x) ( ((x)[0]<<24) + ((x)[1]<<16) + ((x)[2]<<8) + (x)[3])

// a few prototypes
int64_t INT8(unsigned char* x);
double DBL8(unsigned char * x);
void byteswap(void* data,int size,int n);
int fa_mmax(int n,int nsmax,int nmsmax);
void fa_countval_spec(int* nmsmax,int *nsmax,int* sptrunc,int*result);
double fast_pow(double x, int n);
double IBMfloat(unsigned char*x);
