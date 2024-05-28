#include "rfa.h"
// split a stream of bits (data section of a grib message) into values
// actually integer-valued, but cast to double for convenience...
void fa_grib_expand(unsigned char*inbuf,int nbits, int nval, double *fieldout,
    double minval, double maxval, double scale){
  int i,bbits,intbuf;
  int bitmask,pp1,pp2;
  long int buff;
//  double scale;

  bbits = 0;
  pp1 = bitmask = (1<< nbits) -1;
  pp2 = (1<<(nbits-1));
//  if (isnan(scale)) scale=(maxval-minval)/bitmask;
#ifdef DEBUG
  Rprintf("grib-expander:\n");
  Rprintf("nbits=%i,nval=%i\n",nbits,nval);
  Rprintf("minval=%lf, maxval=%lf, scale=%lf\n",minval,maxval,scale);
#endif

  buff = 0;
  for(i=0;i < nval;i++){
    while (bbits < nbits) { buff = (buff<<8) + *inbuf++ ; bbits += 8;}
    bbits -= nbits;
    intbuf = (buff>>bbits) & bitmask;
// the following makes sure that you get exactly minval and maxval
// without any IEEE differences
    if(intbuf < pp2 || isnan(maxval)) fieldout[i] = minval + scale*intbuf;
    else fieldout[i] = maxval - scale*(pp1 - intbuf);
// for non-arpege style (maxval undefined)
//  fieldout[i] = minval + scale * intbuf;
  }
}

void fa_spectral_combine(double* data1,double* data2,double* data,
                         int nsmax,int nmsmax,
                         int sptrunc,int sppow,int* ERR){
// combine the grib-compacted (high wave number) stream data1 with the non-compacted data2
// data1 is scaled by a power of the laplacian
  int i,n,m,mmax,msplit;
  double LAP1,LAP;

  *ERR=1;
  for(n=0;n<= nsmax;n++){
    mmax = fa_mmax(n,nsmax,nmsmax);
    msplit = (n==0) ? mmax : MAX(0,sptrunc - n) ;
    m=0;
    while(m <= msplit) {
      for(i=0;i<4;i++) *(data++) = *(data2++);
      m++;
    }
    while(m <= mmax) {
// LAP(m,n) = (m^2 + n^2)^sppow
// TO DO: check in original code whether sppow>2 is even possible...
// We define using 1./LAP for 1 single reason:
//   it comes out numerically closer to Rfa3.
// We get bit-identical results for sppow=0,1,2!
      if (sppow > 0){
        LAP = LAP1 = (n*n + m*m);
        i=1;
        while(i++ < sppow) LAP *= LAP1;
        LAP = 1./LAP;
        for(i=0;i<4;i++) *(data++) = *(data1++)*LAP;
      }
      else  for(i=0;i<4;i++) *(data++) = *(data1++);
      m++;
    }
  }
  *ERR=0;
}


void fa_spectral_order(double* data,int* nmsmax,int* nsmax,
    int* nx,int* ny, Rcomplex* fftdata){
// re-order the raw data stream into complex FFT components
// Rcomplex is a struct with double i and r 
// Usually compatible with C99 "double complex"

// output is a matrix NX x NY in the 2D-FFT format of R

  int i,m,n,ind1,ind2,ind3,ind4,Mmin,Mmax,Nmin,Nmax,offset,mmax;
  double a0,a1,a2,a3;

  Mmin = floor((*nx-1)/2.);
  Nmin = floor((*ny-1)/2.);
  Mmax = floor((*nx)/2.);
  Nmax = floor((*ny)/2.);
// initialise to 0 !
  for(i=0;i<(*nx * *ny);i++) fftdata[i].r = fftdata[i].i = 0.;

  offset = 0;
  for(n=0;n<= *nsmax;n++){
    mmax = fa_mmax(n,*nsmax,*nmsmax);
    for(m=0;m<= mmax;m++) {
// m>=0,n>=0 (Mmin<0)
//      ind1 = (n-Nmin)*(*nx) + (m-Mmin);
      ind1 = n*(*nx) + m;
      fftdata[ind1].r = (a0=data[offset]) - (a3=data[offset+3]);
      fftdata[ind1].i = (a1=data[offset+1]) + (a2=data[offset+2]);

      if (m>0 && m<=Mmin) {
// m<0,n>=0
        ind2 = n*(*nx) + (*nx - m);
        fftdata[ind2].r = a0 + a3;
        fftdata[ind2].i = a1 - a2;
      }

      if(n>0 && n<=Nmin) {
// m>=0,n<0
        ind3 = (*ny-n)*(*nx) + m;
        fftdata[ind3].r = a0 + a3;
        fftdata[ind3].i = -a1 + a2;
        if (m>0 && m<=Mmin) {
// m<0,n<0
          ind4 = (*ny - n)*(*nx) + (*nx - m);
          fftdata[ind4].r = a0 - a3 ;
          fftdata[ind4].i = -a1 - a2;
        }
      }
      offset+=4;
    }
  }
}


// read a grib-0 message
// GRIB-0 as found in FA files (according to DECOGA)
// values must be allocated with length nval

void fa_grib0(unsigned char* grib,int griblen,int nval,double* values,
              double minval, double maxval,int* ERR){
  int len0,len1,len2,len3,len4,len5;
  int i,npadding,flag1,flag2;
  int nbits;
  unsigned char *sec0,*sec1,*sec4,*sec5;
  double zscale, refval;
  int bitmask, scalefactor;

  *ERR=1; // default error value
  if(griblen < 4+24+11+4) {
    Rprintf("ERROR: data stream is too short to contain regular GRIB-0 message.\n");
    *ERR=-10;
    return;
  }

//section 0
  sec0 = grib;
  len0 = 4;
  if( (sec0[0] != 'G') || (sec0[1] != 'R') || (sec0[2] != 'I') || (sec0[3] != 'B') ){
    Rprintf("ERROR section-0\n");
    *ERR=-10;
    return;
  }

//section 1
  sec1 = sec0 + len0;
  len1 = INT3(sec1);
  if( (len1 != 24) || (sec1[3] != 0) ) {
    Rprintf("ERROR section-1\n");
    *ERR=-1;
    return;
  }
  if (sec1[7] & 128) {
    Rprintf("ERROR: Section 2 not allowed!\n");
    *ERR=-2;
    return;
  }
  if (sec1[7] & 64) {
    Rprintf("ERROR: Section 3 not allowed!\n");
    *ERR=-3;
    return;
  }
#ifdef DEBUG
  Rprintf("GRIB section1:\n");
  Rprintf("[1-3]: length : %i\n[4] grib edition : %i\n",len1,(int)sec1[3]);
  for(i=1;i<20;i++) Rprintf("byte[%i]: %i\n",i+4,(int)sec1[i+3]); //same numbering as DECOGA.F
  Rprintf("LEVEL: type=%d, lev=%d\n",sec1[9],INT2((sec1+10)));
  Rprintf("DATE: %i-%i-%i %i:%i\n",sec1[12],sec1[13],sec1[14],sec1[15],sec1[16]);
#endif

//section 2 : not allowed
  len2 = 0;
//section 3 : not allowed
  len3 = 0;
//section 4 (check)
  sec4 = sec1 + len1;
  len4 = INT3(sec4);
  if (len4 < 11) {
    Rprintf("ERROR: section 4 too short.\n");
    *ERR=-4;
    return;
  }
  if(griblen != (len0 + len1 + len2 + len3 + len4 + 4)){
    Rprintf("ERROR: data stream is shorter than GRIB-0 message.\n");
    Rprintf("message length=%d, expected griblen=%d\n",
        len0 + len1 + len2 + len3 + len4 + 4,griblen);
    Rprintf("len0=%d len1=%d len2=%d len3=%d len4=%d len5=%d\n",
        len0,len1,len2,len3,len4,4);
    *ERR=-10;
    return;
  }
  grib += len4;
  
//section 5
  sec5 = sec4 + len4;
  len5 = 4;
  if( (sec5[0] != '7') || (sec5[1] != '7') || (sec5[2] != '7') || (sec5[3] != '7') ){
    Rprintf("ERROR section-5\n");
    *ERR=-5;
    return;
  }

// start reading section 4: we check a few important bytes
  npadding = sec4[3] & 15;
  flag1 = sec4[3] & 128;// ???
  flag2 = sec4[3] & 64; // complex packing...
  nbits = sec4[10]; // also provided by FA header. We don't check it explicitly.
  refval=IBMfloat(sec4+6);
  scalefactor=INT2((sec4+4)); //first bit is sign bit! Max value is 2^15.
  if (scalefactor > 32768) scalefactor = 32768 - scalefactor;
// these are never used for normal ALADIN/ARPEGE files  
// OOOPS: apparently they are (grib type=1 in stead of 2)
#ifdef DEBUG
  Rprintf("Sector 4:\n");
  Rprintf("len=%i, padding=%i\n",len4, npadding);
  Rprintf("Spher.Harm.=%i,complex packing=%i\n",flag1,flag2);
  Rprintf("scalefactor=%i\n",(int) scalefactor);
  Rprintf("refval=%lf\n",refval);
  Rprintf("nbits=%i\n",nbits);
#endif
  if ( nval * nbits != (len4-11)*8-npadding ) {
    Rprintf("Error in grib message length.\n");
    *ERR=-6;
    return;
  }
// type-1 grib files: get minval and scale from the GRIB section, but maxval not used
// type-2 grib files: minval and maxval are given, scale is calculated from both
  if (isnan(minval)) {
    minval = refval;
    maxval = NAN;
    zscale = fast_pow(2., scalefactor);
    // FIXME: if abs(scalefactor) > 63 then  "<<" can not work!
    // There are some tricks to be a bit faster, maybe, knowing
    //  that scalefactor is an integer. Found some code on internet.
//    if (scalefactor>0) zscale = ( ((long) 1)<<scalefactor);
//    else zscale = 1./( ((long) 1)<< abs(scalefactor));
  }
  else {
    bitmask = (1 << nbits) -1;
    zscale = (maxval-minval)/bitmask;
  }
  fa_grib_expand(sec4+11,nbits,nval,values,minval, maxval, zscale);
  *ERR=0;
  return;
}

// Decode a FA record
void fa_decode(unsigned char* ibuf,int* buflen,double*data,int* ndata,
               int* nsmax,int* nmsmax,int* ERR){

  int spectral,grib,nbits,sptrunc,pow; // possible FA header parameters
  double minval,maxval;
  unsigned char *obuf,*gribmessage;
  int i,j,err1,padding,nval,nval1,nval2,griblen;
  double *part1,*part2; // used for combining spectral data
  int little_endian=( *(uint16_t*)"\0\xff" > 255); // TRUE if little-endian
  unsigned char spbuf[8];
  float sp;

// FIXME: nsmax,nmsmax,are only necessary for calculating NVAL etc.
// ARPEGE files have a different formula...
  *ERR=1; // if anything goes wrong unexpectedly, ERR=1 is returned
  grib = INT8(ibuf);
  spectral = INT8(ibuf+8);
  ibuf += 16;

#ifdef DEBUG
  Rprintf("grib=%i spectral=%i\n",grib,spectral);
#endif
  nval = *ndata;

// no grib compactification -> simply copy the data (but mind the endianness and possibly single precision!)
// NOTE: it's probably better in this case to avoid the C code. Just do it in R.
  if (grib<=0) {
    if( *buflen == 8*nval + 16) {
      obuf = (unsigned char*) data;
      for(i=0; i < nval*8; i++) *(obuf++)=*(ibuf++);
      obuf = (unsigned char*) data;
      if (little_endian) byteswap(obuf,8,nval);
    }
    else if ( *buflen == 4*nval + 16) {
//      Rprintf("WARNING: assuming single precision data.\n");
      // single_precision = 1;
      // THIS IS HARDER, because we have to return double to R
      // If nval is odd, there may be some trailing zeroes.
      // Note that the values are swapped 2 by 2.
      // We have byteswap in groups of 8 bytes.
      for (i=0; i<nval; i+=2) {
        for (j=0; j<8; j++) spbuf[j] = *(ibuf++) ;
        if (little_endian) byteswap(spbuf, 8, 1); // NOT (, 4, 2)
        data[i] = (double) *((float*) spbuf);
        if (i < nval-1) data[i+1] = (double) *((float*) (spbuf+4));
      }
        
      // FIXME: in case of Single Precision: byteswap by 4 or by 8?
      // In echkevo, it is by 8!
      // byteswap(obuf, 8, nval/2);
      // So now the values may still be swapped 2 by 2
    }
    else {
      Rprintf("ERROR: FA length=%d, nval=%d\n",*buflen,nval);
      return;
    }
    *ERR=0;
    return;
  }
// grib compactification
  else {
    nbits=INT8(ibuf);
    ibuf += 8;
#ifdef DEBUG
    Rprintf("NBITS=%i\n",nbits);
#endif
    if (spectral) {
      sptrunc=INT8(ibuf);
      pow=INT8(ibuf+8);
      ibuf += 16;
      if (grib >= 2) {
        minval = DBL8(ibuf);
        maxval = DBL8(ibuf+8);
        ibuf += 16;
      }
      else {
        minval = maxval = NAN;
      }

#ifdef DEBUG
      Rprintf("SPTRUNC=%i, POW=%i\n",sptrunc,pow);
#endif
      fa_countval_spec(nmsmax,nsmax,&sptrunc,&nval1);
//      part1= (double*) malloc(nval1 * sizeof(double));
      part1= (double*) R_alloc(nval1 , sizeof(double));
// total length of the grib message (all sectors): 4 + 24 + (11+data) + 4
// add bytes to data sector to get even lengthto an even number of bytes!
      griblen = 32 + 11 + ceil(nval1*nbits/8.);
      griblen += (griblen%2);
      fa_grib0(ibuf,griblen,nval1,part1,minval,maxval,&err1);
      if(err1){
        return;
      }
//      fa_rescale(part1,nval1,minval,maxval,nbits);
// the grib message is padded to a multiple of 8 bytes in the FA message
      padding = (griblen%8) ? (8-griblen%8) : 0;
// skip to the non-compacted part of spectral data
// which is stored as big-endian doubles
      ibuf += griblen + padding;
      nval2 = nval - nval1;
      if(7*8 + griblen+padding+8*nval2 != *buflen) {
        Rprintf("ERROR: data and message lengths do not correspond.\n");
        Rprintf("7*8 + griblen+padding+8*nval2=7*8 + %i + %i + 8*%i = %i != %i.\n",
               griblen,padding,nval2,7*8+griblen+padding+nval2*8,*buflen);
//        free(part1);
        return;
      }
//      part2= (double*) malloc(nval2 * sizeof(double));
      part2= (double*) R_alloc(nval2 , sizeof(double));
      obuf = (unsigned char*) part2;
      for(i=nval2*8;i;i--) *(obuf++)=*(ibuf++);
      if(little_endian) byteswap(part2,8,nval2);
// combine both parts of spectral data
#ifdef DEBUG
      Rprintf("combining two streams for spectral data.\n");
      Rprintf("nval1=%i , nval2=%i , nval=%i \n",nval1,nval2,nval1+nval2);
#endif
      fa_spectral_combine(part1,part2,data, *nsmax, *nmsmax,sptrunc,pow,&err1);
//      free(part1);
//      free(part2);
    }
    else { // gridpoint data
      if (grib>=2) {
        minval = DBL8(ibuf);
        maxval = DBL8(ibuf+8);
        ibuf += 16;
      }
      else {
        minval = maxval = NAN;
      }
// be careful for rounding errors!
      griblen = 32 + 11 + ceil((nval*nbits)/8.);
//      griblen += (griblen%4)?(4-(griblen%4)):0;
      if(griblen%2) griblen++;
      fa_grib0(ibuf,griblen,nval,data,minval,maxval,&err1);
      if(err1){
        return;
      }
    }
  }
  *ERR=0;
  return;
}


