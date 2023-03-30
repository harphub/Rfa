// FA encoding
// FAR FROM READY !!!
#include "rfa.h"
void INT8w(unsigned char* x,int64_t ll){
  int i;
  for(i = 0;i<8;i++) x[7-i] = (ll >> (i*8) )&255;
}

void DBL8w(unsigned char * x,double val){
  int64_t i;
  i = *(int64_t*)&val;
  INT8w(x,i);
}

// expected length of a complete FA data message with spectral data
void fa_message_length_spectral(int *nmsmax,int *nsmax,int *nbits,int *sptrunc,int* result){
  int header_length, data_length, griblen,nval,nval1,nval2,i;
// header: 2 to 7 doubles
  header_length = 2 + (*nbits>0)*5;
// data length
  i=-1;
  fa_countval_spec(nmsmax,nsmax,&i,&nval);
  if(nbits>0){
      fa_countval_spec(nmsmax,nsmax,sptrunc,&nval1);
// the grib message is padded to a multiple of 8 bytes
      griblen = ceil( (32 + 11 + nval1* *nbits/8)/8); 

      nval2=nval-nval1;
      data_length=griblen + nval2;
    }
    else data_length=nval;
    
  *result = 8*(header_length + data_length);
}


// add a finished FA field to the file
// fa_write_field(
// check available space
// change date
// write name & location is info sectors

// encode a set of integers to a stream 
// the integers are masked by (2^nbits - 1), to be sure we don't run into problems
// but in principle, the values should only be 0...(2^nbits-1) anyway
void fa_grib_squeeze(unsigned char*bitstream,int streamlen,int nbits,
                     double*data, int nval, double minval, double maxval)
{
  int i,bbits,length;
  int64_t bitmask, buff,tval;
  double scale,meanval;

  bitmask = (1<< nbits) -1; // the maximum integer value allowed

  scale = (maxval-minval)/bitmask;
  meanval = (minval+maxval)/2.;

  buff = bbits = 0;
  for(i=0;i < nval;i++){
// rescale by minval & maxval, to get an integer [0,pp1]
    if(data[i] < meanval) tval = round((data[i]-minval)/scale);
    else tval = bitmask - round((maxval-data[i])/scale);
    tval &= bitmask; // just for extreme safety if minval or maxval are wrong

    buff = (buff<<nbits) + tval;
    bbits += nbits;
    while (bbits >= 8) {
      bbits -= 8;
      *bitstream++ = (buff>>bbits) & 255;
    } 
  }
// after the last value, we may have a few bits left (max 7)
// we pad with zeroes
  if (bbits)  *bitstream++ = (buff << (8-bbits)) & 255;
}

void fa_spectral_split(double* data1,double* data2,double* data,
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
      for(i=0;i<4;i++) *(data2++) = *(data++);
      m++;
    }
    while(m <= mmax) {
// LAP(m,n) = (m^2 + n^2)^pow
      if(sppow>0){
        LAP = LAP1 = n*n + m*m;
// faster?: LAP = (LAP1 += 2*m+1);
        i=1;
// usually, pow=0, 1 or 2, so don't bother about a more sophisticated power function
        while(i++ < sppow) LAP *= LAP1;
        for(i=0;i<4;i++) *(data1++) = *(data++)*LAP;
      }
      else for(i=0;i<4;i++) *(data1++) = *(data++);
      m++;
    }
  }
  *ERR=0;
}


void fa_spectral_order_inv(double* data,int* nmsmax,int* nsmax,int* nx,int* ny, Rcomplex* fftdata){
// re-order FTT components from R into the raw data stream quadruplets
// Rcomplex is a struct with double i and r 
// Usually compatible with C99 "double complex"
// nx,ny are always even
// as a result, the components range from [-m/2+1,m/2]
// we may safely assume that nmsmax, nsmax are not higher than these values
  int m,n,ind1,ind2,ind3,ind4,Mmax,Nmax,offset,mmax;
  double a0,a1,a2,a3;

//  Mmin = floor((*nx-1)/2.);
//  Nmin = floor((*ny-1)/2.);
  Mmax = floor((*nx)/2.);
  Nmax = floor((*ny)/2.);
  offset = 0;
  for(n=0;n<= *nsmax;n++){
    mmax = fa_mmax(n,*nsmax,*nmsmax);
    for(m=0;m<= mmax;m++) {
      a0=a1=a2=a3=0;
// m>=0,n>=0 (Mmin<0)
      ind1 = n*(*nx) + m;
//      fftdata[ind1].r = a0 - a3;
//      fftdata[ind1].i = a1 + a2;
      if (m==0 && n==0) {
        a0=fftdata[ind1].r;
        a1=a2=a3=0;
      }
//      else if (n==0) a1=a3=0;
        
      else if (m>0 && m<Mmax) {
// m<0,n>=0
        ind2 = n*(*nx) + (*nx - m);
//        fftdata[ind2].r = a0 + a3;
//        fftdata[ind2].i = a1 - a2;
        a0 =(fftdata[ind1].r + fftdata[ind2].r)/2.;
        a3 =(-fftdata[ind1].r + fftdata[ind2].r)/2.;
        a1 =(fftdata[ind1].i + fftdata[ind2].i)/2.;
        a2 =(fftdata[ind1].i - fftdata[ind2].i)/2.;
      }
      else if(n>0 && n<Nmax) {
// m>=0,n<0
//        fftdata[ind3].r = a0 + a3;
//        fftdata[ind3].i = -a1 + a2;
        ind3 = (*ny-n)*(*nx) + m;
        a0 =(fftdata[ind1].r + fftdata[ind3].r)/2.;
        a3 =(-fftdata[ind1].r + fftdata[ind3].r)/2.;
        a1 =(fftdata[ind1].i - fftdata[ind3].i)/2.;
        a2 =(fftdata[ind1].i + fftdata[ind3].i)/2.;
        
      }
      data[offset]=a0;
      data[offset+1]=a1;
      data[offset+2]=a2;
      data[offset+3]=a3;
      offset+=4;
    }
  }
#ifdef DEBUG
  Rprintf("fa_spectral_reorder_inv wrote to length %d.\n",offset);
#endif
}


// write a grib-0 message
// GRIB-0 as found in FA files (according to DECOGA.F)

void fa_grib0_write(unsigned char* grib,int griblen,double* values,int nval,
                    int nbits, double minval, double maxval, int* ERR){
  int len0,len1,len4,len5;
  int npadding,flag1,flag2;
  int scalefactor;
  unsigned char *sec0,*sec1,*sec4,*sec5;

  *ERR=1; // default error value
#ifdef DEBUG
  Rprintf("grib0_write: griblen=%i,nval=%i,nbits=%i\n",griblen,nval,nbits);
#endif
  len4 = 11 + ceil(nval*nbits/8);
  len4 += (len4%2)?1:0;

  len0 = 4;
  len1 = 24;
  len5 = 4;
  if(griblen != len0 + len1 + len4 + len5) {
    Rprintf("ERROR: data stream is not the right length. Allocated=%i, needed=%i\n",griblen,len0 + len1 + len4 + len5);
    *ERR=-10;
    return;
  }

//section 0
  sec0 = grib;
  sec0[0]='G';
  sec0[1]='R';
  sec0[2]='I';
  sec0[3]='B';

//section 1
// Only a few entries are really important:
  sec1 = sec0 + len0;
  sec1[0]=(len1>>16)&255;//always =0
  sec1[1]=(len1>>8)&255; //=0
  sec1[2]=(len1)&255;    //=24
  sec1[3]=0; 
// some FLAGS that have to be zero:
  sec1[7]=0;

// the other parameters are mostly irrelevant for FA files
// ALADIN puts correct values for some, but e.g. gl ignores this

//origin=MF (or anything)
  sec1[4]=98;
// model and grid identification
// MF sets grid ID to 254. 255 would mean that it is defined in section 2.
  sec1[5]=1; sec1[6]=254;
//parameter -> set to "undefined"
  sec1[8]=255;

// vertical level: this is probably irrelevant in FA, because you get this from the field name.
// ALADIN writes the correct type (e.g. 109 for hybrid levels), but sets level to zero anyway...
// CHECK THAT MODEL RUNS OK IF THIS IS CHANGED.
  sec1[9]=1; // type of level (1=sfc etc)
  sec1[10]=sec1[11]=0; // level 
// date : just write a default date, it's totally irrelevant
  sec1[12]=93;sec1[13]=9;sec1[14]=2;sec1[15]=sec1[16]=0; //YY,MM,DD,HH:MM
// time unit, range & flags
  sec1[17]=sec1[18]=sec1[19]=sec1[20]=0;
// average or accumulated, number missing
// we set this to "product is valid at time T + 0
  sec1[21]=10;sec1[22]=sec1[23]=0;
// date etc -> NEVER used...
  if (grib[7] & 128) {
    Rprintf("ERROR: Section 2 not allowed!\n");
    *ERR=-2;
    return;
  }
  if (grib[7] & 64) {
    Rprintf("ERROR: Section 3 not allowed!\n");
    *ERR=-3;
    return;
  }
  grib += len1;

//section 2 : not allowed
//section 3 : not allowed
//section 4 (check)
  sec4 = sec1 + len1;
  sec4[0] = (len4>>16)&255;
  sec4[1] = (len4>>8)&255;
  sec4[2] = (len4)&255;
// number of padding bits (0-15):
  npadding=8*(len4-11)-nbits*nval;
  sec4[3] = npadding; 
//sec4[3] also contains some flags...
// 1: GP or SP  2: packing 3:float or int 4: more flags in byte14
// these must all be 0, so nothing left to do..

// sec4[4-5]: scale factor
// sec4[6-9]: refvalue in IBM-float
  sec4[4]=sec4[5]=sec4[6]=sec4[7]=sec4[8]=sec4[9]=0;
// nbits is checked by XRD code, so it must be set correctly
  sec4[10]=nbits;

//section 5
  sec5=sec4+len4;
  sec5[0]='7';
  sec5[1]='7';
  sec5[2]='7';
  sec5[3]='7';

// start writing section 4

  fa_grib_squeeze(sec4+11,len4-11,nbits,values,nval,minval,maxval);
  *ERR=0;
  return;
}

// Encode a FA record
// output is to obuf (byte stream), not to a file.
void fa_encode(unsigned char* obuf,int* buflen,double*data,int* ndata,
               int* nbits,int* sptrunc,int *spectral,int* lgrib,int*pow,
               int* nsmax,int* nmsmax,int*ndgl,int*ndlon,int* ERR){

  double minval,maxval;
  int little_endian=( *(uint16_t*)"\0\xff" > 255); // TRUE if little-endian
  unsigned char *ibuf,*gribmessage;
  int i,err1,padding,nval,nval1,nval2,griblen,grib;
  double *part1,*part2; // used for combining spectral data

  *ERR=1; // if anything goes wrong unexpectedly, ERR=1 is returned

#ifdef DEBUG
  Rprintf("nmsmax=%i,nsmax=%i,ndlon=%i,ndgl=%i,nbits=%i,sptrunc=%i,spectral=%i,lgrib=%i,pow=%i\n",
         *nmsmax,*nsmax,*ndlon,*ndgl,*nbits,*sptrunc,*spectral,*lgrib,*pow);
  Rprintf("buflen=%i, \n",*buflen);
#endif
  if(*spectral) {
    i=-1;
    fa_countval_spec(nmsmax,nsmax,&i,&nval);
  }
  else nval = *ndgl * *ndlon;

  if(nval != *ndata) {
    Rprintf("ERROR: expected nval=%d, but input data has length ndata=%d\n",nval,*ndata); 
    return;
  }
#ifdef DEBUG
  Rprintf("Starting to write output buffer.\n");
#endif
  if(*lgrib == 1) grib=2; // 2 signifies a standard ALADIN grib message
  else if (*lgrib==0) grib=0;
  else {
    Rprintf("ERROR: lgrib must be TRUE or FALSE.\n");
    return;
  }
  INT8w(obuf,(int64_t) grib );
  INT8w(obuf+8,(int64_t) *spectral);
  obuf += 16;
// no grib compactification -> simply copy the data (but mind the endianness!)
  if(grib<=0) {
#ifdef DEBUG
    Rprintf("No grib compactification.\n");
#endif
    if( *buflen != 8*nval + 16) {
      Rprintf("ERROR: FA legth=%d, nval=%d\n",*buflen,nval);
      return;
    }
    ibuf = (unsigned char*) data;
    for(i=nval*8;i;i--) *(obuf++)=*(ibuf++);
    obuf -= 8*nval;
    if (little_endian) byteswap(obuf,8,nval);
    *ERR=0;
    return;
  }
// grib compactification
  else {
#ifdef DEBUG
    Rprintf("data has grib compactification.\n");
#endif
    INT8w(obuf,(int64_t) *nbits);
    obuf += 8;
    if(*spectral){
      INT8w(obuf,(int64_t) *sptrunc);
      INT8w(obuf+8,(int64_t) *pow);
      obuf += 16;
// first split the spectral data 
      fa_countval_spec(nmsmax,nsmax,sptrunc,&nval1);
      part1= (double*) malloc(nval1 * sizeof(double));
      nval2 = nval - nval1;
      part2= (double*) malloc(nval2 * sizeof(double));
// TO DO: find the right criterion for sppow
//        it probably depends on min and max of part1
//        or on the decay of data as a function of (m^2+n^2) ...
      fa_spectral_split(part1,part2,data, *nsmax, *nmsmax,*sptrunc,*pow,&err1);

// calculate minimum and maximum of part1
      minval=maxval=part1[0];
      for(i=1;i<nval1;i++){
        if(part1[i]<minval) minval=part1[i];
        if(part1[i]>maxval) maxval=part1[i];
      }
      DBL8w(obuf, minval);
      DBL8w(obuf+8, maxval);
      obuf += 16;
// total length of the grib message (all sectors): 4 + 24 + (11+data) + 4
// and padded to an even number of bytes
      griblen = 32 + 11 + ceil(nval1* *nbits/8.);
      griblen += (griblen%2)?1:0;
// the grib message is padded to a multiple of 8 bytes in the FA message
      if (griblen % 8) {
#ifdef DEBUG
        Rprintf("Padding grib message with %i zeroes.\n",8-griblen%8);
#endif
        padding=8-griblen%8;
      } else padding=0;
      if(griblen + padding + 8*7 + 8*nval2 != *buflen){
        Rprintf("ERROR: the allocated message length %i is not equal to the expected length %i.\n",
               *buflen,griblen+padding+8*7 + 8*nval2);
        return;
      }

      fa_grib0_write(obuf,griblen,part1,nval1,*nbits,minval,maxval,&err1);
      if(err1){
        return;
      }
      obuf += griblen;
// the grib message is padded to a multiple of 8 bytes in the FA message
      if (griblen % 8) for(i=griblen%8;i;i--) *(obuf++)=0;

// now the non-compacted part of spectral data
// which is stored as big-endian doubles
      if(little_endian) byteswap(part2,8,nval2);
      ibuf = (unsigned char*) part2;
      for(i=nval2*8;i;i--) *(obuf++)=*(ibuf++);
      free(part1);
      free(part2);
    }
    else { // gridpoint data
//      fa_rescale_inv(data,nval,minval,maxval,nbits);
      griblen = 32 + 11 + ceil(nval* *nbits/8.);
      griblen += (griblen%2)?1:0;
#ifdef DEBUG
      Rprintf("gridpoint data in GRIB format.\nexpected griblen=%i\n",griblen);
#endif
// the grib message is padded to a multiple of 8 bytes in the FA message
      if (griblen % 8) {
        padding=8-griblen%8;
      } else padding=0;
      if(griblen + padding + 8*5 != *buflen){
        Rprintf("ERROR: the allocated message length %i is not equal to the expected length %i.\n",
               *buflen,griblen+padding+5*8);
        return;
      }
      minval=maxval=data[0];
      for(i=1;i<nval;i++){
        if(data[i]<minval) minval=data[i];
        if(data[i]>maxval) maxval=data[i];
      }
#ifdef DEBUG
      Rprintf("Writing minval=%lf, maxval=%lf\n",minval,maxval);
#endif
      DBL8w(obuf,minval);
      DBL8w(obuf+8,maxval);
      obuf += 16;
#ifdef DEBUG
      Rprintf("Calling fa_grib0_write.\n");
#endif
      fa_grib0_write(obuf,griblen,data,nval,*nbits,minval,maxval,&err1);
      if(err1){
        return;
      }
      obuf += griblen;
#ifdef DEBUG
      Rprintf("Padding grib message with %i zeroes.\n",8-griblen%8);
#endif
      if(padding) for(i=padding;i;i--) *(obuf++)=0;
    }
  }
  *ERR=0;
  return;
}


