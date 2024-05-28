#include "rfa.h"
// parse an FA file. Return field names, byte locations etc.
// nfields is number of actual fields PLUS 7 or 8 (frame)
// the first 7 (8) entries are removed later from the field list
// foffset,flen: fields
// hoffset,hlen: holes
// We assume all vectors have been correctly allocated by the calling routine!!!

// malloc replaced by R_alloc, free no longer required
// Hopefully, that will keep the crashes away...

// We use doubles for byte addresses (easier to pass to R)
// Because integers are 32 bit, and since we multiply by 8 we got problems with huge files.

// fastfind: given a field name, find it's index, byte location & length
// We use DOUBLE for tar_offset and foffset, because 32-bit integers are too limited
void fa_fastfind_name(char **filename, double *tar_offset,char **fnm, char **fname,
                      double *foffset, int *flen, int *findex,int *err){
  FILE* fafile;
  int blocksize, nfields, maxfields,nameblock_size,nlist;
  int i,j,k,lfound;
  int64_t header[22],buff[3],next_byte,name_section_offset;
  char sbuf[17];
  int little_endian=( *(uint16_t*)"a" < 255); // TRUE if little-endian

#ifdef DEBUG
  Rprintf("Fastfind in file %s\n",*filename);
  Rprintf("tar_offset=%i, fname=%s\n",(int)*tar_offset,*fnm);
#endif

  sbuf[16]='\0';
  lfound=0;
  if(!(fafile=fopen(*filename,"r"))) {
    *err=10;
    Rprintf("Couldn't open file!\n");
    return;
  };

  if(*tar_offset) fseek(fafile,(int64_t) *tar_offset,SEEK_SET);

  k=fread(header,8,22,fafile);
  if(little_endian) byteswap(header,8,22);
  blocksize=header[0]*8;
  nfields=header[5]; // this includes the holes
  maxfields=header[12]; // number of fields per sequence 
  nameblock_size=header[19]; // how many "blocks" for info? Almost always 1

// this is looped if there are multiple name sectors
  nlist=ceil((double)(nfields)/maxfields);
#ifdef DEBUG
  if(nlist>1) Rprintf("There are %i name sections!\n",nlist);
  Rprintf("nameblock_size=%i, maxfields=%i\n",nameblock_size,maxfields);
#endif
// go to the name sector (first 7 or 8 fields will be the frame & date specifications)
  name_section_offset=*tar_offset + blocksize;
  fseek(fafile, name_section_offset, SEEK_SET);
  *findex=0;
// start looking from 7 (8)? the fields 0-6 (7) are the FA frame
// Nah, maybe you want to find the frame fields using this function.
// sometimes you may have to skip to next name sector...`
  for(i=0;i<nfields;i++){
    k=fread(sbuf,1,16,fafile);
//    if(!strcmp(sbuf,*fname)) { 
    if(strstr(sbuf,*fnm)) { 
#ifdef DEBUG
      Rprintf("Found field %s at i=%i\n",sbuf,i);
#endif
      *findex=i%maxfields;
      strcpy(*fname,sbuf);
      lfound=1;
      break;
    }
// if you have reached the end of the name section, you must skip to next part
// the address is hidden at the END of the header section
    if((i+1)%maxfields==0){
      j = (i+1)/maxfields;
#ifdef DEBUG
      Rprintf("Reached end of name section. index=%i\n", j);
#endif
      fseek(fafile, blocksize - 8*j, SEEK_SET);
      k=fread(buff,1,8,fafile);
      if(little_endian) byteswap(buff,8,1);
      name_section_offset = (buff[0]-1) * blocksize;

#ifdef DEBUG
      Rprintf("++++++ next name sector expected at %i\n",name_section_offset);
#endif
      fseek(fafile, name_section_offset, SEEK_SET);
    }
  }

  if(lfound){
#ifdef DEBUG
      Rprintf("Found field %s at index=%i\n",*fname,*findex);
#endif
    fseek(fafile,name_section_offset + nameblock_size*blocksize + *findex*16,SEEK_SET);
//    fseek(fafile, nameblock_size*blocksize,SEEK_CUR);
    k=fread(buff,2,8,fafile);
    if(little_endian) byteswap(buff,8,2);
    *flen=buff[0]*8;
    *foffset=(double) (*tar_offset+8*(buff[1]-1));
#ifdef DEBUG
      Rprintf("flen=%i, foffset=%lf\n",*flen,*foffset);
#endif
    *err=0;
  }
  else *err=1;
  fclose(fafile);
}

// Parse a complete file and return the list of fields
// FIXME: this may crash if the header data is incorrect
//        a simple error would be nicer
void fa_parse_file(char** filename,double* tar_offset,
		   int* ninfields,
                   char** fnames,double* foffset,int* flen,int*findex,
                   int* spectral,int* ngrib, int* nbits,int* sptrunc,int* sppow,
                   double* hoffset,int* hlen,int*hindex,int* lparse,int* err){

  FILE* fafile;
  int blocksize, nfields, nholes, nmeta;
  int i,j,k,ccfields,ccholes,ccholes2,ccfields2,*is_hole;
  int maxfields,nameblock_size;
  int nlist,ndata,ll;
  int64_t header[22],buff[3],next_byte,name_section_offset,address_section_offset;
  int little_endian=( *(uint16_t*)"a" < 255); // TRUE if little-endian
  
  char sbuf[17],*empty="                ";

  *err=0;
  if(!(fafile=fopen(*filename,"r"))) {
    *err=10;
    Rprintf("Couldn't open file!\n");
    return;
  };

#ifdef DEBUG
  Rprintf("Opening %s at offset %i\n",*filename,*tar_offset);
  if (little_endian) Rprintf("This is a little_endian machine.\n");
  else Rprintf("This is a big_endian machine.\n");
#endif
  if(*tar_offset) {
#ifdef DEBUG
    Rprintf("Jumping to tar_offset.\n");
#endif
    fseek(fafile, *tar_offset, SEEK_SET);
  }
  k = fread(header,8,22,fafile);
  if (little_endian) byteswap(header,8,22);
#ifdef DEBUG
  Rprintf("HEADER:\n");
  for(i=0 ; i<22 ; i++) Rprintf("header[%i]: %i\n",i+1,(int)header[i]);
#endif
  blocksize=header[0]*8;
  nholes=header[20];
  nfields=header[5]-nholes;
  maxfields=header[12]; // number of fields per sequence 
  nameblock_size=header[19]; // how many "blocks" for info? Almost always 1

#ifdef DEBUG
  Rprintf("Expecting %i fields and %i holes.\n", nfields, nholes);
  Rprintf("field names per section: %i\n", maxfields);
#endif
  if (nfields != *ninfields) {
    Rprintf("Number of expected fields does not correspond with declaration\n");
    fclose(fafile);
    return;
  }
  sbuf[16]='\0';
// # name sections:
// check file size (not OK for tar archive...)
//  maxpos = blocksize *  header[4];
//  if (fseek(fafile, 0, SEEK_END) != maxpos) {
//      Rprintf("ERROR occured. File appears to be corrupted.\n");
//      fclose(fafile);
//      return;
//    }
// this is looped if there are multiple name sectors
  ccfields = ccholes = 0;
  ccholes2 = ccfields2 = 0;
  nlist = ceil((double)(nfields+nholes)/maxfields);
#ifdef DEBUG
  if (nlist>1) Rprintf("There are %i name sections!\n",nlist);
#endif
//  is_hole=malloc(sizeof(int)*(maxfields));
  is_hole = (int*) R_alloc(maxfields, sizeof(int));
  name_section_offset=blocksize;
  for(ll=1 ; ll <= nlist;ll++){
    if (ll < nlist) ndata=maxfields;
    else ndata=(nfields+nholes)-(nlist-1)*maxfields;
#ifdef DEBUG
    Rprintf("===Starting name sector %i\n",ll);
    Rprintf("===Expecting %i fields\n",ndata);
#endif
    fseek(fafile,*tar_offset + name_section_offset,SEEK_SET);
    for (i=0 ; i<ndata ; i++){
      k = fread(sbuf, 1, 16, fafile);
//      Rprintf("article %i: %s\n",i+1,sbuf);
   
// check for holes, which have an empty name field
      if(strcmp(sbuf,empty)) { 
        if(ccfields<nfields) {
          findex[ccfields] = i; // so the first field in every name sector gets index 0
          strcpy(fnames[ccfields++], sbuf);
          is_hole[i]=0;
        }
        else {Rprintf("ERROR: nfields doesn't match\n");fclose(fafile);*err = -1;break;}
      }
      else {
        if(ccholes<nholes) {
          is_hole[i]=1;
          hindex[ccholes++] = i;
        }
        else {Rprintf("ERROR: nholes doesn't match\n");fclose(fafile);*err = -1;break;}
      } 
    }
    if(*err) break;

// read offset and length
// TODO: maybe speed up by reading/byteswapping all numbers at once
    address_section_offset=name_section_offset + nameblock_size*blocksize;
    fseek(fafile,*tar_offset + address_section_offset,SEEK_SET);
    for(i=0 ; i<ndata ; i++){
      k = fread(buff, 2, 8, fafile);
      if (little_endian) byteswap(buff,8,2);
// the file has length and start position in 8-byte words. 
// We change this to bytes and to offset (easier when using "fseek")
// also, we add the tar_offset
      if (is_hole[i]) {
        hoffset[ccholes2] = *tar_offset+8*(buff[1]-1);
        hlen[ccholes2] = buff[0]*8;
        ccholes2++;
      }
      else {
        foffset[ccfields2] = *tar_offset+8*(buff[1]-1);
        flen[ccfields2] = buff[0]*8;
        ccfields2++;
      }
    }
#ifdef DEBUG
    Rprintf("---FIELDS:\n 1. %s at %i, length=%i\n",
        fnames[0],(int) foffset[0],flen[0]);
    Rprintf("...\n %i. %s at %i, length=%i\n",
        ccfields,fnames[ccfields-1],
        (int) foffset[ccfields-1], flen[ccfields-1]);
    Rprintf("ccfields=%i, ccfields2=%i, ccholes=%i,ccholes2=%i\n",
        ccfields,ccfields2,ccholes,ccholes2);
#endif
    if (ll < nlist){
// The location of the "extension" is at the end of the first (header) sector
      fseek(fafile, *tar_offset + blocksize - 8*ll, SEEK_SET);
      k=fread(buff,1,8,fafile);
      if(little_endian) byteswap(buff,8,1);
      name_section_offset = (buff[0]-1) * blocksize;

#ifdef DEBUG
      Rprintf("++++++ next name sector expected at %i\n",name_section_offset);
#endif
    }
//
  }
//  free(is_hole);
  if(*err) {
    Rprintf("ERROR occured.\n");
    fclose(fafile);
    return;
  }
// read data headers (spectral, grib, truncation)
// the first 7 (8) fields are the "frame", so we skip those
  if (!*lparse) {
#ifdef DEBUG
    Rprintf("Not parsing data sector.\n");
#endif
    fclose(fafile);
    return;
  }
#ifdef DEBUG
    Rprintf("Parsing data sector.\n");
#endif
  if (nfields < 7) {
    Rprintf("The file contains no data fields. Parsing is pointless.\n");
    fclose(fafile);
    return;
  }
  else if (nfields < 7) {
    Rprintf("The file contains less than 7 data fields. Not enough for data frame.\n");
  }
  else {
    if (nfields == 7 || strcmp(fnames[7], "DATX-DES-DONNEES")) nmeta=7 ;
    else nmeta=8;
  }
#ifdef DEBUG
    Rprintf("Parsing data sector.\n");
#endif

  for (i=nmeta; i<nfields; i++){
#ifdef DEBUG
    Rprintf("field %i\n", i);
#endif
    if (fseek(fafile,foffset[i],SEEK_SET) != 0) {
      Rprintf("ERROR occured. File appears to be corrupted.\n");
      fclose(fafile);
      return;
    }
    if ( (k=fread(buff,8,2,fafile)) != 2) {
      Rprintf("ERROR occured. File appears to be corrupted.\n");
      fclose(fafile);
      return;
    }
    if (little_endian) byteswap(buff,8,2);
//    Rprintf("%s: grib=%li, spectral=%li\n",fnames[i],buff[0],buff[1]);
//    grib[i]=buff[0]; // GRIB or not
    ngrib[i] = buff[0];
    spectral[i] = buff[1];
    if (buff[0]>=1) { // GRIB compactification
      if ((k=fread(buff,8,1,fafile)) != 1) {
        Rprintf("ERROR occured. File appears to be corrupted.\n");
        fclose(fafile);
        return;
      }
      if(little_endian) byteswap(buff,8,1);
      nbits[i]=buff[0];
      if(spectral[i]){
        if ((k=fread(buff,8,2,fafile)) != 2) {;
          Rprintf("ERROR occured. File appears to be corrupted.\n");
          fclose(fafile);
          return;
        }
        if(little_endian) byteswap(buff,8,2); // sptrunc and lagrangian power
        sptrunc[i]=buff[0];
        sppow[i]=buff[1];
      }
      else sptrunc[i]=sppow[i]=-1;
    }
    else nbits[i]=sptrunc[i]=sppow[i]=-1;
  }

  fclose(fafile);
}

