## R has seperate file pointers for read and write for "r+b" files
## so when calling "seek", make sure you move the right one !!!

#############################################
###    FAenc_msg                          ###
### FA message encoding, no file writing. ###
#############################################
"FAenc_msg" <- function(data,faframe,lspec=NULL,lgrib=NULL,nbits=NULL,sptrunc=NULL,sppow=NULL,quiet=TRUE){
### this only does the actual encoding, not writing to file etc.
### you need to provide a frame if you want to write spectral data etc.
### but it could be just a list with the necessary values...
### all encoding parameters are NULL: default is to use those from the original
### if the field is new, some defaults will be taken
  warning("FAenc is still in an experimental phase. Be careful!",call. =FALSE,immediate. =TRUE)

  if (inherits(faframe,"FAfile")) faframe <- attr(faframe,"frame")

### default values
### TO DO: base FA routines have some criterion to vary sppow...
  if (is.null(lspec)) lspec <- FALSE
  if (is.null(lgrib)) lgrib <- TRUE
  if (is.null(nbits)) nbits <- 16
  if (is.null(sptrunc)) sptrunc <- 10
  if (is.null(sppow)) sppow <- 1

### the full domain C+I+E
  ndlon <- faframe$ndlon
  ndgl  <- faframe$ndgl
### the actual domain C+I
  ndlux <- faframe$ndlux
  ndgux <- faframe$ndgux
### the spectral truncation
  nsmax <- faframe$nsmax
  nmsmax<- faframe$nmsmax

  if (dim(data)[1]!=ndlon | dim(data)[2]!=ndgl) {
    if (dim(data)[1] != ndlux | dim(data)[2]!=ndgux){
      stop("data has wrong dimension.")
    }
    if (!quiet) cat("Getting the data into the right size with biper!\n")
    data <- biper(data, ext=c(ndlon-ndlux, ndgl-ndgux))
  }

## 1. find message length so you can allocate before calling C

  if (!lspec) {
    nval <- ndlon * ndgl
    data <- as.vector(data)
    if (!lgrib) {
      mlen <- 16 + 8*nval
    } else {
      griblen <- 32 + 11 + ceiling(nval*nbits/8) # plus 1 byte if griblen is odd, but we pad to 8n anyway
## padding to multiple of 8 bytes:
      if(griblen%%8) griblen <- griblen + (8-griblen%%8)
      mlen <- 40 + griblen
    }
  } else {
    data.fft <- fft(data)
    data <- FAfft2raw(data.fft, nmsmax, nsmax, ndlon, ndgl, quiet=quiet)/(ndlon*ndgl)
    nval <- length(data)
    if (!lgrib) {
      mlen <- 16 + 8*nval
    } else {
      nval1 <- .C("fa_countval_spec",as.integer(nmsmax),as.integer(nsmax),
                   as.integer(sptrunc),result=integer(1),PACKAGE="Rfa")$result
      nval2 <- nval - nval1
      griblen <- 32 + 11 + ceiling(nval1*nbits/8) 
      if (griblen%%8) griblen <- griblen + (8-griblen%%8)
      mlen <- 56 + griblen + 8*nval2
    }
  }

## 2. call C routine for actual encoding
  if (!quiet) cat("mlen=",mlen,"nval=",nval,"nbits=",nbits,"sptrunc=",sptrunc,
                  "\nlspec=",lspec,"lgrib=",lgrib,"sppow=",sppow,"nsmax=",nsmax,"nmsmax=",nmsmax,"\n")
## avoid NA values in the call to C: you can allow them, but then you can't see errors as well
  if (is.na(sptrunc)) {
    sppow <- 0
    sptrunc <- 0
  }
  encoded <- .C("fa_encode",obuf=raw(mlen),buflen=as.integer(mlen),data=as.vector(data),ndata=as.integer(nval),
                   as.integer(nbits),as.integer(sptrunc),spectral=as.integer(lspec),lgrib=as.integer(lgrib),
                   pow=as.integer(sppow),nsmax=as.integer(nsmax),nmsmax=as.integer(nmsmax),
                   ndgl=as.integer(ndgl),ndlon=as.integer(ndlon),
                   ERR=integer(1))
  if (encoded$ERR !=0) stop("an error seems to have occured")
### now we still have to write the message to the file
# for testing, we stop here
  return(encoded$obuf)
}
#####################################
FArename <- function(fa, field, newname, quiet=TRUE) {
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa,quiet=quiet)
  if (!is.null(attr(fa,"tarfile"))) stop("Manipulation of files in an archive is forbidden. Read only!")
  if (nchar(newname)<5) stop("Field name must have at least 5 characters.")
  newname <- format(newname, width=16)
  if (!is.na(FAfind(fa, newname))) {
    cat("WARNING: field",newname,"already exists. Remove that first.\n")
    return(fa)
  }
  fnr <- FAfind(fa, field)
  if (is.na(fnr)) {
    cat("WARNING: field",field,"doesn't exist.\n")
    return(fa)
  }
  if (!quiet) cat("Changing field",fa$list$name[fnr],"to",newname,".")
  findex <- fa$list$index[fnr]
  header <- attr(fa, "header")
  ff <- file(attr(fa, "filename"), open="r+b")
  on.exit(try(close(ff), silent=TRUE))
## set "open" byte
#  seek(ff,16,rw="write")
#  writeBin(as.integer(0),ff,size=8,endian="big")

  blocksize <- header[1]*8
  if (!quiet) cat("Setting namefield at location",blocksize+16*findex,"\n")
  seek(ff, blocksize+16*findex, rw="write")
  writeChar(newname, ff, eos=NULL)
  close(ff)
  return(FAopen(attr(fa,"filename")))
}

#######################################
###            FAremove             ###
### remove a field from an FA file  ###
#######################################
FAremove <- function(fa,field,quiet=TRUE){
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa,quiet=quiet)
  if (!is.null(attr(fa,"tarfile"))) stop("Manipulation of files in an archive is forbidden. Read only!")
  fnr <- FAfind(fa,field)
  if (is.na(fnr)) {
    cat("WARNING: field",field,"doesn't exist.\n")
    return(fa)
  }
### the "real" field index takes the 7 (8) frame fields and holes into account:
  findex <- fa$list$index[fnr]
  pos <- fa$list$offset[fnr]
  len <- fa$list$length[fnr]
  if (!quiet) cat("fnr=",fnr,"findex=",findex,"pos=",pos,"len=",len,"\n")
  header <- attr(fa,"header")
  ff <- file(attr(fa,"filename"),open="r+b")
  on.exit(try(close(ff),silent=TRUE))
## set "open" byte
  seek(ff,16,rw="write")
  writeBin(as.integer(0),ff,size=8,endian="big")

### 1. set the name to "empty"
  blocksize <- header[1]*8
  if (!quiet) cat("Setting namefield to empty at location",blocksize+16*findex,"\n")
  seek(ff,blocksize+16*findex,rw="write")
  empty <- "                "
  writeChar(empty,ff,eos=NULL)

### 2. go to address section
# in fact we don't need to change anything
# except if there already was a "small" hole...
# OR if we are erasing the last field: then we remove everything
# BUG: if findex > header[13]: continued after data section
# BUG 2: if we remove the lst field and want to "erase completely", it should not be empty but "**FIN D'INDEX**"
  nnamesect <- header[7] # #of name sectors, usually 1
  seek(ff,(1+nnamesect)*blocksize+16*findex,rw="write")
  if (findex == header[6]){  # it's the last entry in the data sector, so we can remove all signs of it.
    writeBin(as.integer(c(0,0)),ff,size=8,endian="big")
  } else {
## we recalculate len because it may be shortened by previous re-writes (very unlikely...)
## the next sector may be a hole, so you must somehow use findex
    PPP <- rbind(fa$list[,c("offset","index")],fa$holes[,c("offset","index")])
    len <- PPP$offset[which(PPP$index==findex+1)] - fa$list$offset[fnr]
    if(!quiet) cat("Given length:",fa$list$length[fnr],"calculated from offsets:",len,"\n")
    writeBin(as.integer(len/8),ff,size=8,endian="big")
  }

### 3. set data section to zero
  if (!quiet) cat("Setting data sector to zero.\n")
  seek(ff,pos,rw="write")
  writeBin(as.integer(rep(0,len/8)),ff,size=8,endian="big")
### 4. adapt the header
# 1 more hole, except if it was the last field
# in that case we just remove all traces of it
  if (findex < header[6]){
    header[21] <- header[21]+1
  } else {
    header[6] <- header[6]-1
  }
# data length is decreased
  header[9] <- header[9] - fa$list$length[fnr]/8

# if we removed the longest field, we must change maxlength in the header
# BUGLET: this should also look at frame sectors, even if they are probably shorter
  if (!quiet) cat("max length:",header[8]*8,"this field:",len,"other fields:",
     max(fa$list$length[-fnr]),"\n")
  if (len==header[8]*8 & max(fa$list$length[-fnr])<header[8]*8) {
    if(!quiet) cat("We removed the largest field. Resetting maxsize in header.\n")
    header[8] <- max(fa$list$length[-fnr])/8
  }
# set modification date (not important at all...)
  header <- FAheader.date(header)
  attr(fa,"header") <- header
  FAwrite_header(ff, header)
  close(ff)

## finish by re-parsing the file to get field list etc. up to date
  return(FAopen(attr(fa,"filename")))
}
##################################################################
###    FAenc                                                   ###
### the basic interface for writing a field to file            ###
### if data is already encoded, it is a "raw" byte vector      ###
##################################################################
FAenc <- function(fa,fieldname,data,lspec=NULL,lgrib=NULL,nbits=NULL,
                  sptrunc=NULL,sppow=NULL,overwrite=TRUE,quiet=TRUE){
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa,quiet=quiet)
  if (!is.null(attr(fa,"tarfile"))) 
    stop("Manipulation of files in an archive is forbidden. Read only!")
  fnr <- FAfind(fa,fieldname)

  if (length(fnr)==0 || is.na(fnr)) {
    if (!quiet) cat("This is a new data field, not a replacement.\n")
    if (!is.raw(data)) data <- FAenc_msg(data,attr(fa,"frame"),lspec=lspec,lgrib=lgrib,nbits=nbits,
                       sptrunc=sptrunc,sppow=sppow,quiet=quiet)
    fa <- FAadd(fa, fieldname, data, quiet=quiet)
    return(invisible(fa))
  }

  if (!quiet) cat("This is an existing data field (fnr=",fnr,"), data will be replaced.\n")

# in this case, by default you use exactly the same encoding as the original 
  if (!is.raw(data)){
    if (is.null(lspec)) lspec <- fa$list$spectral[fnr]
    if (is.null(lgrib)) lgrib <- !is.na(fa$list$nbits[fnr])
    if (is.null(nbits)) nbits <- fa$list$nbits[fnr]
    if (is.null(sptrunc)) sptrunc <- fa$list$sptrunc[fnr]
    if (is.null(sppow)) sppow <- fa$list$sppow[fnr]
    if (!quiet) cat("Encoding options: lspec=",lspec,"lgrib=",lgrib,"\n",
                   "nbits=",nbits,"sptrunc=",sptrunc,"sppow=",sppow,"\n")
    data <- FAenc_msg(data,attr(fa,"frame"),lspec=lspec,lgrib=lgrib,
                  nbits=nbits,sptrunc=sptrunc,sppow=sppow,quiet=quiet)
  }

  mlen <- length(data)
  fieldname <- fa$list$name[fnr]
  findex <- fa$list$index[fnr]
  fpos <- fa$list$offset[fnr]
  flen <- fa$list$length[fnr]
  header <- attr(fa,"header")

  if (!quiet) cat("mlen=",mlen,"flen=",flen,"fieldname=",fieldname,"findex=",findex,"\n")

  if (flen < mlen) {
    if (!quiet) cat("The replacement data message is longer than the original.\n")
###special cases:
### - if its the last field, you don't have to leave a hole (but that works OK with FAremove)
### - if the current flen is already leaving  a small hole
### But in both cases, you can just remove the field and add the replacement.
### It may end up at a different place, but who cares?
    fa <- FAremove(fa,fieldname,quiet=quiet)
    fa <- FAadd(fa,fieldname,data,quiet=quiet)
# set header correctly ([12] +1)
    header <- FAheader.date(attr(fa,"header"))
    header[12] <- header[12]+1
    attr(fa,"header") <- header
    FAwrite_header(fa)
    return(invisible(fa))
  }

#-- write to file
  ff <- file(attr(fa, "filename"), open="r+b")
  blocksize <- header[1] * 8

  if (!quiet) cat("Writing data mlen=", mlen, "at position", fpos, "\n")
  seek(ff, fpos, rw="write")
  writeBin(as.raw(data),ff)
  if (flen > mlen) {
    if (!quiet) cat("New data is shorter. This will leave a small hole.\n")
## add zeroes to fill up the hole
    writeBin(as.raw(rep(0,flen-mlen)),ff)
## change the data length in the length sector and in the header
    seek(ff,2*blocksize+16*findex,rw="write")
    writeBin(as.integer(length(data)),ff,size=8,endian="big")
### adapt the header
    header[11] <- header[11]+1   # "shorter rewrite" counter
    header[9] <- header[9] - (flen-mlen)/8 # total data length
#-- if this was the longest data set, check for new maximum
    if (flen/8 == header[8]) {
#-- BUGLET: ?do holes count also? ?we should also consider the 7 (8) metafields?
      header[8] <- max(c(fa$list$length[-fnr],fa$holes$length))/8
    } 
  }

  header <- FAheader.date(header) 
  attr(fa,"header") <- header
  FAwrite_header(ff, header)
  close(ff)

  return(invisible(fa))
}

######################################
###    FAadd                       ###
### add a new field to an FA file  ###
######################################
FAadd <- function(fa, fieldname, data, quiet=TRUE, ...){
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa,quiet=quiet)
  if (!is.null(attr(fa,"tarfile"))) stop("Manipulation of files in an archive is forbidden. Read only!")
  filename <- attr(fa,"filename")

  if (nchar(fieldname) < 5) stop("Field name must have at least 5 characters.")
  fieldname <- format(fieldname, width=16)

  if (!is.raw(data)) data <- FAenc_msg(data, fa, quiet=quiet, ...)
  mlen <- length(data)

  ff <- file(filename, open="r+b")
  on.exit(try(close(ff), silent=TRUE))

  header <- attr(fa, "header")
  blocksize <- header[1]*8

#-- is there a hole big enough?
  hh <- fa$holes[which(fa$holes$length >= mlen),]
  if (dim(hh)[1] > 0) {
    hh <- hh[order(hh$length),]
    ## if hole length larger than mlen -> "shorter rewrite" ???
    if (!quiet) cat("filling hole!\n")
    header[21] <- header[21] - 1 # remove 1 hole from counter
    pp <- hh$offset[1]
    ll <- hh$length[1] # we overwrite this with mlen anyway
    findex <- hh$index[1]
  } else {
## TODO: this assumes exactly 3 meta-sectors
## multiple name sectors etc not yet OK
    if (!quiet) cat("Appending at the end of file.\n")
    filelen <- header[5] * blocksize
    header[6] <- header[6]+1 # 1 field more
    lastf <- dim(fa$list)[1]
    if (lastf == 0) { #no fields yet in file, only frame -> offset ???
      nmeta <- attr(fa, "nmeta")
      findex <- nmeta # "index" starts from 0, so nmeta means the first field after frame
      pp <- fa$meta$offset[nmeta] + fa$meta$length[nmeta]
    } else {
      findex <- fa$list$index[lastf]+1 # this may be more than lastf (frame & holes)
      pp <- fa$list$offset[lastf] + fa$list$length[lastf]
#      cat("lastf", lastf, "pp", pp, "mlen", mlen,"filelen", filelen,"\n",sep=" ")
    }
    if (!quiet) cat("filelen=", filelen, "pp=", pp, "mlen=", mlen, "\n");
    if (pp + mlen > filelen){
      addsector <- ceiling((pp+mlen-filelen)/blocksize)
      if (!quiet) {
        cat("message does not fit into file length",filelen,".\n")
        cat("adding", addsector, "sectors to the file.\n")
      }
      seek(ff, filelen, rw="write")
      writeBin(as.integer(rep(0, addsector*blocksize/8)), ff, size=8, endian="big")
      header[5]  <- header[5]  + addsector
      header[22] <- header[22] + addsector
    }
  }
## write field name
  if (!quiet) cat("Writing fieldname", fieldname, "at position", blocksize+16*findex, "\n")
  seek(ff, blocksize+16*findex, rw="write")
  writeChar(fieldname, ff, nchars=16, eos=NULL)
  seek(ff, blocksize+16*findex, rw="read")
## write field length and location
  if (!quiet) cat("Writing len", mlen/8, "and offset", pp/8, "at position", 2*blocksize+16*findex, "\n")
  seek(ff, 2*blocksize+16*findex, rw="write")
## add 1 to the offset: FA files internally start from 1, not 0
  writeBin(as.integer(c(mlen/8, pp/8+1)), ff, size=8, endian="big")
## write data itself
  seek(ff, pp, rw="write")
  writeBin(data,ff)
## update header ([max|tot] data length, edit time)
  if (mlen/8 < header[7]) header[7] <- mlen/8 # new "longest field" value?
  if (mlen/8 > header[8]) header[8] <- mlen/8 # new "longest field" value?
  header[9] <- header[9] + mlen/8
  header <- FAheader.date(header)
  attr(fa, "header") <- header
  FAwrite_header(ff, header)
  close(ff)
  invisible(FAopen(filename))
}


FAfft2raw <- function(fftdata,nmsmax,nsmax,ndlon,ndgl,quiet=TRUE){
### C routine, but R expects FFT components ordered c(0:m,(-m):(-1)), not (-m):m
  if (!quiet) cat("nmsmax=",nmsmax,"nsmax=",nsmax,"ndlon=",ndlon,"ndgl=",ndgl,"dim=",dim(data),"\n")
  nval <- .C("fa_countval_spec",as.integer(nmsmax),as.integer(nsmax),
                                as.integer(-1),nval=integer(1))$nval
  if (!quiet) cat("Expecting",nval,"values (",nval/4,"quadruplets).\n")
  .C("fa_spectral_order_inv",rawdata=double(nval),
                      nmsmax=as.integer(nmsmax),nsmax=as.integer(nsmax),
                      nx=as.integer(ndlon),ny=as.integer(ndgl),
                      fftdata=as.vector(fftdata))$rawdata
}


