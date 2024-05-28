## decode a single FA message
## the message is in inbuf (as raw bytestream data)
FAdec_msg <- function(inbuf, faframe, clip=TRUE, outform="G", quiet=TRUE){
## is it spectral data or gridpoint?
  ngrib <- readBin(inbuf[1:8], what="int", size=8, n=1, endian="big")
  lspec <- readBin(inbuf[9:16], what="int", size=8, n=1, endian="big")
  if (ngrib > 0) {
    nbits <- readBin(inbuf[17:24], what="int", size=8, n=1, endian="big")
  }

  if (is.null(faframe$nmsmax) | is.null(faframe$nsmax) | 
     is.null(faframe$ndgl) | is.null(faframe$ndlon) ) stop("Incorrect FA frame.") 

  if (ngrib < 100) {
    # CLASSIC FA
    if (lspec) {
      ndata <- .C("fa_countval_spec",as.integer(faframe$nmsmax),as.integer(faframe$nsmax),
                        as.integer(-1),result=integer(1),PACKAGE="Rfa")$result
    } else {
      ndata <- faframe$ndgl * faframe$ndlon
    }
    if (!quiet) {
      cat("FA message length:",length(inbuf),"bytes\n")
      cat("Allocating vectors for ndata=",ndata,"values.\n")
    }
    if (ngrib <= 0 & !lspec) {
      # non-compacted grid point data. COULD be from surfex!
      if (ndata == (length(inbuf) - 16)/8 ) {
        data <- readBin(inbuf[-(1:16)], what = "numeric", size=8, n=ndata, endian="big")
        sfx_missing <- 1.E20
      } else if (ndata == (length(inbuf) - 16)/4 ) {
        # single precision: missing data value in surfex is "rounded":
        sfx_missing <- 100000002004087734272
        data <- readBin(inbuf[-(1:16)], what = "numeric", size=4, n=ndata, endian="big")
        # And, as you may expect, the 32-bit values are SWAPPED 2 BY 2
        # like in echkevo, because of the use of 64-bit big-endian integers
        swap <- TRUE
        if (swap) {
          ttt1 <- data[c(TRUE, FALSE)]
          data[c(TRUE, FALSE)] <- data[c(FALSE, TRUE)]
          data[c(FALSE, TRUE)] <- ttt1
        }
      } else {
        stop("Unexpected data length ",length(inbuf) - 16, " for ", ndata, " data points.")
      }
      # We assume that the exact value sfx_missing will not occur in the wild...
      data[data == sfx_missing] <- NA

    } else {
      cdata <- .C("fa_decode",ibuf=as.raw(inbuf),buflen=as.integer(length(inbuf)),
                         data=numeric(ndata),ndata=as.integer(ndata),
                         nsmax=as.integer(faframe$nsmax),nmsmax=as.integer(faframe$nmsmax),
                         err=integer(1),PACKAGE="Rfa")
      if (cdata$err>0) stop("An error occured.")
      data <- cdata$data
    }
  } else {
    # NEW (eccodes) FA
    if (!require(Rgrib2)) {
      stop("This FA file uses eccodes for GRIB encoding. Requires Rgrib2 package to read.")
    }
    msg <- inbuf[-(1:24)]
    if (all(msg[5:8] == charToRaw("BIRG"))) {
      msg <- .C("fa_byteswap", data=msg, size = as.integer(8), n=as.integer(length(msg)/8))$data
    } else if (any(msg[1:4] != charToRaw("GRIB"))) {
      stop("The grib sector seems corrupted?")
    }
    # Use Rgrib2 to decode
    # NOTE: sometimes, the units in GRIB are different from the default in FA
    #       so we may need to re-scale...
    #       That's why we must call Ghandle, not just Gdec.
    gh <- Rgrib2::Ghandle(msg)
    data <- Rgrib2::Gdec(gh)
    scaling <- Rgrib2::Ginfo(gh, IntPar=c("FMULTM", "FMULTE"))
    zmult <- scaling$FMULTM * 10.^scaling$FMULTE
    data <- data / zmult
    # TODO: re-order the spectral components? Or is it OK?

  }
  if (!lspec) {
    if (is.element(outform,c("R", "S"))) warning("FA data is in grid format, can not convert.")
    FAdata <- matrix(data,ncol=faframe$ndgl)
  } else {
    if (outform=="R") return(data)  # just the raw spectral data
#-- re-arrange in a matrix of (complex) Fourier components
    FAdata <- FAraw2fft(data,faframe$nmsmax,faframe$nsmax,
                             faframe$ndlon,faframe$ndgl)
    if (outform=="S") return(FAdata)
    FAdata <- Re(fft(FAdata,inverse=TRUE))
  }
  if (clip) {
    if (!quiet) cat("Clipping data to physical domain",faframe$ndlun,":",faframe$ndlux,",",
                    faframe$ndgun,":",faframe$ndgux,"\n")
    FAdata <- FAdata[faframe$ndlun:faframe$ndlux,
                            faframe$ndgun:faframe$ndgux]
  }
#-- 5. Add meta-data as attributes? No, let FAopen do that
  return(FAdata)
}

 
### the basic function for FA decoding
FAdec <- function(fa, ...){
  UseMethod("FAdec")
}

### decode from a filename without FAopen: slightly faster for a single field
### but slower for multiple fields from the same file.
### This function exists for one context: you need to extract a single field from MANY files
###   because you can then provide faframe and avoid reading the metadata every time.
##FAdec.FAfile <- function(fa,field,clip=TRUE,outform="G",quiet=TRUE,...)

FAdec.character <- function(fa, field , clip=TRUE, outform="G", archname=NULL, tar.offset=NULL,
                  faframe=NULL, fatime=NULL, quiet=TRUE, ...){
# allow fa to be a character string with attributes (just like FAfile objects)
# if the attributes are undefined, it just gives NULL
  if (is.null(archname)) archname <- attr(fa, "tarfile")
  if (is.null(tar.offset)) tar.offset <- attr(fa, "tar.offset")

  if (is.null(archname)){
    tar.offset <- 0
    filename <- fa
    if (!file.exists(filename)) stop(paste("File",filename,"not found."))
  } else {
    filename <- archname
    if (!file.exists(filename)) stop(paste("File",filename,"not found."))
    if (is.null(tar.offset)) tar.offset <- FindInTar(archname, fa, quiet=quiet)
    if (is.na(tar.offset)) stop(sprintf("File %s not found in archive %s\n", fa, archname))
  }

## we no longer use "fa_fastfind_index". Cumbersome to debug and maintain.
  if (is.numeric(field) || length(field)>1) {
    fa1 <- FAopen(fa, archname=archname, tar.offset=tar.offset, lparse=FALSE, quiet=quiet)
    return(FAdec.FAfile(fa1, field=field, clip=clip, outform=outform, quiet=quiet, ...))
  }

## if we use strstr() to search, we don't need to format to 16 characters
  fnm <- field
  if (!quiet) cat(filename,tar.offset,fnm,"\n")
  # NOTE: for a large tar file, 32bit integers are too small to contain tar.offset!
  # so we use DOUBLE (because R doesn't support 64bit integer)
  fastfind <- .C("fa_fastfind_name",as.character(path.expand(filename)),
                  tar.offset=as.numeric(tar.offset),
                  fnm=fnm,fname="1234567890123456",foffset=numeric(1),
                  flen=integer(1),findex=integer(1),err=integer(1))
  fname <- fastfind$fname
  if (fastfind$err>0) {
    if (fastfind$err==1) stop("Field ", field, " not found.")
    else stop()
  }

  fpos <- fastfind$foffset
  flen <- fastfind$flen
  findex <- fastfind$findex
# if a frame is provided, don't read it from the file
# this is mainly useful when reading a single field from very many different files (climate data etc)
# FEATURE: if you provide a faframe, fatime is NOT read
  if (is.null(faframe)) {
    metadata <- FAread_meta(fa, archname=archname)
    faframe <- FAframe(metadata)
    fatime <- FAtime(metadata)
  }
#-- 2. read message from file
  inbuf <- FAread_msg(filename, fpos, flen)
#-- 3. send FA record to C code
  result <- FAdec_msg(inbuf, faframe, clip=clip, outform=outform, quiet=quiet)
#-- 4. Add meta-data as attributes
  if (outform=="G"){
    info <- c(FAdescribe(fname), time=list(fatime),
              origin = if (is.null(archname)) fa else sprintf("%s in %s", fa, archname))
    result <- meteogrid::as.geofield(result, domain=FAdomain(faframe), info=info)
  }
  result
}

FAdec.FAfile <- function(fa, field, clip=TRUE, outform="G", quiet=TRUE, drop_missing=TRUE, ...){
  faframe <- attr(fa,"frame")
  if (faframe$FAtype!="aladin") {
    cat("WARNING: this is an arpege file. Decoding is not supported.\nAnything may happen!\n")
  }
  #-- 1. fix the filename
  filename <- if (is.null(attr(fa, "tarfile"))) attr(fa, "filename") else attr(fa, "tarfile")

  #-- 2. fix index of the field
  fnr <- vapply(field, function(f) FAfind(fa, f)[1], 0.)

  #
  if (any(is.na(fnr))) { 
    missing <- which(is.na(fnr))
    warning("Fields ", paste(field[missing], collapse=", "), " not found.")
    if (drop_missing) fnr <- fnr[!is.na(fnr)]
    if (length(fnr)==0) stop("Nothing to decode.")
  }
  #
  myfun <- function(ff) {
    if (is.na(ff)) {
      NA 
    } else {
      fname <- fa$list$name[ff]
      fpos <- fa$list$offset[ff]
      flen <- fa$list$length[ff]
      inbuf <- FAread_msg(filename, fpos, flen)
      FAdec_msg(inbuf, faframe, clip=clip, outform=outform, quiet=quiet)
    }
  }
  #-- 4. Add meta-data as attributes
  if (!clip) dims <- c(faframe$ndlon, faframe$ndgl) 
  else dims <- c(faframe$ndlux - faframe$ndlun + 1,faframe$ndgux - faframe$ndgun + 1)

  if (outform=="G"){
    ### TODO: add unit information in info$units?
    ### could be parallellised (mc...) but probably not worth it
    if (is.null(attr(fa,"tarfile"))) origin <- attr(fa,"filename")
    else origin <- sprintf( "%s in %s",attr(fa,"filename"),attr(fa,"tarfile"))
    if (length(fnr) > 1) {
      result <- vapply(fnr, myfun, FUN.VALUE=array(1., dim=dims))
      result <- meteogrid::as.geofield(result, domain=fa, 
                          info = list("origin" = origin, "time" = attr(fa, "time")), 
			  extra_dim = list("prm"=fa$list$name[fnr]))
    } else {
      info <- c(FAdescribe(fa$list$name[fnr]), time=list(attr(fa, "time")),
                origin=origin)
      result <- meteogrid::as.geofield(myfun(fnr), domain=fa, info=info)
    }
#    dim(result) <- c("x"=NULL, "y"=NULL, "name"=fa$list$name[fnr])
#    names(dim(result)) <- c("x", "y", "name")
#    dimnames(result) <- list(NULL, NULL, fa$list$name[fnr])
  } else {
    if (length(fnr) == 1) {
      result <- myfun(fnr)
    } else if (outform == "M") {
      result <- vapply(fnr, myfun, FUN.VALUE=array(1., dim=dims))
    } else {
    # this is a #LIST#, no longer a matrix!!!
      result <- lapply(fnr, myfun)
    }
  }
  result
}

FAread_msg <- function(fa, fpos, flen) {
  # if "fa" is a connection, just jump to the location
  # if it is a FA file, check whether it is part of a tar archive
  # BUT: we assume the offset is allready correct, because that is standard in FAopenTar
  if (!inherits(fa, "connection")) {
    on.exit(try(close(fa), silent=TRUE))
    if (inherits(fa, "character")) filename <- fa
    else if (inherits(fa, "FAfile")) {
      if (!is.null(attr(fa, "tarfile"))) {
        filename <- attr(fa, "tarfile")
        if (fpos < attr(fa, "tar.offset")) stop("FAread_msg: field position smaller than tar offset.")
      } else filename <- attr(fa, "filename")
    } else stop("FAread_msg error: bad fa")

    fa <- file(filename, open="rb")
  }
  #--  jump to data location & read message
  seek(fa, fpos)
  inbuf <- readBin(fa, "raw", n=flen)
  if (length(inbuf) != flen){
    stop("Unable to read full FA message. Broken file?")
  }
  inbuf
}

FAraw2fft <- function(rawdata,nmsmax,nsmax,ndlon,ndgl){
### C routine, but R expects FFT components ordered c(0:m,(-m):(-1)), not (-m):m
  data <- .C("fa_spectral_order",data=rawdata,nmsmax=as.integer(nmsmax),nsmax=as.integer(nsmax),
                      nx=as.integer(ndlon),ny=as.integer(ndgl),
                      fftdata=complex(ndlon*ndgl,0,0))$fftdata
  matrix(data,ncol=ndgl,nrow=ndlon,byrow=F)
}


