######################################
### Rfa: FA file format API for R  ###
### v4.0: Alex Deckmyn, 2015       ###
######################################
### Rfa without the XRD and FA libraries !

## 1. Some methods (cfr. meteogrid package)

print.FAfile <- function(x, ...){
    cat(paste("File", attr(x,"filename"), ":\n"))
    if(!is.null(attr(x,"tarfile"))) cat(paste("part of tar archive",attr(x,"tarfile")),"\n")
    cat(paste("Valid at",
            as.character(attr(x,"validdate"),format="%Y/%m/%d %H:%M:%S")),"\n")
    cat(paste("Contains",attr(x,"nfields"), "fields.\n"))
}

print.FAframe <- function(x, ...){
    cat("FA frame",x$name,"\n")
    cat("physical grid dimension:",x$ndlux - x$ndlun +1,"x",
                                   x$ndgux - x$ndgun + 1,
        "nlev=",x$nlev,"\n")
    cat("nsmax=",x$nsmax,"nmsmax=",x$nmsmax,"\n")
    cat("ndlon=",x$ndlon,"ndgl=",x$ndgl,"\n")
}

################

## 2. The main function for examining a FA file

FAopen <- function(filename, archname=NULL, tar.offset=NULL, lparse=TRUE, quiet=TRUE){
  if (inherits(filename, "FAfile")) return(filename)
  if (!is.character(filename)) stop("Not a regular filename")
### this routine understands the "large-scale" structure of the FA file.
### It is an LFI file with some special records.
### The first 3 parts of the FA or LFI file have dimension blocksize
### tar.offset is used if the file is part of a tar archive
### lparse=TRUE by default. It is sometimes useful to have an overview of
### spectral, truncation, nbits in grib compactification etc.
  if (is.null(archname)) archname <- attr(filename, "tarfile")
  if (is.null(tar.offset)) tar.offset <- attr(filename, "tar.offset")

  if (is.null(archname)) {
    tar.offset <- 0
    fnam <- path.expand(filename)
  } else {
    if (is.null(tar.offset)) tar.offset <- FindInTar(archname,filename,quiet=quiet)
    if (!quiet) cat("Tarfile offset:",tar.offset,"\n")
    fnam <- archname
  }
  if (!file.exists(fnam)) stop(sprintf("File %s does not exist.",fnam))
  ff <- file(fnam,"rb")
  seek(ff, tar.offset)
  FAheader <- readBin(ff, what="integer", n=22, size=8, endian="big")
  close(ff)
  if (!quiet) print(FAheader)

  #-- some first-order checking:
  if (FAheader[4] != 22 || FAheader[2] != 16) stop("This probably isn't a regular FA file.")

  #-- there may be "holes" in the file
  #-- these are just the same as a logical article, but have a blanc name
  nholes <- FAheader[21]
  nfields <- FAheader[6] - nholes ### the first 7(8) fields are the "frame"!
  if (!quiet) cat("nfields=", nfields, "nholes=", nholes, "\n")
  #-- if the file path contains a ~, the C code will crash dramatically!
  #-- so we use path.expand for safety
  #-- offset values are passed as double, because we need 64bits, especially in tar archives,
  #-- to overcome 4GB file size limit.
  # FIXME: for safety, we *must* also pass nfields so it can be checked!
  faparse <- .C("fa_parse_file",filename=path.expand(fnam),tar_offset=as.numeric(tar.offset),
                               nfields=as.integer(nfields),
                               fnames=rep("                ",nfields),foffset=numeric(nfields),
                               flen=integer(nfields),findex=integer(nfields),
                               spectral=integer(nfields),nbits=integer(nfields),
                               sptrunc=integer(nfields),sppow=integer(nfields),
                               hoffset=numeric(nholes+1),hlen=integer(nholes+1),
                               hindex=integer(nholes+1),lparse=as.integer(lparse),err=integer(1))
  metadata <- FAread_meta(filename, archname=archname)
#-- recent addition: frame may have 8th field DATX-DES-DONNEES
#-- so we can no longer hard-code the number "7"!
  nmeta <- length(metadata)
  falist <- data.frame(name=faparse$fnames,offset=faparse$foffset,length=faparse$flen,
                 spectral=faparse$spectral,nbits=faparse$nbits,
                 sptrunc=faparse$sptrunc,sppow=faparse$sppow,index=faparse$findex,
                 stringsAsFactors=FALSE)[-(1:nmeta),]
  metalist <- data.frame(name=faparse$fnames, offset=faparse$foffset, length=faparse$flen,
                 index=faparse$findex, stringsAsFactors=FALSE)[(1:nmeta),]

  rownames(falist) <- NULL; # remove the nmeta frame sectors from the list
  if (lparse) {
    falist$spectral <- (falist$spectral==1)
    falist$nbits[which(falist$nbits== -1)] <- NA
    falist$sptrunc[!falist$spectral | is.na(falist$nbits)] <- NA
    falist$sppow[!falist$spectral | is.na(falist$nbits)] <- NA
  } else {
    falist$spectral <- NULL
    falist$nbits <- NULL
    falist$sptrunc <- NULL
    falist$sppow <- NULL
  }
  if (nholes > 0) {
    faholes <- data.frame(offset=faparse$hoffset,length=faparse$hlen,index=faparse$hindex)[1:nholes,]
  } else {
    faholes <- data.frame(offset=numeric(0),length=numeric(0),index=numeric(0))
  }
  result <- list(list=falist, holes=faholes, meta=metalist)
  class(result) <- c("FAfile",class(result))
  attr(result, "header") <- FAheader
  attr(result, "filename") <- filename
  attr(result, "nmeta") <- nmeta
  attr(result, "nfields") <- nfields-nmeta
  attr(result, "nholes") <- nholes
  attr(result, "frame") <- FAframe(metadata)
  attr(result, "domain") <- FAdomain(attr(result,"frame"))
  attr(result, "time") <- FAtime(metadata)
  attr(result, "tarfile") <- archname
  attr(result, "tar.offset") <- tar.offset
  result
}

FAfind <- function(fa, field, as_name=FALSE){
  # numeric fields: just check for maximum
  if (is.numeric(field)) {
    fnr <- if (field <= attr(fa, "nfields")) field else  NA
  } else {
    # now try a perfect match
    fnr <- match(sprintf("%-16.16s", field), fa$list$name)
    if (is.na(fnr)) {
      # maybe the name was incomplete ("CLSTEMP")? 
      # not recommended, but I do it myself...
      fnr <- grep(field, fa$list$name, fixed=TRUE)
#      if (length(fnr)==0)
    }
  }
  if (as_name) fa$list$name[fnr] else fnr
}

FAgrep <- function(fa, field, as_name=TRUE) {
  if (!inherits(fa, "FAfield")) fa <- FAopen(fa)
  grep(field, fa$list$name, value=as_name)
}  


"FAdescribe" <- function(fname) {
  if (substr(fname,1,5)=="TROPO") {
    levelinfo <- "Tropopause"
    varname <- substr(fname,6,16)
  } else if (substr(fname, 1, 4)=="SURF") {
    levelinfo <- "Surface"
    varname <- substr(fname, 5, 16)
  } else if (substr(fname,1,4)=="PROF") {
    levelinfo <- "Prof"
    varname <- substr(fname,5,16)
  } else if (substr(fname,1,3)=="JET") {
    levelinfo <- "Max wind speed level"
    varname <- substr(fname,4,16)
  }  else if (substr(fname,1,3)=="MER" |substr(fname,1,3)=="MSL" ) {
    levelinfo <- "MSL"
    varname <- substr(fname,4,16) 
  } else if (substr(fname,1,1)=="S" & regexpr("[0-9]{3}",substr(fname,2,4))==1) {
    levelinfo <- paste("Hybrid Level",as.numeric(substr(fname,2,4)))
    varname <- substr(fname,5,16)
  } else if (substr(fname,1,1)=="P" & regexpr("[0-9]{5}",substr(fname,2,6))==1 ) {
    levelinfo <- paste(as.numeric(substr(fname,2,6))/100,"hPa")
    varname <- substr(fname,7,16)
  } else if (substr(fname,1,1)=="H" &  regexpr("[0-9]{5}",substr(fname,2,6))==1) {
    levelinfo <- paste(as.numeric(substr(fname,2,6)),"m")
    varname <- substr(fname,7,16)
  } else if (substr(fname,1,1)=="V" &  regexpr("[0-9]{3}",substr(fname,2,4))==1) {
    levelinfo <- paste("Iso Pot. Vort",as.numeric(substr(fname,2,4)))
    varname <- substr(fname,5,16)
  } else if (substr(fname,1,1)=="T" &  regexpr("[0-9]{3}",substr(fname,2,4))==1) {
    levelinfo <- paste(as.numeric(substr(fname,2,4)),"K")
    varname <- substr(fname,5,16)
  } else {
    levelinfo <- ""
    varname <- fname
  }

  list(name=fname,variable=varname,level=levelinfo)
}

