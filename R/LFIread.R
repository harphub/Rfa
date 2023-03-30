### Rlfi : a set of routines for decoding LFI files from SURFEX
### Uses only R routines for reading, but requires meteogrid for projections
### Alex Deckmyn, KMI
### Original version: 2011-08-10

### Only validated for SURFEX!!!
### This has NEVER been validated for other data (e.g. Arome) 

### to change a series of integers into ASCII characters
### useful for reading some character-valued records

int2ascii <- function(x) rawToChar(as.raw(x))

print.LFIfile <- function(x,...){
    print.noquote(paste("File", attr(x,"filename"), ":"))
    print.noquote(paste("Valid at",
            as.character(attr(x,"validdate"),format="%Y/%m/%d %H:%M:%S")))
    print.noquote(paste("Contains",attr(x,"nfields"), "fields."))
}

LFIopen <- function(filename, quiet=TRUE){
  on.exit(try(close(lfi), silent=TRUE))
  lfi <- file(filename, open="rb")
  LFIheader <- readBin(lfi, what="integer", n=22, size=8, endian="big")
  if (!quiet) message("LFI header:", LFIheader, "\n")

  nfields <- LFIheader[6]
  blocksize <- LFIheader[1]*8
### skip to the second "record"
  if (nfields==0) {
    warning("This appears to be an empty LFI file!")
    result <- list()
  } else {
    if (!quiet) message("Scanning data names (sector 2)\n")
    seek(lfi, blocksize)
#  allnames <- readChar(lfi,16*nfields)  ### should be OK if only ASCII code
### ATTENTION: Name length is now only 12 characters
### but still "16" in the header.
### Also, the stored length is 16 characters, so the last 4 can be random values
### so set the last for characters to "space"
### they will probably fix this in the future...
    rawnames <- readBin(lfi, "raw", n=16*nfields)
    rawnames[c(rep(FALSE,12), rep(TRUE, 4))] <- as.raw(32)
    allnames <- rawToChar(rawnames)
    LFInames <- substring(allnames, seq(1, 16*nfields, by=16), seq(16, 16*nfields, by=16))
### skip to list of offsets and lengths starting at 12000h
    if(!quiet) cat("Scanning offset and data length (sector 3)\n")
    seek(lfi,2*blocksize)
    sizepos <- readBin(lfi, what="integer", size=8, n=2*nfields, endian="big")
    LFIpos <- sizepos[c(FALSE, TRUE)]
    LFIlen <- sizepos[c(TRUE, FALSE)]
    result <- list(list=data.frame(name=LFInames,
                                   offset=8*(LFIpos-1),
                                   length=8*LFIlen, stringsAsFactors=FALSE))
  }
  attr(result, "header") <- LFIheader
  attr(result, "nfields") <- nfields
  class(result) <- c("LFIfile", class(result))
  attr(result, "filename") <- filename
### by encapsulating with try(), this will not crash on e.g. non-surfex LFI files
  try(attr(result, "domain") <- LFIdomain(result,quiet=quiet))
  try(attr(result, "time") <- list(validdate = LFItime(result)))
  result
}

LFIdomain <- function(lfi,quiet=TRUE){
### we assume Lambert projection!
  LFIproj <- meteogrid::trim(int2ascii(LFIread(lfi,"GRID_TYPE",type="integer")))
  if(!quiet) cat("Projection=",LFIproj,"\n")
  if(LFIproj=="CONF PROJ"){
    nx <- LFIread(lfi,"IMAX","integer")
    ny <- LFIread(lfi,"JMAX","integer")
    SW <- c(LFIread(lfi,"LONOR","numeric"), LFIread(lfi,"LATOR","numeric"))
    Lat0 <- LFIread(lfi,"LAT0","numeric")
    Lon0 <- LFIread(lfi,"LON0","numeric")
#  rpk <-  LFIread(lfi,"RPK","numeric")  ### not necessary
#  beta <- LFIread(lfi,"BETA","numeric") ### assumed 0
    xhat <- LFIread(lfi,"XHAT","numeric")
    yhat <- LFIread(lfi,"YHAT","numeric")
    dx <- xhat[2] - xhat[1]
    dy <- yhat[2] - yhat[1]

    projection=list(proj="lcc",lon_0=Lon0,lat_1=Lat0,lat_2=Lat0,R=6371229)
### BUGFIX: the SW given above is in fact the lat/lon of the non-physical point!!!
### so we must move it by (dx,dy) to get the domain right!
    pxy1 <- meteogrid::project(SW,proj=projection,inv=FALSE)
    px1 <- pxy1[1] + dx
    py1 <- pxy1[2] + dy
    xy1 <- meteogrid::project(c(px1, py1), proj = projection, inv = TRUE)
    SW <- c(xy1$x, xy1$y)
    px2 <- px1 + (nx - 1) * dx
    py2 <- py1 + (ny - 1) * dy
    xy2 <- meteogrid::project(c(px2,py2),proj=projection,inv=TRUE)
    NE <- c(xy2$x,xy2$y)
    result <- list(projection=projection,nx=nx,ny=ny,SW=SW,NE=NE,dx=dx,dy=dy)
    class(result) <- c("geodomain",class(result))
    return(result)
  } else if ( LFIproj=="LONLAT REG" ) {
    projection <- list(proj="latlong")
    nx <- LFIread(lfi,"IMAX","integer")
    ny <- LFIread(lfi,"JMAX","integer")
    SW <- c(LFIread(lfi,"LONOR","numeric"), LFIread(lfi,"LATOR","numeric"))
    xhat <- LFIread(lfi,"XHAT","numeric")
    yhat <- LFIread(lfi,"YHAT","numeric")
    dx <- xhat[2] - xhat[1]
    dy <- yhat[2] - yhat[1]
    SW <- SW + c(dx,dy)
    NE <- SW + c((nx-1)*dx,(ny-1)*dy)
    result <- list(projection=projection,nx=nx,ny=ny,SW=SW,NE=NE,dx=dx,dy=dy)
    return(result)
  } else if ( LFIproj=="GAUSS" ) {
    warning("Grid type GAUSS not supported.")
    return(NULL)
  } else if (LFIproj=="CARTESIAN") {
    warning("Grid type CARTESIAN: no geographic information.")
    return(NULL)
  } else if(LFIproj=="NONE"){
    warning("Grid type is NONE.")
    return(NULL)
  } else {
    warning(paste("Unknown grid type",LFIproj,"."))
    return(NULL)
  }
}

LFItime <- function(lfi) {
### ONLY VALID DATE & TIME !!!
### date as 3 integers (y,m,d):
  LFIdat <- LFIread(lfi, "DTCUR%TDATE", type="integer")
### time as 1 real (seconds)
  LFItim <- LFIread(lfi, "DTCUR%TIME", type="numeric")
  validdate <- ISOdatetime(LFIdat[1], LFIdat[2], LFIdat[3], 0, 0, 0) + LFItim
  validdate
### if the file is called AROMOUT_.00xx.lfi we know it is a +xx forecast

}

LFIdec <- function(lfi, field, ...){
  if (!inherits(lfi, "LFIfile")) lfi <- LFIopen(lfi)
  data <- LFIread(lfi, field, add.attr=TRUE, ...)
  if (is.null(data)) return(NULL)
  gridtype <- attr(data, "gridtype")
  nx <- attr(lfi, "domain")$nx
  ny <- attr(lfi, "domain")$ny
  if (gridtype==4 & length(data)==(nx+2)*(ny+2) ){
    result <- matrix(data,ncol=ny+2)[2:(nx+1),2:(ny+1)]
    class(result) <- c("geofield")
    attr(result, "domain") <- attr(lfi,"domain")
    attr(result, "info") <- list(name = attr(data, "name"),
                                 time = attr(lfi, "time"))
    return(result)
  }
  else stop("Not a X_Y grid?")
}

### a plain reading of data field (no interpretation)
LFIread <- function(lfi, field,type="numeric", missing=1.0E+20,
                    add.attr=FALSE, quiet=TRUE){
  if (!inherits(lfi,"LFIfile")) lfi <- LFIopen(lfi)
  if (is.numeric(field)) {
    findex <- field
  }
  else if ( is.na(findex <- match(meteogrid::trim(field),meteogrid::trim(lfi$list$name))) ) {
    grepfield <- grep(field,lfi$list$name,fixed=TRUE)
    if(length(grepfield)==0) {
      warning("non-existant field")
      return(NULL)
    }
    if(length(grepfield)>1) {
      cat("WARNING: non-unique field:\n ",
          paste(lfi$list$name[grepfield],collapse="\n  "),"\n")
      return(NULL)
    }
    findex <- grepfield
  }
  field <- lfi$list$name[findex]
  if (!quiet) cat("Reading ",lfi$list$name[findex],"(",findex,") from ",attr(lfi,"filename"),"\n")
  mpos <- lfi$list$offset[findex]
  mlen <- lfi$list$length[findex]
  a <- file(attr(lfi,"filename"),open="rb")
  on.exit(try(close(a)))
### position is given in words of 8 octets
  seek(a,mpos)
  gridtype <- readBin(a,what="integer",size=8,n=1,endian="big")
  commentlen <- readBin(a,what="integer",size=8,n=1,endian="big")
  if (!quiet) cat("gridtype=",gridtype,"   commentlen=",commentlen,"\n")
  if(commentlen > 0)  {
    icomm <- readBin(a,what="integer",size=8,n=commentlen,endian="big")
    mcomment <- int2ascii(icomm)
    if(!quiet) print(mcomment)
  } else {
    mcomment <- " "
  }
### we assume the data is NOT compressed
  zzz <- readBin(a,what=type,n=mlen/8-2-commentlen,size=8,endian="big")
  if(add.attr) {
    attr(zzz, "mcomment") <- mcomment
    attr(zzz, "gridtype") <- gridtype
    attr(zzz, "name") <- field
  }
### sometimes, missing == 999
### in that case, dx could be > missing.
### so the comparison has to be "==missing"
### Hopefully, this works in all cases.
  zzz[which(zzz==missing)] <- NA
  zzz
}

### a plain replacement of data field (no interpretation)
### just a quick bit of coding.
LFIreplace <- function(lfi, field, data, missing=1.0E+20, quiet=TRUE){
  if (is.character(lfi)) lfi <- LFIopen(lfi)
  if (is.numeric(field)) {
    findex <- field
  } else if ( is.na(findex <- match(meteogrid::trim(field),meteogrid::trim(lfi$list$name))) ) {
    grepfield <- grep(field,lfi$list$name)
    if (length(grepfield)==0) {
      warning("non-existant field")
      return(NULL)
    }
    if (length(grepfield)>1) {
      cat("WARNING: non-unique field:\n ",
          paste(lfi$list[grepfield],collapse="\n  "),"\n")
      return(NULL)
    }
    findex <- grepfield
  }
  field <- lfi$list$name[findex]
  if (!quiet) cat("locating ",lfi$list$name[findex],"(",findex,") from ",attr(lfi,"filename"),"\n")
  pos <- lfi$list$offset[findex]
  len <- lfi$list$length[findex]
  if (!quiet) cat("pos=",pos,"  len=",len,"\n")
  a <- file(attr(lfi,"filename"),open="r+b")
  on.exit(try(close(a)))
### position is given in words of 8 octets
  seek(a, where=pos, rw="read")
  if (!quiet) cat("rpos=",seek(a,rw="read"),"\n")

  gridtype <- readBin(a,what="integer", size=8, n=1, endian="big")
  commentlen <- readBin(a,what="integer", size=8, n=1, endian="big")

  datalen <- len/8 - 2 - commentlen
  if (!quiet) cat("gridtype=",gridtype,"   commentlen=",commentlen,"\n")
  if (!quiet) cat("Expected data length:",datalen,"\n")
  nx <- attr(lfi,"domain")$nx
  ny <- attr(lfi,"domain")$ny

  if (gridtype==4) {
    if (!quiet) cat("Expected data dimensions:",nx," x ",ny,"\n")

    if (dim(data)[1] == nx & dim(data)[2]==ny){
      data2 <- matrix(missing,nrow=nx+2,ncol=ny+2)
      data2[2:(nx+1),2:(ny+1)] <- data
      newdata <- as.vector(data2)
    } else if(dim(newdata)[1] == nx+2 & dim(newdata)[2]==ny+2){
      newdata <- as.vector(data)
    } else stop("New data does not have right dimensions.")
  } else {
    if (length(data) != datalen) stop("New data does not have right length.")
    newdata <- data
  }

  if (commentlen > 0) {
    icomm <- readBin(a,what="integer",size=8,n=commentlen,endian="big")
    mcomment <- int2ascii(icomm)
    if (!quiet) print(mcomment)
  }
### we assume the data is NOT compressed
  if (length(newdata) != datalen) stop("new data has wrong length!")
  data[is.na(newdata)] <- missing
### now we set the "write" pointer to the same place as the "read" pointer
  wpos <- seek(a)
  if (!quiet) cat("Position: ",wpos,"\n")
  seek(a, rw="write", where=wpos)
  writeBin(as.vector(newdata), a, size=8, endian="big")
}
