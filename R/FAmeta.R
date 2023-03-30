### because the frame and time "objects" are so different,
### I decided not to read them with standard FAdec.
### This is even independent from FAopen !
FAread_meta <- function(filename, archname=NULL, quiet=TRUE){
  if (inherits(filename,"FAfile")) { # it's an FA object, not a filename
    archname <- attr(filename,"tarfile")
    tar.offset <- attr(filename,"tar.offset")
    filename <- attr(filename,"filename")
  } else {
    if (is.null(archname)) tar.offset <- 0
    else tar.offset <- FindInTar(archname, filename)
  }
   
  on.exit(try(close(ff)))
  if (is.null(archname)) ff <- file(filename, open="rb") else ff <- file(archname, open="rb")
  seek(ff, tar.offset)

#-- in fact the second value gives the length of this section. But it is always 22.
  FAheader <- readBin(ff, what="integer", n=22, size=8, endian="big")
  blocksize <- FAheader[1]*8
  nfields <- FAheader[6]
  nholes <- FAheader[21]
  if (nfields - nholes < 7) stop("File contains less than 7 data sectors, so no complete frame.")
  seek(ff, tar.offset + blocksize)
### it is possible to have 'holes' in the meta section!
### so we must be prepared for that.
### let's read to 8+nholes (if possible)
  nread <- min(8+nholes, nfields, blocksize/16)
  allnames <- rawToChar(readBin(ff, "raw", n=nread*16))
  FAnames.all <- substring(allnames, seq(1, nread*16, by=16), seq(16, 16*(nread+1), by=16))
  is.hole <- (FAnames.all=="                ")
  FAnames <- FAnames.all[!is.hole]
  if (FAnames[1] != "CADRE-DIMENSIONS") stop("FRAME not at expected location. Not a regular FA file.")
  if (FAnames[2] != "CADRE-FRANKSCHMI") stop("FRAME not at expected location. Not a regular FA file.")
  if (FAnames[3] != "CADRE-REDPOINPOL") stop("FRAME not at expected location. Not a regular FA file.")
  if (FAnames[4] != "CADRE-SINLATITUD") stop("FRAME not at expected location. Not a regular FA file.")
  if (FAnames[5] != "CADRE-FOCOHYBRID") stop("FRAME not at expected location. Not a regular FA file.")
  if (FAnames[7] != "DATE-DES-DONNEES") stop("FRAME not at expected location. Not a regular FA file.")
# check for new DATX_DES_DONNEES field
  if (length(FAnames)>7 && FAnames[8]=="DATX-DES-DONNEES") {
    if (!quiet) cat("Found DATX-DES-DONNEES. Adding to meta data.\n")
    nmeta <- 8
  } else nmeta <- 7
  FAnames <- FAnames[1:nmeta]
  seek(ff, tar.offset+2*blocksize)
  FAlen <- numeric(nmeta)
  FApos <- numeric(nmeta)
  j <- 1
  for (i in 1:nmeta) {
    ## skip any holes
    while( is.hole[j] ) { 
      readBin(ff, what="integer", size=8, n=2, endian="big")
      j <- j+1
    }
    FAlen[i] <- readBin(ff, what="integer", size=8, n=1, endian="big")
    FApos[i] <- readBin(ff, what="integer", size=8, n=1, endian="big")
    j <- j+1
  }
  
  metalist <- list()
  for (i in 1:nmeta) {
    seek(ff,tar.offset+(FApos[i]-1)*8)
    fatype <- ifelse(is.element(i, c(1,3,6,7,8)), "integer", "numeric")
    metalist[[i]] <- readBin(ff, what=fatype, n=FAlen[i], size=8, endian="big")
  }
  names(metalist) <- FAnames ## this also saves the name of the frame, which is in FAnames[6]
  attr(metalist,"pos") <- FApos
  attr(metalist,"len") <- FAlen
  metalist
}

FAread_header <- function(fa) {
  if (!inherits(fa, "connection")) {
    if (inherits(fa,"FAfile")) filename <- attr(fa,"filename")
    else if (is.character(fa)) filename <- fa
    else stop("Not a valid file.")
    fa <- file(filename,open="rb")
    on.exit(close(fa))
  }
  seek(fa, 0, rw="read")
  header <- readBin(fa,what="integer",n=22,size=8,endian="big")
  header
}

FAheader.date <- function(header, quiet=TRUE){
# set modification date (not important at all...)
  now <- Sys.time()
  now.date <- as.integer(format(now,"%Y%m%d"))
  now.time <- as.integer(format(now,"%H%M%S"))
  if (!quiet) cat("Modification time:",now.date,now.time,"\n")
  header[16] <- now.date
  header[17] <- now.time
  if (header[18]==header[14] & header[19]==header[15]){
    header[18] <- now.date
    header[19] <- now.time
  }
  return(header)
}

FAwrite_header <- function(fa, header=NULL){
  if (inherits(fa,"FAfile") && is.null(header)) header <- attr(fa,"header")
  if (is.null(header) || length(header) != 22) stop("No header provided or incorrect length.")
  if (!inherits(fa, "connection")) {
    if (inherits(fa,"FAfile")) filename <- attr(fa,"filename")
    else if (is.character(fa)) filename <- fa
    else stop("Not a valid file.")
    fa <- file(filename, open="wb")
    on.exit(close(fa))
  }
  seek(fa, 0, rw="write")
  writeBin(as.integer(header), fa, size=8, endian="big")
}


### modify e.g. the date information:
### first read all meta files using FAread_meta()
### then write with this function
FAwrite_meta <- function(filename, metadata) {
  if (!is.character(filename)) filename <- attr(filename,"filename")
  oldmeta <- FAread_meta(filename)
  nmeta <- length(oldmeta)

  ff <- file(filename, open="r+b")
  on.exit(try(close(ff)))

  pos <- attr(oldmeta,"pos")
  len <- attr(oldmeta,"len")
  for (nm in names(metadata)) {
    i <- match(nm, names(oldmeta))
    if (is.na(i)) {
      if (length(metadata[[nm]])==1 && ! names(oldmeta)[6] %in% names(metadata)) {
### Change the name of the frame itself (i.e. the name of meta field 6)
        i <- 6
        newname <- format(nm, width=16)
        seek(ff, 0, rw="read")
        header <- readBin(ff,what="integer",n=22,size=8,endian="big")
        blocksize <- header[1]*8
        seek(ff, blocksize+16*5, rw="write")
        writeChar(newname, ff, eos=NULL)
      } else stop(paste(nm,"not found."))
    }
    if (length(metadata[[nm]]) != length(oldmeta[[i]])) stop(paste(nm,": length has changed."))
    seek(ff,(pos[i]-1)*8,rw="write")
    numtype <- if (i %in% c(1,3,6,7,8)) as.integer else as.numeric
    writeBin(numtype(metadata[[nm]]),ff,size=8,endian="big")
  }
}

### DANGEROUS: renaming meta-data (only frame name and DATX-DES-DONNEES are reasonably safe)
FArename.meta <- function(fa, field, newname, quiet=TRUE) {
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa,quiet=quiet)
  if (!is.null(attr(fa,"tarfile"))) stop("Manipulation of files in an archive is forbidden. Read only!")
  oldmeta <- FAread_meta(fa)
# only frame name may be changed by giving 'field' as a number
  if (field==6) field <- names(oldmeta)[6]
  if (field==8) field <- "DATX-DES-DONNEES"
  i <- match(field, names(oldmeta)) 
  if (is.na(i)) stop(paste(field,"not found."))

  if (nchar(newname)<5) stop("Field name must have at least 5 characters.")
  newname <- format(newname, width=16)

  ff <- file(attr(fa, "filename"), open="r+b")
  on.exit(try(close(ff), silent=TRUE))
  header <- attr(fa, "header")
  blocksize <- header[1]*8
  if (!quiet) cat("Setting namefield at location",blocksize+16*(i-1),"\n")
  seek(ff, blocksize+16*(i-1), rw="write")
  writeChar(newname, ff, eos=NULL)
  close(ff)
  return(FAopen(attr(fa,"filename")))
}

##########################
### DECODING META DATA ###
##########################

FAtime <- function(metadata){
  if (is.numeric(metadata)) { # it's a vector of 11 (or 22) integers
    timelist <- metadata[1:11]
    timex <- if (length(metadata)<22) NULL else metadata[12:22]
  } else {
    timelist <- metadata$"DATE-DES-DONNEES"
    timex <- metadata$"DATX-DES-DONNEES"
  }
  if (!is.null(timex) && timex[1]==1) {
    extended <- TRUE
  } else {
    extended <- FALSE
    timex <- rep(0,11)
  }

### timelist and the possible extended timex are vectors of 11 integers each
### if there is extended time data, we must check if it is relevant: if lead time is more exact than that in hours.
  basedate <- ISOdatetime(timelist[1], timelist[2], timelist[3], timelist[4],
                          timelist[5], timex[3], tz="UTC")
  validdate <- basedate
  if (extended) ldt <- timex[4]
  else ldt <- timelist[7]*3600

  if (timelist[9] == 10) {
    if (ldt %% 3600 == 0) forecast <- sprintf("+%i%s", ldt/3600, ifelse(timelist[6]==1,"h","?"))
    else forecast <- sprintf("+%is", ldt)
    validdate <- basedate + ldt
  } else {
    if (timelist[9] == 1) forecast <- "Initialized"
    else forecast <- "Uninitialized"
  }
  result <- list(
#  result <- paste(format(basedate,"%Y/%m/%d z%H:%M"),forecast)
# NOTE: ldt should be in seconds from now on...
               leadtime    = ldt / 3600,
               basedate    = basedate,
               validdate   = validdate,
               forecast    = forecast,
               initialised = (timelist[9] == 1)
             )
  if (extended) result$tstep <- timex[7]
  result
}

FAframe <- function(metadata){
  frame <- list()
  FAdim <- metadata$"CADRE-DIMENSIONS"
  FAfs  <- metadata$"CADRE-FRANKSCHMI"
  FArpp <- metadata$"CADRE-REDPOINPOL"
  FAsll <- metadata$"CADRE-SINLATITUD"
  FAlev <- metadata$"CADRE-FOCOHYBRID"
#-- the sixth field does not have a fixed name
#-- it is the name of the frame, and the content is 1 integer (==1, the FA version)
  FAname <- names(metadata)[6]
  FAver  <- metadata[[6]]
  if (FAver != 1) warning(paste("FA version is",FAver))
  FAtype=ifelse(FAdim[5]<0,"aladin","arpege")
  if (FAtype=="aladin"){
#-- "old style" aladin geometry (up to about 2006)
    if (FAsll[1]>=0){ 
      frame$nroteq <- FAsll[1]
      frame$rpk    <- FAsll[10]
      frame$lonr   <- FAsll[2]*180/pi
      frame$latr   <- FAsll[3]*180/pi
      frame$lon0   <- FAsll[8]*180/pi
      frame$lat0   <- FAsll[9]*180/pi
      frame$lon1   <- FAsll[4]*180/pi
      frame$lat1   <- FAsll[5]*180/pi
      frame$lon2   <- FAsll[6]*180/pi
      frame$lat2   <- FAsll[7]*180/pi
      frame$delx   <- FAsll[15]
      frame$dely   <- FAsll[16]
      frame$exwn   <- FAsll[17]
      frame$eywn   <- FAsll[18]
      frame$lx     <- FAsll[13]
      frame$ly     <- FAsll[14]
      frame$nsotrp <- FAsll[11]
      frame$ngivo  <- FAsll[12]
    } else {
#-- newer version:
      frame$nroteq <- FAsll[1]
      frame$rpk    <- FAsll[2]
      frame$lon0   <- FAsll[3]*180/pi
      frame$lat0   <- FAsll[4]*180/pi
      frame$lonc   <- FAsll[5]*180/pi
      frame$latc   <- FAsll[6]*180/pi
      frame$lon1   <- FAsll[13]*180/pi
      frame$lat1   <- FAsll[14]*180/pi
      frame$lon2   <- FAsll[15]*180/pi
      frame$lat2   <- FAsll[16]*180/pi
      frame$delx   <- FAsll[7]
      frame$dely   <- FAsll[8]
      frame$lx     <- FAsll[9]
      frame$ly     <- FAsll[10]
      frame$xwn   <- FAsll[11]
      frame$ywn   <- FAsll[12]
    }
### for LatLon grid, also delx etc are in radians !
    if (frame$rpk < 0){
      frame$delx <- frame$delx*180/pi
      frame$dely <- frame$dely*180/pi
      frame$lx <- frame$lx*180/pi
      frame$ly <- frame$ly*180/pi
    }
    frame$nsmax  <- FAdim[1]
    frame$ndgl   <- FAdim[2]
    frame$ndlon  <- FAdim[3]
    frame$nlev   <- FAdim[4]
    frame$nmsmax <- -FAdim[5]

    frame$ndlun <- FArpp[3]
    frame$ndlux <- FArpp[4]
    frame$ndgun <- FArpp[5]
    frame$ndgux <- FArpp[6]
#    frame$nbzonl<- FArpp[7] 
#    frame$nbzong<- FArpp[8]
  } else {
    cat("This appears to be an Arpege file.\n")
    frame$nsmax  <- FAdim[1]
    frame$ndgl   <- FAdim[2]
    frame$ndlon  <- FAdim[3]
    frame$nlev   <- FAdim[4]
    frame$nsttyp <- FAdim[5]
## stretching
    frame$latsin <- FAfs[1]
    frame$loncos <- FAfs[2]
    frame$lonsin <- FAfs[3]
    frame$stretch <- FAfs[4]

  }

  frame$FAdim <- FAdim
  frame$FAfs  <- FAfs
  frame$FArpp <- FArpp
  frame$FAsll <- FAsll
  frame$FAlev <- FAlev

  frame$levels <- list(
    refpressure = frame$FAlev[1],
    A = frame$FAlev[2:(frame$nlev+2)],
    B = frame$FAlev[(frame$nlev+3):(2*(frame$nlev+1)+1)]
    )

  frame$FAname <- FAname
  frame$FAver  <- FAver
  frame$FAtype <- FAtype
  class(frame)=c("FAframe",class(frame))
  frame
}

FAdomain <- function(faframe,quiet=TRUE){
### Create the "geodomain" specification for a given FA frame
  if (!inherits(faframe,"FAframe")) stop("Not a FAframe object.")
  if (faframe$FAtype != "aladin") stop("Only aladin type supported for now!")  

### FArpp(2) : -1 means gridpoint on C+I+E, +1 spectral
### For the physical domain: no extension zone !!!
  nx <- faframe$ndlux - faframe$ndlun + 1
  ny <- faframe$ndgux - faframe$ndgun + 1 
  if (!quiet) cat("C+I dim:",nx,"x",ny,"\n")
  if (!quiet) cat("with extension zone:",faframe$ndgl,"x",faframe$ndlon,"\n")
  if (faframe$rpk == -9){
    projtype <- list(proj="latlong")
  } else if (faframe$nroteq == -2 || faframe$lat0 == 0) {
    projtype <- "mercator"
  } else if (faframe$rpk == 1) {
    projtype <- "stere"
  } else {
    projtype <- "lambert"
  }

  domain <- meteogrid::Make.domain(
              projtype = projtype, reflon = faframe$lon0, reflat = faframe$lat0,
              nxny = c(nx,ny), dxdy = c(faframe$delx, faframe$dely),
              clonlat = c(faframe$lonc, faframe$latc), 
              exey = c(faframe$ndlon - nx, faframe$ndgl-ny),
              tilt = faframe$lon0,
              earth = list(R = 6371229))

## for old FA files that don't have clonlat:
## in the future: (re-)calculate clonlat, dxdy (ref DomainExtent) rather than going back to NE, SW
  if (is.null(domain$clonlat)) {
#    xy <- meteogrid::project(list(x = c(faframe$lon1, faframe$lat1), y = c(faframe$lon2,faframe$lat2)),
#                       proj = domain$projection)
    domain$clonlat <- NULL
    domain$SW <- c(faframe$lon1, faframe$lat1)
    domain$NE <- c(faframe$lon2, faframe$lat2)
  }
  domain
}

