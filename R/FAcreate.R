###################################
###     FAcreate                ###
### create a new FA file        ###
###################################

FAmake.header <- function(sector.size=1000, nsector=4){
  now <- Sys.time()
  now.date <- as.integer(format(now,"%Y%m%d"))
  now.time <- as.integer(format(now,"%H%M%S"))
# default header: blocksize 
  header <- rep(0,22)
  header[1]  <- sector.size # sector size (in 8 byte words, not bytes)
  header[2]  <- 16   # word size (always 16)
  header[3]  <- 0    # correct closure (always 0)
  header[4]  <- 22   # header size (always 22)
  header[5]  <- nsector  # nr of sectors (at least 3 for meta-data)
  header[6]  <- 0    # nr of fields (FA frame needs 7!)
  header[7]  <- 0    # shortest field length (almost always 1)
  header[8]  <- 0    # longest field length
  header[9]  <- 0    # data length
  header[10] <- 0    # nr. of rewrites
  header[11] <- 0    # nr. of shorter rewrites
  header[12] <- 0    # nr. of longer rewrites
  header[13] <- header[1]/2    # fields per sector
  header[14] <- now.date # creation
  header[15] <- now.time
  header[16] <- now.date # first mod.
  header[17] <- now.time
  header[18] <- now.date # last mod.
  header[19] <- now.time
  header[20] <- 1    # index sectors (always 1?)
  header[21] <- 0    # nr of holes
  header[22] <- nsector    # sectors? what is different
  header
}

### in the future, extended will be set to TRUE by default
FAmake.time <- function(fcdate="1970010100", ldt=0, unit='h',init=1,
                        previous=0, first=0, extended=FALSE, tstep=1)
{
  # TODO: go for the long date specification (up to seconds) as default?
  # FIXME: this needs a bit of work!
  if (is.character(fcdate) | is.numeric(fcdate)) fcdate <- as.POSIXct(fcdate, format="%Y%m%d%H")

  fatime <- numeric(11)
  fatime[1] <- as.numeric(format(fcdate,'%Y'))
  fatime[2] <- as.numeric(format(fcdate,'%m'))
  fatime[3] <- as.numeric(format(fcdate,'%d'))
  fatime[4] <- as.numeric(format(fcdate,'%H'))
  fatime[5] <- as.numeric(format(fcdate,'%M'))
  fatime[6] <- switch(unit,"h"=1,"m"=2,0)
  fatime[7] <- round(ldt/switch(unit,"h"=3600,"m"=60,1))
  fatime[8] <- 0 # ???
  if (ldt > 0) fatime[9] <- 10
  else fatime[9] <- init
  fatime[10] <- round(previous/switch(unit,"h"=3600,"m"=60,1))
  fatime[11] <- 0 # ???
## since cy40, you can have a second vector with extended time information
  if (extended) {
    fxtime <- numeric(11)
    fxtime[1] <- 1
    fxtime[2] <- 0
    fxtime[3] <- as.numeric(format(fcdate,'%S'))
    fxtime[4] <- ldt
    fxtime[5] <- previous
    fxtime[6] <- first
    fxtime[7] <- tstep
    fxtime[8] <- 0
    fxtime[9] <- 0
    fxtime[10] <- 0
    fxtime[11] <- 0
    fatime <- c(fatime, fxtime)
  }
  as.integer(fatime)
}

# create a new file with the same settings and frame, copy some fields
# if "newfile" exists and overwrite==FALSE we just copy the fields
FAduplicate <- function(fa, newfile, fields=character(0), overwrite=FALSE) {
  if (!inherits(fa, "FAfile")) fa <- FAopen(fa)
  if (inherits(newfile, "FAfile")) newfile <- attr(newfile, "filename")

  header <- FAread_header(fa)
  if (!file.exists(newfile) || overwrite) {
    metadata <- FAread_meta(fa)
#  frame <- FAframe(metadata)
#  time <- FAtime(metadata)
  # create with small number of sectors
    newfa <- FAcreate(newfile, sector.size=header[1], nsector=5, frame=metadata, overwrite=overwrite)
  } else {
    # don't create, just open the existing file
    newfa <- FAopen(newfile)
  }

  # copy the frame

  # copy other fields
  for (ff in fields) {
    # just read the raw binary data and write it to the new file
    # allow for more general "grep" like "^S...TEMPE" FAgrep()
    pos <- FAgrep(fa, ff, as_name=FALSE)
    if (length(pos) > 0 && !is.na(pos)) {
      for (pp in pos) { # there may be several fields matching! (e.g. "SURFPREC"
        x <- FAread_msg(fa, fpos=fa$list$offset[pp], flen=fa$list$length[pp])
        newfa <- FAenc(newfa, fieldname=fa$list$name[pp], data=x)
      }
    } else {
      warning("field ", ff, " not found.")
    }
  }
  invisible(newfa)
}

## TO DO: time format may be as a list of two vectors
FAcreate <- function(filename, frame, time=NULL, sector.size=1000, 
                     nsector=10, overwrite=FALSE)
{
  if (file.exists(filename)) {
    if (!overwrite) stop("File already exists.")
    else {
      warning("Overwriting existing file!",call.=FALSE,immediate.=TRUE)
      file.remove(filename)
    }
  }
#--
#-- 1. fix the header
#--
  header <- FAmake.header(sector.size, nsector)
  blocksize <- sector.size*8

#--
#-- 2. Create the file itself
#--

  ff <- file(filename, "wb")
  on.exit(try(close(ff), silent=TRUE))

  ## header section
  seek(ff, 0, rw="write")
  writeBin(as.integer(header), ff, size=8, endian="big")
  writeBin(as.integer(rep(0, sector.size-22)), ff, size=8, endian="big")
  ## name section:
  writeChar(rep("**FIN D'INDEX **", header[13]), ff, eos=NULL)
  ## address section:
  writeBin(as.integer(rep(0, sector.size)), ff, size=8, endian="big")
  ## data sections
  writeBin(as.integer(rep(0, sector.size*(nsector-3))), ff, size=8, endian="big")

#--
#-- 3. Write the meta data if provided (frame and date/time)
#-- make sure the number of sectors is enough for the frame!
#-- "frame" could also be a metadata list
  if (!is.null(frame)) {
    if (inherits(frame, "FAframe")) {
      if (is.null(time)) time <- FAmake.time()
      if (is.character(time)) time <- FAmake.time(fcdate=time)
      if (length(time)==11) nmeta <- 7 else nmeta <- 8
      framename <- format(frame$name, width=16)
      metadata <- list("CADRE-DIMENSIONS" = frame[["FAdim"]],
                       "CADRE-FRANKSCHMI" = frame[["FAfss"]],
                       "CADRE-REDPOINPOL" = frame[["FArpp"]],
                       "CADRE-SINLATITUD" = frame[["FAsll"]],
                       "CADRE-FOCOHYBRID" = frame[["FAlev"]],
                       "FRAME-NAME      " = 1,
                       "DATE-DES-DONNEES" = time[1:11])
      if (nmeta==8) metadata <- c(metadata, "DATX-DES-DONNEES" = time[12:22])
      names(metadata)[6] <- framename
    } else {
      # we assume "frame" actually contains the full metadata
      metadata <- frame
      nmeta <- length(metadata)
    }
    metalengths <- vapply(metadata, length, 1)
    totlength <- sum(metalengths)
    if (totlength > nsector * sector.size) {
      warning("metadata needs more sectors than available.")
      nsector <- ceil(totlength / sector.size)
    }
## the name sector
    seek(ff, blocksize, rw="write")
    writeChar(names(metadata), ff, eos=NULL)

## the addres/length sector
    len <- vapply(metadata, length, 1)
    pos <- 3*sector.size + c(0, cumsum(len)[1:(nmeta-1)]) + 1
#    print(pos)
#    print(len)
    seek(ff, 2*blocksize, rw="write")
    for(i in 1:nmeta) writeBin(as.integer(c(len[i], pos[i])), ff, size=8, endian="big")

## the data itself
    for(i in 1:nmeta){
      seek(ff, (pos[i]-1)*8, rw="write")
      numtype <- if (i %in% c(1,3,6,7,8)) as.integer else as.numeric
      writeBin(numtype(metadata[[i]]), ff, size=8, endian="big")
    }
## fix the header for field number
    header[6] <- nmeta
    header[7] <- min(metalengths)
    header[8] <- max(metalengths)
    header[9] <- totlength
    seek(ff, 5*8, rw="write")
    writeBin(as.integer(header[6:9]), ff, size=8, endian="big")
  }
#--
#-- 4. wrap-up
#--
  close(ff)
  ## don't try to FAopen() if there is no frame? there is no data yet, anyway
  invisible(FAopen(filename))
}

##############################################
###    FAmake.frame                        ###
### create a new "frame" from domain input ###
##############################################
FAmake.frame = function(domain, extension=c(11,11), relaxation=c(8,8), nsmax=NULL, nmsmax=NULL,
                        lineargrid=TRUE, levels=NULL, name="FA-FRAME        ",
                        sptrunc=10){
  warning("BEWARE: FAmake.frame is new and mostly untested!")
  if (!inherits(domain,"geodomain")) domain <- attributes(domain)$domain
  frame <- list(FAname=name,FAver=1,FAtype="aladin")
 
  fulldim <- c(domain$nx,domain$ny) + extension

# 1. Dimensions ("CADRE-DIMENSIONS")
  ndlon <- fulldim[1]
  ndgl <- fulldim[2]
  if (is.null(levels)) nlev <- 1
  else nlev <- length(levels$A) - 1

  if (is.null(nmsmax)) nmsmax <- ifelse(lineargrid,floor((fulldim[1]-1)/2),floor((fulldim[1]-1)/3))
  if (is.null(nsmax)) nsmax <- ifelse(lineargrid,floor((fulldim[2]-1)/2),floor((fulldim[2]-1)/3))
  FAdim <- c(nsmax,ndgl,ndlon,nlev,-nmsmax)
  frame$nsmax  <- FAdim[1]
  frame$ndgl   <- FAdim[2]
  frame$ndlon  <- FAdim[3]
  frame$nlev   <- FAdim[4]
  frame$nmsmax <- -FAdim[5]
  frame$FAdim <- FAdim

# 2. Vertical levels ("CADRE-FOCOHYBRID")
  if (is.null(levels)) levels <- list(refpressure=102500,A=0,B=1)
  frame$levels <- levels
  frame$FAlev <- c(levels$refpressure,levels$A,levels$B)

# 3. Projection ("CADRE-SINLATITUD")
  proj <- attributes(domain)$projection
  domex <- meteogrid::DomainExtent(domain)

  lonc <- domex$clonlat[1]
  latc <- domex$clonlat[2]
  lon1 <- domex$SW[1]
  lat1 <- domex$SW[2]
  lon2 <- domex$NE[1]
  lat2 <- domex$NE[2]
  if(proj$proj=="lcc"){
    nroteq <- -1
    rpk  <- sin(lat0)
    lon0 <- proj$lon_0/180*pi
    lat0 <- proj$lat_1/180*pi
    delx <- domain$dx
    dely <- domain$dy
  }
  else if (proj$proj=="latlong"){
    nroteq <- -1
    rpk  <- -9
    lon0 <- 0
    lat0 <- 0
    delx <- domain$dx/180*pi
    dely <- domain$dy/180*pi
  }
  else stop("Only lambert and latlong allowed.")

  lx <- delx*ndlon
  ly <- dely*ndgl
  xwn <- 2*pi/lx
  ywn <- 2*pi/ly

## FAsll has length 18
  FAsll <- c(nroteq, rpk, lon0, lat0, lonc, latc, delx, dely,
             lx, ly, xwn, ywn, lon1, lat1, lon2, lat2, 0, 0)
  frame$nroteq <- FAsll[1]
  frame$rpk    <- FAsll[2]
  frame$lon0   <- FAsll[3]*180/pi
  frame$lat0   <- FAsll[4]*180/pi
  frame$lonc   <- FAsll[5]*180/pi
  frame$latc   <- FAsll[6]*180/pi
  frame$delx   <- FAsll[7]
  frame$dely   <- FAsll[8]
  frame$lx     <- FAsll[9]
  frame$ly     <- FAsll[10]
  frame$xwn    <- FAsll[11]
  frame$ywn    <- FAsll[12]
  frame$lon1   <- FAsll[13]*180/pi
  frame$lat1   <- FAsll[14]*180/pi
  frame$lon2   <- FAsll[15]*180/pi
  frame$lat2   <- FAsll[16]*180/pi
  frame$FAsll  <- FAsll

# 4. grid properties ("CADRE-REDPOINPOL")
  ndlun <- 1  # not the most general...
  ndgun <- 1
  ndlux <- ndlun + domain$nx - 1
  ndgux <- ndgun + domain$ny - 1
  nbzonl <- relaxation[1]
  nbzong <- relaxation[2]

  FArpp <- numeric(8 + 2*nsmax + 4)
  if (any(extension>0)) extzone <- 1 else extzone <- 0
## extzone <- -1 also has a meaning, something like "only gridpoint, no spectral" but I never saw such a file..
## maybe typical for fullpos? 
  FArpp[1:8] <- c(sptrunc,extzone,ndlun,ndlux,ndgun,ndgux,nbzonl,nbzong)

## the rest of the vector is something completely different
  nnn <- (0:nsmax)/nsmax
  tt1 <- cumsum(floor(nmsmax*sqrt(abs(1-nnn^2))+1E-10)*4) # this gives half of the numbers, the end points of segments
  tt0 <- c(0,tt1[1:(nsmax-1)]+1)
# now merge the two vectors
  FArpp[seq(9 ,12 + 2*nsmax,by=2)] <- tt0
  FArpp[seq(10,12 + 2*nsmax,by=2)] <- tt1

  frame$ndlun <- FArpp[3]
  frame$ndlux <- FArpp[4]
  frame$ndgun <- FArpp[5]
  frame$ndgux <- FArpp[6]
  frame$nbzonl <- FArpp[7] 
  frame$nbzong <- FArpp[8]
  frame$FArpp <- FArpp

# 5. ("CADRE-FRANKSCHMI")
# hardly ever used for ALADIN: just 4 zeros.
  FAfs <- rep(0,4)
  faframe$FAfs <- FAfs

# wrap-up
  class(frame)=c("FAframe",class(frame))
  frame
}


