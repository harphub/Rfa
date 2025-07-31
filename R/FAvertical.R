#######################################
FApressure <- function(fa, lev, SP=FAdec(fa, "SURFPRESSION    ") ){
## returns a field with the local pressure at a given hybrid level
## needs surface pressure (a geofield or just a vector, matrix... of values)
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa, lparse=FALSE)
  if (attr(fa,"frame")$nlev < lev) stop(paste("Only",attr(fa,"frame")$nlev,"levels available"))
### usually, surface pressure is stored as logarithm
  if (SP[1,1]<50) SP <- exp(SP)
  pref <- attr(fa,"frame")$levels$refpressure
### Half level pressures below and above:
  ph1 <- pref * attr(fa,"frame")$levels$A[(lev+1)] + SP * attr(fa,"frame")$levels$B[(lev+1)]
  ph0 <- pref * attr(fa,"frame")$levels$A[lev] + SP * attr(fa,"frame")$levels$B[lev]

### this formula follows the FullPos 'Scientist guide' documentation (Ryad El Katib, 2002)
### REFERENCE -> GPPREH, GPXYB, GPREF 
### assuming LAPRXPK=.F. (if .T. -> simpler formula without alpha, just (ppl1+ppl2)/2 )
  alpha <- if (lev==1) 1 else 1-ph0/(ph1-ph0)*log(ph1/ph0)
  pp <- ph1 * exp(-alpha)
  if (inherits(SP, "geofield")) {
    attributes(pp)$info$name <- sprintf("S%03iPRESSURE",lev)
    attributes(pp)$info$level <- sprintf("Hybrid level %i",lev)
  }
  pp
}

###########################

FAslice <- function(fa, par="TEMPERATURE", axis="X", n=1, type="S",
                    plot=FALSE, plot.function=contour, ...){
### note that "type" is the kind of INPUT fields you want to read!

  if (!inherits(fa,"FAfile")) fa <- FAopen(fa)

  if (axis=="X") {
    npoints <- attr(fa,"domain")$nx
    get.data <- function(fa,fld) FAdec(fa, fld)[,n]
  } else {
    npoints <- attr(fa,"domain")$ny
    get.data <- function(fa,fld) FAdec(fa, fld)[n,]
 }

  if (type=="S") {
    nlev <- attributes(fa)$frame$nlev
    get.name <- function(lev, par) sprintf("S%03i%s", lev, par)
  } else {
# TO DO: FIX pressures >= 100000
    pvar <- grep(paste("^P[0-9]{5}",substr(par,1,10), sep=""), fa$list$name, value=TRUE)
    plev <- sort(as.numeric(unique(substr(pvar,2,6)))) # may be different levels than 'plevels.out'!
    nlev <- length(plev)
    get.name <- function(lev, par) sprintf("P%05i%s", plev[lev], substr(par,1,10))
  }

# now we fill an array with the 'slice' data
  result <- array(NA,dim=c(npoints, nlev))
  for (ll in 1:nlev){
    field <- get.name(ll,par)
    if (!is.na(FAfind(fa,field))) result[,ll] <- get.data(fa,field)
  }

#  colnames(result) <- if (type=="S") 1:nlev else plev
  mylevels <- if (type=="S") 1:nlev else plev

  if (plot) {
    plot.function(1:npoints,mylevels,result,ylim=rev(range(mylevels)), asp=NULL, ...)
#    title(ylab="level")
#    Axis(mylevels,side=2)
  }
  invisible(result)
}

################

FApressures.local <- function(faframe, surfpressure){
## returns for a single point (identified only by the locat pressure)
## a vector of pressures for the (full) hybrid levels
## the code follows the FullPos 'Scientist guide' documentation (Ryad El Katib, 2002)
  if(!inherits(faframe,"FAframe")) faframe <- attr(faframe,"frame")

  nlev <- faframe$nlev
  A <- faframe$levels$A
  B <- faframe$levels$B
  pref <- faframe$levels$refpressure
  if (surfpressure < 10000) surfpressure <- exp(surfpressure) 
#stop("Surface pressure must be in Pa! Minimum is 100 hPa = 10000 Pa")

  result <- .C("fa_pressures", A=A, B=B, pref=pref, nlev=as.integer(nlev),
                psurf=surfpressure, pressure=numeric(nlev))
  exp(result$pressure)
}

FAsounding <- function(fa, par="TEMPERATURE", lon=NULL, lat=NULL, index=NULL, id=NULL,
                       levtype="S",
                       plevels.out=NULL, method="bilin"){
### either give lon/lat or index=cbind(i,j)
### if there are multiple locations, you should give a code (id)
### to distinguish them in the resulting table
### method can be "closest" "bilin" "bicubic"
### the output vertical co-ordinate is always pressure

### TO: "index" may contain non-integer values, so may also need interpolation

  if (!inherits(fa,"FAfile")) fa <- FAopen(fa)
  npar <- length(par)

  do.interp <- !is.null(lon)
  if (do.interp) {
    weights <- meteogrid::point.interp.init(domain=attr(fa,"domain"), lon=lon, lat=lat,
                                            method=method)
    if (any(is.na(weights))) {
      bad <- unique(which(is.na(weights), arr.ind=TRUE )[,"row"])
      stop("points outside domain:\n",paste(bad,":", lon[bad],",", lat[bad],"\n"))
    }
    npoints <- length(lon)
# a function that does the interpolation or index selection:
    get.data <- function(fa,fld) meteogrid::point.interp(
                    infield=FAdec(fa,fld), weights=weights, method=method)
  } else if (!is.null(index)) {
    if (is.list(index)) index <- cbind(index$i, index$j)
    npoints <- dim(index)[1]
    get.data <- function(fa, fld) FAdec(fa, fld)[index]
  } else {
    stop("Either lon/lat or index must be given!")
  }

  LID <- !is.null(id)
  if (!LID && npoints>1) {
    if (!is.null(index)) {
      id <- apply(index, 1, paste,collapse=".") # use index as label
    } else id <- seq_along(lon) # just number the points
    LID <- TRUE
  }

  if (levtype == "S") { # hybrid model levels
    nlev <- attributes(fa)$frame$nlev
    surfpres <- get.data(fa, "SURFPRESSION")
    if (any(surfpres< 10000)) surfpres <- exp(surfpres)
## TODO: find the minumum and maximum hybrid levels needed for interpolation to plevels.out
##       -> decode only the necessary fields.
## TODO: if npar>1, try to do all interpolations in 1 call 
##       -> only calculate pressures and interpol weights once
    if (!is.null(plevels.out)) {
      result <- data.frame("p"=rep(plevels.out, npoints)/100)
      if (LID) result <- cbind("id"=rep(id,each=length(plevels.out)),result)
    } else {
      pressures <- vapply(X   = surfpres,
                          FUN = FApressures.local,
                          FUN.VALUE = numeric(nlev), faframe=fa)
      result <- data.frame("p"=as.vector(pressures)/100,
                           "model_level"=rep(1:nlev, npoints))
      if (LID) result$id <- rep(id, each=nlev)
    }
  } else if (levtype=="P") { # pressure levels directly from FA file
    # no vertical interpolation: use plevels.out only for subselection
    # TODO: find all the pressure levels P.....XXX in the fa file
    if (is.null(plevels.out)) {
      flist <- unique(substr(grep("^P[[:digit:]]{5}", fa$list$name, value=TRUE), 2, 6))
      plist <- as.numeric(flist)/100
      plist[plist==0] <- 1000
      plevels.out <- sort(plist)
    }
    nlev <- length(plevels.out)
    plist <- sprintf("%05i", plevels.out * 100)
    plist[plevels.out == 1000] <- "00000"
    result <- data.frame("p"=rep(plevels.out, npoints))
  }
  for (pp in 1:npar) {
    sounding <- matrix(NA, nrow=nlev, ncol=npoints) # dim=c(nlev,npoints)
# we will pass this to C as a vector, so make sure that nlev is the first dimension
# so every points is represented by a vector of nlev values
    for (ll in 1:nlev){
      #print(ll)
      if (levtype == "S") field <- sprintf("S%03i%-12.12s", ll, par[pp])
      else if (levtype == "P") field <- sprintf("P%s%-10.10s", plist[ll], par[pp])
      else stop("unknown levtype", levtype)
      fnum <- FAfind(fa, field)
      if (length(fnum) == 1 && !is.na(fnum) ) sounding[ll,] <- get.data(fa, field)
    }
    if (levtype == "S" && !is.null(plevels.out)) {
      npo <- length(plevels.out)
      sounding <- .C("fa_interp2", A=attr(fa, "frame")$levels$A, B=attr(fa, "frame")$levels$B,
                               pref=attr(fa, "frame")$levels$refpressure,
                               psurf=surfpres, nlev=as.integer(nlev),
                               v_in=sounding, npoints=as.integer(npoints),
                               p_out=numeric(npo), n_out=as.integer(npo),
                               v_out=numeric(npo*npoints))$v_out
    }

    result[[ par[pp] ]] <- as.vector(sounding)
  }
  return(result)
}
# TODO: with 3d geofields it is now in fact more logical to
# have a function for reading 3d data and another for interpolating
# that 3d interpol /could/ be in meteogrid, but better not
# usually, vertical interpol will use hybrid co-ordinates
# BUT: by combining, you only need to read the hybrid levels necessary
# for the interpolation, not all -> much faster

FAdec3d <- function(fa, par="TEMPERATURE", levtype="S", plevels.out=NULL){
  if (!inherits(fa,"FAfile")) fa <- FAopen(fa)
  domain <- attr(fa, "domain")
  npoints <- domain$nx *  domain$ny

  surfpres <- FAdec(fa,"SURFPRESSION")
  if (any(surfpres< 10000)) surfpres <- exp(surfpres)

  if (!is.null(plevels.out)) {
    # get only hybrid levels that are needed to interpolate to plevels.out
    pmin <- FApressures.local(fa, min(surfpres))
    pmax <- FApressures.local(fa, max(surfpres))
    minlev <- max(which(pmax < min(plevels.out)*100), 1)
    maxlev <- min(which(pmin > max(plevels.out)*100), attr(fa, "frame")$nlev)
  } else {
    # no interpolation to pressure levels: get *all* hybrid levels
    minlev <- 1
    maxlev <- attr(fa, "frame")$nlev
  }
  nlev <- maxlev - minlev + 1  
#  cat("npoints=",npoints,"\n")
#  cat("minlev=",minlev,"maxlev=",maxlev,"\n")

  fields <- sprintf("S%03i%s", minlev + 1:nlev - 1, par)
  field3d <- FAdec(fa, fields)

  if (is.null(plevels.out) || nlev == 1) {
  # no interpolation to pressure levels, 
  # OR by some incredible luck, the hybrid level is exactly at 1 pressure...
    result <- field3d
    # better just leave original field names as they come from FAdec
    attributes(result)$info$name <- par
    attributes(result)$info$z_type <- "hybrid"
    
    names(dim(result)) <- c("x", "y", "level")
    dimnames(result) <- list("x"=NULL, "y"=NULL, "hybrid"=minlev:maxlev)
  } else {
    npo <- length(plevels.out)
    if (any(plevels.out != sort(plevels.out))) {
      message("Output pressures are not in ascending order! Re-ordering...")
      plevels.out <- sort(plevels.out)
    }
    # we will pass this to C as a vector,
    # make sure that nlev is the first dimension
    # so every points is represented by a vector of nlev values
    field3d <- aperm(field3d, c(3,1,2))
    field3i <- .C("fa_interp2", A=attr(fa, "frame")$levels$A[minlev:(maxlev+1)],
                               B=attr(fa, "frame")$levels$B[minlev:(maxlev+1)],
                               pref=attr(fa, "frame")$levels$refpressure,
                               psurf=surfpres, nlev=as.integer(nlev),
                               v_in=field3d, npoints=as.integer(npoints),
                               p_out=log(plevels.out*100), n_out=as.integer(npo),
                               v_out=numeric(npo*npoints))$v_out
    if (npo==1) {
      info <- list(name=sprintf("P%05i%s",plevels.out*100,par),
                   level=paste(plevels.out,"hPa"), variable=par)
      result <- meteogrid::as.geofield(field3i, domain=domain, 
                                       info=info, time=attr(surfpres,"time"))
    } else {
      # need to shuffle the dimensions for 3d geofield object
      field3i <- array(field3i, dim=c(npo, domain$nx, domain$ny))
      result <- meteogrid::as.geofield(aperm(field3i, c(2,3,1)), domain, 
				       info=list(name=par, z_type="hPa"),
				       time=attr(surfpres,"time"))
      names(dim(result))[3] <- "hPa"
      dimnames(result) <- list("x"=NULL, "y"=NULL, "hPa"=plevels.out)
      # by using "name" dimension, [[.]] will give it as name!
#      names(dim(result))[3] <- "name"
#      dimnames(result)[3] <- sprintf("P%5.5i%s", plevels.out*100, par)

    }
  }

  return(result)
}

