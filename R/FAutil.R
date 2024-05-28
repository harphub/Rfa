### Re-order ALADIN model order of raw spectral data to FA file order
### useful for some rare cases (Data assimilation)
FArawreorder <- function(rawdata,nmsmax,nsmax,inv=FALSE){
### NEVER checked for nsmax != nmsmax !!!
### inv doesn't make a difference if nmsmax=nsmax
#  maxm <- floor(sqrt(1 - (0:nsmax)^2/nsmax^2) * nmsmax + 1e-10) + 1
  warning("THIS FUNCTION IS NOT VALIDATED FOR NMSMAX != NSMAX !!!")
  ll <- .C("fa_countval_spec",as.integer(nmsmax),as.integer(nsmax),
           sptrunc=as.integer(-1),result=integer(1))$result
  if(length(rawdata) != ll)
    stop(paste("Incorrect data. Length should be ",ll))

  .C("fa_rawreorder",fieldin=as.double(rawdata),fieldout=numeric(ll),
     as.integer(nmsmax),as.integer(nsmax))$fieldout
}


############ MAYBE USEFUL?

ecto <- function(data,NMSMAX=maxM,NSMAX=maxN){
### spectral energy (like the FAmous ecto routine)
### not a very efficient implementation, mainly to check my FFT representation
  maxM <- dim(data)[1] %/% 2
  minM <- (-dim(data)[1] ) %/%2+1
  maxN <- dim(data)[2] %/% 2
  minN <- (-dim(data)[2] ) %/%2+1

### construct matrices that give the m and n position in the fft matrix
  M <-  matrix(c(0:maxM,minM:(-1)),ncol=ncol(data),nrow=nrow(data),byrow=FALSE)
  N <-  matrix(c(0:maxN,minN:(-1)),ncol=ncol(data),nrow=nrow(data),byrow=TRUE)

  ect <- rep(0,max(NMSMAX,NSMAX)+1)

  K <- sqrt(M^2+N^2* NMSMAX^2/NSMAX^2)
  frac <- K - floor(K + 1e-10)
### there's a chance of rounding error if k is int - 10e-12 i.s.o int
### frac would then be almost 1 i.s.o. =0
### better take floor(K+1e-10)

  energy <- abs(data)^2

  ect[1] <- energy[1,1]
### and now the *really* slow part:
  for(k in 2:length(ect)){
    ect[k] <- sum( (energy*(1-frac))[(K>=k-1) & (K < k)] ) +
                     sum( (energy*frac)[(K > k-2) & (K < k-1)] )
  }
  ect
}

ectoplot <- function(data,add=FALSE,...){
  ect <- ecto(data)
  maxMN <- length(ect)-1

  if (add) points(1:maxMN,ect[-1],type="o",...)
  else plot(1:maxMN,ect[-1],log="xy",type="o",...)
}
################################################################

### bi-periodise a field using cubic splines
biper <- function(data, ext = c(11,11), bc=0){
  if (length(ext) == 1) ext <- rep(ext, 2)
  if (any(dim(data) < 3 )) stop("grid dimensions must be >= 3")
  realdim <- dim(data) 
  newdim <- dim(data) + ext
  data2 <- matrix(0, nrow=newdim[1],ncol=newdim[2])
  data2[1:realdim[1], 1:realdim[2]] <- data
  # call C-spline code
  res <- .C("biper", data=as.numeric(as.vector(data2)),
              as.integer(realdim[1]), as.integer(realdim[2]),
              as.integer(newdim[1]), as.integer(newdim[2]), as.integer(bc))$data
  # smoothing
  res <- .C("smooth_extension", data=res,
              as.integer(realdim[1]), as.integer(realdim[2]),
              as.integer(newdim[1]), as.integer(newdim[2]), as.integer(bc))$data
  matrix(res, nrow=newdim[1], ncol=newdim[2])

}

add_ext <- function(data, ext=c(11,11), bc=0) {
  realdim <- dim(data) 
  newdim <- dim(data) + ext
  data2 <- matrix(0, nrow=newdim[1],ncol=newdim[2])
  data2[1:realdim[1], 1:realdim[2]] <- data

  matrix(res,nrow=newdim[1],ncol=newdim[2])
}

smooth_ext <- function(data, ext=c(11,11)){
#  stop("biperiodicisation not available.")
  newdim <- dim(data) 
  realdim <- dim(data) - ext
  res <- .C("smooth_extension", data=as.numeric(as.vector(data)),
              as.integer(realdim[1]), as.integer(realdim[2]),
              as.integer(newdim[1]), as.integer(newdim[2]))$data
  matrix(res,nrow=newdim[1],ncol=newdim[2])
}

#######################################
### all possible sizes for FA domains (powers of 2,3,5)
### there must be at least 1 power of 2
FAsizes <- function(nmin,nmax){
  max2 <- floor(log(nmax)/log(2))
  max3 <- floor(log(nmax)/log(3))
  max5 <- floor(log(nmax)/log(5))

  allpow <- array(NA,dim=c(max2+1,max3+1,max5+1))
  for(i in 1:max2) for(j in 0:max3) for(k in 0:max5)
    allpow[(i+1),(j+1),(k+1)] <- 2^i*3^j*5^k
  allpow[allpow<nmin | allpow>nmax] <- NA
  sort(as.vector(allpow))
}
#######################################
