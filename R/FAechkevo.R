### basic support for echkevo files
### these are like FA files: the frame is identical
### but the rest is completely different, and even a bit buggy
### if the ECHK file *comes from* a little-endian machine, some
### parameters are swapped. YOU CAN'T CHECK THIS, because it doesn't
### depend on where you are running Rfa, only on the origin of the file.

FAechk.open <- function(filename,lswap=TRUE){
## set lswap=FALSE if the echkevo file was PRODUCED by a big-endian HPC
  fa <- FAopen(filename,lparse=FALSE)
  ff <- file(attr(fa,"filename"),open="rb")
  on.exit(try(close(ff),silent=TRUE))
### first meta field:
  i1 <- FAfind(fa,"REDOCU0000000000")
  p1 <- fa$list$offset[i1]
  l1 <- fa$list$length[i1]
### l1 = 8*(nval), but the data is originally 32bit (F**K echkevo)
### they are stored as byteswapped 64 bit
### so we need to do some weird stuff to get the numbers right!
### we get val2 val1 val4 val3 ...
  nval1 <- l1/8 # THIS should be the number of fields
  seek(ff,p1)
### make sure you read an EVEN number of numbers!
  if(lswap){
    nvh <- ceiling(nval1/2)
    val1 <- readBin(ff,"integer",size=4,n=nvh*2,endian="big")
    reord <- rep(2*(1:nvh),each=2) + c(0,-1) #=2,1,4,3,6,5,...
    val1 <- val1[reord][1:nval1]
  } else {
    val1 <- readBin(ff,"integer",size=4,n=nval1,endian="big")
  }

  nfrq <- val1[1]      ## frequency (in timesteps) of export, often ==1
  nfld <- val1[2]      ## nr of fields (par * lev) exported
  npoints <- val1[3]   ## nr of points exported
  nlev <- val1[4]      ## number of levels in the model run (NFLEV)
  nthr <- val1[5]      ## number of thermodynamic variables, often ==1, but 3 if NH dynamics
  npas <- val1[6]      ## other ("passive") fields, currently always 0
  lnsp <- as.logical(val1[7]) ## TRUE means pressure is still ln(pres) ???

  nstep <- attr(fa,"nfields") -3 # the 7 frame fields have already been extracted! 
  fld <- val1[8:(7+nfld)]
## TODO: check how general this is:
## IF NH-DYN: 2 more fields
  if (max(fld) > (4 + nthr)*nlev + 1) stop("Inconsistency in NTHR vs fields.")

  if (nthr==0) thrnames <- NULL
  else if (nthr==1) thrnames <- "temp"
  else thrnames <- c("temp",paste("THR",2:nthr,sep=""))

  allfields <- c(sprintf("%s%03i",rep(c("vor","div","u","v",thrnames),each=nlev),1:nlev),
                   "pressure")
  names(fld) <- allfields[fld]
  index <- as.data.frame(matrix(val1[(nfld+8):(nfld+7+2*npoints)],ncol=2,byrow=TRUE))
  names(index) <- c("i","j")
  info <- list(npoints=npoints,nfld=nfld,nlev=nlev,nstep=nstep,fld=fld,nfrq=nfrq,
               nthr=nthr,npas=npas,lnsp=lnsp)

###################
  i2 <- FAfind(fa,"RESTEP0000000000")
  p2 <- fa$list$offset[i2]
  l2 <- fa$list$length[i2]
  seek(ff,p2)
  info$tstep <- readBin(ff,"double",n=1,size=8,endian="big")

###################
  i3 <- FAfind(fa,"REGEOM0000000000")
  p3 <- fa$list$offset[i3]
  l3 <- fa$list$length[i3]
  seek(ff,p3)
  mgeo <- readBin(ff,"double",n=5*npoints+1,size=8,endian="big")

  geo <- data.frame( i=index$i,j=index$j,
                     lon=mgeo[1:npoints]*180/pi,
                     lat=mgeo[(npoints+1):(2*npoints)]*180/pi,
                     mapfactor=mgeo[(2*npoints+1):(3*npoints)],
                     CorX=mgeo[(3*npoints+1):(4*npoints)],
                     CorY=mgeo[(4*npoints+1):(5*npoints)])

  result <- list(fa=fa,info=info,geo=geo)
  result
}

FAechk.read <- function(fe,tstep) {
### read 1 *timestep*
  ff <- file(attr(fe$fa,"filename"),open="rb")
  on.exit(try(close(ff),silent=TRUE))
### first meta field:
  fname <- sprintf("RE%4i0000000000",tstep)
  i1 <- FAfind(fe$fa,fname)
  p1 <- fe$fa$list$offset[i1]
  l1 <- fe$fa$list$length[i1]
  nval1 <- l1/8
  seek(ff,p1)
  data <- readBin(ff,"double",size=8,endian="big",n=nval1)
### this data is organised
  nfld <- fe$info$nfld
  npoints <- fe$info$npoints
#  cat("nval1=",nval1,"nfld=",nfld,"npoints=",npoints,"\n")
  result <- matrix(data[(2*nfld+1):nval1],nrow=nfld,ncol=npoints,byrow=TRUE)
  result
}

FAechkevo <- function(filename,lswap=TRUE){
  fe <- FAechk.open(filename,lswap=lswap)
  data <- array(NA,dim=c(fe$info$nstep,fe$info$nfld,fe$info$npoints))

  for(i in 1:fe$info$nstep) {
    data[i,,] <- FAechk.read(fe,i-1)
  }
  dimnames(data) <- list(
                         time=fe$info$tstep*fe$info$nfrq*(0:(fe$info$nstep-1)),
                         field=names(fe$info$fld),
                         point=1:fe$info$npoints)
  c(fe,data=list(data))
} 
