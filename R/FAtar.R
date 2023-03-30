### for tar archives containing multiple FA files
### especially different lead times from a single run
### if the tar file contains subdirectories, these also show up as separate entries in the listing
### So to cope with this, FAopen must be "forgiving"

#' @export
FindInTar <- function(archname, filename, quiet=TRUE){
  # open archive
  if (!file.exists(archname)) stop("File ", archname, " not found")
  tf <- file(archname, open='rb')
  on.exit(close(tf))
  # offset of first entry is zero
  offset <- 0
  while (TRUE) {
    # goto beginning of entry
    seek(tf, offset)
    # read file name
    fName <- readBin(tf,what="char",n=1,size=100)
    if (nchar(fName)==0) stop(sprintf('file %s not found in archive\n',filename))
    # read size
    seek(tf, offset+124, origin='start')
# file size is encoded as a length 12 octal string, 
# with the last character being '\0' (so 11 actual characters)
    sz <- readChar(tf, nchars=11) 
    # convert string to number of bytes
    sz <- sum(as.numeric(strsplit(sz,'')[[1]])*8^(10:0))
    if (!quiet) cat(sprintf('entry %s, %i bytes\n', fName, sz))
    if (fName==filename) {
      break
    } else {
      # goto the next message
      offset <- offset+512*(ceiling(sz/512)+1)
    }
  }
  # skip entry header
  offset <- offset+512
  return(offset)
}

# in a gzipped file, using seek() to skip to a byte location is unsafe
# so we do it by reading raw bytes
# just to avoid having to read multiple GB of data, we do it in steps
# default step is 10MB
zip_skip <- function(zfile, skip, by=10^7) {
  while (skip > by) {
    readBin(zfile, "raw", by)
    skip <- skip - by
  }
  readBin(zfile, "raw", skip)
  zfile
}

#' Parse a (gzipped) tar file and make a list of the files
#' @param archname The name of a (gzipped) tar file
#' @param gzip Set to TRUE if the archive is a gzipped file
#' @return A named list of the files contained in the archive. 
#'     Every list element has 3 attributes: filename, (byte) location and file size.
#' @export
ParseTar <- function(archname, gzip=FALSE) {
  blocksize <- 512
  if (!file.exists(archname)) stop("File ", archname, " not found")
  # open archive
  if (gzip) tf <- gzfile(archname, open='rb')
  else tf <- file(archname, open='rb')
  on.exit(close(tf))
  # offset of first entry is zero
  fnames <- list()

  offset <- 0
  nfile <- 0
  while (TRUE) {
    # goto beginning of entry
    # in gzipped archives and on windows: avoid using seek()!!!
    if (!gzip) seek(tf, offset)
    ### but for up to ~10MB : no problem, I guess
    else if (offset > seek(tf)) {
      ## large file: split readBin into blocks of e.g. 1E7 bytes (~10MB)?
      # readBin(tf, what = "raw", n = offset - seek(tf))
      tf <- zip_skip(tf, offset - seek(tf))
    }
    # read file name
    # readBin(..., what="char") can give errors in gzipped file
    header <- readBin(tf, what="raw", n=blocksize)
    if (length(header) < blocksize) break
    # a tar archive usually ends with two 0-filled records
    if (all(header == 0)) break
    # check UStar format
    # it /should/ be a \0 terminated string "ustar"
    # but older GNU-tar files have "ustar  "
    # so don't be too strict...
    check <- rawToChar(header[258:262])
    if (check != "ustar") stop("Probably not a correct tar archive.")
    # size
    # file size is encoded as a length 12 octal string, 
    # with the last character being '\0' (so 11 actual characters)
    sz <- rawToChar(header[125:136])
    # convert string to number of bytes
    fsize <- sum(as.numeric(strsplit(sz,'')[[1]])*8^(10:0))
    # is it a real file? (not a link, directory, PAX header...)
    typeflag <- rawToChar(header[157])
    if (typeflag %in% c("", "0")) {
      nfile <- nfile + 1
      # name
      fName <- rawToChar(header[1:100])
      if (nchar(fName)==0) break
      # prefix?
      prefix <- rawToChar(header[346:500])
      if (!prefix=="") fName <- paste(prefix, fName, sep="/")

      fnames <- c(fnames, fName)
      attr(fnames[[nfile]], "tar.offset") <- offset + blocksize 
      attr(fnames[[nfile]], "tarfile") <- archname
      attr(fnames[[nfile]], "size") <- fsize
      # cat(sprintf('entry %s, %i bytes (type %s)\n', fName, fsize, typeflag))
    } else {
      ### "x": PAX header, "g": global PAX header, "5"=(sub)directory
      # cat(sprintf('entry %s, %i bytes, TYPE: %s\n', rawToChar(header[1:100]), fsize, typeflag))
    }
    # goto the next (header) message
    # skip actual file
    offset <- offset + blocksize*(ceiling(fsize/blocksize) + 1)
  }
# return a named list of characters strings with attributes?
  names(fnames) <- fnames
  return(fnames)
}

FAopenTar <- function(archname, lparse=FALSE, quiet=TRUE, gzip=FALSE){
### parse a tar file and create a list of FAfile objects
### by default don't parse the complete files
###   (so you don't see which fields are spectral etc)
### only the field names and data locations.
### this is faster and you probably don't really need the extra info anyway.
### it may be interesting to index simply by lead time of the file
### because often tar files are archives of a single run

  if (!file.exists(archname)) stop("File ", archname, " not found")
  archlist <- ParseTar(archname, gzip=gzip)
  if (!quiet) print(archlist)
# archlist now has tarfile and tar.offset as attributes for every entry
# no need to give them explicitely (but you can if you want...)

# the error-catching is very ugly. There must be a cleaner way.
  result <- lapply(archlist, function(ff) try(FAopen(ff, lparse=lparse, quiet=quiet), silent=TRUE))  
  names(result) <- archlist

# any errors? if yes, take only the entries that were OK.
  NOK <- vapply(result, function(x) inherits(x,"try-error"), FUN.VALUE=TRUE)

  if (any(NOK)) result <- result[!NOK]
  if (all(NOK)) stop("Archive does not contain any valid FA files.") 
  result
}  

