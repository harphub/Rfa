\name{LFIreplace}
\alias{LFIreplace}
\title{Encode an LFI field (SURFEX). Only replacement of an existing field is possible.}
\description{
Encodes a matrix to the given field in a LFI file.
}

\usage{
LFIreplace(lfi,field,data,missing=1.0E+20,quiet=TRUE)
}

\arguments{
\item{lfi}{An object as returned by \code{LFIopen}. \code{lfi}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{LFIopen}). }
\item{field}{Can be either an integer, in which case it indicates the
  position of the field in the fieldlist, or else it can be a character
  string. If \code{field} is a string, it may any part of the complete
  name of the field, as long as it is unique.}
\item{data}{A matrix with the new values. It must have the right dimensions of the original field. LFI fields have a "skin" of missing values that surround the data. Thus, \code{data} must have either dimensions (nx,ny), in which case NA values are added at the border, or (nx+2,ny+2). }
\item{missing}{Data value that should be interpreted as missing value. NA values in \code{data} are replaced by this value before encoding. }
\item{quiet}{FALSE turns on the diagnostic messages.}
}
\value{
The return value is a \code{LFIfile} class object with the new contents.
}
\details{
The field values are replaced in the LFI file by the given data.
}
\seealso{\code{\link{LFIopen}},\code{\link{LFIdec}}
}
\keyword{file}
