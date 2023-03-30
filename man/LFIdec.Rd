\name{LFIopen}
\alias{LFIopen}
\alias{LFIdec}
\alias{LFIread}
\alias{LFItime}
\alias{LFIdomain}
\title{Decode an LFI field (SURFEX)}
\description{
Read and/or decodes fields of a LFI file.
}
\usage{
LFIopen(filename,quiet=TRUE)
LFIdec(lfi,field,...)
LFIread(lfi,field,type="numeric",missing=1.0E+20,add.attr=FALSE,quiet=TRUE)
LFItime(lfi)
LFIdomain(lfi,quiet=TRUE)
}

\arguments{
\item{filename}{Guess what...}
\item{lfi}{An object as returned by \code{LFIopen}. \code{lfi}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{LFIopen}). }
\item{field}{Can be either an integer, in which case it indicates the
  position of the field in the fieldlist, or else it can be a character
  string. If \code{field} is a string, it may any part of the complete
  name of the field, as long as it is unique.}
\item{type}{Expected data type.}
\item{missing}{Data value that should be interpreted as missing value. This is replaced by NA in the output. }
\item{add.attr}{If TRUE, extra attributes (\code{gridtype}, \code{mcomment}) are added to the decoded field.}
\item{quiet}{FALSE turns on the diagnostic messages.}
\item{...}{options for LFIread.}
}
\value{
An \code{LFIfile} object is a list with two components:
\itemize{
\item{list}{A list of all available fields: name and (byte) loxcation in the file. The columns contain byte address etc.}
\item{holes}{Not implemented for LFI. NULL. }
}
A \code{LFIfile} object also has a list of attributes:
\itemize{
\item{header}{The LFI header: a vector of 22 integers.}
\item{filename}{Full name of the file.}
\item{nfields}{Number of fields contained in the file.}
\item{domain}{All geographical meta data: domain size, projection, vertical levels etc.}
\item{validdate}{Time information. LFI files contain no information about leadtime!}
}

\code{LFIdec} returns a \code{geofield}, but note that the date information is minimal, because LFI files don't even store forecast lead time separately.


}
\details{
As a user, you will probably only need to use \code{LFIopen} and \code{LFIdec}. The other functions are called implicitely.
LFI files don't store all date information. Only the valid date is available (forecast date + range). So there is no way to know the forecast date and the lead time separately, except e.g. from the file name, but that is up to the user.
}
\seealso{\code{\link{LFIreplace}},\code{\link{FAopen}},\code{\link{FAdec}}
}
\keyword{file}
