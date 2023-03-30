\name{FAread_meta}
\alias{FAread_meta}
\alias{FAread_header}
\alias{FAwrite_meta}
\alias{FAwrite_header}
\title{Read/write FA meta data}
\description{
Functions to read/write the 'meta-data' encoded in the file header and 7 or 8 'frame' and 'time' articles.
}
\usage{
FAread_meta(filename, archname=NULL, quiet=TRUE)
FAread_header(fa)
FAwrite_meta(filename, metadata)
FAwrite_header(fa, header=NULL)
}

\arguments{
\item{filename}{Guess what.}
\item{fa}{A \code{FAfile}, a file name or an open file connection.}
\item{archname}{The file may be part of a tar archive. Archives are always considered read-only!}
\item{metadata}{Usually a list containing the 6 'frame' and 1 or 2 'time' articles as numerical vectors, as output by \code{FAread.meta}. For \code{FAframe()} the time components may be missing, and for \code{FAtime} only the 1 or 2 time vectors need to be provided. In \code{FAtime} it may even be simply a vector of length 11 or 22. }
\item{header}{A vector of length 22. If missing, the header is taken from the FAfile object.}
\item{quiet}{Set to FALSE for (a lot of) diagnostic information.}
}
\value{
  \code{FAdomain} returns a \code{geofield} object, i.e. a
  matrix containing numerical data. with following attributes:
  \item{domain}{A list of all grid and projection parameters, as returned
    by \code{FAgrid}.}
  \item{time}{A list of analysis and forecast time parameters, as
    returned by \code{FAtime}.}
  \item{info}{A description of the field, as returned by \code{FAdescribe}.}
}
\details{
\code{FAwrite_meta} rewrites the metadata fields of an FA file (the list may be partial: only the provided new fields will be overwritten). This is mainly useful for modifying date/time information. The new data must have exactly the same dimensions as the original.
}
\seealso{\code{\link{FAopen}}
}
\keyword{file}
