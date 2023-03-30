\name{FAframe}
\alias{FAframe}
\alias{FAtime}
\alias{FAdomain}
\title{FA meta data ('frame', 'time')}
\description{
Functions to decode the 'meta-data' encoded in the first 7 or 8 articles of an FA file.
}
\usage{
FAframe(metadata)
FAtime(metadata)
FAdomain(faframe, quiet=TRUE)
}

\arguments{
\item{metadata}{Usually a list containing the 6 'frame' and 1 or 2 'time' articles as numerical vectors, as output by \code{FAread.meta}. For \code{FAframe()} the time components may be missing, and for \code{FAtime} only the 1 or 2 time vectors need to be provided. In \code{FAtime} it may even be simply a vector of length 11 or 22. }
\item{faframe}{A \code{FAframe} class object as returned by \code{FAframe}.}
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
\seealso{\code{\link{FAopen}},\code{\link{FAread.meta}}
}
\keyword{file}
