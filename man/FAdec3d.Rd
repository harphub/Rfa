\name{FAdec3d}
\alias{FAdec3d}
\title{Extract 3d data}
\description{
Extracts a 3d parameter from an FA file.
}

\usage{
FAdec3d(fa, par="TEMPERATURE", plevels.out=NULL, gf=TRUE)
}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. \code{fa}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{FAopen}). }
\item{par}{A character string denoting the field to extract. Default is 'TEMPERATURE'.}
\item{plevels.out}{A vector of pressure levels (in hPa). If it is NULL, the model levels are returned. If it is defined, the pseudo-sounding is interpolated (using linear interpolation on log(p)).}
\item{gf}{If TRUE, the result is a list of geofields. If FALSE, the result is a plain 3d array. If only a single level is returned, it is always a geofield.}
}
\value{
  Either a geofield (if length(plevels.out)==1), a list of geofields (gf=TRUE) or an array with dimensions \code{c(length(plevels.out), nx, ny)}. nx and ny are the dimensions of the domain. If \code{plevels.out=NULL}, the first dimension is the number of hybrid levels.

  Note that if the lowest model level has a pressure that is lower than a requested pressure level (e.g. on top of a mountain there may be no 1000hPa) the local value is set to NA. There is no extrapolation.
}
\seealso{\code{\link{FApressure}}, \code{\link[meteogrid]{point.interp}},
         \code{\link{FAsounding}}
}
\examples{
\dontrun{
  bel1 <- FAopen('ICMSHABOF+0001')
  P35 <- FAdec3d(bel1,"TEMPERATURE", plevels.out=seq(100,1000,by=100))
 ### temperature inetrpolated to pressures 100, 200, ... , 1000hPa
}}
\keyword{file}
