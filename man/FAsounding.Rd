\name{FAsounding}
\alias{FAsounding}
\title{Extract pseudo-soundings}
\description{
Extracts vertical samples from a FA file for 1 or more points and parameters.
}

\usage{
FAsounding(fa,par="TEMPERATURE", lon=NULL, lat=NULL, index=NULL, 
           id=NULL, levtype="S", plevels.out=NULL, method="bilin")
}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. \code{fa}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{FAopen}). }
\item{par}{A character string (or a vector of strings) denoting the fields to extract. Default is 'TEMPERATURE'.}
\item{lon,lat}{Co-ordinates at which to extract the pseudo-sounding. The fields are interpolated to these locations using 'method'.}
\item{index}{In stead of longitudes and latitudes, one may also give a set of inideces, a matrix with two columns, e.g. cbind(i,j).}
\item{id}{If multiple points are extracted, they must be distinguished by an ID. This may be the WMO code, but if it is not provided, the default is to either simply number from 1 to N, or use the index values.}
i\item{levtype}{Type of vertical levels in the FA file. "S" for hybrid levels (default), "P" for pressure levels. In the latter case, \code{plevels.out} is used for subselecting pressure levels, not for interpolation.}
\item{method}{Interpolation method. May be 'bilin', 'closets' or 'bicubic'.}
\item{plevels.out}{A vector of pressure levels (in hPa). If it is NULL, the model levels are returned. If it is defined, the pseudo-sounding is interpolated (using linear interpolation on log(p)).}
}
\value{
  A \code{data.frame} containing columns "model_level","pressure" and one column for every entry in \code{par}. If \code{pressure.levels} is given, there is no column "model_level". }
\seealso{\code{\link{FApressure}}, \code{\link[meteogrid]{point.interp}},
         \code{\link{FAdec3d}}
}
\examples{
\dontrun{
  bel1 <- FAopen('ICMSHABOF+0001')
  P35 <- FApressure(bel1,35) ### pressure field at hybrid model level 35
}}
\keyword{file}
