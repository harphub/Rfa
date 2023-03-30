\name{FApressure}
\alias{FApressure}
\alias{FApressures.local}
\title{Calculate the pressure at a hybrid model level of a FA field}
\description{
Calculates pressure field at model levels.
}

\usage{
FApressure(fa,lev,SP = FAdec(fa, "SURFPRESSION    "))
FApressures.local(faframe,surfpressure)
}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. \code{fa}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{FAopen}).}
\item{faframe}{A \code{FAframe}, but it may also be any class that has a FAframe attribute.}
\item{lev}{An integer indicating the model level ("S00n...") to calculate pressure field.}
\item{SP}{Surface pressure field. By default it is taken from the FA file. If it is not available there, SP must be given explicitely.}
\item{surfpressure}{A single numerical value.}
}
\value{
  \code{FApressure} a \code{geofield} containg the pressure field at the requested hybrid model level. \code{FApressures.local} returns a vector of pressures for all the hybrid levels, based on the given surface pressure.}
\details{The pressure calculation uses the formula P = A*P_ref + B*P_surf for half-levels. Then an interpolation is done to full levels.
}
\seealso{\code{\link{FAopen}},\code{\link{FAdec}},\code{\link{iview}}
}
\examples{
\dontrun{
  bel1 <- FAopen('ICMSHABOF+0001')
  P35 <- FApressure(bel1,35) ### pressure field at hybrid model level 35
}}
\keyword{file}
