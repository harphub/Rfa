\name{FAmake.time}
\alias{FAmake.time}
\title{Create a list of time entries for a FA file.}
\description{
Encode or remove one field into an existing FA file opened by
\code{FAopen}. These functions are quite new and should be used with caution!
}

\usage{
FAmake.time(fcdate="1970010100", ldt=0, unit='h', init=1, previous=0, first=0, extended=FALSE, tstep=1)
}
\arguments{
\item{fcdate}{Forecast data, a POSIXct class object. A string or simple number YYYYMMDDHH is also possible.}
\item{ldt}{Forecast range in seconds.}
\item{unit}{Unit used for \code{range}. Usually 'h'. But ay be 'm' for minutes.}
\item{init}{Should the time vector be labelled as 'initialised' or not?}
\item{previous, first}{Range (in seconds!) of 'previous' file and 'start file' for accumulations etc.}
\item{extended}{Should we add a second set of 11 integers for more time meta data?}
\item{tstep}{Only used when \code{extended=TRUE}: time step in seconds.}
}
\value{Returns a vector of length 11 for FA file time sector.
}
\seealso{\code{\link{FAcreate}}, \code{\link{FAmake.header}}, \code{\link{FAmake.frame}}
}
\keyword{file}
