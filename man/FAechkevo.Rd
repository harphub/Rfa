\name{FAechkevo}
\alias{FAechk.open}
\alias{FAechk.read}
\alias{FAechkevo}
\title{Read FA echkevo files (experimental)}
\description{
Echkevo files are also partially FA files (frame), but the other fields are different.
}

\usage{
FAechk.open(filename, lswap=TRUE)
FAechk.read(fe, tstep)
FAechkevo(filename, lswap=TRUE)
}

\arguments{
\item{filename}{A character string with complete path.}
\item{lswap}{Set to FALSE if the echkevo file was PRODUCED by a big-endian HPC. This covers a bug in echkevo.}
\item{fe}{Output from \code{FAechk.open}.}
\item{tstep}{An integer. Returns the output (all points and fields) for a single timestep.}
}

\value{
\code{FAechkevo} returns a list with various elements:
\item{fa}{The output from \code{FAopen}. This contains all frame information etc.}
\item{info}{Again a list, with information about which fields and grid points are in the echkevo output.}
\item{docu}{The numerical values in the "docu" field. Better use the interpreted version in \code{info}.}
\item{tstep}{Timestep in seconds.}
\item{geo}{}
\item{data}{The actual data. This is a 3 dimensional array with dimensions(timesteps, fields, points).}
}
\seealso{\code{\link{FAdec}}
}
\keyword{file}
