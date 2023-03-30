\name{utilities}
\alias{FAdescribe}
\alias{FAsizes}
\alias{FAfind}
\alias{ParseTar}
\title{Simple utilities}
\description{
Various functions that may make life easier.
}

\usage{
FAdescribe(fname)
FAsizes(nmin, nmax)
FAfind(fa, field)
ParseTar(archname, quiet = TRUE)
}

\arguments{
\item{fname}{The name of a field in a FA file. A string of length 16}
\item{nmin,nmax}{Minimum and maximum of the vector of admissible FA grid sizes.}
\item{fa}{A FAfile}
\item{field}{An integer or the name of a field.}
\item{archive}{a .tar file.}
\item{quiet}{set FALSE for more diagnostic output.}
}
\value{
\itemize{
\item{FAsizes returns a vector of admissible FA grid sizes (powers of 2,3,5).}
\item{FAdescribe returns a 'description' by splitting a FA field name into level information and parameter name.}
\item{FAfind searches for a field name in the list of fields in a FAfile. This uses grep, so giving less than 16 characters may work, but it may return multiple matches.}
\item{ParseTar returns a list with the byte location of all files in a tar archive.}
}}
\keyword{file}
