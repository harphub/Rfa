\name{ecto}
\alias{ecto}
\title{Spectral energy of a field.}
\description{This function implements the basic behaviour of the
  \emph{ectoplasm} program. It calculates the spectral
  energy as a function of wave number, for a given field.}
\usage{ecto(data,NMSMAX=maxM,NSMAX=maxN)}
\arguments{
  \item{data}{A geofield object or a data matrix .}
\item{NMSMAX,NSMAX}{Highest wave numbers in X and Y direction. If not
  specified, they are derived from the matrix dimensions of
  \code{data}, which may give a slightly different result.}
}
\value{
  A vector of length \code{max(NMSAX,NSMAX)+1} containing the spectral
  energy per wave number.}
