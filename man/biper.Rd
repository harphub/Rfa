\name{biper}
\alias{biper}
\title{Biperiodic continuation of a field.}
\description{
This is essentially a wrapper to the ETIBIHI routine for biperiodicisation
of a 2D field using isotropic splines.
}
\usage{biper(data,newdim=dim(data)+11)}
\arguments{
  \item{data}{A geofield object or a data matrix .}
\item{newdim}{The dimensions of the final field. By default, the
  extension zone is 11 points wide.}
}
