\name{internal}
\alias{FArawreorder}
\alias{FAraw2fft}
\alias{FAfft2raw}
\alias{FAheader.date}
\title{Internal routines}
\description{
You will probably never want to use these functions. But who knows...}
\usage{
FArawreorder(rawdata, nmsmax, nsmax, inv = FALSE)
FAraw2fft(rawdata,nmsmax,nsmax,ndlon,ndgl)
FAfft2raw(fftdata,nmsmax,nsmax,ndlon,ndgl,quiet=TRUE)
FAheader.date(header,quiet=TRUE)
}

\arguments{
\item{rawdata}{A vector of 'raw' (quadruplet) spectral components.}
\item{nmsmax, nsmax}{Spectral truncations.}
\item{fftdata}{A matrix of complex spectral components.}
\item{ndlon,ndgl}{Grid dimensions.}
\item{inv}{For the inverse transformation, of course.}
\item{header}{A vector of length 22, the header of a FA/LFI file.}
\item{quiet}{Set to FALSE for (a lot of) debugging information.}
}
\value{
\itemize{
   \item{FArawreorder changes the order of the quadruplets. This is of importance when reading some GRIB files for data assimilation B-matrix calculations, because there the quadruplets are ordered differently!}
  \item{FAraw2fft and FAfft2raw convert from a vector of quadruplets to a complex matrix of FFT components and back.}
  \item{FAheader.date updates the 'last change' date in the header.}
}}
\keyword{file}
