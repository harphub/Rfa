\name{FAslice}
\alias{FAslice}
\title{Extract vertical cross section}
\description{
Extracts a vertical cross section from a FA file for 1 parameter, along any line of the grid.
}

\usage{
FAslice(fa,par="TEMPERATURE", axis="X", n=1, type="S", plot=FALSE, 
        plot.function=contour, ...)
}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. \code{fa}
  may also be a character string, in which case it is interpreted as a
  filename (and automatically opened using \code{FAopen}). }
\item{par}{A character string denoting the field to extract. Default is 'TEMPERATURE'.}
\item{type}{The level type to extract. Can be 'S' for hybrid levels or 'P' for pressure levels.}
\item{axis}{Can be 'X' or 'Y' and denotes in which direction the slice is taken.}
\item{n}{The index (along the other axis) where to take the cross section.}
\item{plot}{If TRUE, the vertical cross-section is plotted.}
\item{plot.function, ...}{The function that should draw the plot, and any options that can be passed.}
}
\value{
  A matrix with a cross section of the chosen parameter.}
\seealso{\code{\link{FApressure}}, \code{\link{FAsounding}}
}
\examples{
\dontrun{
  bel1 <- FAopen('ICMSHABOF+0001')
  P35 <- FAslice(bel1,axis="X",15)
  P35 <- FAslice(bel1,axis="X",15, plot=TRUE, plot.function=limage,
                 color.palette=irainbow)
}}
\keyword{file}
