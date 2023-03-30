\name{FAmake.header}
\alias{FAmake.header}
\title{Create a basic header (22 reals) for a FA file}
\description{
Reurns a basic header for a FAfile. The only option is the sector size.}

\usage{
FAmake.header(sector.size=1000)
}
\arguments{
\item{sector.size}{The size of the data sectors for a FA file in 16bit words.}
}
\value{A header is a vector of 22 reals.}
\seealso{\code{\link{FAcreate}}, \code{\link{FAmake.frame}}, \code{\link{FAmake.time}}
}
\keyword{file}
