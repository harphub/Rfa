\name{FAmake.frame}
\alias{FAmake.frame}
\title{Create a new faframe class object}
\description{
Encode or remove one field into an existing FA file opened by
\code{FAopen}. These functions are quite new and should be used with caution!
}

\usage{
FAmake.frame(domain,extension=c(11,11),relaxation=c(8,8),nsmax=NULL,nmsmax=NULL,
             lineargrid=TRUE,levels=NULL,
             name="FA-FRAME        ", sptrunc=10)
}
\arguments{
\item{domain}{A \code{geodomain} or having a domain as attribute.}
\item{extension}{Dimension of the extension zone in X and Y directions.}
\item{relaxation}{Size of teh relaxation zone for Davis coupling. NOT USED?}
\item{nmsmax, nsmax}{Spectral truncation. Default value depends on the choice of linergrid.}
\item{lineargrid}{If TRUE, the default values of nmsmax and nsmax are given by floor((nx-1)/2). If FALSE, quadtratic grid is taken: floor((nx-1)/3). If nsmax and nmsmax are given explicitely, this option is ignored.}
\item{levels}{Hybrid level definition. Default is 1 level.}
\item{name}{Name for the frame as written in the FA file.}
\item{sptrunc}{The default wavelength for truncating the GRIB compacted part. Best to leave this at 10.}
}
\value{Returns a FAframe class object. Only 'aladin' type.
}
\seealso{\code{\link{FAcreate}}, \code{\link{FAmake.header}}, \code{\link{FAmake.time}}
}
\keyword{file}
