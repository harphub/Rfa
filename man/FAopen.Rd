\name{FAopen}
\alias{FAopen}
\alias{FAopenTar}
\alias{FAopen}
\title{Open an FA file}
\description{
Opens an FA file and returns a \code{FAlist} object. The data itself is
not decoded, but domain specifications, date, list of fields etc are returned.
}

\usage{
FAopen(filename,archname=NULL,tar.offset=NULL,lparse=TRUE,quiet=TRUE)
FAopenTar(archname,lparse=FALSE,quiet=TRUE)
}

\arguments{
\item{filename}{A character string with complete path.}
\item{archname}{If the FA file is part of a tar archive. A character string with complete path.}
\item{tar.offset}{Byte location of the FA file within the tar archive. If NULL, the file is searched for.}
\item{lparse}{If FALSE, the FA file(s) are only searched for all available fields and their bit adresses. If TRUE, the FA file(s) are parsed completely and the main characteristics off every field are also returned (spectral, grib encoding etc.). This is quit fast and can be useful (it is not vital), but may be slow for large tar archives over a slow conection.}
\item{quiet}{If FALSE, several diagnostic messages are printed.}
}

\value{
An \code{FAfile} object is a list with two components:

\item{list}{A list of all available fields (excluding the first 7). The columns contain byte address etc. The number of columns depends on the value of \code{lparse}}
\item{holes}{a list of "holes" in the file resulting from rewrites and removals of fields. Only important for encoding routines.}

A \code{FAfile} object also has a list of attributes:

\item{filename}{Full name of the file.}
\item{tarfile,tar.offset}{NULL and 0 for individual files.}
\item{nfields,nholes}{Number of fields (holes) contained in the file.}
\item{frame}{All geographical meta data: domain size, projection, vertical levels etc.}
\item{time}{Time information (analysis and forecast).}

\code{FAopenTar} returns a list of such objects.
}
\seealso{\code{\link{FAdec}}
}
\keyword{file}
