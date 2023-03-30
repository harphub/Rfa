\name{FAcreate}
\alias{FAcreate}
\title{Create a new FA file}
\description{
Creates a new FA file based on domain and time information.
}
\usage{
FAcreate(filename, frame, time=FAmake.time(), sector.size=1000, 
         nsector=10, overwrite=FALSE)
}

\arguments{
\item{filename}{Guess what? A string with the file name.}
\item{frame}{An \code{FAframe} class object, e.g. created by \code{FAmake.frame}.}
\item{time}{A vector of 11 integers representing date and forecast range. Default is to cal \code{FAmake.time()}, creating a default date.}
\item{sector.size}{Size of the data sectors (in 8 byte words). Since a field name is 16 bytes, ths is double the number of fields that can be encoded without having to do some sector magic.}
\item{nsector}{How many data sectors of the given size do you want? Absolute minimum is 4. Putting more just gives a larger file to begin with. When you write new fields to the file, data sectors are added on the fly when necessary.}
\item{overwrite}{Just a simple safety switch.}
}
\value{
This function creates a new file (with only the 7 meta-information fields) and returns a \code{FAfile} class object for the newly created file.
}
\seealso{\code{\link{FAopen}},\code{\link{FAmake.frame}}
}
\keyword{file}
