\name{FAenc}
\alias{FAenc}
\alias{FAenc_msg}
\alias{FAadd}
\alias{FAremove}
\alias{FArename}
\alias{FArename_meta}
\title{Encode/remove/rename a field in an existing FA file.}
\description{
Encode, remove or rename one field in an existing FA file opened by
\code{FAopen}. These functions are quite new and should be used with caution!
}

\usage{
FAenc(fa, fieldname, data, lspec=NULL, lgrib=NULL,
      nbits=NULL, sptrunc=NULL, sppow=NULL,
      overwrite=TRUE, quiet=TRUE)
FAenc_msg(data, faframe, lspec=NULL, lgrib=NULL, nbits=NULL,
          sptrunc=NULL, sppow=NULL, quiet=TRUE)
FAadd(fa, fieldname, data, quiet=TRUE, ...)
FAremove(fa, field, quiet=TRUE)
FArename(fa, field, newname, quiet=TRUE)
FArename_meta(fa, field, newname, quiet=TRUE)

}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. If
  fa is a character string, it is interpreted as a file name and
  opened using FAopen.}
\item{data}{The data field to be encoded. If the field dimension are
  smaller than required by the dataframe, the domain is extended.}
\item{fieldname, newname}{A name for the new field. It must be at least 5 characters
  long and may not yet exist in the FAfile. If the name already exists,
  FAenc will replace it by the new values (unless \code{overwrite}==FALSE). }
\item{field}{Name or number of the field to be removed/renamed.}
\item{lspec}{FALSE stores the data in gridpoint format, TRUE in spectral
  format.  }
\item{lgrib}{Grib encoding.}
\item{nbits}{Number of bits per value in the GRIB encoding.}
\item{sptrunc}{Spectral truncation.}
\item{sppow}{Power of the laplacian used for rescaling.}
\item{overwrite}{ If FALSE, the function exits if the field name
  already exists in the given FA file. If \code{overwrite} is
  \code{TRUE} an existing field with the same name will be
  overwritten.}
\item{faframe}{An FAframe.}
\item{quiet}{FALSE turns on a lot of diagnostic messages.}
\item{...}{Options to FAens.msg}
}
\value{
  All these functions return a modified \code{FAfile} object resulting from opening the modified file. See the example for typical way of use.
}
\details{
  \code{FAenc} is the main function for encoding data into FA format, and probably the only one you should use. The other functions are rather internal functions that may be useful in special cases.

To write a data matrix, it must have the correct dimensions. This
  often means that a biperiodicisation zone must be added. This can be
  done using the function \code{biper}. If the data matrix doesn't have a
  extension zone yet, \code{FAenc} will call \code{biper} automatically.

If the file is not large enough to contain the new data, it is extended. For this reason, these routines will not accept files that are part of a tar archive.

Be very careful with \code{FArename.meta}: only the frame name (6) or 'DATX-DES-DONNEES' (8) may be renamed without corrupting the file completely.

}
\seealso{\code{\link{FAopen}},\code{\link{FAdec}},\code{\link{biper}}
}
\examples{
\dontrun{
  fa1 <- FAopen("someFAfile")
  fa1 <- FAremove(fa1,"someField")
}}
\keyword{file}
