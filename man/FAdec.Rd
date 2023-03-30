\name{FAdec}
\alias{FAdec}
\alias{FAmsg.dec}
\alias{FAdec.character}
\alias{FAdec.FAfile}
\title{Decode an FA field}
\description{
Decodes one field of an FA file. Calls different methods depending on the first argument.
}

\usage{
FAdec(fa,...)
\method{FAdec}{character}(fa,field,clip=TRUE,outform="G",archname=NULL,
                tar.offset=NULL,faframe=NULL,fatime=NULL,quiet=TRUE,...)
\method{FAdec}{FAfile}(fa,field,clip=TRUE,outform="G",quiet=TRUE,...)
FAmsg.dec(inbuf,faframe,clip=TRUE,outform="G",quiet=TRUE)
}

\arguments{
\item{fa}{A \code{FAfile} object as returned by \code{FAopen}. \code{fa}
  may also be a character string, in which case it is interpreted as a
  filename. If multiple fields 
  are decoded from a single file, it is considerably faster to \code{FAopen} it
  only once! For a single field, }
\item{field}{Can be either an integer, in which case it indicates the
  position of the field in the fieldlist (the frame is not counetd!), or else it can be a character
  string. \code{field} may be any part of the complete
  name of the field, as long as it is unique. If various names fit the given string, an error will thrown only if \code{fa} is an \code{FAfile} object. If \code{fa} is a character string (and the file has thus not been parsed completely the first fit will be chosen, even if there are more.}
\item{clip}{If TRUE, the extension zone is not included in the result.}
\item{outform}{If "G", the result is a meteogrid object. "M" returns the same data matrix but without all the meteogrid attributes. This is a bit faster and may be useful when decoding many fields into a larger structure. Other values are
  "S" for returning the spectral components and "R" for the raw data
  encoded in the file. These alternatives are only applicable if the
  field is stored on spectral format, and the result does not include
  geographical information. }
\item{archname,tar.offset}{If \code{fa} is a file name, it may be part of a tar archive.} 
\item{faframe,fatime}{If \code{fa} is a file name, these variables can be added. This speeds up the decoding (since the frame is not read from the file). If \code{faframe} is provided, \code{fatime} will also not be read from the file, even if it is NULL. The rationale is that we need the frame for decoding the object, but the time info is not vital.}
\item{inbuf}{A bitstream containing an FA article in binary form.} 
\item{quiet}{Set to TRUE for (a lot of) debugging information.}
\item{...}{Not used.}
}
\value{
  A \code{geofield} object, i.e. a
  matrix containing numerical data. with following attributes:
  \item{domain}{A list of all grid and projection parameters, as returned
    by \code{FAgrid}.}
  \item{time}{A list of analysis and forecast time parameters, as
    returned by \code{FAtime}.}
  \item{info}{A description of the field, as returned by \code{FAdescribe}.}
}
\details{\code{FAdec} is a function with two different methods for FAfile and character objects. Some of the arguments are only accepted if \code{fa} is a character string. \code{FAmsg.dec} is an internal function, which does the actual decoding of a raw bitstream.
}
\seealso{\code{\link{FAopen}},\code{\link{FAenc}},\code{\link[meteogrid]{iview}}
}
\examples{
\dontrun{
  bel1 <- FAopen('ICMSHABOF+0001')
  oro <- FAdec(bel1,'SPECSURFGEOP')
}}
\keyword{file}
