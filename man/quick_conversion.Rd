\name{A4-quick-conversion}
\alias{A4-quick-conversion}
\alias{qDF}
\alias{qDT}
\alias{qM}
\alias{mctl}
\alias{mrtl}
\alias{qF}
\alias{qG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quick Data Conversion}
\description{
Convert common data types quickly, without method dispatch and extensive checks:
\itemize{
\item \code{qDF} and \code{qDT} convert vectors, matrices, higher-dimensional arrays and suitable lists to \code{data.frame} and \code{data.table} respectively.
\item \code{qM} converts vectors, higher-dimensional arrays, data.frames and suitable lists to matrix.
\item \code{mctl} and \code{mrtl} are exported C++ functions that column- or row-wise convert a matrix to list, \code{data.frame} or \code{data.table}. They are used internally by \code{qDF} and \code{qDT}, \code{\link{dapply}}, \code{\link{BY}}, etc...
\item \code{qF} converts atomic vectors to factor.
\item \code{qG} is a programmers function that converts atomic vectors to a quick-group (class attribute 'qG') - a kind of factor light which is just an integer vector with an attribute 'N.groups' indicating the number of groups. (Thus saving time consumed by converting the levels of a numeric vector to character).
}
}
\usage{
qDF(X, row.names.col = FALSE)
qDT(X, row.names.col = FALSE)
qM(X)
mctl(X, names = FALSE, ret = 0L)
mrtl(X, names = FALSE, ret = 0L)
qF(x, ordered = TRUE, na.exclude = TRUE)
qG(x, ordered = TRUE, na.exclude = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix, higher-dimensional array, data.frame or list. \code{mctl} and \code{mrtl} only take matrices.}
  \item{x}{a atomic vector or factor.}
  \item{row.names.col}{should a column capturing names or row.names be added? i.e. when converting atomic objects to data.frame or data.frame to data.table. Can be logical \code{TRUE}, which will add a column \code{"row.names"} in front, or can supply a name for the column i.e. \code{"column1"}.}
  \item{ordered}{logical. \code{TRUE} sorts the levels and returns an ordered factor.}
  \item{names}{logical. Should the list be named?}
  \item{ret}{integer. 3 return options: 0L - return a list, 1L - return a data.frame, 2L - return a data.table.}
  \item{na.exclude}{logical. \code{TRUE} preserves missing values in character or numeric vectors (i.e. no underlying integer is generated to represent \code{NA}, and functions like \code{na_rm} will remove missing values from the factor). \code{FALSE} creates an integer for \code{NA} and also attaches a class 'na.included', which will make computations with \code{\link{collapse}} functions faster (since no \code{NA} checking is required). \emph{Note} that in both cases a \code{'NA'} level is attached, so this option only concerns the underlying integer vector and the response to functions dealing with missing values (like \code{na_rm}).}
}
\value{
\code{qDF} - returns a data.frame\cr
\code{qDT} - returns a data.table\cr
\code{qM} - returns a matrix\cr
\code{mctl}, \code{mrtl} - return a list, data.frame or data.table\cr
\code{qF} - returns a factor\cr
\code{qG} - returns a quick-group (= integer vector)
}
\seealso{
\code{\link{GRP}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
mtcarsM <- qM(mtcars)      # Matrix from data.frame
mtcarsDT <- qDT(mtcarsM)   # data.table from matrix columns
mrtl(mtcarsM, TRUE, 2L)    # data.table from matrix rows, etc...
qDF(mtcarsM, "cars")       # Adding a row.names column when converting from matrix
qDT(mtcars, "cars")        # Saving row.names when converting data.frame to data.table

cylF <- qF(mtcars$cyl)     # Factor from atomic vector
cylG <- qG(mtcars$cyl)     # Quick-group from atomic vector
cylG                       # See the simple structure of this object
}
\keyword{manip}
\keyword{documentation}
