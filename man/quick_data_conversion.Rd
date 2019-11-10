\name{quick_data_conversion}
\alias{quick_data_conversion}
\alias{qDF}
\alias{qDT}
\alias{qM}
\alias{qF}
\alias{qG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quick Data Conversion}
\description{
Convert common data types quickly, without method dispatch and extensive checks:
\itemize{
\item \code{qDF} and \code{qDT} convert vectors, matrices, higher-dimensional arrays and suitable lists to \code{data.frame} and \code{data.table} respectively.
\item \code{qM} converts vectors, higher-dimensional arrays, data.frames and suitable lists to matrix.
\item \code{mctl} and \code{mrtl} are exported C++ functions that column -or row-wise convert a matrix to list, data.frame or data.table.
\item \code{qF} converts atomic vectors to factor.
\item \code{qG} converts atomic vectors to a quick-group (class attribute 'qG') - a kind of factor light which is just an integer vector with an attribute 'N.groups' indicating the number of groups. (This is convenient especially for numeric vectors because it saves converting the levels to character which can require a lot of time).
}
}
\usage{
qDF(X)
qDT(X)
qM(X)
mctl(X, names = FALSE, ret = 0L)
mrtl(X, names = FALSE, ret = 0L)
qF(x, ordered = TRUE)
qG(x, ordered = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A vector, matrix, higher-dimensional array, data.frame or list. \code{mctl} and \code{mrtl} only take matrices.}
  \item{x}{A atomic vector.}
  \item{ordered}{should groups be ordered, so that 1 represents the first group, 2 the second etc.?}
  \item{names}{should the list be named?}
  \item{ret}{3 return options: 0L - return a list, 1L - return a data.frame, 2L - return a data.table.}
}
\value{
\code{qDF} - returns a data.frame\cr
\code{qDT} - returns a data.table\cr
\code{qM} - returns a matrix\cr
\code{mctl}, \code{mrtl} - return a list, data.frame or data.table\cr
\code{qF} - returns a factor\cr
\code{qG} - returns a quick-group (= integer vector)
}
\seealso{\code{\link{GRP}}}
\examples{
mtcarsM <- qM(mtcars)      # Matrix from data.frame
mtcarsDT <- qDT(mtcarsmat) # data.table from matrix columns, etc...
mrtl(mtcarsM, TRUE, 2L)    # data.table from matrix rows, etc...

cylF <- qF(mtcars$cyl)     # Factor from atomic vector
cylG <- qG(mtcars$cyl)     # Quick-group from atomic vector
cylG                       # See the simple structure of this object
}
