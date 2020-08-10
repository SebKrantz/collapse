\name{A4-quick-conversion}
\alias{A4-quick-conversion}
\alias{qDF}
\alias{qDT}
\alias{qM}
\alias{mctl}
\alias{mrtl}
\alias{as.numeric_factor}
\alias{as.character_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quick Data Conversion}
\description{
Convert common data objects quickly, without method dispatch and extensive checks:
\itemize{
\item \code{qDF} and \code{qDT} convert vectors, matrices, higher-dimensional arrays and suitable lists to data frame and \emph{data.table}, respectively.
\item \code{qM} converts vectors, higher-dimensional arrays, data frames and suitable lists to matrix.
\item \code{mctl} and \code{mrtl} column- or row-wise convert a matrix to list, data frame or \emph{data.table}. They are used internally by \code{qDF} and \code{qDT}, \code{\link{dapply}}, \code{\link{BY}}, etc\dots
\item \code{\link{qF}} converts atomic vectors to factor (documented on a separate page).
\item \code{as.numeric_factor} and \code{as.character_factor} convert factors, or all factor columns in a data frame / list, to numeric or character (by converting the levels).
}
}
\usage{
qDF(X, row.names.col = FALSE)
qDT(X, row.names.col = FALSE)
qM(X)
mctl(X, names = FALSE, return = "list")
mrtl(X, names = FALSE, return = "list")
as.numeric_factor(X, keep.attr = TRUE)
as.character_factor(X, keep.attr = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, factor, matrix, higher-dimensional array, data frame or list. \code{mctl} and \code{mrtl} only accept matrices, \code{as.numeric_factor} and \code{as.character_factor} only accept factors, data frames or lists.}
  \item{row.names.col}{should a column capturing names or row.names be added? i.e. when converting atomic objects to data frame or data frame to \emph{data.table}. Can be logical \code{TRUE}, which will add a column \code{"row.names"} in front, or can supply a name for the column i.e. \code{"column1"}.}
  \item{names}{logical. Should the list be named using row/column names from the matrix?}
    \item{return}{an integer or string specifying what to return. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "list"   \tab\tab returns a plain list \cr
                 2 \tab\tab "data.frame" \tab\tab returns a data frame \cr
                 3 \tab\tab "data.table"   \tab\tab returns a \emph{data.table} \cr
  }
  }
  \item{keep.attr}{logical. \code{TRUE} keeps all attributes of factor variables apart from the levels and class attributes (such as variable labels etc.).}
}
\value{
\code{qDF} - returns a data.frame\cr
\code{qDT} - returns a data.table\cr
\code{qM} - returns a matrix\cr
\code{mctl}, \code{mrtl} - return a list, data frame or \emph{data.table} \cr
\code{qF} - returns a factor\cr
\code{as.numeric_factor} - returns X with factors converted to numeric variables\cr
\code{as.character_factor} - returns X with factors converted to character variables
}
\seealso{
\code{\link{qF}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
mtcarsM <- qM(mtcars)                   # Matrix from data.frame
mtcarsDT <- qDT(mtcarsM)                # data.table from matrix columns
head(mrtl(mtcarsM, TRUE, "data.frame")) # data.frame from matrix rows, etc..
head(qDF(mtcarsM, "cars"))              # Adding a row.names column when converting from matrix
head(qDT(mtcars, "cars"))               # Saving row.names when converting data frame to data.table

cylF <- qF(mtcars$cyl)                  # Factor from atomic vector
cylF

## Factor to numeric conversions
identical(mtcars,  as.numeric_factor(dapply(mtcars, qF)))
}
\keyword{manip}
\keyword{documentation}
