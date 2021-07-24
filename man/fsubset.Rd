\name{fsubset}
\alias{fsubset}
\alias{sbt}
\alias{ss}
\alias{fsubset.default}
\alias{fsubset.matrix}
\alias{fsubset.data.frame}
\title{Fast Subsetting Matrix-Like Objects}
\description{
 \code{fsubset} returns subsets of vectors, matrices or data frames which meet conditions. It is programmed very efficiently and uses C source code from the \emph{data.table} package. Especially for data frames it is significantly (4-5 times) faster than \code{\link{subset}} or \code{dplyr::filter}. The methods also provide enhanced functionality compared to \code{\link{subset}}. The function \code{ss} provides an (internal generic) programmers alternative to \code{[} that does not drop dimensions and is significantly faster than \code{[.data.frame}.
}
\usage{
fsubset(x, \dots)
sbt(x, \dots)     # Shortcut for fsubset

\method{fsubset}{default}(x, subset, \dots)

\method{fsubset}{matrix}(x, subset, \dots, drop = FALSE)

\method{fsubset}{data.frame}(x, subset, \dots)

# Fast subsetting (replaces `[`, programmers choice)
ss(x, i, j)
}

\arguments{
  \item{x}{object to be subsetted.}
  \item{subset}{logical expression indicating elements or rows to keep:
    missing values are taken as \code{FALSE}. The default and matrix methods only support logical vectors or row-indices (or a character vector of rownames if the matrix has rownames; the data frame method also supports logical vectors or row-indices). }
  \item{\dots}{For the matrix or data frame method: multiple comma-separated expressions indicating columns to select. Otherwise: further arguments to be passed to or from other methods.}
  \item{drop}{passed on to \code{[} indexing operator. Only available for the matrix method.}
  \item{i}{positive or negative row-indices or a logical vector to subset the rows of \code{x}.}
  \item{j}{a vector of column names, positive or negative indices or a suitable logical vector to subset the columns of \code{x}. \emph{Note:} Negative indices are converted to positive ones using \code{j <- seq_along(x)[j]}.}
}
\details{
  \code{fsubset} is a generic function, with methods supplied for vectors, matrices, and data
  frames (including lists). It represents an improvement in both speed and functionality over \code{\link{subset}}. The function \code{ss} is an improvement of \code{[} to subset (vectors) matrices and data frames without dropping dimensions. It is significantly faster than \code{[.data.frame}. For subsetting columns alone, please see \link[=fselect]{selecting and replacing columns}.

  For ordinary vectors, the result is \code{.Call(C_subsetVector, x, subset)}, where \code{C_subsetVector} is an internal function in the \emph{data.table} package. The \code{subset} can be integer or logical. Appropriate errors are delivered for wrong use.

  For matrices the implementation is all base-R but slightly more efficient and more versatile than \code{\link{subset.matrix}}. Thus it is possible to \code{subset} matrix rows using logical or integer vectors, or character vectors matching rownames. The \code{drop} argument is passed on to the indexing method for matrices.

  For both matrices and data frames, the \code{\dots} argument can be used to subset columns, and is evaluated in a non-standard way. Thus it can support vectors of column names, indices or logical vectors, but also multiple comma separated column names passed without quotes, each of which may also be replaced by a sequence of columns i.e. \code{col1:coln}, and new column names may be assigned e.g. \code{fsubset(data, col1 > 20, newname =  col2, col3:col6)} (see examples).

  For data frames, the \code{subset} argument is also evaluated in a non-standard way. Thus next to vector of row-indices or logical vectors, it supports logical expressions of the form \code{col2 > 5 & col2 < col3} etc. (see examples). The data frame method uses \code{C_subsetDT}, an internal C function from the \emph{data.table} package to subset data frames, hence it is significantly faster than \code{\link{subset.data.frame}}. If fast data frame subsetting is required but no non-standard evaluation, the function \code{ss} is slightly simpler and faster.

Factors may have empty levels after subsetting; unused levels are not automatically removed.  See \code{\link{fdroplevels}} for a way to drop all unused levels from a data frame.

}
\value{
  An object similar to \code{x} containing just the selected elements (for
  a vector), rows and columns (for a matrix or data frame).
}
\note{
No replacement method \code{fsubset<-} or \code{ss<-} is offered in \emph{collapse}. For efficient subset relacement (without copying) use \code{data.table::set}, which can also be used with data frames and tibbles.
}
\seealso{
  \code{\link{fselect}},
  \code{\link{get_vars}},
  \code{\link{ftransform}},
\link[=A3-fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
fsubset(airquality, Temp > 90, Ozone, Temp)
fsubset(airquality, Temp > 90, OZ = Ozone, Temp) # With renaming
fsubset(airquality, Day == 1, -Temp)
fsubset(airquality, Day == 1, -(Day:Temp))
fsubset(airquality, Day == 1, Ozone:Wind)
fsubset(airquality, Day == 1 & !is.na(Ozone), Ozone:Wind, Month)

ss(airquality, 1:10, 2:3)         # Significantly faster than airquality[1:10, 2:3]
fsubset(airquality, 1:10, 2:3)    # This is possible but not advised
}
\keyword{manip}
