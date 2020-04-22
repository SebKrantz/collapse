\name{fsubset}
\alias{fsubset}
\alias{ss}
\alias{fsubset.default}
\alias{fsubset.matrix}
\alias{fsubset.data.frame}
\title{Fast Subsetting Vectors, Matrices and Data Frames}
\description{
  Return subsets of vectors, matrices or data frames which meet conditions. \code{fsubset} is programmed very efficiently and uses C source code from the \emph{data.table} package. Especially for data.frame's it is significantly (4-5 times) faster than \code{base::subset} (or \code{dplyr::filter}). The methods also provide more functionality compared to \code{base::subset}.
}
\usage{
fsubset(x, \dots)

\method{fsubset}{default}(x, subset, \dots)

\method{fsubset}{matrix}(x, subset, \dots, drop = FALSE)

\method{fsubset}{data.frame}(x, subset, \dots)

# Convenience function for fast subsetting data.frames without lazy-eval
ss(data, subset, select)
}

\arguments{
  \item{x}{object to be subsetted.}
  \item{data}{a data.frame.}
  \item{subset}{logical expression indicating elements or rows to keep:
    missing values are taken as false. The default and matrix methods only support logical vectors or row-indices (or a character vector of rownames if the matrix has rownames). The data.frame method also supports logical vectors or row-indices. The function \code{\link{ss}} is slightly faster for subsetting data.frame with indices or logical vectors. }
  \item{\dots}{For the matrix data.frame method: multiple comma-separated expressions indicating columns to select. Otherwise: further arguments to be passed to or from other methods.}
  \item{select}{a vector or column names, indices or a suitable logical vector.}
  \item{drop}{passed on to \code{[} indexing operator. Only available for the matrix method.}
}
\details{
  This is a generic function, with methods supplied for matrices, data
  frames and vectors (including lists).

  For ordinary vectors, the result is \code{.Call(C_subsetVector, x, subset)}, where \code{C_subsetVector} is an internal function in the \emph{data.table} package. The \code{subset} can be integer or logical. Appropriate errors are delivered for wrong use.

  For matrices the implementation is all base-R but slightly more efficient and more versatile than \code{base::subset.matrix}. Thus it is possible to \code{subset} matrix rows using logical or integer vectors, or character vectors matching rownames. The \code{drop} argument is passed on to the indexing method for matrices.

  For both matrices and data frames, the \code{...} argument can be used to subset columns, and is evaluated in a non-standard way. Thus it can support vectors of column names, indices or logical vectors, but also multiple comma separated column names passed without quotes, each of which may also be replaced by a sequence of columns i.e. \code{col1:coln} (see examples).

  For data frames, the \code{subset} argument is also evaluated in a non-standard way. Thus next to vector of row-indices or logical vectors, it supports logical expressions of the form \code{col2 > 5 & col2 < col3} etc. (see examples). The data frame method uses \code{C_subsetDT}, an internal C function from the \emph{data.table} package to subset data.frames, hence it is significantly faster than \code{base::subset.data.frame}. If fast data frame subsetting is required but no non-standard evaluation, the convenience function \code{\link{ss}} is slightly simpler and faster.

Factors may have empty levels after subsetting; unused levels are not automatically removed.  See \code{\link{droplevels}} for a way to drop all unused levels from a data frame.

}
\value{
  An object similar to \code{x} containing just the selected elements (for
  a vector), rows and columns (for a matrix or data frame).
}
\seealso{
  \code{\link{[}}, % = ./Extract.Rd
  \code{\link{fselect}},
  \code{\link{get_vars}},
  \code{\link{ftransform}},
  \code{\link{droplevels}}
}
\examples{
fsubset(airquality, Temp > 80, Ozone, Temp)
fsubset(airquality, Day == 1, -Temp)
fsubset(airquality, Day == 1, -(Day:Temp))
fsubset(airquality, Day == 1, Ozone:Wind)
fsubset(airquality, Day == 1 & !is.na(Ozone), Ozone:Wind, Month)

fsubset(airquality, 1:10, 2:3)
ss(airquality, 1:10, 2:3)         # Slightly faster !
}
\keyword{manip}
