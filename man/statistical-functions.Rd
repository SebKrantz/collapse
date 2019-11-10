\name{fast statistical functions}
\alias{fstats}
\alias{fast-statistical-functions}
\alias{fast statistical functions}
\title{Fast (Grouped, Weighted) Statistical Functions for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
With \code{fsum}, \code{fprod}, \code{fmean}, \code{fmedian}, \code{fmode}, \code{fvar}, \code{fsd}, \code{fmin}, \code{fmax}, \code{ffirst}, \code{flast}, \code{fNobs} and \code{fNdistinct}, (and the vector-valued \code{fscale}, \code{B} and \code{W}) \code{collapse} presents a coherent set of extremely fast and flexible statistcal functions (S3 generics) to perform column-wise, grouped and weighted computuations on matrix-like objects (atomic vectors, matrices and data.frames, with special support for \code{dplyr} grouped tibbles and data.table's). 

}
\usage{
## All functions (FUN) follow a common syntax in 4 methods:

\method{FUN}{default}(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,] 
    use.g.names = TRUE, ...)

\method{FUN}{matrix}(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,] 
    drop = TRUE, use.g.names = TRUE, ...)

\method{FUN}{data.frame}(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,] 
    drop = TRUE, use.g.names = TRUE, ...)

\method{FUN}{grouped_df}(x, [w = NULL,] TRA = FALSE, [na.rm = TRUE,] 
    drop.groups = FALSE, [drop.w = TRUE,] ...)
}
\arguments{
\item{x}{a vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values. Supported by \code{fmean, fvar, fsd} and \code{fmode}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to TRUE in all functions and implemented at very little computational cost. Not available for \code{fNobs}.}

\item{drop}{\emph{matrix and data.frame method:} Drop dimensions and return an atomic vector, if \code{g = NULL} and \code{TRA = FALSE}.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop.groups}{\emph{grouped_df method:} Logical. Remove grouping variables after computation.}

\item{drop.w}{\emph{grouped_df method:} Logical. Remove weighting variable after computation (if contained in \code{grouped_df}).}

\item{...}{extra arguments to some functions, i.e. the algorithm used.}
}
\details{
Please see the documentation of individual functions. 
}
\value{
The computed value of \code{x}. 
}
\examples{
## default vector method
mpg <- mtcars$mpg
fsum(mpg)                         # Simple sum
fsum(mpg, TRA = "\%")              # Simple transformation: obtain percentages of mpg
fsum(mpg, mtcars$cyl)             # Grouped sum
fsum(mpg, mtcars$cyl, TRA = "\%")  # Percentages by group

## data.frame Method
fsum(mtcars)
fsum(mtcars, TRA = "\%")
fsum(mtcars, mtcars[c(2,8:9)])           # Grouped column sum
g <- GRP(mtcars, ~ cyl + vs + am)        # Here precomputing the groups!
fsum(mtcars, g)                          # Faster !!
fsum(mtcars, g, TRA = "\%")

## matrix method
m <- qM(mtcars)
fsum(m)
fsum(m, g) # ...

## method for grouped tibbles - for use with dplyr:
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg,carb) \%>\% fsum
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsum("\%")

}

