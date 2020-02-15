\name{ffirst-flast}
\alias{ffirst}
\alias{ffirst.default}
\alias{ffirst.matrix}
\alias{ffirst.data.frame}
\alias{ffirst.grouped_df}
\alias{flast}
\alias{flast.default}
\alias{flast.matrix}
\alias{flast.data.frame}
\alias{flast.grouped_df}
\title{Fast (Grouped) First and Last Value for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{ffirst} and \code{flast} are S3 generic functions that (column-wise) returns the first and last values in \code{x}, (optionally) grouped by \code{g}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (groupwise) first and last values.
}

\usage{
ffirst(x, ...)
flast(x, ...)
%## Default S3 methods:
\method{ffirst}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, ...)
\method{flast}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, ...)
%## S3 methods for class 'matrix'
\method{ffirst}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, drop = TRUE, ...)
\method{flast}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ...)
%## S3 methods for class 'data.frame'
\method{ffirst}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, drop = TRUE, ...)
\method{flast}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ...)
%## S3 methods for class 'grouped_df'
\method{ffirst}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
       use.g.names = FALSE, keep.group_vars = TRUE, ...)
\method{flast}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
      use.g.names = FALSE, keep.group_vars = TRUE, ...)
}
\arguments{
\item{x}{a vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values and choose the first / last non-missing value i.e. if the first (1) / last (n) value is NA, take the second (2) / second-to-last (n-1) value etc...}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{...}{arguments to be passed to or from other methods.}

}
\value{
\code{ffirst} returns the first value in \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its first value, grouped by \code{g}. Similarly \code{flast} returns the last value in \code{x}, ...
}
\seealso{
\link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
ffirst(airquality$Ozone)                   # Simple first value
ffirst(airquality$Ozone, airquality$Month) # Grouped first value
ffirst(airquality$Ozone, airquality$Month,
       na.rm = FALSE)                      # Grouped first, but without skipping initial NA's

## data.frame method
ffirst(airquality)
ffirst(airquality, airquality$Month)
ffirst(airquality, airquality$Month, na.rm = FALSE) # Again first Ozone measurement in month 6 is NA

## matrix method
aqm <- qM(airquality)
ffirst(aqm)
ffirst(aqm, airquality$Month) # etc...

## method for grouped tibbles - for use with dplyr
library(dplyr)
airquality \%>\% group_by(Month) \%>\% ffirst
airquality \%>\% group_by(Month) \%>\% select(Ozone) \%>\% ffirst(na.rm = FALSE)

# Note: All examples generalize to flast!
}
\keyword{univar}
\keyword{manip}
