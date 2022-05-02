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
ffirst(x, \dots)
flast(x, \dots)

\method{ffirst}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, \dots)
\method{flast}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, \dots)

\method{ffirst}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, drop = TRUE, \dots)
\method{flast}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, \dots)

\method{ffirst}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
       use.g.names = TRUE, drop = TRUE, \dots)
\method{flast}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, \dots)

\method{ffirst}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
       use.g.names = FALSE, keep.group_vars = TRUE, \dots)
\method{flast}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
      use.g.names = FALSE, keep.group_vars = TRUE, \dots)
}
\arguments{
\item{x}{a vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. \code{TRUE} skips missing values and returns the first / last non-missing value i.e. if the first (1) / last (n) value is \code{NA}, take the second (2) / second-to-last (n-1) value etc..}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\value{
\code{ffirst} returns the first value in \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its first value, grouped by \code{g}. Similarly \code{flast} returns the last value in \code{x}, \dots
}
\note{
Both functions are significantly faster if \code{na.rm = FALSE}, particularly \code{ffirst} which can take direct advantage of the 'group.starts' elements in \code{\link{GRP}} objects.
}
\seealso{
\link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
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
ffirst(aqm, airquality$Month) # etc..
\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
airquality \%>\% group_by(Month) \%>\% ffirst()
airquality \%>\% group_by(Month) \%>\% select(Ozone) \%>\% ffirst(na.rm = FALSE)
}
# Note: All examples generalize to flast.
}
\keyword{univar}
\keyword{manip}
