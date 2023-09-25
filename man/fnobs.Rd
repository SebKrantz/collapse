\name{fnobs}
\alias{fnobs}
\alias{fnobs.default}
\alias{fnobs.matrix}
\alias{fnobs.data.frame}
\alias{fnobs.grouped_df}
\title{Fast (Grouped) Observation Count for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fnobs} is a generic function that (column-wise) computes the number of non-missing values in \code{x}, (optionally) grouped by \code{g}. It is much faster than \code{sum(!is.na(x))}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped) observation count.
}

\usage{
fnobs(x, \dots)

\method{fnobs}{default}(x, g = NULL, TRA = NULL, use.g.names = TRUE, \dots)

\method{fnobs}{matrix}(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, \dots)

\method{fnobs}{data.frame}(x, g = NULL, TRA = NULL, use.g.names = TRUE, drop = TRUE, \dots)

\method{fnobs}{grouped_df}(x, TRA = NULL, use.g.names = FALSE, keep.group_vars = TRUE, \dots)
}
\arguments{
\item{x}{a vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "NA"     |     1 - "fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
\code{fnobs} preserves all attributes of non-classed vectors / columns, and only the 'label' attribute (if available) of classed vectors / columns (i.e. dates or factors). When applied to data frames and matrices, the row-names are adjusted as necessary.
}
\value{
Integer. The number of non-missing observations in \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its number of non-missing observations, grouped by \code{g}.
}
\seealso{
\code{\link{fndistinct}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
fnobs(airquality$Solar.R)                   # Simple Nobs
fnobs(airquality$Solar.R, airquality$Month) # Grouped Nobs

## data.frame method
fnobs(airquality)
fnobs(airquality, airquality$Month)
fnobs(wlddev)                               # Works with data of all types!
head(fnobs(wlddev, wlddev$iso3c))

## matrix method
aqm <- qM(airquality)
fnobs(aqm)                                  # Also works for character or logical matrices
fnobs(aqm, airquality$Month)
\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
airquality \%>\% group_by(Month) \%>\% fnobs()
wlddev \%>\% group_by(country) \%>\%
           select(PCGDP,LIFEEX,GINI,ODA) \%>\% fnobs()
}
}
\keyword{univar}
\keyword{manip}
