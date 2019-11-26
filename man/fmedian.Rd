\name{fmedian}
\alias{fmedian}
\alias{fmedian.default}
\alias{fmedian.matrix}
\alias{fmedian.data.frame}
\alias{fmedian.grouped_df}
\title{Fast (Grouped) Median Value for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmedian} is a generic function that computes the (column-wise) median value of all values in \code{x}, (optionally) grouped by \code{g}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped) median value.
}
\usage{
fmedian(x, ...)

\method{fmedian}{default}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
        use.g.names = TRUE, ...)

\method{fmedian}{matrix}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
        use.g.names = TRUE, drop = TRUE, ...)

\method{fmedian}{data.frame}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
        use.g.names = TRUE, drop = TRUE, ...)

\method{fmedian}{grouped_df}(x, TRA = FALSE, na.rm = TRUE,
        use.g.names = FALSE, keep.group_vars = TRUE, ...)
}
\arguments{
\item{x}{a numeric vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = FALSE}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{...}{arguments to be passed to or from other methods.}
}
\details{
Median value estimation is done using \code{std::nth_element} in C++, which is an efficient partial sorting algorithm. A downside of this is that vectors need to be copied first and then partially sorted, thus \code{fmedian} currently requires additional memory equal to the size of the object (\code{x}).

Grouped computations are currently performed by mapping the data to a sparse-array directed by \code{g} and then partially sorting each row (group) of that array. For reasons I don't fully understand this requires less memory than a full deep copy which is done with no groups.

When applied to data frame's with groups or \code{drop = FALSE}, \code{fmedian} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed objects. The attributes of the \code{data.frame} itself are also perserved.

}
\value{
The median value of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its median value, grouped by \code{g}.
}
\seealso{
\code{\link{fmean}}, \link{Fast Statistical Functions}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fmedian(mpg)                         # Simple median value
fmedian(mpg, TRA = "-")              # Simple transformation: Subtract median value
fmedian(mpg, mtcars$cyl)             # Grouped median value
fmedian(mpg, mtcars[c(2,8:9)])       # More groups...
g <- GRP(mtcars, ~ cyl + vs + am)    # Precomputing groups gives more speed !!
fmedian(mpg, g)
fmedian(mpg, g, TRA = "-")           # Groupwise subtract median value

## data.frame method
fmedian(mtcars)
fmedian(mtcars, TRA = "-")
fmedian(mtcars, g)
fmedian(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fmedian(m)
fmedian(m, TRA = "-")
fmedian(m, g) # etc...

## method for grouped tibbles - for use with dplyr
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian("-")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fmedian
}

