\name{fmode}
\alias{fmode}
\alias{fmode.default}
\alias{fmode.matrix}
\alias{fmode.data.frame}
\alias{fmode.grouped_df}
\title{Fast (Grouped, Weighted) Statistical Mode for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmode} is a generic function and returns the (column-wise) statistical mode i.e. the most frequent value of \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w}.
The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) mode.
}
\usage{
fmode(x, ...)

\method{fmode}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, ...)

\method{fmode}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ...)

\method{fmode}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ...)

\method{fmode}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE, ...)
}
\arguments{
\item{x}{a vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE}, \code{NA} is treated as any other value.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain \code{max} of weighting variable after computation (if contained in \code{grouped_df}).}

\item{...}{arguments to be passed to or from other methods.}

}
\details{
\code{fmode} implements a pretty fast algorithm to find the statistical mode utilizing index- hashing implemented in the \code{Rcpp::sugar::IndexHash} class.

If all values are distinct, the first value is returned. If there are multiple distinct values having the top frequency, the first value established as having the top frequency when passing through the data from element 1 to element n is returned. If \code{na.rm = FALSE}, \code{NA} is not removed but treated as any other value (i.e. it's frequency is counted). If all values are \code{NA}, \code{NA} is always returned.

The weighted mode is computed by summing up the weights for all distinct values and choosing the value with the largest sum. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

This all seamlessly generalizes to grouped computations, which are currently performed by mapping the data to a sparse-array directed by \code{g} and then going group-by group.

\code{fmode} preserves all the attributes of the objects it is applied to (apart from names or row-names which are adjusted as necessary). If a data frame is passed to \code{fmode} and \code{drop = TRUE}, \code{base::unlist} will be called on the result, which might or might not be sensible depending on the data at hand.
}
\value{
The statistical mode of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its mode, grouped by \code{g}. See also Details.
}
\seealso{
\code{\link{fmean}}, \code{\link{fmedian}}, \link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## World Development Data
attach(wlddev)
## default vector method
fmode(PCGDP)                    # Numeric mode
fmode(PCGDP, iso3c)             # Grouped numeric mode
fmode(PCGDP, iso3c, LIFEEX)     # Grouped and weighted numeric mode
fmode(region)                   # Factor mode
fmode(date)                     # Date mode (defaults to first value since panel is balanced)
fmode(country)                  # Character mode (also defaults to first value)
fmode(OECD)                     # Logical mode
                                # ...all the above can also be performed grouped and weighted
## matrix method
m <- qM(airquality)
fmode(m)
fmode(m, na.rm = FALSE)         # NA frequency is also counted
fmode(m, airquality$Month)      # Groupwise
fmode(m, w = airquality$Day)    # Weighted: Later days in the month are given more weight
fmode(m>50, airquality$Month)   # Groupwise logical mode
                                # etc ...
## data.frame method
fmode(wlddev)                   # Gives one row
fmode(wlddev, drop = TRUE)      # calling unlist -> coerce to character vector
fmode(wlddev, iso3c)            # Grouped mode
fmode(wlddev, iso3c, LIFEEX)    # Grouped and weighted mode

detach(wlddev)
}
\keyword{univar}
\keyword{manip}
