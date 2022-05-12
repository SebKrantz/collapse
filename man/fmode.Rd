\name{fmode}
\alias{fmode}
\alias{fmode.default}
\alias{fmode.matrix}
\alias{fmode.data.frame}
\alias{fmode.grouped_df}
\title{Fast (Grouped, Weighted) Statistical Mode for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmode} is a generic function and returns the (column-wise) statistical mode i.e. the most frequent value of \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w}.
The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) mode. Ties between multiple possible modes can be resolved by taking the minimum, maximum, (default) first or last occurring mode.
}
\usage{
fmode(x, \dots)

\method{fmode}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, ties = "first", nthreads = 1L, \dots)

\method{fmode}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ties = "first", nthreads = 1L, \dots)

\method{fmode}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, ties = "first", nthreads = 1L, \dots)

\method{fmode}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
      ties = "first", nthreads = 1L, \dots)
}
\arguments{
\item{x}{a vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE}, \code{NA} is treated as any other value.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{ties}{an integer or character string specifying the method to resolve ties between multiple possible modes i.e. multiple values with the maximum frequency or sum of weights:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "first"   \tab\tab take the first occurring mode. \cr
                 2 \tab\tab "min" \tab\tab take the smallest of the possible modes. \cr
                 3 \tab\tab "max"   \tab\tab take the largest of the possible modes. \cr
                 4 \tab\tab "last"   \tab\tab take the last occurring mode. \cr
                }
        \emph{Note:} \code{"min"/"max"} don't work with character data. % For logical data \code{TRUE} will be chosen unless \code{ties = "min"}.
        See also Details.
  }

  \item{nthreads}{integer. The number of threads to utilize. Parallelism is across groups for grouped computations and at the column-level otherwise. }

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain \code{sum} of weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
\code{fmode} implements a pretty fast C-level hashing algorithm inspired by the \emph{kit} package to find the statistical mode. % utilizing index- hashing implemented in the \code{Rcpp::sugar::IndexHash} class.

%If all values are distinct, the first value is returned. If there are multiple distinct values having the top frequency, the first value established as having the top frequency when passing through the data from element 1 to element n is returned.
If \code{na.rm = FALSE}, \code{NA} is not removed but treated as any other value (i.e. it's frequency is counted). If all values are \code{NA}, \code{NA} is always returned.

The weighted mode is computed by summing up the weights for all distinct values and choosing the value with the largest sum. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

It is possible that multiple values have the same mode (the maximum frequency or sum of weights). Typical cases are simply when all values are either all the same or all distinct. In such cases, the default option \code{ties = "first"} returns the first occurring value in the data reaching the maximum frequency count or sum of weights. For example in a sample \code{x = c(1, 3, 2, 2, 4, 4, 1, 7)}, the first mode is 2 as \code{fmode} goes through the data from left to right. \code{ties = "last"} on the other hand gives 1. It is also possible to take the minimum or maximum mode, i.e. \code{fmode(x, ties = "min")} returns 1, and \code{fmode(x, ties = "max")} returns 4. It should be noted that options \code{ties = "min"} and \code{ties = "max"} give unintuitive results for character data (no strict alphabetic sorting, similar to using \code{<} and \code{>} to compare character values in R). These options are also best avoided if missing values are counted (\code{na.rm = FALSE}) since no proper logical comparison with missing values is possible: With numeric data it depends, since in C++ any comparison with \code{NA_real_} evaluates to \code{FALSE}, \code{NA_real_} is chosen as the min or max mode only if it is also the first mode, and never otherwise. For integer data, \code{NA_integer_} is stored as the smallest integer in C++, so it will always be chosen as the min mode and never as the max mode. For character data, \code{NA_character_} is stored as the string \code{"NA"} in C++ and thus the behavior depends on the other character content. \code{fmode} also implements a fast method for logical values which does not support the options \code{"first"/"last"} i.e. \code{TRUE} is returned unless \code{ties = "min"}.

% This all seamlessly generalizes to grouped computations, which are performed by mapping the data to a sparse-array (except for logical values) and then going group-by group.

\code{fmode} preserves all the attributes of the objects it is applied to (apart from names or row-names which are adjusted as necessary in grouped operations). If a data frame is passed to \code{fmode} and \code{drop = TRUE} (the default), \code{\link{unlist}} will be called on the result, which might not be sensible depending on the data at hand.
}
\value{
The (\code{w} weighted) statistical mode of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighed) mode. %See also Details.
}
\seealso{
\code{\link{fmean}}, \code{\link{fmedian}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
x <- c(1, 3, 2, 2, 4, 4, 1, 7, NA, NA, NA)
fmode(x)                            # Default is ties = "first"
fmode(x, ties = "last")
fmode(x, ties = "min")
fmode(x, ties = "max")
fmode(x, na.rm = FALSE)             # Here NA is the mode, regardless of ties option
fmode(x[-length(x)], na.rm = FALSE) # Not anymore..

## World Development Data
attach(wlddev)
## default vector method
fmode(PCGDP)                      # Numeric mode
head(fmode(PCGDP, iso3c))         # Grouped numeric mode
head(fmode(PCGDP, iso3c, LIFEEX)) # Grouped and weighted numeric mode
fmode(region)                     # Factor mode
fmode(date)                       # Date mode (defaults to first value since panel is balanced)
fmode(country)                    # Character mode (also defaults to first value)
fmode(OECD)                       # Logical mode
                                  # ..all the above can also be performed grouped and weighted
## matrix method
m <- qM(airquality)
fmode(m)
fmode(m, na.rm = FALSE)         # NA frequency is also counted
fmode(m, airquality$Month)      # Groupwise
fmode(m, w = airquality$Day)    # Weighted: Later days in the month are given more weight
fmode(m>50, airquality$Month)   # Groupwise logical mode
                                # etc..
## data.frame method
fmode(wlddev)                      # Calling unlist -> coerce to character vector
fmode(wlddev, drop = FALSE)        # Gives one row
head(fmode(wlddev, iso3c))         # Grouped mode
head(fmode(wlddev, iso3c, LIFEEX)) # Grouped and weighted mode

detach(wlddev)
}
\keyword{univar}
\keyword{manip}
