\name{fmedian}
\alias{fmedian}
\alias{fmedian.default}
\alias{fmedian.matrix}
\alias{fmedian.data.frame}
\alias{fmedian.grouped_df}
\title{Fast (Grouped, Weighted) Median Value for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmedian} is a generic function that computes the (column-wise) median value of all values in \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) median value.

}
\usage{
fmedian(x, \dots)

\method{fmedian}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
        use.g.names = TRUE, nthreads = 1L, \dots)

\method{fmedian}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
        use.g.names = TRUE, drop = TRUE, nthreads = 1L, \dots)

\method{fmedian}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
        use.g.names = TRUE, drop = TRUE, nthreads = 1L, \dots)

\method{fmedian}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
        use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
        nthreads = 1L, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values, but only if \code{x} is also missing. }

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{nthreads}{integer. The number of threads to utilize. Parallelism is across groups for grouped computations and at the column-level otherwise. No parallelism is available for weighted computations. }

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain summed weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}
}
\details{
Median value estimation is done using \code{std::nth_element} in C++, which is an efficient partial sorting algorithm. A downside of this is that vectors need to be copied first and then partially sorted, thus \code{fmedian} currently requires additional memory equal to the size of the vector (\code{x} or a column of \code{x}).

Grouped computations are currently performed by mapping the data to a sparse-array and then partially sorting each row (group) of that array. Because of compiler optimizations this requires less memory than a full deep copy done with no groups.  % For reasons I don't fully understand this requires less memory than a full deep copy which is done with no groups.

The weighted median is defined as the element \code{k} from a set of sorted elements, such that the sum of weights of all elements larger and all elements smaller than k is \code{<= sum(w)/2}. If the half-sum of weights (\code{sum(w)/2}) is reached exactly for some element k, then (summing from the lower end) both k and k+1 would qualify as the weighted median (and some possible additional elements with zero weights following k would also qualify). \code{fmedian} solves these ties by taking a simple arithmetic mean of all elements qualifying as the weighted median.

The weighted median is computed using \code{\link{radixorder}} to first obtain an ordering of all elements, so it is considerably more computationally expensive than the unweighted version. With groups, the entire vector is also ordered, and the weighted median is computed in a single ordered pass through the data (after group-summing the weights, skipping weights for which \code{x} is missing).

If \code{x} is a matrix or data frame, these computations are performed independently for each column. When applied to data frames with groups or \code{drop = FALSE}, \code{fmedian} preserves all column attributes. The attributes of the data frame itself are also preserved.


%smaller to the left and to the right of
%With weights,
%a lower weighted median is calculated using a fast radix-sort of the values in each column. The lower weighted median is defined as the smallest element \code{sort(x)[i]} for which \code{sum(w[order(x)][-seq_len(i)]) <= sum(w)/2}.

}
\value{
The (\code{w} weighted) median value of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) median.
}
\seealso{
\code{\link{fnth}}, \code{\link{fmean}}, \code{\link{fmode}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fmedian(mpg)                         # Simple median value
fmedian(mpg, w = mtcars$hp)          # Weighted median: Weighted by hp
fmedian(mpg, TRA = "-")              # Simple transformation: Subtract median value
fmedian(mpg, mtcars$cyl)             # Grouped median value
fmedian(mpg, mtcars[c(2,8:9)])       # More groups..
g <- GRP(mtcars, ~ cyl + vs + am)    # Precomputing groups gives more speed !
fmedian(mpg, g)
fmedian(mpg, g, mtcars$hp)           # Grouped weighted median
fmedian(mpg, g, TRA = "-")           # Groupwise subtract median value
fmedian(mpg, g, mtcars$hp, "-")      # Groupwise subtract weighted median value

## data.frame method
fmedian(mtcars)
head(fmedian(mtcars, TRA = "-"))
fmedian(mtcars, g)
fmedian(fgroup_by(mtcars, cyl, vs, am))   # Another way of doing it..
fmedian(mtcars, g, use.g.names = FALSE)   # No row-names generated

## matrix method
m <- qM(mtcars)
fmedian(m)
head(fmedian(m, TRA = "-"))
fmedian(m, g) # etc..
%\donttest{ % No code relying on suggested package
library(dplyr)
# grouped_df method
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian()
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian(hp)             # Weighted
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fmedian()              # Faster grouping!
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fmedian(TRA = "-")     # De-median
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg, hp) \%>\%    # Faster selecting
      fmedian(hp, "-")  # Weighted de-median mpg, using hp as weights
%}
}
\keyword{univar}
\keyword{manip}
