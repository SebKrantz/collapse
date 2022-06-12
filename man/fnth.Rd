\name{fnth}
\alias{fnth}
\alias{fnth.default}
\alias{fnth.matrix}
\alias{fnth.data.frame}
\alias{fnth.grouped_df}
\title{
Fast (Grouped, Weighted) N'th Element/Quantile for Matrix-Like Objects
}
\description{
\code{fnth} (column-wise) returns the n'th smallest element from a set of unsorted elements \code{x} corresponding to an integer index (\code{n}), or to a probability between 0 and 1. If \code{n} is passed as a probability, ties can be resolved using the lower, upper, or (default) average of the possible elements. These are discontinuous and fast methods to estimate a sample quantile.
}
\usage{
fnth(x, n = 0.5, \dots)

\method{fnth}{default}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, ties = "mean", nthreads = 1L, \dots)

\method{fnth}{matrix}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, ties = "mean", nthreads = 1L, \dots)

\method{fnth}{data.frame}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, ties = "mean", nthreads = 1L, \dots)

\method{fnth}{grouped_df}(x, n = 0.5, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
     ties = "mean", nthreads = 1L, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{n}{the element to return using a single integer index such that \code{1 < n < NROW(x)}, or a probability \code{0 < n < 1}. See Details. }

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{nthreads}{integer. The number of threads to utilize. Parallelism is across groups for grouped computations and at the column-level otherwise. No parallelism is available for weighted computations. }

\item{ties}{an integer or character string specifying the method to resolve ties between adjacent qualifying elements:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "mean"   \tab\tab take the arithmetic mean of all qualifying elements. \cr
                 2 \tab\tab "min" \tab\tab take the smallest of the elements. \cr
                 3 \tab\tab "max"   \tab\tab take the largest of the elements. \cr
                }
  }

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain \code{sum} of weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
This is an R port to \code{std::nth_element}, an efficient partial sorting algorithm in C++. It is also used to calculated the median (in fact the default \code{fnth(x, n = 0.5)} is identical to \code{fmedian(x)}, so see also the details for \code{\link{fmedian}}).

\code{fnth} generalizes the principles of median value calculation to find arbitrary elements. It offers considerable flexibility by providing both simple order statistics and simple discontinuous quantile estimation. Regarding the former, setting \code{n} to an index between 1 and \code{NROW(x)} will return the n'th smallest element of \code{x}, about 2x faster than \code{sort(x, partial = n)[n]}. As to the latter, setting \code{n} to a probability between 0 and 1 will return the corresponding element of \code{x}, and resolve ties between multiple qualifying elements (such as when \code{n = 0.5} and \code{x} is even) using the arithmetic average \code{ties = "mean"}, or the smallest \code{ties = "min"} or largest \code{ties = "max"} of those elements.

If \code{n > 1} is used and \code{x} contains missing values (and \code{na.rm = TRUE}, otherwise \code{NA} is returned), \code{n} is internally converted to a probability using \code{p = (n-1)/(NROW(x)-1)}, and that probability is applied to the set of complete elements (of each column if \code{x} is a matrix or data frame) to find the \code{as.integer(p*(fnobs(x)-1))+1L}'th element (which corresponds to option \code{ties = "min"}). Note that it is necessary to subtract and add 1 so that \code{n = 1} corresponds to \code{p = 0} and \code{n = NROW(x)} to \code{p = 1}. %So if \code{n > 1} is used in the presence of missing values, and the default \code{ties = "mean"} is enabled, the resulting element could be the average of two elements.

When using grouped computations (supplying a vector or list to \code{g} subdividing \code{x}) and \code{n > 1} is used, it is transformed to a probability \code{p = (n-1)/(NROW(x)/ng-1)} (where \code{ng} contains the number of unique groups in \code{g}) and \code{ties = "min"} is used to sort out clashes. This could be useful for example to return the n'th smallest element of each group in a balanced panel, but with unequal group sizes it more intuitive to pass a probability to \code{n}.

If weights are used, the same principles apply as for weighted median calculation: A target partial sum of weights \code{p*sum(w)} is calculated, and the weighted n'th element is the element k such that all elements smaller than k have a sum of weights \code{<= p*sum(w)}, and all elements larger than k have a sum of weights \code{<= (1 - p)*sum(w)}. If the partial-sum of weights (\code{p*sum(w)}) is reached exactly for some element k, then (summing from the lower end) both k and k+1 would qualify as the weighted n'th element (and some possible additional elements with zero weights following k would also qualify). If \code{n > 1}, the lowest of those elements is chosen (congruent with the unweighted behavior), %(ensuring that \code{fnth(x, n)}) and \code{fnth(x, n, w = rep(1, NROW(x)))}, always provide the same outcome)
but if \code{0 < n < 1}, the \code{ties} option regulates how to resolve such conflicts, yielding lower-weighted, upper-weighted or (default) average weighted n'th elements.

The weighted n'th element is computed using \code{\link{radixorder}} to first obtain an ordering of all elements, so it is considerably more computationally expensive than the unweighted version. With groups, the entire vector is also ordered, and the weighted n'th element is computed in a single ordered pass through the data (after calculating partial-group sums of the weights, skipping weights for which \code{x} is missing).

If \code{x} is a matrix or data frame, these computations are performed independently for each column. Column-attributes and overall attributes of a data frame are preserved (if \code{g} is used or \code{drop = FALSE}).

}
\value{
The (\code{w} weighted) n'th element of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) n'th element.

}
\seealso{
\code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fnth(mpg)                         # Simple nth element: Median (same as fmedian(mpg))
fnth(mpg, 5)                      # 5th smallest element
sort(mpg, partial = 5)[5]         # Same using base R, fnth is 2x faster.
fnth(mpg, 0.75)                   # Third quartile
fnth(mpg, 0.75, w = mtcars$hp)    # Weighted third quartile: Weighted by hp
fnth(mpg, 0.75, TRA = "-")        # Simple transformation: Subtract third quartile
fnth(mpg, 0.75, mtcars$cyl)             # Grouped third quartile
fnth(mpg, 0.75, mtcars[c(2,8:9)])       # More groups..
g <- GRP(mtcars, ~ cyl + vs + am)       # Precomputing groups gives more speed !
fnth(mpg, 0.75, g)
fnth(mpg, 0.75, g, mtcars$hp)           # Grouped weighted third quartile
fnth(mpg, 0.75, g, TRA = "-")           # Groupwise subtract third quartile
fnth(mpg, 0.75, g, mtcars$hp, "-")      # Groupwise subtract weighted third quartile

## data.frame method
fnth(mtcars, 0.75)
head(fnth(mtcars, 0.75, TRA = "-"))
fnth(mtcars, 0.75, g)
fnth(fgroup_by(mtcars, cyl, vs, am), 0.75)   # Another way of doing it..
fnth(mtcars, 0.75, g, use.g.names = FALSE)   # No row-names generated

## matrix method
m <- qM(mtcars)
fnth(m, 0.75)
head(fnth(m, 0.75, TRA = "-"))
fnth(m, 0.75, g) # etc..
\donttest{ % No code relying on suggested package
library(dplyr)
## grouped_df method
mtcars \%>\% group_by(cyl,vs,am) \%>\% fnth(0.75)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fnth(0.75, hp)           # Weighted
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fnth(0.75)              # Faster grouping!
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fnth(0.75, TRA = "/")   # Divide by third quartile
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg, hp) \%>\%     # Faster selecting
      fnth(0.75, hp, "/")  # Divide mpg by its third weighted group-quartile, using hp as weights
}
}
\keyword{univar}
\keyword{manip}
