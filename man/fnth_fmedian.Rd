\name{fnth-fmedian}
\alias{fnth}
\alias{fnth.default}
\alias{fnth.matrix}
\alias{fnth.data.frame}
\alias{fnth.grouped_df}
\alias{fmedian}
\alias{fmedian.default}
\alias{fmedian.matrix}
\alias{fmedian.data.frame}
\alias{fmedian.grouped_df}
\title{
Fast (Grouped, Weighted) N'th Element/Quantile for Matrix-Like Objects
}
\description{
\code{fnth} (column-wise) returns the n'th smallest element from a set of unsorted elements \code{x} corresponding to an integer index (\code{n}), or to a probability between 0 and 1. If \code{n} is passed as a probability, ties can be resolved using the lower, upper, or average of the possible elements, or, since v1.9.0, continuous quantile estimation. The new default is quantile type 7 (as in \code{\link{quantile}}). For \code{n > 1}, the lower element is always returned (as in \code{sort(x, partial = n)[n]}). See Details.

\code{fmedian} is a simple wrapper around \code{fnth}, which fixes \code{n = 0.5} and (default) \code{ties = "mean"} i.e. it averages eligible elements. See Details. %Users may prefer a quantile based definition of the weighted median.
}
\usage{
fnth(x, n = 0.5, \dots)
fmedian(x, \dots)

\method{fnth}{default}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
     use.g.names = TRUE, ties = "q7", nthreads = .op[["nthreads"]],
     o = NULL, check.o = is.null(attr(o, "sorted")), \dots)
\method{fmedian}{default}(x, \dots, ties = "mean")

\method{fnth}{matrix}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
     use.g.names = TRUE, drop = TRUE, ties = "q7", nthreads = .op[["nthreads"]], \dots)
\method{fmedian}{matrix}(x, \dots, ties = "mean")

\method{fnth}{data.frame}(x, n = 0.5, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
     use.g.names = TRUE, drop = TRUE, ties = "q7", nthreads = .op[["nthreads"]], \dots)
\method{fmedian}{data.frame}(x, \dots, ties = "mean")

\method{fnth}{grouped_df}(x, n = 0.5, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
     use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
     ties = "q7", nthreads = .op[["nthreads"]], \dots)
\method{fmedian}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
        use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
        ties = "mean", nthreads = .op[["nthreads"]], \dots)

}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{n}{the element to return using a single integer index such that \code{1 < n < NROW(x)}, or a probability \code{0 < n < 1}. See Details. }

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values only where \code{x} is also missing.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{ties}{an integer or character string specifying the method to resolve ties between adjacent qualifying elements:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "mean"   \tab\tab take the arithmetic mean of all qualifying elements. \cr
                 2 \tab\tab "min" \tab\tab take the smallest of the elements. \cr
                 3 \tab\tab "max"   \tab\tab take the largest of the elements. \cr
                 5-9 \tab\tab "qn" \tab\tab continuous quantile types 5-9, see \code{\link{fquantile}}. \cr
                }
}

\item{nthreads}{integer. The number of threads to utilize. Parallelism is across groups for grouped computations on vectors and data frames, and at the column-level otherwise. See Details. }

\item{o}{integer. A valid ordering of \code{x}, e.g. \code{radixorder(x)}. With groups, the grouping needs to be accounted e.g. \code{radixorder(g, x)}.}

\item{check.o}{logical. \code{TRUE} checks that each element of \code{o} is within \code{[1, length(x)]}. The default uses the fact that orderings from \code{\link{radixorder}} have a \code{"sorted"} attribute which let's \code{fnth} infer that the ordering is valid. The length and data type of \code{o} is always checked, regardless of \code{check.o}.}

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain \code{sum} of weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{for \code{fmedian}: further arguments passed to \code{fnth} (apart from \code{n}). If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
For v1.9.0 \code{fnth} was completely rewritten in C and offers significantly enhanced speed and functionality. It uses a combination of quickselect, quicksort, and radixsort algorithms, combined with several (weighted) quantile estimation methods and, where possible, OpenMP multithreading. This synthesis can be summarised as follows:

\itemize{
\item without weights, quickselect is used to determine a (lower) order statistic. If \code{ties \%!in\% c("min", "max")} a second order statistic is found by taking the max of the upper part of the partitioned array, and the two statistics are averaged using a simple mean (\code{ties = "mean"}), or weighted average according to a \code{\link{quantile}} method (\code{ties = "q5"-"q9"}). For \code{n = 0.5}, all supported quantile methods give the sample median. With matrices, multithreading is always across columns, for vectors and data frames it is across groups unless \code{is.null(g)} for data frames.

\item with weights and no groups (\code{is.null(g)}), \code{\link{radixorder}} is called internally (on each column of \code{x}). The ordering is used to sum the weights in order of \code{x} and determine weighted order statistics or quantiles. See details below. Multithreading is disabled as \code{\link{radixorder}} cannot be called concurrently on the same memory stack.

\item with weights and groups (\code{!is.null(g)}), R's quicksort algorithm is used to sort the data in each group and return an index which can be used to sum the weights in order and proceed as before. This is multithreaded across columns for matrices, and across groups otherwise.

\item in \code{fnth.default}, an ordering of \code{x} can be supplied to '\code{o}' e.g. \code{fnth(x, 0.75, o = radixorder(x))}. This dramatically speeds up the estimation both with and without weights, and is useful if \code{fnth} is to be invoked repeatedly on the same data. With groups, \code{o} needs to also account for the grouping e.g. \code{fnth(x, 0.75, g, o = radixorder(g, x))}. Multithreading is possible across groups. See Examples.
}

%This is an R port to \code{std::nth_element}, an efficient partial sorting algorithm in C++. It is also used to calculated the median (in fact the default \code{fnth(x, n = 0.5)} is identical to \code{fmedian(x)}, so see also the details for \code{\link{fmedian}}).

% \code{fnth} generalizes the principles of median value calculation to find arbitrary elements. It offers considerable flexibility by providing both simple order statistics and simple discontinuous quantile estimation. Regarding the former, setting \code{n} to an index between 1 and \code{NROW(x)} will return the n'th smallest element of \code{x}, about 2x faster than \code{sort(x, partial = n)[n]}. As to the latter, setting \code{n} to a probability between 0 and 1 will return the corresponding element of \code{x}, and resolve ties between multiple qualifying elements (such as when \code{n = 0.5} and \code{x} is even) using the arithmetic average \code{ties = "mean"}, or the smallest \code{ties = "min"} or largest \code{ties = "max"} of those elements.

If \code{n > 1}, the result is equivalent to (column-wise) \code{sort(x, partial = n)[n]}. Internally, \code{n} is converted to a probability using \code{p = (n-1)/(NROW(x)-1)}, and that probability is applied to the set of non-missing elements to find the \code{as.integer(p*(fnobs(x)-1))+1L}'th element (which corresponds to option \code{ties = "min"}). % Note that it is necessary to subtract and add 1 so that \code{n = 1} corresponds to \code{p = 0} and \code{n = NROW(x)} to \code{p = 1}. %So if \code{n > 1} is used in the presence of missing values, and the default \code{ties = "mean"} is enabled, the resulting element could be the average of two elements.
When using grouped computations with \code{n > 1}, \code{n} is transformed to a probability \code{p = (n-1)/(NROW(x)/ng-1)} (where \code{ng} contains the number of unique groups in \code{g}).

If weights are used and \code{ties = "q5"-"q9"}, weighted continuous quantile estimation is done as described in \code{\link{fquantile}}.

For \code{ties \%in\% c("mean", "min", "max")}, a target partial sum of weights \code{p*sum(w)} is calculated, and the weighted n'th element is the element k such that all elements smaller than k have a sum of weights \code{<= p*sum(w)}, and all elements larger than k have a sum of weights \code{<= (1 - p)*sum(w)}. If the partial-sum of weights (\code{p*sum(w)}) is reached exactly for some element k, then (summing from the lower end) both k and k+1 would qualify as the weighted n'th element. If the weight of element k+1 is zero, k, k+1 and k+2 would qualify... . If \code{n > 1}, k is chosen (consistent with the unweighted behavior). %(ensuring that \code{fnth(x, n)}) and \code{fnth(x, n, w = rep(1, NROW(x)))}, always provide the same outcome)
If \code{0 < n < 1}, the \code{ties} option regulates how to resolve such conflicts, yielding lower (\code{ties = "min"}: k), upper (\code{ties = "max"}: k+2) or average weighted (\code{ties = "mean"}: mean(k, k+1, k+2)) n'th elements.

Thus, in the presence of zero weights, the weighted median (default \code{ties = "mean"}) can be an arithmetic average of >2 qualifying elements. Users may prefer a quantile based weighted median by setting \code{ties = "q5"-"q9"}, which is a continuous function of \code{p} and ignores elements with zero weights.

For data frames, column-attributes and overall attributes are preserved if \code{g} is used or \code{drop = FALSE}.

}
\value{
The (\code{w} weighted) n'th element/quantile of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) n'th element/quantile.

}
\seealso{
\code{\link{fquantile}}, \code{\link{fmean}}, \code{\link{fmode}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
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
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg, hp) \%>\%    # Faster selecting
      fnth(0.75, hp, "/")  # Divide mpg by its third weighted group-quartile, using hp as weights

# Efficient grouped estimation of multiple quantiles
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\%
    fmutate(o = radixorder(GRPid(), mpg)) \%>\%
    fsummarise(mpg_Q1 = fnth(mpg, 0.25, o = o),
               mpg_median = fmedian(mpg, o = o),
               mpg_Q3 = fnth(mpg, 0.75, o = o))
}

## fmedian()
fmedian(mpg)                         # Simple median value
fmedian(mpg, w = mtcars$hp)          # Weighted median: Weighted by hp
fmedian(mpg, TRA = "-")              # Simple transformation: Subtract median value
fmedian(mpg, mtcars$cyl)             # Grouped median value
fmedian(mpg, mtcars[c(2,8:9)])       # More groups..
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
fmedian(m)
head(fmedian(m, TRA = "-"))
fmedian(m, g) # etc..
\donttest{ % No code relying on suggested package
# grouped_df method
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian()
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmedian(hp)             # Weighted
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fmedian()              # Faster grouping!
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fmedian(TRA = "-")     # De-median
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg, hp) \%>\%   # Faster selecting
      fmedian(hp, "-")  # Weighted de-median mpg, using hp as weights
}
}
\keyword{univar}
\keyword{manip}
