\name{fsum}
\alias{fsum}
\alias{fsum.default}
\alias{fsum.matrix}
\alias{fsum.data.frame}
\alias{fsum.grouped_df}
\title{Fast (Grouped, Weighted) Sum for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fsum} is a generic function that computes the (column-wise) sum of all values in \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w} (e.g. to calculate survey totals). The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) sum.
}
\usage{
fsum(x, \dots)

\method{fsum}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, fill = FALSE, nthreads = 1L, \dots)

\method{fsum}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, fill = FALSE, nthreads = 1L, \dots)

\method{fsum}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, fill = FALSE, nthreads = 1L, \dots)

\method{fsum}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE,
     keep.w = TRUE, fill = FALSE, nthreads = 1L, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{fill}{logical. Initialize result with \code{0} instead of \code{NA} when \code{na.rm = TRUE} e.g. \code{fsum(NA, fill = TRUE)} returns \code{0} instead of \code{NA}. }

\item{nthreads}{integer. The number of threads to utilize. See Details. }

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain summed weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
% Non-grouped sum computations internally utilize long-doubles in C++, for additional numeric precision.

Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping them in the computation (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{\link{sum}} which just runs through without any checks).

The weighted sum (e.g. survey total) is computed as \code{sum(x * w)}, but in one pass and about twice as efficient. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and are therefore extremely fast. See Benchmark and Examples below.

When applied to data frames with groups or \code{drop = FALSE}, \code{fsum} preserves all column attributes. The attributes of the data frame itself are also preserved.

Since v1.6.0 \code{fsum} explicitly supports integers. Integers are summed using the long long type in C which is bounded at +-9,223,372,036,854,775,807 (so ~4.3 billion times greater than the minimum/maximum R integer bounded at +-2,147,483,647). If the value of the sum is outside +-2,147,483,647, a double containing the result is returned, otherwise an integer is returned. With groups, an integer overflow error is provided if the sum in any group is outside +-2,147,483,647. Data should be coerced to double beforehand in such cases.

Multithreading, added in v1.8.0, applies at the column-level unless \code{g = NULL} and \code{nthreads > NCOL(x)}. Parallelism over groups is not available because sums are computed simultaneously within each group. \code{nthreads = 1L} uses a serial version of the code, not parallel code running on one thread. This serial code is always used with less than 100,000 obs (\code{length(x) < 100000} for vectors and matrices), because parallel execution itself has some overhead.

}
\value{
The (\code{w} weighted) sum of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) sum.
}
\section{See Also}{
\code{\link{fprod}}, \code{\link{fmean}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\section{Examples}{\preformatted{
## default vector method
mpg <- mtcars$mpg
fsum(mpg)                         # Simple sum
fsum(mpg, w = mtcars$hp)          # Weighted sum (total): Weighted by hp
fsum(mpg, TRA = "\%")              # Simple transformation: obtain percentages of mpg
fsum(mpg, mtcars$cyl)             # Grouped sum
fsum(mpg, mtcars$cyl, mtcars$hp)  # Weighted grouped sum (total)
fsum(mpg, mtcars[c(2,8:9)])       # More groups..
g <- GRP(mtcars, ~ cyl + vs + am) # Precomputing groups gives more speed !
fsum(mpg, g)
fmean(mpg, g) == fsum(mpg, g) / fnobs(mpg, g)
fsum(mpg, g, TRA = "\%")           # Percentages by group

## data.frame method
fsum(mtcars)
fsum(mtcars, TRA = "\%")
fsum(mtcars, g)
fsum(mtcars, g, TRA = "\%")

## matrix method
m <- qM(mtcars)
fsum(m)
fsum(m, TRA = "\%")
fsum(m, g)
fsum(m, g, TRA = "\%")
% \donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsum(hp)  # Weighted grouped sum (total)
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fsum(hp) # Equivalent and faster !!
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fsum(TRA = "\%")
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg) \%>\% fsum()
% }
}
}
\section{Benchmark}{\preformatted{
## This compares fsum with data.table (2 threads) and base::rowsum
# Starting with small data
mtcDT <- qDT(mtcars)
f <- qF(mtcars$cyl)

library(microbenchmark)
microbenchmark(mtcDT[, lapply(.SD, sum), by = f],
               rowsum(mtcDT, f, reorder = FALSE),
               fsum(mtcDT, f, na.rm = FALSE), unit = "relative")

                              expr        min         lq      mean    median        uq       max neval cld
 mtcDT[, lapply(.SD, sum), by = f] 145.436928 123.542134 88.681111 98.336378 71.880479 85.217726   100   c
 rowsum(mtcDT, f, reorder = FALSE)   2.833333   2.798203  2.489064  2.937889  2.425724  2.181173   100  b
     fsum(mtcDT, f, na.rm = FALSE)   1.000000   1.000000  1.000000  1.000000  1.000000  1.000000   100 a

# Now larger data
tdata <- qDT(replicate(100, rnorm(1e5), simplify = FALSE)) # 100 columns with 100.000 obs
f <- qF(sample.int(1e4, 1e5, TRUE))                        # A factor with 10.000 groups

microbenchmark(tdata[, lapply(.SD, sum), by = f],
               rowsum(tdata, f, reorder = FALSE),
               fsum(tdata, f, na.rm = FALSE), unit = "relative")

                              expr      min       lq     mean   median       uq       max neval cld
 tdata[, lapply(.SD, sum), by = f] 2.646992 2.975489 2.834771 3.081313 3.120070 1.2766475   100   c
 rowsum(tdata, f, reorder = FALSE) 1.747567 1.753313 1.629036 1.758043 1.839348 0.2720937   100  b
     fsum(tdata, f, na.rm = FALSE) 1.000000 1.000000 1.000000 1.000000 1.000000 1.0000000   100 a
}}
\keyword{univar}
\keyword{manip}
