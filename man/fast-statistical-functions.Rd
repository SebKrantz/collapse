\name{fast-statistical-functions}
\alias{A1-fast-statistical-functions}
\alias{fast-statistical-functions}
\alias{.FAST_STAT_FUN}
\alias{.FAST_FUN}
\title{Fast (Grouped, Weighted) Statistical Functions for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
With \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{fnth}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fnobs}} and \code{\link{fndistinct}}, \emph{collapse} presents a coherent set of extremely fast and flexible statistical functions (S3 generics) to perform column-wise, grouped and weighted computations on vectors, matrices and data frames, with special support for grouped data frames / tibbles (\emph{dplyr}) and \emph{data.table}'s.
}
\section{Usage}{\if{html}{\out{<div class="sourceCode r">}}\preformatted{
## All functions (FUN) follow a common syntax in 4 methods:
FUN(x, ...)

## Default S3 method:
FUN(x, g = NULL, [w = NULL,] TRA = NULL, [na.rm = TRUE,]
    use.g.names = TRUE, [nthreads = 1L,] ...)

## S3 method for class 'matrix'
FUN(x, g = NULL, [w = NULL,] TRA = NULL, [na.rm = TRUE,]
    use.g.names = TRUE, drop = TRUE, [nthreads = 1L,] ...)

## S3 method for class 'data.frame'
FUN(x, g = NULL, [w = NULL,] TRA = NULL, [na.rm = TRUE,]
    use.g.names = TRUE, drop = TRUE, [nthreads = 1L,] ...)

## S3 method for class 'grouped_df'
FUN(x, [w = NULL,] TRA = NULL, [na.rm = TRUE,]
    use.g.names = FALSE, keep.group_vars = TRUE,
    [keep.w = TRUE,] [stub = TRUE,] [nthreads = 1L,] ...)
}\if{html}{\out{</div>}}
}
\section{Arguments}{
\tabular{lll}{
\code{x} \tab \tab a vector, matrix, data frame or grouped data frame (class 'grouped_df'). \cr

\code{g} \tab \tab a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}. \cr

\code{w} \tab \tab a numeric vector of (non-negative) weights, may contain missing values. Supported by \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fnth}}, \code{\link{fvar}}, \code{\link{fsd}} and \code{\link{fmode}}. \cr

\code{TRA} \tab \tab an integer or quoted operator indicating the transformation to perform:
0 - "na"     |     1 - "fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}. \cr

\code{na.rm} \tab \tab logical. Skip missing values in \code{x}. Defaults to \code{TRUE} in all functions and implemented at very little computational cost. Not available for \code{\link{fnobs}}. \cr

\code{use.g.names} \tab \tab logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s. \cr

\code{nthreads} \tab \tab integer. The number of threads to utilize. Supported by \code{\link{fsum}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fnth}}, \code{\link{fmode}} and \code{\link{fndistinct}}. \cr

\code{drop} \tab \tab \emph{matrix and data.frame methods:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}. \cr

\code{keep.group_vars}  \tab \tab \emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation. By default grouping variables are added, even if not present in the grouped_df. \cr

\code{keep.w} \tab \tab \emph{grouped_df method:} Logical. \code{TRUE} (default) also aggregates weights and saves them in a column, \code{FALSE} removes weighting variable after computation (if contained in \code{grouped_df}). \cr

\code{stub} \tab \tab \emph{grouped_df method:} Character. If \code{keep.w = TRUE} and \code{stub = TRUE} (default), the aggregated weights column is prefixed by the name of the aggregation function (mostly \code{"sum."}). Users can specify a different prefix through this argument, or set it to \code{FALSE} to avoid prefixing.\cr

\code{\dots} \tab \tab arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly (except for the grouped_df method which always returns visible output). \cr
}
}
\section{Details}{
Please see the documentation of individual functions.
}
\section{Value}{
 \code{x} suitably aggregated or transformed. Data frame column-attributes and overall attributes are generally preserved if the output is of the same data type.
}
\section{Related Functionality}{
\itemize{
\item Functions \code{\link{fquantile}} and \code{\link{frange}} are for atomic vectors.

\item Panel-decomposed (i.e. between and within) statistics as well as grouped and weighted skewness and kurtosis are implemented in \code{\link{qsu}}.

\item The vector-valued functions and operators \code{\link{fcumsum}}, \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fhdbetween]{fhdbetween/HDB}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=fhdwithin]{fhdwithin/HDW}}, \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}} and \code{\link[=fgrowth]{fgrowth/G}} are grouped under \link[=data-transformations]{Data Transformations} and \link[=time-series-panel-series]{Time Series and Panel Series}. These functions also support \link[=indexing]{indexed data} (\emph{plm}).
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=data-transformations]{Data Transformations}, \link[=time-series-panel-series]{Time Series and Panel Series}
}
\section{Examples}{\if{html}{\out{<div class="sourceCode r">}}\preformatted{
## default vector method
mpg <- mtcars$mpg
fsum(mpg)                         # Simple sum
fsum(mpg, TRA = "/")              # Simple transformation: divide all values by the sum
fsum(mpg, mtcars$cyl)             # Grouped sum
fmean(mpg, mtcars$cyl)            # Grouped mean
fmean(mpg, w = mtcars$hp)         # Weighted mean, weighted by hp
fmean(mpg, mtcars$cyl, mtcars$hp) # Grouped mean, weighted by hp
fsum(mpg, mtcars$cyl, TRA = "/")  # Proportions / division by group sums
fmean(mpg, mtcars$cyl, mtcars$hp, # Subtract weighted group means, see also ?fwithin
      TRA = "-")

## data.frame method
fsum(mtcars)
fsum(mtcars, TRA = "\%")                  # This computes percentages
fsum(mtcars, mtcars[c(2,8:9)])           # Grouped column sum
g <- GRP(mtcars, ~ cyl + vs + am)        # Here precomputing the groups!
fsum(mtcars, g)                          # Faster !!
fmean(mtcars, g, mtcars$hp)
fmean(mtcars, g, mtcars$hp, "-")         # Demeaning by weighted group means..
fmean(fgroup_by(mtcars, cyl, vs, am), hp, "-")  # Another way of doing it..


fmode(wlddev, drop = FALSE)              # Compute statistical modes of variables in this data
fmode(wlddev, wlddev$income)             # Grouped statistical modes ..

## matrix method
m <- qM(mtcars)
fsum(m)
fsum(m, g) # ..

## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars |> group_by(cyl,vs,am) |> select(mpg,carb) |> fsum()
mtcars |> fgroup_by(cyl,vs,am) |> fselect(mpg,carb) |> fsum() # equivalent and faster !!
mtcars |> fgroup_by(cyl,vs,am) |> fsum(TRA = "\%")
mtcars |> fgroup_by(cyl,vs,am) |> fmean(hp)         # weighted grouped mean, save sum of weights
mtcars |> fgroup_by(cyl,vs,am) |> fmean(hp, keep.group_vars = FALSE)
}\if{html}{\out{</div>}}
}
\section{Benchmark}{\if{html}{\out{<div class="sourceCode r">}}\preformatted{
## This compares fsum with data.table (2 threads) and base::rowsum
# Starting with small data
mtcDT <- qDT(mtcars)
f <- qF(mtcars$cyl)

library(microbenchmark)
microbenchmark(mtcDT[, lapply(.SD, sum), by = f],
               rowsum(mtcDT, f, reorder = FALSE),
               fsum(mtcDT, f, na.rm = FALSE), unit = "relative")

#                              expr        min         lq      mean    median        uq       max neval cld
# mtcDT[, lapply(.SD, sum), by = f] 145.436928 123.542134 88.681111 98.336378 71.880479 85.217726   100   c
# rowsum(mtcDT, f, reorder = FALSE)   2.833333   2.798203  2.489064  2.937889  2.425724  2.181173   100  b
#     fsum(mtcDT, f, na.rm = FALSE)   1.000000   1.000000  1.000000  1.000000  1.000000  1.000000   100 a

# Now larger data
tdata <- qDT(replicate(100, rnorm(1e5), simplify = FALSE)) # 100 columns with 100.000 obs
f <- qF(sample.int(1e4, 1e5, TRUE))                        # A factor with 10.000 groups

microbenchmark(tdata[, lapply(.SD, sum), by = f],
               rowsum(tdata, f, reorder = FALSE),
               fsum(tdata, f, na.rm = FALSE), unit = "relative")

#                              expr      min       lq     mean   median       uq       max neval cld
# tdata[, lapply(.SD, sum), by = f] 2.646992 2.975489 2.834771 3.081313 3.120070 1.2766475   100   c
# rowsum(tdata, f, reorder = FALSE) 1.747567 1.753313 1.629036 1.758043 1.839348 0.2720937   100  b
#     fsum(tdata, f, na.rm = FALSE) 1.000000 1.000000 1.000000 1.000000 1.000000 1.0000000   100 a
}\if{html}{\out{</div>}}
}

\keyword{univar}
\keyword{manip}
\keyword{documentation}
