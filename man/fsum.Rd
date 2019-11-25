\name{fsum}
\alias{fsum}
\alias{fsum.default}
\alias{fsum.matrix}
\alias{fsum.data.frame}
\alias{fsum.grouped_df}
\title{Fast (Grouped) Sum for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fsum} is a generic function that computes the (column-wise) sum of all values in \code{x}, (optionally) grouped by \code{g}. The \code{TRA} argument can further be used to transform \code{x} using its (grouped) sum.
}
\usage{
\method{fsum}{default}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
     use.g.names = TRUE, ...)

\method{fsum}{matrix}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, ...)

\method{fsum}{data.frame}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, ...)

\method{fsum}{grouped_df}(x, TRA = FALSE, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, ...)
}
\arguments{
\item{x}{a numeric vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = FALSE}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{...}{arguments to be passed to or from other methods.}

}
\details{
Non-grouped sum computations internally utilize long-doubles in C++, for additional numeric precision.

Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping them in the computation (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{base::sum} which just runs through without any checks).

This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast. See Benchmark and Examples below.

When applied to data frame's with groups or \code{drop = FALSE}, \code{fsum} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed objects. The attributes of the \code{data.frame} itself are also perserved.

}
\value{
The sum of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its sum, grouped by \code{g}.
}
\section{See Also}{
\link{Fast Statistical Functions}
}
\section{Examples}{\preformatted{
## default vector method
mpg <- mtcars$mpg
fsum(mpg)                         # Simple sum
fsum(mpg, TRA = "\%")              # Simple transformation: obtain percentages of mpg
fsum(mpg, mtcars$cyl)             # Grouped sum
fsum(mpg, mtcars[c(2,8:9)])       # More groups...
g <- GRP(mtcars, ~ cyl + vs + am) # Precomputing groups gives more speed !!
fsum(mpg, g)
fmean(mpg, g) == fsum(mpg, g) / fNobs(mpg, g)
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

## method for grouped tibbles - for use with dplyr
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsum
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsum("\%")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fsum
}}
\section{Benchmark}{\preformatted{
## Let's run some benchmarks and compare fsum against data.table and base::rowsum
# Starting with small data
mtcDT <- qDT(mtcars)
f <- qF(mtcars$cyl)

library(microbenchmark)
microbenchmark(mtcDT[, lapply(.SD, sum), by = f],
               rowsum(mtcDT, f, reorder = FALSE),
               fsum(mtcDT, f, na.rm = FALSE), unit = "relative")
# My results:
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
# My results:
                              expr      min       lq     mean   median       uq       max neval cld
 tdata[, lapply(.SD, sum), by = f] 2.646992 2.975489 2.834771 3.081313 3.120070 1.2766475   100   c
 rowsum(tdata, f, reorder = FALSE) 1.747567 1.753313 1.629036 1.758043 1.839348 0.2720937   100  b
     fsum(tdata, f, na.rm = FALSE) 1.000000 1.000000 1.000000 1.000000 1.000000 1.0000000   100 a
}}

