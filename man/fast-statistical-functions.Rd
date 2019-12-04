\name{fast-statistical-functions}
\alias{fast-statistical-functions}
\alias{.FAST_STAT_FUN}
\title{Fast (Grouped, Weighted) Statistical Functions for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
With \code{\link{fsum}}, \code{\link{fprod}}, \code{\link{fmean}}, \code{\link{fmedian}}, \code{\link{fmode}}, \code{\link{fvar}}, \code{\link{fsd}}, \code{\link{fmin}}, \code{\link{fmax}}, \code{\link{ffirst}}, \code{\link{flast}}, \code{\link{fNobs}} and \code{\link{fNdistinct}}, \code{collapse} presents a coherent set of extremely fast and flexible statistcal functions (S3 generics) to perform column-wise, grouped and weighted computuations on atomic vectors, matrices and data.frames, with special support for \code{dplyr} grouped tibbles and data.table's.

(\emph{Note}: The vector-valued functions and operators \code{\link[=fscale]{fscale/STD}}, \code{\link[=fbetween]{fbetween/B}}, \code{\link[=fHDbetween]{fHDbetween/HDB}}, \code{\link[=fwithin]{fwithin/W}}, \code{\link[=fHDwithin]{fHDwithin/HDW}}, \code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D}} and \code{\link[=fgrowth]{fgrowth/G}} are documented under \link[=data-transformations]{Data Transformations} and \link[=time-series-panel-series]{Time-Series and Panel-Series}. These functions also support \code{plm::pseries} and \code{plm::pdata.frame}'s.)

}
\section{Usage}{\preformatted{
## All functions (FUN) follow a common syntax in 4 methods:
FUN(x, ...)

## Default S3 method:
FUN(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,]
    use.g.names = TRUE, ...)

## S3 method for class 'matrix'
FUN(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,]
    use.g.names = TRUE, drop = TRUE, ...)

## S3 method for class 'data.frame'
FUN(x, g = NULL, [w = NULL,] TRA = FALSE, [na.rm = TRUE,]
    use.g.names = TRUE, drop = TRUE, ...)

## S3 method for class 'grouped_df'
FUN(x, [w = NULL,] TRA = FALSE, [na.rm = TRUE,]
    use.g.names = FALSE, keep.group_vars = TRUE, [keep.w = TRUE,] ...)
}
}
\section{Arguments}{
\tabular{lll}{
\code{x} \tab \tab a vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}). \cr \cr \cr

\code{g} \tab\tab a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}. \cr \cr \cr

\code{w} \tab\tab a numeric vector of (non-negative) weights, may contain missing values. Supported by \code{\link{fmean}}, \code{\link{fvar}}, \code{\link{fsd}} and \code{\link{fmode}}. \cr \cr \cr

\code{TRA} \tab\tab an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}. \cr \cr \cr

\code{na.rm} \tab\tab logical. Skip missing values in \code{x}. Defaults to \code{TRUE} in all functions and implemented at very little computational cost. Not available for \code{\link{fNobs}}. \cr \cr \cr

\code{use.g.names} \tab\tab make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles. \cr \cr \cr

\code{drop} \tab\tab \emph{matrix and data.frame methods:} Drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = FALSE}. \cr \cr \cr

\code{keep.group_vars}  \tab\tab \emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation. \cr \cr \cr

\code{keep.w} \tab\tab \emph{grouped_df method:} Logical. \code{TRUE} also aggregates weights and saves them in a column, \code{FALSE} removes weighting variable after computation (if contained in \code{grouped_df}). \cr \cr \cr

\code{...} \tab\tab arguments to be passed to or from other methods, and extra arguments to some functions, i.e. the algorithm used to compute variances etc. \cr \cr \cr
}
}
\details{
Please see the documentation of individual functions.
}
\value{
 \code{x} aggregated. With data.frame's column-attributes and overall attributes are preserved.
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=data-transformations]{Data Transformations}, \link[=time-series-panel-series]{Time-Series and Panel-Series}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fsum(mpg)                         # Simple sum
fsum(mpg, TRA = "\%")              # Simple transformation: obtain percentages of mpg
fsum(mpg, mtcars$cyl)             # Grouped sum
fmean(mpg, mtcars$cyl)            # Grouped mean
fmean(mpg, w = mtcars$hp)         # Weighted mean, weighted by hp
fmean(mpg, mtcars$cyl, mtcars$hp) # Grouped mean, weighted by hp
fsum(mpg, mtcars$cyl, TRA = "\%")  # Percentages by group
fmean(mpg, mtcars$cyl, mtcars$hp, # Replace vector elements with their weighted group-mean
      TRA = "replace")

## data.frame method
fsum(mtcars)
fsum(mtcars, TRA = "\%")
fsum(mtcars, mtcars[c(2,8:9)])           # Grouped column sum
g <- GRP(mtcars, ~ cyl + vs + am)        # Here precomputing the groups!
fsum(mtcars, g)                          # Faster !!
fmean(mtcars, g, mtcars$hp)
fmean(mtcars, g, mtcars$hp, "-")         # demeaning by weighted group means... see also ?W

fmode(wlddev, drop = FALSE)              # Compute statistical modes of variables in this data
fmode(wlddev, wlddev$income)             # grouped statistical modes ..

## matrix method
m <- qM(mtcars)
fsum(m)
fsum(m, g) # ...

## method for grouped tibbles - for use with dplyr
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg,carb) \%>\% fsum
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsum("\%")
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean(hp)         # weighted grouped mean, save sum of weights
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean(hp, keep.group_vars = FALSE)
}

