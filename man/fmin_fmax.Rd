\name{fmax-fmin}
\alias{fmax}
\alias{fmax.default}
\alias{fmax.matrix}
\alias{fmax.data.frame}
\alias{fmax.grouped_df}
\alias{fmin}
\alias{fmin.default}
\alias{fmin.matrix}
\alias{fmin.data.frame}
\alias{fmin.grouped_df}
\title{Fast (Grouped) Maxima and Minima for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmax} and \code{fmin} are generic functions that compute the (column-wise) maximum and minimum value of all values in \code{x}, (optionally) grouped by \code{g}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped) maximum or minimum value.
}
\usage{
fmax(x, \dots)
fmin(x, \dots)
%## Default S3 methods:
\method{fmax}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, \dots)
\method{fmin}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, \dots)
%## S3 methods for class 'matrix'
\method{fmax}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
\method{fmin}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
%## S3 methods for class 'data.frame'
\method{fmax}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
\method{fmin}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
%## S3 methods for class 'grouped_df'
\method{fmax}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, \dots)
\method{fmin}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
Missing-value removal as controlled by the \code{na.rm} argument is done at no extra cost since in C++ any logical comparison involving \code{NA} or \code{NaN} evaluates to \code{FALSE}. Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{base::max} and \code{base::min} which just run through without any checks).

This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

When applied to data frame's with groups or \code{drop = FALSE}, \code{fmax} and \code{fmin} preserve all column attributes (such as variable labels) but do not distinguish between classed and unclassed objects. The attributes of the \code{data.frame} itself are also perserved.
}
\value{
\code{fmax} returns the maximum value of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its maximum value, grouped by \code{g}. Analogous, \code{fmin} returns the minimum value ...
}
\seealso{
\link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fmax(mpg)                         # maximum value
fmin(mpg)                         # minimum value (all examples below use fmax but apply to fmin)
fmax(mpg, TRA = "\%")              # Simple transformation: Take percentage of maximum value
fmax(mpg, mtcars$cyl)             # Grouped maximum value
fmax(mpg, mtcars[c(2,8:9)])       # More groups...
g <- GRP(mtcars, ~ cyl + vs + am) # Precomputing groups gives more speed !!
fmax(mpg, g)
fmax(mpg, g, TRA = "\%")           # Groupwise percentage of maximum value
fmax(mpg, g, TRA = "replace")     # Groupwise replace by maximum value

## data.frame method
fmax(mtcars)
fmax(mtcars, TRA = "\%")
fmax(mtcars, g)
fmax(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fmax(m)
fmax(m, TRA = "\%")
fmax(m, g) # etc...

## method for grouped tibbles - for use with dplyr
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax("\%")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fmax
}

