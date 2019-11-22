\name{Maxima & Minima}
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
\code{fmax} and \code{fmin} are generic functions that compute the (column-wise) maximum and minimum value of all values in \code{x}, (optionally) grouped by \code{g}. The \code{TRA} argument can further be used to transform \code{x} using its (grouped) maximum or minimum value.
}
\usage{
## Default S3 method:
fmax(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...)
fmin(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, ...)

## S3 method for class 'matrix'
fmax(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
fmin(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)

## S3 method for class 'data.frame'
fmax(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)
fmin(x, g = NULL, TRA = FALSE, na.rm = TRUE, use.g.names = TRUE, drop = TRUE, ...)

## S3 method for class 'grouped_df'
fmax(x, TRA = FALSE, na.rm = TRUE, keep.group_vars = TRUE, ...)
fmin(x, TRA = FALSE, na.rm = TRUE, keep.group_vars = TRUE, ...)
}
\arguments{
\item{x}{a numeric vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to TRUE and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector, if \code{g = NULL} and \code{TRA = FALSE}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}
}
\details{
Missing-value removal as controlled by the \code{na.rm} argument is done at no extra cost since in C++ any logical comarison involving \code{NA} or \code{NaN} evaluates to \code{FALSE}. Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{base::max} and \code{base::min} which just run through without any checks).

This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

When applied to data frame's with groups or \code{drop = FALSE}, \code{fmax} and \code{fmin} preserve all column attributes (such as variable labels) but do not distinguish between classed and unclassed objects. The attributes of the \code{data.frame} itself are also perserved.
}
\value{
\code{fmax} returns the maximum value of \code{x}, grouped by \code{g}, or (if \code{TRA} is used) \code{x} transformed by its maximum value, grouped by \code{g}. Analogous, \code{fmin} returns the minimum value ...
}
\section{Examples}{\preformatted{
## default vector method
mpg <- mtcars$mpg
fmax(mpg)                         # maximum value
fmin(mpg)                         # minimum value (all further examples use fmax but also apply to fmin)
fmax(mpg, TRA = "\%")              # Simple transformation: Take percentage of maximum value
fmax(mpg, mtcars$cyl)             # Grouped maximum value
fmax(mpg, mtcars[c(2,8:9)])       # More groups...
g <- GRP(mtcars, ~ cyl + vs + am) # Precomputing groups gives more speed !!
fmax(mpg, g)
fmax(mpg, g, TRA = "\%")           # Groupwise percentage of maximum value
fmax(mpg, g, TRA = "replace")     # Groupwise replace by maximum value

## data.frame Method
fmax(mtcars)
fmax(mtcars, TRA = "\%")
fmax(mtcars, g)
fmax(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fmax(m)
fmax(m, TRA = "\%")
fmax(m, g) # etc...

## method for grouped tibbles - for use with dplyr:
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax("\%")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fmax
}}

