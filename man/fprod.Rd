\name{fprod}
\alias{fprod}
\alias{fprod.default}
\alias{fprod.matrix}
\alias{fprod.data.frame}
\alias{fprod.grouped_df}
\title{Fast (Grouped) Product for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fprod} is a generic function that computes the (column-wise) product of all values in \code{x}, (optionally) grouped by \code{g}. The \code{TRA} argument can further be used to transform \code{x} using its (grouped) product.
}
\usage{
\method{fprod}{default}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
      use.g.names = TRUE, ...)

\method{fprod}{matrix}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
      drop = TRUE, use.g.names = TRUE, ...)

\method{fprod}{data.frame}(x, g = NULL, TRA = FALSE, na.rm = TRUE,
      drop = TRUE, use.g.names = TRUE, ...)

\method{fprod}{grouped_df}(x, TRA = FALSE, na.rm = TRUE,
      drop.groups = FALSE, ...)
}
\arguments{
\item{x}{a numeric vector, matrix, data.frame or grouped tibble (\code{dplyr::grouped_df}).}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*". See \code{\link{TRA}}.}

\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to TRUE and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{drop}{\emph{matrix and data.frame method:} drop dimensions and return an atomic vector, if \code{g = NULL} and \code{TRA = FALSE}.}

\item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}

\item{drop.groups}{\emph{grouped_df method:} Logical. Remove grouping variables after computation.}
}
\details{
Non-grouped product computations internally utilize long-doubles in C++, for additional numeric precision.

Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping them in the computation (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{base::prod} which just runs through without any checks).

This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

When applied to data frame's with groups or \code{drop = FALSE}, \code{fprod} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed objects. The attributes of the \code{data.frame} itself are also perserved.

}
\value{
The product of \code{x}, grouped by \code{g}, or (if \code{TRA} is used) \code{x} transformed by its product, grouped by \code{g}.
}
\section{Examples}{\preformatted{
## default vector method
mpg <- mtcars$mpg
fprod(mpg)                         # Simple product
fprod(mpg, TRA = "/")              # Simple transformation: Divide by product
fprod(mpg, mtcars$cyl)             # Grouped product
fprod(mpg, mtcars[c(2,8:9)])       # More groups...
g <- GRP(mtcars, ~ cyl + vs + am)  # Precomputing groups gives more speed !!
fprod(mpg, g)
fprod(mpg, g, TRA = "/")           # Groupwise divide by product

## data.frame Method
fprod(mtcars)
fprod(mtcars, TRA = "/")
fprod(mtcars, g)
fprod(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fprod(m)
fprod(m, TRA = "/")
fprod(m, g) # etc...

## method for grouped tibbles - for use with dplyr:
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fprod
mtcars \%>\% group_by(cyl,vs,am) \%>\% fprod("/")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fprod
}}

