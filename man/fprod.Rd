\name{fprod}
\alias{fprod}
\alias{fprod.default}
\alias{fprod.matrix}
\alias{fprod.data.frame}
\alias{fprod.grouped_df}
\title{Fast (Grouped, Weighted) Product for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fprod} is a generic function that computes the (column-wise) product of all values in \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) product.
}
\usage{
fprod(x, \dots)

\method{fprod}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, \dots)

\method{fprod}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, \dots)

\method{fprod}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = TRUE, drop = TRUE, \dots)

\method{fprod}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
      use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain product of weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}
}
\details{
Non-grouped product computations internally utilize long-doubles in C++, for additional numeric precision.

%Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping them in the computation (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{\link{prod}} which just runs through without any checks).

%This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

The weighted product is computed as \code{prod(x * w)}. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

%When applied to data frames with groups or \code{drop = FALSE}, \code{fprod} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed objects. The attributes of the data frame itself are also preserved.

For further computational details see \code{\link{fsum}}, which works equivalently.

}
\value{
The (\code{w} weighted) product of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) product.
}
\seealso{
\code{\link{fsum}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fprod(mpg)                         # Simple product
fprod(mpg, w = mtcars$hp)          # Weighted product
fprod(mpg, TRA = "/")              # Simple transformation: Divide by product
fprod(mpg, mtcars$cyl)             # Grouped product
fprod(mpg, mtcars$cyl, mtcars$hp)  # Weighted grouped product
fprod(mpg, mtcars[c(2,8:9)])       # More groups..
g <- GRP(mtcars, ~ cyl + vs + am)  # Precomputing groups gives more speed !
fprod(mpg, g)
fprod(mpg, g, TRA = "/")           # Groupwise divide by product

## data.frame method
fprod(mtcars)
head(fprod(mtcars, TRA = "/"))
fprod(mtcars, g)
fprod(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fprod(m)
head(fprod(m, TRA = "/"))
fprod(m, g) # etc..
%\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fprod(hp)   # Weighted grouped product
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fprod(hp)  # Equivalent and faster
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fprod(TRA = "/")
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fselect(mpg) \%>\% fprod()
%}
}
\keyword{univar}
\keyword{manip}
