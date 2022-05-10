\name{fmin-fmax}
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

\method{fmax}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, \dots)
\method{fmin}{default}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, \dots)

\method{fmax}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
\method{fmin}{matrix}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)

\method{fmax}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)
\method{fmin}{data.frame}(x, g = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, \dots)

\method{fmax}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, \dots)
\method{fmin}{grouped_df}(x, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}
}
\details{
Missing-value removal as controlled by the \code{na.rm} argument is done at no extra cost since in C++ any logical comparison involving \code{NA} or \code{NaN} evaluates to \code{FALSE}. Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{\link{max}} and \code{\link{min}} which just run through without any checks).

%This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

%When applied to data frames with groups or \code{drop = FALSE}, \code{fmax} and \code{fmin} preserve all column attributes (such as variable labels) but do not distinguish between classed and unclassed objects. The attributes of the data frame itself are also preserved.

For further computational details see \code{\link{fsum}}.

}
\value{
\code{fmax} returns the maximum value of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped) maximum value. Analogous, \code{fmin} returns the minimum value \dots
}
\seealso{
\link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fmax(mpg)                         # Maximum value
fmin(mpg)                         # Minimum value (all examples below use fmax but apply to fmin)
fmax(mpg, TRA = "\%")              # Simple transformation: Take percentage of maximum value
fmax(mpg, mtcars$cyl)             # Grouped maximum value
fmax(mpg, mtcars[c(2,8:9)])       # More groups..
g <- GRP(mtcars, ~ cyl + vs + am) # Precomputing groups gives more speed !
fmax(mpg, g)
fmax(mpg, g, TRA = "\%")           # Groupwise percentage of maximum value
fmax(mpg, g, TRA = "replace")     # Groupwise replace by maximum value

## data.frame method
fmax(mtcars)
head(fmax(mtcars, TRA = "\%"))
fmax(mtcars, g)
fmax(mtcars, g, use.g.names = FALSE) # No row-names generated

## matrix method
m <- qM(mtcars)
fmax(m)
head(fmax(m, TRA = "\%"))
fmax(m, g) # etc..
%\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax()
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmax("\%")
mtcars \%>\% group_by(cyl,vs,am) \%>\% select(mpg) \%>\% fmax()
%}
}
\keyword{univar}
\keyword{manip}
