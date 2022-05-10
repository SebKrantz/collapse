\name{fvar-fsd}
\alias{fvar}
\alias{fvar.default}
\alias{fvar.matrix}
\alias{fvar.data.frame}
\alias{fvar.grouped_df}
\alias{fsd}
\alias{fsd.default}
\alias{fsd.matrix}
\alias{fsd.data.frame}
\alias{fsd.grouped_df}
\title{Fast (Grouped, Weighted) Variance and Standard Deviation for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns
\description{
\code{fvar} and \code{fsd} are generic functions that compute the (column-wise) variance and standard deviation of \code{x}, (optionally) grouped by \code{g} and/or frequency-weighted by \code{w}. The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) variance/sd.
}
\usage{
fvar(x, \dots)
fsd(x, \dots)

\method{fvar}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, stable.algo = TRUE, \dots)
\method{fsd}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
    use.g.names = TRUE, stable.algo = TRUE, \dots)

\method{fvar}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, \dots)
\method{fsd}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
    use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, \dots)

\method{fvar}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, \dots)
\method{fsd}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = TRUE,
    use.g.names = TRUE, drop = TRUE, stable.algo = TRUE, \dots)

\method{fvar}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
     use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
     stable.algo = TRUE, \dots)
\method{fsd}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = TRUE,
    use.g.names = FALSE, keep.group_vars = TRUE, keep.w = TRUE,
    stable.algo = TRUE, \dots)
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

\item{keep.w}{\emph{grouped_df method:} Logical. Retain summed weighting variable after computation (if contained in \code{grouped_df}).}

\item{stable.algo}{logical. \code{TRUE} (default) use Welford's numerically stable online algorithm. \code{FALSE} implements a faster but numerically unstable one-pass method. See Details. }

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
\emph{Welford's online algorithm} used by default to compute the variance is well described \href{https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance}{here} (the section \emph{Weighted incremental algorithm} also shows how the weighted variance is obtained by this algorithm).

If \code{stable.algo = FALSE}, the variance is computed in one-pass as \code{(sum(x^2)-n*mean(x)^2)/(n-1)}, where \code{sum(x^2)} is the sum of squares from which the expected sum of squares \code{n*mean(x)^2} is subtracted, normalized by \code{n-1} (Bessel's correction). This is numerically unstable if \code{sum(x^2)} and \code{n*mean(x)^2} are large numbers very close together, which will be the case for large \code{n}, large \code{x}-values and small variances (catastrophic cancellation occurs, leading to a loss of numeric precision). Numeric precision is however still maximized through the internal use of long doubles in C++, and the fast algorithm can be up to 4-times faster compared to Welford's method.

The weighted variance is computed with frequency weights as \code{(sum(x^2*w)-sum(w)*weighted.mean(x,w)^2)/(sum(w)-1)}. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

%Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping the values (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned.

%This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

%When applied to data frames with groups or \code{drop = FALSE}, \code{fvar/fsd} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed object (thus applying \code{fvar/fsd} to a factor column will give a 'malformed factor' error). The attributes of the data frame itself are also preserved.
For further computational detail see \code{\link{fsum}}.
}
\value{
\code{fvar} returns the (\code{w} weighted) variance of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) variance. \code{fsd} computes the standard deviation of \code{x} in like manor.
}
\references{
Welford, B. P. (1962). Note on a method for calculating corrected sums of squares and products. \emph{Technometrics}. 4 (3): 419-420. doi:10.2307/1266577.
}
\seealso{
\link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
fvar(mtcars$mpg)                            # Simple variance (all examples also hold for fvar!)
fsd(mtcars$mpg)                             # Simple standard deviation
fsd(mtcars$mpg, w = mtcars$hp)              # Weighted sd: Weighted by hp
fsd(mtcars$mpg, TRA = "/")                  # Simple transformation: scaling (See also ?fscale)
fsd(mtcars$mpg, mtcars$cyl)                 # Grouped sd
fsd(mtcars$mpg, mtcars$cyl, mtcars$hp)      # Grouped weighted sd
fsd(mtcars$mpg, mtcars$cyl, TRA = "/")      # Scaling by group
fsd(mtcars$mpg, mtcars$cyl, mtcars$hp, "/") # Group-scaling using weighted group sds

## data.frame method
fsd(iris)                           # This works, although 'Species' is a factor variable
fsd(mtcars, drop = FALSE)           # This works, all columns are numeric variables
fsd(iris[-5], iris[5])              # By Species: iris[5] is still a list, and thus passed to GRP()
fsd(iris[-5], iris[[5]])            # Same thing much faster: fsd recognizes 'Species' is a factor
head(fsd(iris[-5], iris[[5]], TRA = "/")) # Data scaled by species (see also fscale)

## matrix method
m <- qM(mtcars)
fsd(m)
fsd(m, mtcars$cyl) # etc..
%\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsd()
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsd(keep.group_vars = FALSE) # Remove grouping columns
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsd(hp)      # Weighted by hp
mtcars \%>\% group_by(cyl,vs,am) \%>\% fsd(hp, "/") # Weighted scaling transformation
%}
}
\keyword{univar}
\keyword{manip}
