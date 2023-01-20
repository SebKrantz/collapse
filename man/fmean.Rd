\name{fmean}
\alias{fmean}
\alias{fmean.default}
\alias{fmean.matrix}
\alias{fmean.data.frame}
\alias{fmean.grouped_df}
\title{Fast (Grouped, Weighted) Mean for Matrix-Like Objects}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{fmean} is a generic function that computes the (column-wise) mean of \code{x}, (optionally) grouped by \code{g} and/or weighted by \code{w}.
The \code{\link{TRA}} argument can further be used to transform \code{x} using its (grouped, weighted) mean.
}
\usage{
fmean(x, \dots)

\method{fmean}{default}(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
      use.g.names = TRUE, nthreads = .op[["nthreads"]], \dots)

\method{fmean}{matrix}(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
      use.g.names = TRUE, drop = TRUE, nthreads = .op[["nthreads"]], \dots)

\method{fmean}{data.frame}(x, g = NULL, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
      use.g.names = TRUE, drop = TRUE, nthreads = .op[["nthreads"]], \dots)

\method{fmean}{grouped_df}(x, w = NULL, TRA = NULL, na.rm = .op[["na.rm"]],
      use.g.names = FALSE, keep.group_vars = TRUE,
      keep.w = TRUE, nthreads = .op[["nthreads"]], \dots)
}
\arguments{
\item{x}{a numeric vector, matrix, data frame or grouped data frame (class 'grouped_df').}

\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}

\item{w}{a numeric vector of (non-negative) weights, may contain missing values.}

\item{TRA}{an integer or quoted operator indicating the transformation to perform:
0 - "replace_NA"     |     1 - "replace_fill"     |     2 - "replace"     |     3 - "-"     |     4 - "-+"     |     5 - "/"     |     6 - "\%"     |     7 - "+"     |     8 - "*"     |     9 - "\%\%"     |     10 - "-\%\%". See \code{\link{TRA}}.}


\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost. If \code{na.rm = FALSE} a \code{NA} is returned when encountered.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{nthreads}{integer. The number of threads to utilize. See Details of \code{\link{fsum}}. }

\item{drop}{\emph{matrix and data.frame method:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if \code{g = NULL} and \code{TRA = NULL}.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{keep.w}{\emph{grouped_df method:} Logical. Retain summed weighting variable after computation (if contained in \code{grouped_df}).}

\item{\dots}{arguments to be passed to or from other methods. If \code{TRA} is used, passing \code{set = TRUE} will transform data by reference and return the result invisibly.}

}
\details{
% Non-grouped mean computations internally utilize long-doubles in C++, for additional numeric precision.

% Missing-value removal as controlled by the \code{na.rm} argument is done very efficiently by simply skipping them in the computation (thus setting \code{na.rm = FALSE} on data with no missing values doesn't give extra speed). Large performance gains can nevertheless be achieved in the presence of missing values if \code{na.rm = FALSE}, since then the corresponding computation is terminated once a \code{NA} is encountered and \code{NA} is returned (unlike \code{\link{mean}} which just runs through without any checks).

The weighted mean is computed as \code{sum(x * w) / sum(w)}. If \code{na.rm = TRUE}, missing values will be removed from both \code{x} and \code{w} i.e. utilizing only \code{x[complete.cases(x,w)]} and \code{w[complete.cases(x,w)]}.

%This all seamlessly generalizes to grouped computations, which are performed in a single pass (without splitting the data) and therefore extremely fast.

%When applied to data frames with groups or \code{drop = FALSE}, \code{fmean} preserves all column attributes (such as variable labels) but does not distinguish between classed and unclassed object (thus applying \code{fmean} to a factor column will give a 'malformed factor' error). The attributes of the data frame itself are also preserved.

For further computational details see \code{\link{fsum}}, which works equivalently.

}
\value{
The (\code{w} weighted) mean of \code{x}, grouped by \code{g}, or (if \code{\link{TRA}} is used) \code{x} transformed by its (grouped, weighted) mean.
}
\seealso{
\code{\link{fmedian}}, \code{\link{fmode}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default vector method
mpg <- mtcars$mpg
fmean(mpg)                         # Simple mean
fmean(mpg, w = mtcars$hp)          # Weighted mean: Weighted by hp
fmean(mpg, TRA = "-")              # Simple transformation: demeaning (See also ?W)
fmean(mpg, mtcars$cyl)             # Grouped mean
fmean(mpg, mtcars[8:9])            # another grouped mean.
g <- GRP(mtcars[c(2,8:9)])
fmean(mpg, g)                      # Pre-computing groups speeds up the computation
fmean(mpg, g, mtcars$hp)           # Grouped weighted mean
fmean(mpg, g, TRA = "-")           # Demeaning by group
fmean(mpg, g, mtcars$hp, "-")      # Group-demeaning using weighted group means

## data.frame method
fmean(mtcars)
fmean(mtcars, g)
fmean(fgroup_by(mtcars, cyl, vs, am))  # Another way of doing it..
head(fmean(mtcars, g, TRA = "-"))      # etc..

## matrix method
m <- qM(mtcars)
fmean(m)
fmean(m, g)
head(fmean(m, g, TRA = "-")) # etc..
\donttest{ % No code relying on suggested package
## method for grouped data frames - created with dplyr::group_by or fgroup_by
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean()         # Ordinary
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean(hp)       # Weighted
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean(hp, "-")  # Weighted Transform
mtcars \%>\% group_by(cyl,vs,am) \%>\%
           select(mpg,hp) \%>\% fmean(hp, "-")      # Only mpg

mtcars \%>\% fgroup_by(cyl,vs,am) \%>\%                # Equivalent and faster !
           fselect(mpg,hp) \%>\% fmean(hp, "-")
}
}
\keyword{univar}
\keyword{manip}
