\name{funique}
\alias{funique}
\alias{funique.default}
\alias{funique.data.frame}
\alias{funique.sf}
\alias{funique.pseries}
\alias{funique.pdata.frame}
\alias{fnunique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Unique Elements / Rows
}
\description{
\code{funique} is an efficient alternative to \code{\link{unique}} (or \code{unique.data.table, kit::funique, dplyr::distinct}).

\code{fnunique} is an alternative to \code{NROW(unique(x))} (or \code{data.table::uniqueN, kit::uniqLen, dplyr::n_distinct}).

The \emph{collapse} versions are versatile and highly competitive. % on data frames.
}
\usage{
funique(x, \dots)

\method{funique}{default}(x, sort = FALSE, method = "auto", \dots)

\method{funique}{data.frame}(x, cols = NULL, sort = FALSE, method = "auto", \dots)

\method{funique}{sf}(x, cols = NULL, sort = FALSE, method = "auto", \dots)

# Methods for indexed data / compatibility with plm:

\method{funique}{pseries}(x, sort = FALSE, method = "auto", drop.index.levels = "id", \dots)

\method{funique}{pdata.frame}(x, cols = NULL, sort = FALSE, method = "auto", drop.index.levels = "id", \dots)


fnunique(x) # Fast NROW(unique(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a atomic vector or data frame / list of equal-length columns. }
  \item{sort}{logical. \code{TRUE} orders the unique elements / rows. \code{FALSE} returns unique values in order of first occurrence. }
   \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto"   \tab\tab automatic selection: hash if \code{sort = FALSE} else radix. \cr
                 2 \tab\tab "radix" \tab\tab use radix ordering to determine unique values. Supports \code{sort = FALSE} but only for character data. \cr
                 3 \tab\tab "hash"   \tab\tab use index hashing to determine unique values. Supports \code{sort = TRUE} but only for atomic vectors (default method). \cr
  }
  }
\item{cols}{compute unique rows according to a subset of columns. Columns can be selected using column names, indices, a logical vector or a selector function (e.g. \code{is.character}). \emph{Note:} All columns are returned. }
\item{\dots}{arguments passed to \code{\link{radixorderv}}, e.g. \code{decreasing} or \code{na.last}. Only applicable if \code{method = "radix"}.}
\item{drop.index.levels}{character. Either \code{"id"}, \code{"time"}, \code{"all"} or \code{"none"}. See \link{indexing}.}
}
\details{
If all values/rows are already unique, then \code{x} is returned. Otherwise a copy of \code{x} with duplicate rows removed is returned.  See \code{\link{group}} for some additional computational details.

The \emph{sf} method simply ignores the geometry column when determining unique values.

Methods for indexed data also subset the index accordingly.
}
\note{
No matrix method is provided. Please see \code{kit::funique} for a fast alternative with a matrix method. % The \emph{kit} version is also often faster for vectors.
}
\value{
\code{x} with duplicate elements/rows removed. % sorted in ascending order if \code{sort = TRUE}, and in order of first occurrence if \code{sort = FALSE}.
}
\seealso{
\code{\link{fndistinct}}, \code{\link{group}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}.
}
\examples{
funique(mtcars$cyl)
funique(gv(mtcars, c(2,8,9)))
funique(mtcars, cols = c(2,8,9))
fnunique(gv(mtcars, c(2,8,9)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
