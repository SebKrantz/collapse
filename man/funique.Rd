\name{funique}
\alias{funique}
\alias{funique.default}
\alias{funique.data.frame}
\alias{funique.sf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Unique Elements / Rows
}
\description{
\code{funique} is a substantially faster alternative to \code{\link{unique}}. It is generic with a default vector method and a data frame method.
}
\usage{
funique(x, \dots)

\method{funique}{default}(x, sort = FALSE, method = c("auto", "radix", "hash"), \dots)

\method{funique}{data.frame}(x, cols = NULL, sort = FALSE, \dots)

\method{funique}{sf}(x, cols = NULL, sort = FALSE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a atomic vector or data frame / list of equal-length columns. }
  \item{sort}{logical. \code{TRUE} orders the unique elements / rows. \code{FALSE} returns unique values in order of first occurrence, however numeric values are only provided in order of first occurrence when \code{method = "hash"} i.e. the radix method always sorts numeric data. }
   \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto"   \tab\tab automatic selection: hash for character, logical or if \code{length(x) < 500}, else radix. \cr
                 2 \tab\tab "radix" \tab\tab use radix ordering to determine unique values. \cr
                 3 \tab\tab "hash"   \tab\tab use index hashing to determine unique values. \cr
  }
  }
\item{cols}{compute unique rows according to a subset of columns. Columns can be selected using column names, indices, a logical vector or a selector function (i.e. \code{is.character}). \emph{Note:} All columns are returned. }
\item{\dots}{arguments passed to \code{\link{radixorderv}}, e.g. \code{decreasing} or \code{na.last}. Only applicable in \code{funique.default} if \code{method = "radix"}.}
}
\details{
If \code{x} is a data frame / list and all rows are already unique, then \code{x} is returned. Otherwise a copy of \code{x} with duplicate rows removed is returned.  See \code{\link{qF}} for some additional computational details.

The \emph{sf} method simply ignores the geometry column when determining unique values.
}
\value{
\code{x} with duplicate elements/rows removed. % sorted in ascending order if \code{sort = TRUE}, and in order of first occurrence if \code{sort = FALSE}.
}
\seealso{
\code{\link{qF}}, \link[=A2-fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}.
}
\examples{
funique(mtcars$cyl)
funique(gv(mtcars, c(2,8,9)))
funique(mtcars, cols = c(2,8,9))
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
