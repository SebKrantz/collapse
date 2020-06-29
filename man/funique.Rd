\name{funique}
\alias{funique}
\alias{funique.default}
\alias{funique.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Unique Elements / Rows
}
\description{
\code{funique} is a substantially faster alternative to \code{\link{unique}}. It is generic with a default vector method and a data.frame method.
}
\usage{
funique(x, \dots)

\method{funique}{default}(x, sort = TRUE, method = c("auto", "hash", "radix"), \dots)

\method{funique}{data.frame}(x, cols = NULL, sort = TRUE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a atomic vector or data.frame. }
  \item{sort}{logical. \code{TRUE} (default) orders the uniqe elements / rows. \code{FALSE} returns unique character values in order of first occurrence. Numeric values are only provided in order of first occurrence when \code{method = "hash"}. }
   \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto"   \tab\tab automatic selection: hash for character, logical or if \code{length(x) < 500}, else radix. \cr
                 2 \tab\tab "radix" \tab\tab use radix ordering to determine unique values. \cr
                 3 \tab\tab "hash"   \tab\tab use index hashing to determine unique values. \cr
  }
  }
\item{cols}{compute unique rows according to a subset of columns. Columns can be selected using column names, indices or a selector function (i.e. \code{is.numeric}). \emph{Note:} All columns are returned. }
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
See \code{\link{qF}} for some additional computational details.
}
\value{
Same as \code{\link{unique}}.
}

\seealso{
\code{\link{qF}}, \link[=A2-fast-grouping]{Fast (Ordered) Grouping}, \link[=collapse-documentation]{Collapse Overview}.
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
