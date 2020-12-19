\name{is.regular, is.unlistable}
\alias{is.unlistable}
\alias{is.regular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regular Objects and Unlistable Lists
}
\description{
In \emph{collapse} a regular R object is defined as an R object that is either atomic or a list - checked with \code{is.regular}. A (nested) list composed of regular objects at each level is unlistable - checked with \code{is.unlistable}.
}
\usage{
is.regular(x)
is.unlistable(l, DF.as.list = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an R object.}
  \item{l}{a list.}
    \item{DF.as.list}{logical. \code{TRUE} treats data frames like (sub-)lists; \code{FALSE} like atomic elements.}
}
\details{
\code{is.regular} is simply defined as \code{is.atomic(x) || is.list(x)}. \code{is.unlistable} with \code{DF.as.list = TRUE} is defined as \code{all(rapply(l, is.atomic))}, whereas \code{DF.as.list = FALSE} yields checking using \code{all(unlist(rapply2d(l, is.regular), use.names = FALSE))}.  If \code{l} contains data frames, the latter can be a lot faster than applying \code{is.atomic} to every data frame column.
}
\value{
\code{logical(1)} - \code{TRUE} or \code{FALSE}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{ldepth}}, \code{\link{has_elem}}, \link[=A8-list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
is.regular(list(1,2))
is.regular(2)
is.regular(a ~ c)
l <- list(1, 2, list(3, 4, "b", FALSE))
is.regular(l)
is.unlistable(l)
l <- list(1, 2, list(3, 4, "b", FALSE, e ~ b))
is.regular(l)
is.unlistable(l)

}
\keyword{list}
\keyword{utilities}
