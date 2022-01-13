\name{is_unlistable}
\alias{is_unlistable}
\title{
Unlistable Lists
}
\description{
A (nested) list with atomic objects in all final nodes of the list-tree is unlistable - checked with \code{is_unlistable}.
}
\usage{
is_unlistable(l, DF.as.list = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 % \item{x}{an R object.}
    \item{l}{a list.}
    \item{DF.as.list}{logical. \code{TRUE} treats data frames like (sub-)lists; \code{FALSE} like atomic elements.}
}
\details{
\code{is_unlistable} with \code{DF.as.list = TRUE} is defined as \code{all(rapply(l, is.atomic))}, whereas \code{DF.as.list = FALSE} yields checking using \code{all(unlist(rapply2d(l, function(x) is.atomic(x) || is.list(x)), use.names = FALSE))}, assuming that data frames are lists composed of atomic elements.  If \code{l} contains data frames, the latter can be a lot faster than applying \code{is.atomic} to every data frame column.
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
\code{\link{ldepth}}, \code{\link{has_elem}}, \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
l <- list(1, 2, list(3, 4, "b", FALSE))
is_unlistable(l)
l <- list(1, 2, list(3, 4, "b", FALSE, e ~ b))
is_unlistable(l)

}
\keyword{list}
\keyword{utilities}
