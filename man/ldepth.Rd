\name{ldepth}
\alias{ldepth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the Depth / Level of Nesting of a List
}
\description{
\code{ldepth} provides the depth of a list or list-like structure.
}
\usage{
ldepth(l, DF.as.list = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{DF.as.list}{logical. \code{TRUE} treats data frames like (sub-)lists; \code{FALSE} like atomic elements.}
}
\details{
The depth or level or nesting of a list or list-like structure (e.g. a model object) is found by recursing down to the bottom of the list and adding an integer count of 1 for each level passed. For example the depth of a data frame is 1. If a data frame has list-columns, the depth is 2. However for reasons of efficiency, if \code{l} is not a data frame and \code{DF.as.list = FALSE}, data frames found inside \code{l} will not be checked for list column's but assumed to have a depth of 1.
}
\value{
A single integer indicating the depth of the list.
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
\code{\link{is_unlistable}}, \code{\link{has_elem}}, \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
l <- list(1, 2)
ldepth(l)
l <- list(1, 2, mtcars)
ldepth(l)
ldepth(l, DF.as.list = FALSE)
l <- list(1, 2, list(4, 5, list(6, mtcars)))
ldepth(l)
ldepth(l, DF.as.list = FALSE)
}
\keyword{list}
\keyword{utilities}
