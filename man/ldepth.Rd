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
ldepth(l, DF.as.list = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{DF.as.list}{treat data.frame's as sub-lists?}
}
\details{
The depth or level or nesting of a list or list-like structure (i.e. a classed object) is found by recursing down to the bottom of the list and adding an integer count of 1 for each level passed. For example the depth of a \code{data.frame} is 1. If a \code{data.frame} has list-columns, the depth is 2. However for reasons of efficiency, if \code{l} is not a \code{data.frame} and \code{DF.as.list = TRUE}, \code{data.frame}'s found inside \code{l} will not be checked for list column's but assumed to have a depth of 1.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{is.unlistable}}, \code{\link{has_elem}}, \link[=list-processing]{List Processing}
}
\examples{
l = list(1, 2)
ldepth(l)
l = list(1, 2, mtcars)
ldepth(l)
ldepth(l, DF.as.list = FALSE)
l = list(1, 2, list(4, 5, list(6, mtcars)))
ldepth(l)
ldepth(l, DF.as.list = FALSE)
}
