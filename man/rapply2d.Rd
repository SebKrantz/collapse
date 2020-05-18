\name{rapply2d}
\alias{rapply2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recursively Apply a Function to a List of Data Objects
}
\description{
\code{rapply2d} is a recursive version of \code{lapply} with two key differences to \code{\link{rapply}}: (1) Data frames are considered as final objects, not as (sub-)lists, and (2) the result is never simplified / unlisted.
}
\usage{
rapply2d(l, FUN, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{FUN}{a function that can be applied to all elements in l.}
  \item{...}{additional elements passed to FUN.}
}
\value{
A list of the same structure as \code{l}, where \code{FUN} was applied to all elements.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{unlist2d}}, \link[=A8-list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
l <- list(mtcars, list(mtcars, as.matrix(mtcars)))
rapply2d(l, fmean)
unlist2d(rapply2d(l, fmean))
}
\keyword{manip}
\keyword{list}
