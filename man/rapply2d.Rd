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
rapply2d(l, FUN, \dots, classes = "data.frame")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{FUN}{a function that can be applied to all elements in l.}
  \item{\dots}{additional elements passed to FUN.}
  \item{classes}{character. These are classes of list-based objects inside \code{l} which \code{FUN} should be applied to. Note that \code{FUN} is also applied to all non-list elements in \code{l}. It is thus quite different from the \code{classes} argument to \code{\link{rapply}}.}
}
\value{
A list of the same structure as \code{l}, where \code{FUN} was applied to all final (atomic) elements and list-based objects of a class included in \code{classes}.
}
%% ~Make other sections like Warning with \section{Warning }{\dots.} ~
\seealso{
\code{\link{rsplit}}, \code{\link{unlist2d}}, \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
l <- list(mtcars, list(mtcars, as.matrix(mtcars)))
rapply2d(l, fmean)
unlist2d(rapply2d(l, fmean))
}
\keyword{manip}
\keyword{list}
