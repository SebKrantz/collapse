\name{rapply2d}
\alias{rapply2d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recursively Apply a Function to a List of Data Objects
}
\description{
\code{rapply2d} is a recursive version of \code{lapply} with three differences to \code{\link{rapply}}:

\enumerate{
\item data frames (or other list-based objects specified in \code{classes}) are considered as atomic, not as (sub-)lists
\item \code{FUN} is applied to all 'atomic' objects in the nested list
\item the result is not simplified / unlisted.
}

}
\usage{
rapply2d(l, FUN, \dots, classes = "data.frame")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{FUN}{a function that can be applied to all 'atomic' elements in l.}
  \item{\dots}{additional elements passed to FUN.}
  \item{classes}{character. Classes of list-based objects inside \code{l} that should be considered as atomic. }
}
\value{
A list of the same structure as \code{l}, where \code{FUN} was applied to all atomic elements and list-based objects of a class included in \code{classes}.
}
\note{
The main reason \code{rapply2d} exists is to have a recursive function that out-of-the-box applies a function to a nested list of data frames.

For most other purposes \code{\link{rapply}}, or by extension the excellent \href{https://cran.r-project.org/package=rrapply}{rrapply} function / package, provide more advanced functionality and greater performance.
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
