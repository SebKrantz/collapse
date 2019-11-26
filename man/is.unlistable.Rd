\name{is.regular-is.unlistable}
\alias{is.unlistable}
\alias{is.regular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regular Objects and Unlistable Lists
}
\description{
A regular R object is an R object that is either atomic or a list - checked with \code{is.regular}. A (nested) list composed only of regular objects is unlistable - checked with \code{is.unlistable}.
}
\usage{
is.regular(x)
is.unlistable(l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a R object.}
  \item{l}{a list.}
}
\details{
\code{is.regular} is simply defined as \code{is.atomic(x) || is.list(x)}. \code{is.unlistable} is defined as \code{all(unlist(rapply2d(l, is.regular), use.names = FALSE))}. It could of course also be defined as \code{all(rapply(l, is.atomic))}, but the above is a lot more efficient if \code{l} contains \code{data.frame}'s.
}
\value{
Logical(1) - TRUE or FALSE.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
