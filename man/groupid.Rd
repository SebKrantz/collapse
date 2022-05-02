\name{groupid}
\alias{groupid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Run-Length Type Group-Id
}
\description{
\code{groupid} is an enhanced version of \code{data.table::rleid} for atomic vectors. It generates a run-length type group-id where consecutive identical values are assigned the same integer. It is a generalization as it can be applied to unordered vectors, generate group id's starting from an arbitrary value, and skip missing values.
}
\usage{
groupid(x, o = NULL, start = 1L, na.skip = FALSE, check.o = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an atomic vector of any type. Attributes are not considered.}

  \item{o}{an (optional) integer ordering vector specifying the order by which to pass through \code{x}.}

  \item{start}{integer. The starting value of the resulting group-id. Default is starting from 1.} %For C++ programmers, starting from 0 could be a better choice. }

  \item{na.skip}{logical. Skip missing values i.e. if \code{TRUE} something like \code{groupid(c("a", NA, "a"))} gives \code{c(1, NA, 1)} whereas \code{FALSE} gives \code{c(1, 2, 3)}.}

 \item{check.o}{logical. Programmers option: \code{FALSE} prevents checking that each element of \code{o} is in the range \code{[1, length(x)]}, it only checks the length of \code{o}. This gives some extra speed, but will terminate R if any element of \code{o} is too large or too small. }

}
\value{
An integer vector of class 'qG'. See \code{\link{qG}}.
}

\seealso{
\code{\link{seqid}}, \code{\link{timeid}}, \code{\link{qG}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
groupid(airquality$Month)
groupid(airquality$Month, start = 0)
groupid(wlddev$country)[1:100]

## Same thing since country is alphabetically ordered: (groupid is faster..)
all.equal(groupid(wlddev$country), qG(wlddev$country, na.exclude = FALSE))

## When data is unordered, group-id can be generated through an ordering..
uo <- order(rnorm(fnrow(airquality)))
monthuo <- airquality$Month[uo]
o <- order(monthuo)
groupid(monthuo, o)
identical(groupid(monthuo, o)[o], unattrib(groupid(airquality$Month)))
}
\keyword{manip}
