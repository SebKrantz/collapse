\name{fdroplevels}
\alias{fdroplevels}
\alias{fdroplevels.factor}
\alias{fdroplevels.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Removal of Unused Factor Levels
}
\description{
A substantially faster replacement for \code{\link{droplevels}}.
}
\usage{
fdroplevels(x, ...)

\method{fdroplevels}{factor}(x, ...)

\method{fdroplevels}{data.frame}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a factor, or data frame / list containing one or more factors.}
  \item{\dots}{not used.}
}
\details{
\code{\link{droplevels}} passes a factor from which levels are to be dropped to \code{\link{factor}}, which first calls \code{\link{unique}} and then \code{\link{match}} to drop unused levels. Both functions internally use a hash table, which is highly inefficient. \code{fdroplevels} does not require mapping values at all, but uses a super fast boolean vector method to determine which levels are unused and remove those levels. In addition, if no unused levels are found, \code{x} is simply returned. Any missing values found in \code{x} are efficiently skipped in the process of checking and replacing levels. All other attributes of \code{x} are preserved.
}
\value{
\code{x} with unused factor levels removed.
}
\note{
If \code{x} is malformed e.g. has too few levels, this function can cause a segmentation fault terminating the R session, thus only use with ordinary / proper factors.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qF}}, \code{\link{funique}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
f <- iris$Species[1:100]
fdroplevels(f)
identical(fdroplevels(f), droplevels(f))

fNA <- na_insert(f)
fdroplevels(fNA)
identical(fdroplevels(fNA), droplevels(fNA))

identical(fdroplevels(ss(iris, 1:100)), droplevels(ss(iris, 1:100)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
