\name{t_list}
\alias{t_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Efficient List Transpose
}
\description{
\code{t_list} turns a list of lists inside-out. The performance is quite efficient regardless of the size of the list.
}
\usage{
t_list(l)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list of lists. Elements inside the sublists can be heterogeneous, including further lists.
}
}
\value{
\code{l} transposed such that the second layer of the list becomes the top layer and the top layer the second layer. See Examples.
}

\note{
To transpose a data frame / list of atomic vectors see \code{data.table::transpose()}.
}

\seealso{
\code{\link{rsplit}}, \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# Homogenous list of lists
l <- list(a = list(c = 1, d = 2), b = list(c = 3, d = 4))
str(l)
str(t_list(l))

# Heterogenous case
l2 <- list(a = list(c = 1, d = letters), b = list(c = 3:10, d = list(4, e = 5)))
attr(l2, "bla") <- "abc"  # Attributes other than names are preserved
str(l2)
str(t_list(l2))

rm(l, l2)
}
\keyword{list}
\keyword{manip}
\keyword{utilities}

% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
