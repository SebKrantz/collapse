\name{radixorder}
\alias{radixorder}
\alias{radixorderv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Radix-Based Ordering
}
\description{
A slight modification of \code{\link[=order]{order(..., method = "radix")}} that is more programmer friendly and, importantly, provides features for ordered grouping of data (similar to \code{data.table:::forderv} from which it descended). % \code{radixorderv} is a programmers version directly supporting vector and list input. % Apart from added grouping features, the source code and standard functionality is identical to \code{\link{order(\dots, method = "radix")}.
}
\usage{
radixorder(\dots, na.last = TRUE, decreasing = FALSE, starts = FALSE,
           group.sizes = FALSE, sort = TRUE)

radixorderv(x, na.last = TRUE, decreasing = FALSE, starts = FALSE,
            group.sizes = FALSE, sort = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{comma-separated atomic vectors to order.

}
  \item{x}{
an atomic vector or list of atomic vectors such as a data frame.
}
  \item{na.last}{logical. for controlling the treatment of \code{NA}'s. If \code{TRUE}, missing values in the data are put last; if \code{FALSE}, they are put first; if NA, they are removed.
}
  \item{decreasing}{
logical. Should the sort order be increasing or decreasing? Can be a vector of length equal to the number of arguments in \code{\dots} / \code{x}.
}
  \item{starts}{logical. \code{TRUE} returns an attribute 'starts' containing the first element of each new group i.e. the row denoting the start of each new group if the data were sorted using the computed ordering vector. See Examples.
%%     ~~Describe \code{starts} here~~
}
  \item{group.sizes}{logical. \code{TRUE} returns an attribute 'group.sizes' containing sizes of each group in the same order as groups are encountered if the data were sorted using the computed ordering vector. See Examples.
}
  \item{sort}{logical. This argument only affects character vectors / columns passed. If \code{FALSE}, these are not ordered but simply grouped in the order of first appearance of unique elements. This provides a slight performance gain if only grouping but not alphabetic ordering is required. See also \code{\link{group}}.
%%     ~~Describe \code{sort} here~~
}
}
% \details{
% \code{radixorder} works just like \code{\link[=order]{order(\dots, method = "radix")}}, the source code is the same. However if \code{starts = TRUE}, and attribute
% }
%}
\value{
An integer ordering vector with attributes: Unless \code{na.last = NA} an attribute \code{"sorted"} indicating whether the input data was already sorted is attached. If \code{starts = TRUE}, \code{"starts"} giving a vector of group starts in the ordered data, and if \code{group.sizes = TRUE}, \code{"group.sizes"} giving the vector of group sizes are attached. In either case an attribute \code{"maxgrpn"} providing the size of the largest group is also attached.
}

\author{
The C code was taken - with slight modifications - from \href{https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/main/radixsort.c}{base R source code}, and is originally due to \emph{data.table} authors Matt Dowle and Arun Srinivasan.
}


\seealso{
\link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
radixorder(mtcars$mpg)
head(mtcars[radixorder(mtcars$mpg), ])
radixorder(mtcars$cyl, mtcars$vs)

o <- radixorder(mtcars$cyl, mtcars$vs, starts = TRUE)
st <- attr(o, "starts")
head(mtcars[o, ])
mtcars[o[st], c("cyl", "vs")]  # Unique groups

# Note that if attr(o, "sorted") == TRUE, then all(o[st] == st)
radixorder(rep(1:3, each = 3), starts = TRUE)

# Group sizes
radixorder(mtcars$cyl, mtcars$vs, group.sizes = TRUE)

# Both
radixorder(mtcars$cyl, mtcars$vs, starts = TRUE, group.sizes = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
