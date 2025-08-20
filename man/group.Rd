\name{group}
\alias{group}
\alias{groupv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Hash-Based Grouping
}
\description{
\code{group()} scans the rows of a data frame (or atomic vector / list of atomic vectors), assigning to each unique row an integer id - starting with 1 and proceeding in first-appearance order of the rows. The function is written in C and optimized for R's data structures. It is the workhorse behind functions like \code{\link{GRP}} / \code{\link{fgroup_by}}, \code{\link{collap}}, \code{\link{qF}}, \code{\link{qG}}, \code{\link{finteraction}} and \code{\link{funique}}, when called with argument \code{sort = FALSE}.
}
\usage{
group(\dots, starts = FALSE, group.sizes = FALSE)

groupv(x, starts = FALSE, group.sizes = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{comma separated atomic vectors to group. Also supports a single list of vectors for backward compatibility.}
  \item{x}{an atomic vector or data frame / list of equal-length atomic vectors.}
  \item{starts}{logical. If \code{TRUE}, an additional attribute \code{"starts"} is attached giving a vector of group starts (= index of first-occurrence of unique rows).
}
  \item{group.sizes}{
logical. If \code{TRUE}, an additional attribute \code{"group.sizes"} is attached giving the size of each group.
}
}
\details{
A data frame is grouped on a column-by-column basis, starting from the leftmost column. For each new column the grouping vector obtained after the previous column is also fed back into the hash function so that unique values are determined on a running basis. The algorithm terminates as soon as the number of unique rows reaches the size of the data frame. Missing values are also grouped just like any other values. Invoking arguments \code{starts} and/or \code{group.sizes} requires an additional pass through the final grouping vector.
}
\value{
An object is of class 'qG' see \code{\link{qG}}.
}
\author{
The Hash Function and inspiration was taken from the excellent \emph{kit} package by Morgan Jacob, the algorithm was developed by Sebastian Krantz.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{radixorder}}, \code{\link{GRPid}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# Let's replicate what funique does
g <- groupv(wlddev, starts = TRUE)
if(attr(g, "N.groups") == fnrow(wlddev)) wlddev else
   ss(wlddev, attr(g, "starts"))

}
