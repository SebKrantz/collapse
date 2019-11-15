\name{TRA}
\alias{TRA}
\alias{TRA.default}
\alias{TRA.matrix}
\alias{TRA.data.frame}
\alias{TRA.grouped_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Data by (Groupwise) Replacing or Sweeping out Statistics
}
\description{
\code{TRA} is an S3 generic that efficiently transforms data by either (column-wise) replacing data values with supplied statistics or sweeping the statistics out of the data. \code{TRA} supports grouped operations and/on data frame's, and is thus a generalization of \code{base::sweep}.
}
\usage{
\method{TRA}{default}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{matrix}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{data.frame}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{grouped_df}(x, STATS, FUN = "-", keep.group_keys = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a atomic vector, matrix, data frame or grouped tibble (\code{dplyr::grouped_df}).}
  \item{STATS}{a matching set of summary statistics computed on \code{x}. If \code{g = NULL} (no groups), all methods support an atomic vector of statistics of length \code{NCOL(x)}. The matrix and data.frame methods also support a 1-row matrix or 1-row data.frame/list, repsectively. If groups are supplied to \code{g}, \code{STATS} needs to be of the same type as \code{x} and of appropriate dimensions (such that \code{NCOL(x) == NCOL(STATS)} and \code{NROW(STATS)} matches the number of groups supplied to \code{g} i.e. the number of levels if \code{g} is a factor, with the first row of \code{STATS} corresponding to the first level of \code{g} etc...)}
  \item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}
  \item{FUN}{a function, can be scalar- or vector-valued.}
  \item{...}{further arguments to \code{FUN}.}
  \item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}
  \item{sort}{logical. Sort the groups?. Internally passed to \code{GRP()} or \code{qF()}, and only effective if \code{g} is not already a factor or \code{GRP} object.}
  \item{expand.wide}{logical. If \code{FUN} is a vector-valued function that applied to a vector returns a vector of a fixed length > 1 (such as the \code{quantile} function), \code{expand.wide} can be used to return the result in a wider format (instead of stacking the resulting vectors of fixed length above each other in each output column).}
  \item{parallel}{logical. \code{TRUE} implements simple parallel execution by internally calling \code{parallel::mclapply} instead of \code{base::lapply}.}
  \item{mc.cores}{integer. Argument to \code{parallel::mclapply} indicating the number of cores to use for parallel execution. Can use \code{parallel::detectCores()} to select all available cores. See also \code{?parallel::mclapply}.}
    \item{return}{an integer or string indicating the type of object to return. The default \code{1 - "same"} returns the same object type (i.e. passing a matrix returns a matrix and passig a data frame returns a data frame). \code{2 - "matrix"} always returns the output as matrix, \code{3 - "data.frame"} always returns a data frame and \code{4 - "list"} returns the raw (uncombined) output. \emph{Note}: \code{4 - "list"} works together with \code{expand.wide} to return a list of matrices.}
%    \item{simplify}{logical. Simplify the result to return an object of the same class and with the same attributes. If \code{FALSE}, the raw computation retult in form of a (nested) list is returned.}
}
\details{
\code{TRA} is an efficient and frugal reimplementation of the Split-Apply-Combine computing paradigm. It is faster than \code{base::tapply}, \code{base::by}, \code{base::aggregate} and \code{plyr}, and preserves data attributes just like \code{dapply}.

I note at this point that the philosphy of \code{collapse} is to move beyond this rather slow computing paradigm, which is why the \code{\link{Fast Statistical Functions}} were implemented. However sometimes tasks need to be performed that involve more comlex and customized operations on data, and for these cases \code{TRA} represents a very efficient solution.

\code{TRA} is built principally as a wrapper around \code{lapply(split(x, g), FUN, ...)}, but strongly optimizes on attribute cheking compared to \code{base} R. For more details examine the code yourself or look at the documentation for \code{\link{dapply}} which works very similar (the only difference really is the splitting performed in \code{TRA}).

\code{TRA} is used internally in \code{collap} (\code{collapse}'s main aggregation command) for functions that are not \code{\link{Fast Statistical Functions}}.
}
\value{
\code{X} where \code{FUN} was applied to every column split by \code{g}.
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
v <- iris$Sepal.Length   # A numeric vector
f <- iris$Species        # A factor. Vectors and lists passed to g will internally be converted to factor

## default vector method
TRA(v, f, sum)                          # Sum by species
TRA(v, f, scale)                        # Scale by species (this is only a demo, please use fscale instead)
TRA(v, f, scale, use.g.names = FALSE)   # Omitting auto-generated names
TRA(v, f, quantile)                     # Species quantiles: by default stacked
TRA(v, f, quantile, expand.wide = TRUE) # Wide format

## data.frame method
TRA(num_vars(iris), f, sum)             # Also returns a data.fram
TRA(num_vars(iris), f, sum, return = 2) # Return as matrix ... also works for all the other computations below
TRA(num_vars(iris), f, scale)
TRA(num_vars(iris), f, scale, use.g.names = FALSE)
TRA(num_vars(iris), f, quantile)
TRA(num_vars(iris), f, quantile, expand.wide = TRUE)
TRA(num_vars(iris), f, quantile,        # Return as list of matrices
   expand.wide = TRUE, return = "list")

## matrix method
m <- qM(num_vars(iris))
TRA(m, f, sum)                          # Also return as matrix
TRA(m, f, sum, return = "data.frame")   # Return as data.frame ... also works for all the other computations below
TRA(m, f, scale)
TRA(m, f, scale, use.g.names = FALSE)
TRA(m, f, quantile)
TRA(m, f, quantile, expand.wide = TRUE)
TRA(m, f, quantile, expand.wide = TRUE, # Return as list of matrices
   return = "list")

# grouped tibble method
library(dplyr)
giris <- group_by(iris, Species)
giris \%>\% TRA(sum)                     # Compute sum
giris \%>\% TRA(sum) \%>\% ungroup()       # drop groups attribute
giris \%>\% TRA(sum, use.g.names = TRUE, # Use row.names and
             keep.group_keys = FALSE) # remove 'Species' and groups attribute
giris \%>\% TRA(sum, return = "matrix")  # Return matrix
giris \%>\% TRA(sum, return = "matrix",  # Matrix without groups attribute and row.names
             use.g.names = TRUE) \%>\%
             `attr<-`("groups", NULL)
giris \%>\% TRA(log)                     # Take logs
giris \%>\% TRA(log, use.g.names = TRUE, # Use row.names and
             keep.group_keys = FALSE) # remove 'Species' and groups attribute
giris \%>\% TRA(quantile)                # Compute quantiles (output is stacked)
giris \%>\% TRA(quantile,                # Much better, also keeps 'Species'
             expand.wide = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
