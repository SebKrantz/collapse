\name{BY}
\alias{BY}
\alias{BY.default}
\alias{BY.matrix}
\alias{BY.data.frame}
\alias{BY.grouped_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split-Apply-Combine Computing % (Efficient)
}
\description{
\code{BY} is an S3 generic that efficiently applies functions over vectors or matrix- and data.frame columns by groups, and returns various output formats. A simple parallelism is also available.
}
\usage{
BY(X, ...)

\method{BY}{default}(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same","list"))

\method{BY}{matrix}(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same","matrix","data.frame","list"))

\method{BY}{data.frame}(X, g, FUN, ..., use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same","matrix","data.frame","list"))

\method{BY}{grouped_df}(X, FUN, ..., use.g.names = FALSE, keep.group_vars = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same","matrix","data.frame","list"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a atomic vector, matrix or data frame.}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{FUN}{a function, can be scalar- or vector-valued.}
  \item{...}{further arguments to \code{FUN}.}
  \item{use.g.names}{make group-names and add to the result as names (vector method) or row-names (matrix and data.frame method). No row-names are generated for data.tables and grouped tibbles.}
  \item{sort}{logical. Sort the groups? Internally passed to \code{\link{GRP}} or \code{\link{qF}}, and only effective if \code{g} is not already a factor or \code{\link{GRP}} object.}
  \item{expand.wide}{logical. If \code{FUN} is a vector-valued function returning a vector of fixed length > 1 (such as the \code{\link{quantile}} function), \code{expand.wide} can be used to return the result in a wider format (instead of stacking the resulting vectors of fixed length above each other in each output column).}
\item{parallel}{logical. \code{TRUE} implements simple parallel execution by internally calling \code{parallel::mclapply} instead of \code{base::lapply}.}
  \item{mc.cores}{integer. Argument to \code{parallel::mclapply} indicating the number of cores to use for parallel execution. Can use \code{parallel::detectCores()} to select all available cores. See also \code{?parallel::mclapply}.}
\item{return}{an integer or string indicating the type of object to return. The default \code{1 - "same"} returns the same object type (i.e. passing a matrix returns a matrix and passig a data frame returns a data frame). \code{2 - "matrix"} always returns the output as matrix, \code{3 - "data.frame"} always returns a data frame and \code{4 - "list"} returns the raw (uncombined) output. \emph{Note}: \code{4 - "list"} works together with \code{expand.wide} to return a list of matrices.}
\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}
%    \item{simplify}{logical. Simplify the result to return an object of the same class and with the same attributes. If \code{FALSE}, the raw computation retult in form of a (nested) list is returned.}
}
\details{
\code{BY} is a frugal reimplementation of the Split-Apply-Combine computing paradigm. It is faster than \code{base::tapply}, \code{base::by}, \code{base::aggregate} and \code{plyr}, and preserves data attributes just like \code{\link{dapply}}. % and more versatile though not faster than \code{dplyr}

I note at this point that the philosphy of \code{collapse} is to move beyond this rather slow computing paradigm, which is why the \link[=fast-statistical-functions]{Fast Statistical Functions} were implemented. However sometimes tasks need to be performed that involve more comlex and customized operations on data, and for these cases \code{BY} is a good solution.

\code{BY} is built principally as a wrapper around \code{lapply(split(x, g), FUN, ...)}, but strongly optimizes on attribute cheking compared to \code{base} R. For more details examine the code yourself or look at the documentation for \code{\link{dapply}} which works very similar (the only difference really is the splitting performed in \code{BY}).

\code{BY} is used internally in \code{\link{collap}} (\code{collapse}'s main aggregation command) for functions that are not \link[=fast-statistical-functions]{Fast Statistical Functions}.
}
\value{
\code{X} where \code{FUN} was applied to every column split by \code{g}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dapply}}, \code{\link{collap}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}
}
\examples{
v <- iris$Sepal.Length   # A numeric vector
f <- iris$Species        # A factor. Vectors/lists will internally be converted to factor

## default vector method
BY(v, f, sum)                          # Sum by species
BY(v, f, scale)                        # Scale by species (please use fscale instead)
BY(v, f, scale, use.g.names = FALSE)   # Omitting auto-generated names
BY(v, f, quantile)                     # Species quantiles: by default stacked
BY(v, f, quantile, expand.wide = TRUE) # Wide format

## matrix method
m <- qM(num_vars(iris))
BY(m, f, sum)                          # Also return as matrix
BY(m, f, sum, return = "data.frame")   # Return as data.frame ... also works for computations below
BY(m, f, scale)
BY(m, f, scale, use.g.names = FALSE)
BY(m, f, quantile)
BY(m, f, quantile, expand.wide = TRUE)
BY(m, f, quantile, expand.wide = TRUE, # Return as list of matrices
   return = "list")

## data.frame method
BY(num_vars(iris), f, sum)             # Also returns a data.fram
BY(num_vars(iris), f, sum, return = 2) # Return as matrix ... also works for computations below
BY(num_vars(iris), f, scale)
BY(num_vars(iris), f, scale, use.g.names = FALSE)
BY(num_vars(iris), f, quantile)
BY(num_vars(iris), f, quantile, expand.wide = TRUE)
BY(num_vars(iris), f, quantile,        # Return as list of matrices
   expand.wide = TRUE, return = "list")

## grouped tibble method
library(dplyr)
giris <- group_by(iris, Species)
giris \%>\% BY(sum)                     # Compute sum
giris \%>\% BY(sum, use.g.names = TRUE, # Use row.names and
             keep.group_vars = FALSE) # remove 'Species' and groups attribute
giris \%>\% BY(sum, return = "matrix")  # Return matrix
giris \%>\% BY(sum, return = "matrix",  # Matrix without groups attribute and row.names
             use.g.names = TRUE) \%>\%
             `attr<-`("groups", NULL)
giris \%>\% BY(log)                     # Take logs
giris \%>\% BY(log, use.g.names = TRUE, # Use row.names and
             keep.group_vars = FALSE) # remove 'Species' and groups attribute
giris \%>\% BY(quantile)                # Compute quantiles (output is stacked)
giris \%>\% BY(quantile,                # Much better, also keeps 'Species'
             expand.wide = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line
% \keyword{programming}% use one of  RShowDoc("KEYWORDS")
