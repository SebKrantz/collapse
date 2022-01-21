\name{BY}
\alias{BY}
\alias{BY.default}
\alias{BY.matrix}
\alias{BY.data.frame}
\alias{BY.grouped_df}
\title{
Split-Apply-Combine Computing % (Efficient)
}
\description{
\code{BY} is an S3 generic that efficiently applies functions over vectors or matrix- and data frame columns by groups. Similar to \code{\link{dapply}} it seeks to retain the structure and attributes of the data, but can also output to various standard formats. A simple parallelism is also available.
}
\usage{
BY(x, \dots)

\method{BY}{default}(x, g, FUN, \dots, use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "vector", "list"))

\method{BY}{matrix}(x, g, FUN, \dots, use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "matrix", "data.frame", "list"))

\method{BY}{data.frame}(x, g, FUN, \dots, use.g.names = TRUE, sort = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "matrix", "data.frame", "list"))

\method{BY}{grouped_df}(x, FUN, \dots, keep.group_vars = TRUE, use.g.names = FALSE)
}
\arguments{
  \item{x}{a atomic vector, matrix, data frame or alike object.}
  \item{g}{a \code{\link{GRP}} object, or a factor / atomic vector / list of atomic vectors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{FUN}{a function, can be scalar- or vector-valued. For vector valued functions see \code{expand.wide} and the Note.}
  \item{\dots}{further arguments to \code{FUN}, or to \code{BY.data.frame} for the 'grouped_df' method.}
  \item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}
  \item{sort}{logical. Sort the groups? Internally passed to \code{\link{GRP}}, and only effective if \code{g} is not already a factor or \code{\link{GRP}} object.}
  \item{expand.wide}{logical. If \code{FUN} is a vector-valued function returning a vector of fixed length > 1 (such as the \code{\link{quantile}} function), \code{expand.wide} can be used to return the result in a wider format (instead of stacking the resulting vectors of fixed length above each other in each output column).}
\item{parallel}{logical. \code{TRUE} implements simple parallel execution by internally calling \code{\link{mclapply}} instead of \code{\link{lapply}}.}
  \item{mc.cores}{integer. Argument to \code{\link{mclapply}} indicating the number of cores to use for parallel execution. Can use \code{\link[=detectCores]{detectCores()}} to select all available cores.}
\item{return}{an integer or string indicating the type of object to return. The default \code{1 - "same"} returns the same object type (i.e. class and other attributes are retained if the underlying data type is the same, just the names for the dimensions are adjusted). \code{2 - "matrix"} always returns the output as matrix, \code{3 - "data.frame"} always returns a data frame and \code{4 - "list"} returns the raw (uncombined) output. \emph{Note}: \code{4 - "list"} works together with \code{expand.wide} to return a list of matrices.}
\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation. See also the Note. }
%    \item{simplify}{logical. Simplify the result to return an object of the same class and with the same attributes. If \code{FALSE}, the raw computation retult in form of a (nested) list is returned.}
}
\details{
\code{BY} is a frugal re-implementation of the Split-Apply-Combine computing paradigm. It is faster than \code{\link{tapply}}, \code{\link{by}}, \code{\link{aggregate}} and \emph{plyr}, and preserves data attributes just like \code{\link{dapply}}.  % and more versatile though not faster than \emph{dplyr}

%I note at this point that the philosophy of \emph{collapse} is to move beyond this rather slow computing paradigm, which is why the \link[=fast-statistical-functions]{Fast Statistical Functions} were implemented. However sometimes tasks need to be performed that involve more complex and customized operations on data, and for these cases \code{BY} is a good solution.

It is principally a wrapper around \code{lapply(gsplit(x, g), FUN, \dots)}, that uses \code{\link{gsplit}} for optimized splitting and also strongly optimizes on the internal code compared to \emph{base} R functions. For more details look at the documentation for \code{\link{dapply}} which works very similar (apart from the splitting performed in \code{BY}). The function is intended for simple usage involving data aggregation or flexible computation of summary statistics across groups.
For larger tasks requiring split-apply-combine computing on data frames, the \link[=fast-statistical-functions]{Fast Statistical Functions} and the \emph{data.table} package are more appropriate tools.
}
\note{
\code{BY} can be used with vector-valued functions preserving the length of the data, note however that data are recombined in the order of the groups, not in the order of the original data. It is thus advisable to sort the data by the grouping variable before using \code{BY} with such a function.
}
\value{
\code{X} where \code{FUN} was applied to every column split by \code{g}.
}
\seealso{
\code{\link{dapply}}, \code{\link{collap}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
v <- iris$Sepal.Length   # A numeric vector
f <- GRP(iris$Species)   # A grouping

## default vector method
BY(v, f, sum)                                # Sum by species
head(BY(v, f, scale))                        # Scale by species (please use fscale instead)
head(BY(v, f, scale, use.g.names = FALSE))   # Omitting auto-generated names
BY(v, f, quantile)                           # Species quantiles: by default stacked
BY(v, f, quantile, expand.wide = TRUE)       # Wide format

## matrix method
m <- qM(num_vars(iris))
BY(m, f, sum)                          # Also return as matrix
BY(m, f, sum, return = "data.frame")   # Return as data.frame.. also works for computations below
head(BY(m, f, scale))
head(BY(m, f, scale, use.g.names = FALSE))
BY(m, f, quantile)
BY(m, f, quantile, expand.wide = TRUE)
BY(m, f, quantile, expand.wide = TRUE, # Return as list of matrices
   return = "list")

## data.frame method
BY(num_vars(iris), f, sum)             # Also returns a data.fram
BY(num_vars(iris), f, sum, return = 2) # Return as matrix.. also works for computations below
head(BY(num_vars(iris), f, scale))
head(BY(num_vars(iris), f, scale, use.g.names = FALSE))
BY(num_vars(iris), f, quantile)
BY(num_vars(iris), f, quantile, expand.wide = TRUE)
BY(num_vars(iris), f, quantile,        # Return as list of matrices
   expand.wide = TRUE, return = "list")
 % No code relying on suggested package
## grouped data frame method
library(magrittr) # Note: Used because |> is not available on older R versions
giris <- fgroup_by(iris, Species)
giris \%>\% BY(sum)                      # Compute sum
giris \%>\% BY(sum, use.g.names = TRUE,  # Use row.names and
             keep.group_vars = FALSE) # remove 'Species' and groups attribute
giris \%>\% BY(sum, return = "matrix")   # Return matrix
giris \%>\% BY(sum, return = "matrix",   # Matrix with row.names
             use.g.names = TRUE)
giris \%>\% BY(quantile)                 # Compute quantiles (output is stacked)
giris \%>\% BY(quantile,                 # Much better, also keeps 'Species'
             expand.wide = TRUE)
}
\keyword{manip}
