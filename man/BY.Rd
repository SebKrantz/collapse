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

\method{BY}{default}(x, g, FUN, \dots, use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "vector", "list"))

\method{BY}{matrix}(x, g, FUN, \dots, use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "matrix", "data.frame", "list"))

\method{BY}{data.frame}(x, g, FUN, \dots, use.g.names = TRUE, sort = .op[["sort"]], reorder = TRUE,
   expand.wide = FALSE, parallel = FALSE, mc.cores = 1L,
   return = c("same", "matrix", "data.frame", "list"))

\method{BY}{grouped_df}(x, FUN, \dots, reorder = TRUE, keep.group_vars = TRUE, use.g.names = FALSE)
}
\arguments{
  \item{x}{a vector, matrix, data frame or alike object.}
  \item{g}{a \code{\link{GRP}} object, or a factor / atomic vector / list of atomic vectors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{FUN}{a function, can be scalar- or vector-valued. For vector valued functions see also \code{reorder} and \code{expand.wide}.}
  \item{\dots}{further arguments to \code{FUN}, or to \code{BY.data.frame} for the 'grouped_df' method. Since v1.9.0 data length arguments are also split by groups.}
  \item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). For vector-valued functions (row-)names are only generated if the function itself creates names for the statistics e.g. \code{quantile()} adds names, \code{range()} or \code{log()} don't. No row-names are generated on \emph{data.table}'s.}
  \item{sort}{logical. Sort the groups? Internally passed to \code{\link{GRP}}, and only effective if \code{g} is not already a factor or \code{\link{GRP}} object.}
  \item{reorder}{logical. If a vector-valued function is passed that preserves the data length, \code{TRUE} will reorder the result such that the elements/rows match the original data. \code{FALSE} just combines the data in order of the groups (i.e. all elements of the first group in first-appearance order followed by all elements in the second group etc..). \emph{Note} that if \code{reorder = FALSE}, grouping variables, names or rownames are only retained if the grouping is on sorted data, see \code{\link{GRP}}. }
  \item{expand.wide}{logical. If \code{FUN} is a vector-valued function returning a vector of fixed length > 1 (such as the \code{\link{quantile}} function), \code{expand.wide} can be used to return the result in a wider format (instead of stacking the resulting vectors of fixed length above each other in each output column).}
\item{parallel}{logical. \code{TRUE} implements simple parallel execution by internally calling \code{\link{mclapply}} instead of \code{\link{lapply}}. Parallelism is across columns, except for the default method.}
  \item{mc.cores}{integer. Argument to \code{\link{mclapply}} indicating the number of cores to use for parallel execution. Can use \code{\link[=detectCores]{detectCores()}} to select all available cores.}
\item{return}{an integer or string indicating the type of object to return. The default \code{1 - "same"} returns the same object type (i.e. class and other attributes are retained if the underlying data type is the same, just the names for the dimensions are adjusted). \code{2 - "matrix"} always returns the output as matrix, \code{3 - "data.frame"} always returns a data frame and \code{4 - "list"} returns the raw (uncombined) output. \emph{Note}: \code{4 - "list"} works together with \code{expand.wide} to return a list of matrices.}
\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation. See also the Note. }
%    \item{simplify}{logical. Simplify the result to return an object of the same class and with the same attributes. If \code{FALSE}, the raw computation retult in form of a (nested) list is returned.}
}
\details{
\code{BY} is a re-implementation of the Split-Apply-Combine computing paradigm. It is faster than \code{\link{tapply}}, \code{\link{by}}, \code{\link{aggregate}} and \emph{(d)plyr}, and preserves data attributes just like \code{\link{dapply}}.  % and more versatile though not faster than \emph{dplyr}

%I note at this point that the philosophy of \emph{collapse} is to move beyond this rather slow computing paradigm, which is why the \link[=fast-statistical-functions]{Fast Statistical Functions} were implemented. However sometimes tasks need to be performed that involve more complex and customized operations on data, and for these cases \code{BY} is a good solution.

It is principally a wrapper around \code{lapply(gsplit(x, g), FUN, \dots)}, that uses \code{\link{gsplit}} for optimized splitting and also strongly optimizes on the internal code compared to \emph{base} R functions. For more details look at the documentation for \code{\link{dapply}} which works very similar (apart from the splitting performed in \code{BY}). The function is intended for simple cases involving flexible computation of statistics across groups using a single function e.g. \code{iris |> gby(Species) |> BY(IQR)} is simpler than \code{iris |> gby(Species) |> smr(acr(.fns = IQR))} etc..

% For larger tasks, the \link[=fast-statistical-functions]{Fast Statistical Functions} or the \emph{data.table} package are more appropriate tools.
}
%\note{
%\code{BY} can be used with vector-valued functions preserving the length of the data, note however that, unlike \code{\link{fmutate}}, data is recombined in the order of the groups, not in the order of the original data. It is thus advisable to sort the data by the grouping variable before using \code{BY} with such a function. In particular, in such cases the 'grouped_df' method only keeps grouping columns if data was grouped with \code{fgroup_by(data, ..., sort = TRUE)}, and the grouping algorithm detected that the data is already sorted in the order of the groups (i.e. if \code{attr(with(data, radixorder(...)), "sorted")} is \code{TRUE}), even if \code{keep.group_vars = TRUE}. The same holds for preservation names / rownames in the default, matrix or data frame methods. Basically, \code{BY} is kept as simple as possible without running danger of returning something wrong.
%}
\value{
\code{X} where \code{FUN} was applied to every column split by \code{g}.
}
\seealso{
\code{\link{dapply}}, \code{\link{collap}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
v <- iris$Sepal.Length   # A numeric vector
g <- GRP(iris$Species)   # A grouping

## default vector method
BY(v, g, sum)                                # Sum by species
head(BY(v, g, scale))                        # Scale by species (please use fscale instead)
BY(v, g, fquantile)                          # Species quantiles: by default stacked
BY(v, g, fquantile, expand.wide = TRUE)      # Wide format

## matrix method
m <- qM(num_vars(iris))
BY(m, g, sum)                          # Also return as matrix
BY(m, g, sum, return = "data.frame")   # Return as data.frame.. also works for computations below
head(BY(m, g, scale))
BY(m, g, fquantile)
BY(m, g, fquantile, expand.wide = TRUE)
ml <- BY(m, g, fquantile, expand.wide = TRUE, # Return as list of matrices
         return = "list")
ml
# Unlisting to Data Frame
unlist2d(ml, idcols = "Variable", row.names = "Species")

## data.frame method
BY(num_vars(iris), g, sum)             # Also returns a data.fram
BY(num_vars(iris), g, sum, return = 2) # Return as matrix.. also works for computations below
head(BY(num_vars(iris), g, scale))
BY(num_vars(iris), g, fquantile)
BY(num_vars(iris), g, fquantile, expand.wide = TRUE)
BY(num_vars(iris), g, fquantile,       # Return as list of matrices
   expand.wide = TRUE, return = "list")

## grouped data frame method
giris <- fgroup_by(iris, Species)
giris |> BY(sum)                      # Compute sum
giris |> BY(sum, use.g.names = TRUE,  # Use row.names and
             keep.group_vars = FALSE)  # remove 'Species' and groups attribute
giris |> BY(sum, return = "matrix")   # Return matrix
giris |> BY(sum, return = "matrix",   # Matrix with row.names
             use.g.names = TRUE)
giris |> BY(.quantile)                # Compute quantiles (output is stacked)
giris |> BY(.quantile, names = TRUE,  # Wide output
               expand.wide = TRUE)
}
\keyword{manip}
