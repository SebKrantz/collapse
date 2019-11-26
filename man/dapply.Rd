\name{dapply}
\alias{dapply}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Apply
}
\description{
\code{dapply} efficiently applies functions to columns or rows of matrices and data frame's and (default) returns an object of the same type and with the same attributes, or converts to the other type. A simple parallelism is also available.
}
\usage{
dapply(X, FUN, ..., MARGIN = 2, parallel = FALSE, mc.cores = 1L,
       return = c("same","matrix","data.frame"), drop = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or data frame.}
  \item{FUN}{a function, can be scalar- or vector-valued.}
  \item{...}{further arguments to \code{FUN}.}
  \item{MARGIN}{integer. The margin which \code{FUN} will be applied over. Default \code{2} indicates columns while \code{1} indicates rows. See also Details. }
  \item{parallel}{logical. \code{TRUE} implements simple parallel execution by internally calling \code{parallel::mclapply} instead of \code{base::lapply}.}
  \item{mc.cores}{integer. Argument to \code{parallel::mclapply} indicating the number of cores to use for parallel execution. Can use \code{parallel::detectCores()} to select all available cores. See also \code{?parallel::mclapply}.}
  \item{return}{an integer or string indicating the type of object to return. The default \code{1 - "same"} returns the same object type (i.e. passing a matrix returns a matrix and passig a data frame returns a data frame). \code{2 - "matrix"} always returns the output as matrix and \code{3 - "data.frame"} always returns a data frame.}
  \item{drop}{logical. If the result has only one row or one column, \code{drop = TRUE} will drop dimensions and return a (named) atomic vector.}
}
\details{
\code{dapply} is an efficient command to apply functions to rows or columns of data without loosing information (attributes) about the data or changing the classes or format of the data. It is principally an efficient wrapper around \code{base::lapply} and works as follows:

\itemize{
\item Save the attributes of \code{X}.
\item If \code{MARGIN = 2} (columns), convert matrices to plain lists of columns using \code{\link{mctl}} and remove all attributes from data frames.
\item If \code{MARGIN = 1} (rows), convert matrices to plain lists of rows using \code{\link{mrtl}}. For data frames remove all attributes, efficiently convert to matrix using \code{do.call(rbind, X)} and also convert to list of rows using \code{\link{mrtl}}.
\item Call \code{base::lapply} or \code{parallel::mclapply} on these plain lists (which is faster than calling \code{lapply} on an object with attributes).

\item depending on the requested output type, use \code{base::matrix}, \code{base::unlist} or \code{do.call(cbind, ...)} to convert the result back to a matrix or list of columns.

\item modify the relevant attributes accordingly and efficiently attach to the object again (no further checks, non essential attributes are kept and added at the end of the attribute list).
}

This performance gain from working with plain lists makes \code{dapply} not much slower than calling \code{lapply} itself on a data frame. Because of the conversions involved, row-operations require some memory, but are still faster than \code{base::apply}.
}
\value{
\code{X} where \code{FUN} was applied to every row or column.
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
\code{\link{BY}}, \code{\link{collap}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}
}
\examples{
dapply(mtcars, log)                      # Take natural log of each variable
dapply(mtcars, log, return = "matrix")   # Return as matrix
m <- as.matrix(mtcars)
dapply(m, log)                           # Same thing
dapply(m, log, return = "data.frame")    # Return data frame from matrix
dapply(mtcars, sum); dapply(m, sum)      # Computing sum of each column, return as vector
dapply(mtcars, sum, drop = FALSE)        # This returns a data.frame of 1 row
dapply(mtcars, sum, MARGIN = 1)          # Compute row-sum of each column, return as vector. More efficient than rowSums(mtcars).
dapply(m, sum, MARGIN = 1)               # Same thing for matrices, faster than apply(m, 1, sum)
dapply(m, sum, MARGIN = 1, drop = FALSE) # Gives matrix with one column
dapply(m, quantile, MARGIN = 1)          # Compute row-quantiles
dapply(m, quantile)                      # Column-quantiles
dapply(mtcars, quantile, MARGIN = 1)     # Same for data frames, output is also a data.frame
dapply(mtcars, quantile)

# Let's now take a more complex classed object, like a dplyr grouped tibble
library(dplyr)
gmtcars <- group_by(mtcars,cyl,vs,am)
dapply(gmtcars, log)                     # We can still take logs of rows or columns and get a grouped tibble back
dapply(gmtcars, log, MARGIN = 1)
dapply(gmtcars, quantile, MARGIN = 1)    # Also works for quantiles
dapply(gmtcars, log, return = "matrix")  # Output as matrix also works, groups are attached!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
