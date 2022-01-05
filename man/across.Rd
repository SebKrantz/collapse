\name{across}
\alias{across}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply Functions Across Muliple Columns
}
\description{
\code{across()} can be used inside \code{\link{fmutate}} and \code{\link{fsummarise}} to apply one or more functions to a selection of columns. It is overall very similar to \code{dplyr::across}, but does not support some \code{rlang} specific features. It has some additional features (arguments), and is optimized to work with \emph{collapses}'s, \code{\link{.FAST_FUN}}, yielding much faster computations.
}
\usage{
across(.cols = NULL, .fns, ..., .names = NULL,
       .apply = "auto", .transpose = "auto")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.cols}{select columns using column names and expressions (e.g. \code{a:b} or \code{c(a, b, c:f)}), column indices, logical vectors, or functions yielding a logical value e.g. \code{is.numeric}. \code{NULL} applies functions to all columns except for grouping columns.}
  \item{.fns}{A fuction, character vector of functions or list of functions. Vectors / lists can be named to yield alternative names in the result (see \code{.names}). This argument is evaluated inside \code{substitute()}, and the content (not the names of vectors/lists) is checked against \code{.FAST_FUN} and \code{.OPERATOR_FUN}. Matching functions receive vectorized execution, other functions are applied to the data in a standard way.}
  \item{\dots}{further arguments to \code{.fns}. Arguments are evaluated in the data environment and split by groups as well (for non-vectorized functions if of the same length as the data).}
  \item{.names}{controlls the naming of computed columns. \code{NULL} generates names of the form \code{coli_funj} if multiple functions are used. \code{.names = TRUE} enables this for a single function, \code{.names = FALSE} disables it for multiple functions (sensible for functions such as \code{.OPERATOR_FUN} that rename columns (if \code{.apply = FALSE})). It is also possible to supply a function with two arguments for column and function names e.g. \code{function(c, f) paste0(f, "_", c)}. Finally, you can supply a custom vector of names which must match \code{length(.cols) * length(.fns)}.}
  \item{.apply}{controls whether functions are applied column-by-column (\code{TRUE}) or to multiple columns at once (\code{FALSE}). The default, \code{"auto"}, does the latter for vectorized functions, which have an efficient data frame method. It can also be sensible to use \code{.apply = FALSE} for non-vectorized functions, especially multivariate functions like \code{\link{lm}} or \code{\link{pwcor}}, or functions renaming the data. See Examples.  }
  \item{.transpose}{with multiple \code{.fns}, \code{.transpose} controls whether the result is ordered first by column, then by function (\code{TRUE}), or vice-versa (\code{FALSE}). \code{"auto"} does the former if all functions yield results of the same dimensions (dimensions may differ if \code{.apply = FALSE}). See Examples.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (.cols = NULL, .fns, ..., .names = NULL, .apply = "auto",
    .transpose = "auto")
{
    stop("across() can only work inside fmutate() and fsummarise()")
  }
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
