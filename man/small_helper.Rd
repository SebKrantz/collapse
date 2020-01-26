\name{AA2-small-helpers}
\alias{AA2-small-helpers}
\alias{vlabels}
\alias{vlabels<-}
\alias{namlab}
\alias{add_stub}
\alias{seq_row}
\alias{seq_col}
\alias{setRownames}
\alias{setColnames}
\alias{setDimnames}
\alias{pwcor}
\alias{pwcov}
\alias{all_identical}
\alias{is.categorical}
\alias{is.Date}
\alias{\%!in\%}
\alias{na_rm}
\alias{na_insert}
\alias{fnlevels}
\alias{funique}
\alias{as.numeric_factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{collapse} Small (Helper) Functions
}
\description{
Elementary exported convenience functions in the \code{collapse} package, that help to deal with variable names, labels, missing values, matching and object checking etc... For recoding and replacing values see \code{\link{Recode}}.
}
\usage{
vlabels(X, attrn = "label")   # Get labels of variables in X, in attr(X[[i]], attrn)
vlabels(X, attrn = "label") <- value      # Set labels of variables in X
namlab(X, class = FALSE, attrn = "label") # Return data.frame of names, labels and classes
add_stub(X, stub, pre = TRUE) # Add a stub (i.e. prefix or postfix) to variables
x \%!in\% table                 # The opposite of \%in\%
fnlevels(x)                   # Faster version of nlevels(x) (for factors)
funique(x, ordered = TRUE)    # Faster unique(x) and sort(unique(x)) for vectors
na_rm(x)                      # Remove missing values from vector and return vector
na_insert(X, prop = 0.1)      # Insert missing values at random in vectors, matrices DF's
all_identical(...)            # Check exact equality of multiple objects or list-elements
pwcor(X, ...)                 # Pairwise-correlations
pwcov(X, ...)                 # Pairwise-covariances
seq_row(X)                    # Integer sequences along rows of X
seq_col(X)                    # Integer sequences along columns of X
setRownames(object = nm,      # Set rownames of object and return object
        nm = seq_row(object))
setColnames(object = nm, nm)  # Set colnames of object and return object
setDimnames(object = dn, dn)  # Set dimension names of object and return object
is.categorical(x)             # The opposite of is.numeric
is.Date(x)                    # Check if object is of class "Date", "POSIXlt" or "POSIXct"
as.numeric_factor(X)          # Coerce a factor, or all factor columns in X, to numeric
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a (numeric) matrix or data.frame.}
  \item{object}{an array or data.frame.}
  \item{x, table}{(numeric) vectors.}
  \item{attrn}{Character. Name of attribute to store labels or retrieve labels from.}
  \item{value}{a matching character vector of variable labels.}
  \item{class}{logical. Also show the classes of variables in X in a column?}
  \item{stub}{a single character stub, i.e. "log.", which by default will be pre-applied to all variable or column names in X.}
  \item{pre}{logical. \code{FALSE} will post-apply \code{stub}.}
  \item{nm}{a suitable vector of row- or column-names.}
  \item{dn}{a suitable list of dimension names.}
  \item{ordered}{logical. \code{TRUE} sorts the output, \code{FALSE} is slightly faster.}
  \item{prop}{specify the proportion of observations randomly replaced with \code{NA}.}
  \item{...}{for \code{pwcor} and \code{pwcov}: Other arguments passed to \code{stats::cor} and \code{stats::cov}. For \code{all_identical}, either multiple comma-separated objects or a single list of objects.}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
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
\code{\link{Recode}}, \link[=collapse-documentation]{Collapse Overview} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Variable labels
namlab(wlddev, class = TRUE)
vlabels(wlddev)
vlabels(wlddev) <- vlabels(wlddev)

## Checking exact equality of multiple objects
all_identical(iris, iris, iris, iris)
l <- replicate(100, fmean(num_vars(iris), iris$Species), simplify = FALSE)
all_identical(l)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{attribute}
\keyword{misc}
\keyword{documentation}
