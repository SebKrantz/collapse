\name{Small Helper}
\alias{small_helper}
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
\alias{all.identical}
\alias{is.categorical}
\alias{is.Date}
\alias{\%!in\%}
\alias{na.rm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\code{collapse} Small (Helper) Functions
}
\description{
Elementary exported convenience functions in the \code{collapse} package, that help to deal with variable names, labels, missing values, matching and object checking etc...
}
\usage{
vlabels(X)                    # Get labels of variables in X, stored in attr(X[[i]], "label")
vlabels(X) <- value           # Set labels of variables in X
namlab(X, class = FALSE)      # Return a data.frame of variable names, labels and (optionally) classes
add_stub(X, stub, pre = TRUE) # Add a stub (i.e. prefix or postfix) to the variables / column names of X and return X
seq_row(X)                    # Integer sequences along rows of X
seq_col(X)                    # Integer sequences along columns of X
setRownames(object = nm,      # Set rownames of object and return object. If nm is left unspecified, set an integer sequence
        nm = seq_row(object))
setColnames(object = nm, nm)  # Set colnames of object and return object
setDimnames(object = dn, dn)  # Set dimension names of object and return object
pwcor(X, ...)                 # Pairwise-correlations
pwcov(X, ...)                 # Pairwise-covariances
all.identical(...)            # Check the exact equality of multiple objects or of all elements in a single list
is.categorical(x)             # The opposite of is.numeric
is.Date(x)                    # Check if an object inherits a class "Date", "POSIXlt" or "POSIXct"
x \%!in\% y                     # The opposite of \%in\%
na.rm(x)                      # Remove missing values from vector and return vector
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a (numeric) matrix or data.frame.}
  \item{object}{an array or data.frame.}
  \item{x, y}{(numeric) vectors.}
  \item{class}{logical. Also show the classes of variables in X in a column?}
  \item{stub}{a single character stub, i.e. "log.", which by default will be pre-applied to all variable or column names in X.}
  \item{pre}{logical. \code{FALSE} will post-apply \code{stub}.}
  \item{nm}{a suitable vector of row- or column-names.}
  \item{dn}{a suitable list of dimension names.}
  \item{...}{for \code{pwcor} and \code{pwcov}: Other arguments passed to \code{stats::cor} and \code{stats::cov}. For \code{all.identical}, either multiple comma-separeted objects or a single list of objects.}
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
## Variable labels
namlab(wlddev, class = TRUE)
vlabels(wlddev)
vlabels(wlddev) <- vlabels(wlddev)

## Checking exact equality of multiple objects
all.identical(iris, iris, iris, iris)
l <- replicate(100, fmean(num_vars(iris), iris$Species), simplify = FALSE)
all.identical(l)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
