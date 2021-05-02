\name{AA2-small-helpers}
\alias{AA2-small-helpers}
\alias{.c}
\alias{cinv}
\alias{vlabels}
\alias{vlabels<-}
\alias{vclasses}
\alias{vtypes}
\alias{namlab}
\alias{add_stub}
\alias{rm_stub}
\alias{\%!in\%}
\alias{ckmatch}
\alias{fnlevels}
\alias{fnrow}
\alias{fncol}
\alias{fdim}
\alias{allNA}
\alias{missing_cases}
\alias{na_rm}
\alias{na_omit}
\alias{na_insert}
\alias{seq_row}
\alias{seq_col}
\alias{all_identical}
\alias{all_obj_equal}
\alias{setRownames}
\alias{setColnames}
\alias{setDimnames}
\alias{unattrib}
\alias{setAttrib}
\alias{copyAttrib}
\alias{copyMostAttrib}
\alias{is.categorical}
\alias{is.Date}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Small (Helper) Functions % \emph{collapse}
}
\description{
Convenience functions in the \emph{collapse} package that help to deal with attributes such as variable names and labels, missing values, matching and object checking etc.. Some functions are performance improved replacements for base R functions. % For recoding and replacing values see \code{\link{Recode}}. % for pairwise correlations and covariances see \code{\link{pwcor}}, for summary statistics see \code{\link{qsu}}.
}
\usage{
.c(\dots)                       # Non-standard concatenation i.e. .c(a, b) == c("a", "b")
vlabels(X, attrn = "label")   # Get labels of variables in X, in attr(X[[i]], attrn)
vlabels(X, attrn = "label") <- value      # Set labels of variables in X
vclasses(X)                   # Get classes of variables in X
vtypes(X)                     # Get data storage types of variables in X (calling typeof)
namlab(X, class = FALSE,      # Return data frame of names, labels and classes
          attrn = "label")
add_stub(X, stub, pre = TRUE) # Add a stub (i.e. prefix or postfix) to column names
rm_stub(X, stub, pre = TRUE)  # Remove stub from column names
x \%!in\% table                 # The opposite of \%in\%
ckmatch(x, table,             # Check-match: throws an informative error if non-matched
     e = "Unknown columns:")
fnlevels(x)                   # Faster version of nlevels(x) (for factors)
fnrow(X)                      # Faster nrow for data frames (not faster for matrices)
fncol(X)                      # Faster ncol for data frames (not faster for matrices)
fdim(X)                       # Faster dim for data frames (not faster for matrices)
allNA(x)                      # Check for all missing cases in a vector
missing_cases(X, cols = NULL) # Find missing cases in a matrix / data framme
na_rm(x)                      # Remove missing values from vector and return vector
na_omit(X, cols = NULL,       # Faster na.omit for matrices and data frames
        na.attr = FALSE)
na_insert(X, prop = 0.1)      # Insert missing values at random in vectors, matrices DF's
cinv(X)                       # Choleski (fast) inverse of symmetric, positive def. matrix
all_identical(\dots)            # Check exact equality of multiple objects or list-elements
all_obj_equal(\dots)            # Check near equality of multiple objects or list-elements
seq_row(X)                    # Fast integer sequences along rows of X
seq_col(X)                    # Fast integer sequences along columns of X
setRownames(object,
 nm = if(is.atomic(object))   # Set rownames of object and return object
 seq_row(object) else NULL)
setColnames(object, nm)       # Set colnames of object and return object
setDimnames(object, dn,
            which = NULL)     # Set dimension names of object and return object
unattrib(object)              # Remove all attributes from object
setAttrib(object, a)          # Replace all attributes with list of attributes 'a'
copyAttrib(to, from)          # Copy all attributes from object 'from' to object 'to'
copyMostAttrib(to, from)      # Copy most attributes from object 'from' to object 'to'
is.categorical(x)             # The opposite of is.numeric
is.Date(x)                    # Check if object is of class "Date", "POSIXlt" or "POSIXct"

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or data frame (some functions also support vectors and arrays although that is less common).}
  \item{x, table}{a atomic vector.}
  \item{object, to, from}{a suitable R object.}
  \item{a}{a suitable list of attributes.}
  \item{attrn}{character. Name of attribute to store labels or retrieve labels from.}
  \item{value}{a matching character vector of variable labels.}
  \item{class}{logical. Also show the classes of variables in X in a column?}
  \item{stub}{a single character stub, i.e. "log.", which by default will be pre-applied to all variables or column names in X.}
  \item{pre}{logical. \code{FALSE} will post-apply \code{stub}.}
  \item{cols}{only removes rows with missing values on these columns. Columns can be selected using column names, indices, a logical vector or a selector function (e.g. \code{is.numeric}).}
  \item{na.attr}{logical. \code{TRUE} adds an attribute containing the removed cases. For compatibility reasons this is exactly the same format as \code{na.omit} i.e. the attribute is called "na.action" and of class "omit".}
  \item{nm}{a suitable vector of row- or column-names.}
  \item{dn}{a suitable vector or list of names for dimension(s).}
  \item{which}{integer. If \code{NULL}, \code{dn} has to be a list fully specifying the dimension names of the object. Alternatively, a vector or list of names for dimensions \code{which} can be supplied. See Examples. }
  \item{prop}{double. Specify the proportion of observations randomly replaced with \code{NA}.}
  \item{e}{the error message thrown by \code{ckmatch} for non-matched elements. The message is followed by the comma-separated non-matched elements.}
 \item{\dots}{for \code{.c}: Comma-separated expressions. For \code{all_identical / all_obj_equal}: Either multiple comma-separated objects or a single list of objects in which all elements will be checked for exact / numeric equality.}
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
\link[=collapse-documentation]{Collapse Overview} %% \code{\link{pwcor}}, \code{\link{qsu}}
}
\examples{
## Non-standard concatenation
.c(a, b, "c d", e == f)

## Variable labels
namlab(wlddev, class = TRUE)
vlabels(wlddev)
vlabels(wlddev) <- vlabels(wlddev)

## Stub-renaming
log_mtc <- add_stub(log(mtcars), "log.")
head(log_mtc)
head(rm_stub(log_mtc, "log."))
rm(log_mtc)

## Setting dimension names of an object
head(setRownames(mtcars))
ar <- array(1:9, c(3,3,3))
setRownames(ar)
setColnames(ar, c("a","b","c"))
setDimnames(ar, c("a","b","c"), which = 3)
setDimnames(ar, list(c("d","e","f"), c("a","b","c")), which = 2:3)
setDimnames(ar, list(c("g","h","i"), c("d","e","f"), c("a","b","c")))

## Checking exact equality of multiple objects
all_identical(iris, iris, iris, iris)
l <- replicate(100, fmean(num_vars(iris), iris$Species), simplify = FALSE)
all_identical(l)
rm(l)

## Missing values
mtc_na <- na_insert(mtcars, 0.15)    # Set 15\% of values missing at random
fNobs(mtc_na)                        # See observation count
na_omit(mtc_na)                      # 12x faster than na.omit(mtc_na)
na_omit(mtc_na, na.attr = TRUE)      # Adds attribute with removed cases, like na.omit
na_omit(mtc_na, cols = c("vs","am")) # Removes only cases missing vs or am
na_omit(qM(mtc_na))                  # Also works for matrices
na_omit(mtc_na$vs, na.attr = TRUE)   # Also works with vectors
na_rm(mtc_na$vs)                     # For vectors na_rm is faster ...
rm(mtc_na)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{attribute}
\keyword{misc}
\keyword{documentation}
