\name{small-helpers}
\alias{AA3-small-helpers}
\alias{small-helpers}
\alias{.c}
\alias{vlabels}
\alias{vlabels<-}
\alias{setLabels}
\alias{vclasses}
\alias{namlab}
\alias{add_stub}
\alias{rm_stub}
\alias{\%!in\%}
\alias{massign}
\alias{\%=\%}
\alias{ckmatch}
\alias{all_identical}
\alias{all_obj_equal}
\alias{setRownames}
\alias{setColnames}
\alias{setDimnames}
\alias{unattrib}
\alias{setAttrib}
\alias{copyAttrib}
\alias{copyMostAttrib}
\alias{is_categorical}
\alias{is_date}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Small (Helper) Functions % \emph{collapse}
}
\description{
Convenience functions in the \emph{collapse} package that help to deal with object attributes such as variable names and labels, matching and object checking, and that improve the workflow. % Some functions are performance improved replacements for base R functions. % For recoding and replacing values see \code{\link{Recode}}. % for pairwise correlations and covariances see \code{\link{pwcor}}, for summary statistics see \code{\link{qsu}}.
}
\usage{
.c(\dots)                       # Non-standard concatenation i.e. .c(a, b) == c("a", "b")
nam \%=\% values                # Multiple-assignment e.g. .c(x, y) \%=\% c(1, 2),
massign(nam, values,          # can also assign to different environment.
        envir = parent.frame())
vlabels(X, attrn = "label",   # Get labels of variables in X, in attr(X[[i]], attrn)
        use.names = TRUE)
vlabels(X, attrn = "label") <- value    # Set labels of variables in X (by reference)
setLabels(X, value, attrn = "label",    # Set labels of variables in X (by reference)
          cols = NULL)                  # and return X
vclasses(X, use.names = TRUE) # Get classes of variables in X
namlab(X, class = FALSE,      # Return data frame of names and labels,
  attrn = "label", N = FALSE, # and (optionally) classes, number of observations
  Ndistinct = FALSE)          # and number of non-missing distinct values
add_stub(X, stub, pre = TRUE, # Add a stub (i.e. prefix or postfix) to column names
         cols = NULL)
rm_stub(X, stub, pre = TRUE,  # Remove stub from column names, also supports general
        regex = FALSE,        # regex matching and removing of characters
        cols = NULL, ...)
x \%!in\% table                 # The opposite of \%in\%
ckmatch(x, table,             # Check-match: throws an informative error if non-matched
     e = "Unknown columns:")
all_identical(\dots)            # Check exact equality of multiple objects or list-elements
all_obj_equal(\dots)            # Check near equality of multiple objects or list-elements
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
is_categorical(x)             # The opposite of is.numeric
is_date(x)                    # Check if object is of class "Date", "POSIXlt" or "POSIXct"

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or data frame (some functions also support vectors and arrays although that is less common).}
  \item{x, table}{a (atomic) vector.}
  \item{object, to, from}{a suitable R object.}
  \item{a}{a suitable list of attributes.}
  \item{attrn}{character. Name of attribute to store labels or retrieve labels from.}
  \item{N, Ndistinct}{logical. Options to display the number of observations or number of distinct non-missing values.}
  \item{value}{for \code{whichv} and \code{alloc}: a single value of any vector type. For \code{vlabels<-} and \code{setLabels}: a matching character vector or list of variable labels. }
  \item{use.names}{logical. Preserve names if \code{X} is a list. }
  \item{cols}{integer. (optional) indices of columns to apply the operation to. Note that for these small functions this needs to be integer, whereas for other functions in the package this argument is more flexible. }
  \item{class}{logical. Also show the classes of variables in X in a column?}
  \item{stub}{a single character stub, i.e. "log.", which by default will be pre-applied to all variables or column names in X.}
  \item{pre}{logical. \code{FALSE} will post-apply \code{stub}.}
  \item{regex}{logical. Match pattern anywhere in names using a regular expression and remove it with \code{\link{gsub}}.}
  \item{nm}{a suitable vector of row- or column-names.}
  \item{dn}{a suitable vector or list of names for dimension(s).}
  \item{which}{integer. If \code{NULL}, \code{dn} has to be a list fully specifying the dimension names of the object. Alternatively, a vector or list of names for dimensions \code{which} can be supplied. See Examples. }
  \item{e}{the error message thrown by \code{ckmatch} for non-matched elements. The message is followed by the comma-separated non-matched elements.}
  \item{nam}{character. A vector of object names.}
  \item{values}{a matching atomic vector or list of objects.}
  \item{envir}{the environment to assign into.}
 \item{\dots}{for \code{.c}: Comma-separated expressions. For \code{all_identical / all_obj_equal}: Either multiple comma-separated objects or a single list of objects in which all elements will be checked for exact / numeric equality. For \code{rm_stub}: further arguments passed to \code{\link{gsub}}.}
}
\details{
  \code{copyAttrib} and \code{copyMostAttrib} take a shallow copy of the attribute list, i.e. they don't duplicate in memory the attributes themselves. They also, along with \code{setAttrib}, take a shallow copy of lists passed to the \code{to} argument, so that lists are not modified by reference. Atomic \code{to} arguments are however modified by reference.

\code{copyMostAttrib} copies all attributes except for \code{"names"}, \code{"dim"} and \code{"dimnames"} (like the corresponding C-API function), and further only copies the \code{"row.names"} attribute of data frames if known to be valid. Thus it is a suitable choice if objects should be of the same type but are not of equal dimensions.
}
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
\link[=efficient-programming]{Efficient Programming}, \link[=collapse-documentation]{Collapse Overview} %% \code{\link{pwcor}}, \code{\link{qsu}}
}
\examples{
## Non-standard concatenation
.c(a, b, "c d", e == f)

## Multiple assignment
.c(a, b) \%=\% list(1, 2)
.c(T, N) \%=\% dim(EuStockMarkets)
names(iris) \%=\% iris
list2env(iris)          # Same thing
rm(list = c("a", "b", "T", "N", names(iris)))

## Variable labels
namlab(wlddev)
namlab(wlddev, class = TRUE, N = TRUE, Ndistinct = TRUE)
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{attribute}
\keyword{misc}
\keyword{documentation}
