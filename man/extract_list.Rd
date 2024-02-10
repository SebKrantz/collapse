\name{get_elem}
% \alias{extract_list}
\alias{atomic_elem}
\alias{atomic_elem<-}
\alias{list_elem}
\alias{list_elem<-}
\alias{reg_elem}
\alias{irreg_elem}
\alias{has_elem}
\alias{get_elem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find and Extract / Subset List Elements}
\description{
A suite of functions to subset or extract from (potentially complex) lists and list-like structures. Subsetting may occur according to certain data types, using identifier functions, element names or regular expressions to search the list for certain objects.

\itemize{
\item \code{atomic_elem} and \code{list_elem} are non-recursive functions to extract and replace the atomic and sub-list elements at the top-level of the list tree.
\item \code{reg_elem} is the recursive equivalent of \code{atomic_elem} and returns the 'regular' part of the list - with atomic elements in the final nodes. \code{irreg_elem} returns all the non-regular elements (i.e. call and terms objects, formulas, etc\dots). See Examples.
\item \code{get_elem} returns the part of the list responding to either an identifier function, regular expression, exact element names or indices applied to all final objects. \code{has_elem} checks for the existence of an element and returns \code{TRUE} if a match is found. See Examples.
}
}
\usage{
## Non-recursive (top-level) subsetting and replacing
atomic_elem(l, return = "sublist", keep.class = FALSE)
atomic_elem(l) <- value
list_elem(l, return = "sublist", keep.class = FALSE)
list_elem(l) <- value

## Recursive separation of regular (atomic) and irregular (non-atomic) parts
reg_elem(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE)
irreg_elem(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE)

## Extract elements / subset list tree
get_elem(l, elem, recursive = TRUE, DF.as.list = FALSE, keep.tree = FALSE,
         keep.class = FALSE, regex = FALSE, invert = FALSE, \dots)

## Check for the existence of elements
has_elem(l, elem, recursive = TRUE, DF.as.list = FALSE, regex = FALSE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{value}{a list of the same length as the extracted subset of \code{l}.}
  \item{elem}{a function returning \code{TRUE} or \code{FALSE} when applied to elements of \code{l}, or a character vector of element names or regular expressions (if \code{regex = TRUE}). \code{get_elem} also supports a vector or indices which will be used to subset all final objects.}
  \item{return}{an integer or string specifying what the selector function should return. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "sublist"   \tab\tab subset of list (default) \cr
                 2 \tab\tab "names" \tab\tab column names \cr
                 3 \tab\tab "indices"   \tab\tab column indices \cr
                 4 \tab\tab "named_indices"  \tab\tab named column indices \cr
                 5 \tab\tab "logical"   \tab\tab logical selection vector \cr
                 6 \tab\tab "named_logical"     \tab\tab named logical vector \cr
  }
  \emph{Note}: replacement functions only replace data, names are replaced together with the data.
  }


  \item{recursive}{logical. Should the list search be recursive (i.e. go though all the elements), or just at the top-level?}
  \item{DF.as.list}{logical. \code{TRUE} treats data frames like (sub-)lists; \code{FALSE} like atomic elements.}
  \item{keep.tree}{logical. \code{TRUE} always returns the entire list tree leading up to all matched results, while \code{FALSE} drops the top-level part of the tree if possible.}
  \item{keep.class}{logical. For list-based objects: should the class be retained? This only works if these objects have a \code{[} method that retains the class.}
  \item{regex}{logical. Should regular expression search be used on the list names, or only exact matches?}
  \item{invert}{logical. Invert search i.e. exclude matched elements from the list?}
  \item{\dots}{further arguments to \code{grep} (if \code{regex = TRUE}).}
}
\details{
For a lack of better terminology, \emph{collapse} defines 'regular' R objects as objects that are either atomic or a list. \code{reg_elem} with \code{recursive = TRUE} extracts the subset of the list tree leading up to atomic elements in the final nodes. This part of the list tree is unlistable - calling \code{is_unlistable(reg_elem(l))} will be \code{TRUE} for all lists \code{l}. Conversely, all elements left behind by \code{reg_elem} will be picked up be \code{irreg_elem}. Thus \code{is_unlistable(irreg_elem(l))} is always \code{FALSE} for lists with irregular elements (otherwise \code{irreg_elem} returns an empty list). \cr

If \code{keep.tree = TRUE}, \code{reg_elem}, \code{irreg_elem} and \code{get_elem} always return the entire list tree, but cut off all of the branches not leading to the desired result. If \code{keep.tree = FALSE}, top-level parts of the tree are omitted as far as possible. For example in a nested list with three levels and one data-matrix in one of the final branches, \code{get_elem(l, is.matrix, keep.tree = TRUE)} will return a list (\code{lres}) of depth 3, from which the matrix can be accessed as \code{lres[[1]][[1]][[1]]}. This however does not make much sense. \code{get_elem(l, is.matrix, keep.tree = FALSE)} will therefore figgure out that it can drop the entire tree and return just the matrix. \code{keep.tree = FALSE} makes additional optimizations if matching elements are at far-apart corners in a nested structure, by only preserving the hierarchy if elements are above each other on the same branch. Thus for a list \code{l <- list(list(2,list("a",1)),list(1,list("b",2)))} calling \code{get_elem(l, is.character)} will just return \code{list("a","b")}.
}
% \value{

% }
\seealso{
 \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
m <- qM(mtcars)
get_elem(list(list(list(m))), is.matrix)
get_elem(list(list(list(m))), is.matrix, keep.tree = TRUE)

l <- list(list(2,list("a",1)),list(1,list("b",2)))
has_elem(l, is.logical)
has_elem(l, is.numeric)
get_elem(l, is.character)
get_elem(l, is.character, keep.tree = TRUE)

l <- lm(mpg ~ cyl + vs, data = mtcars)
str(reg_elem(l))
str(irreg_elem(l))
get_elem(l, is.matrix)
get_elem(l, "residuals")
get_elem(l, "fit", regex = TRUE)
has_elem(l, "tol")
get_elem(l, "tol")

}
\keyword{list}
\keyword{manip}
