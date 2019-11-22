\name{extract_list}
\alias{extract_list}
\alias{atomic_elem}
\alias{list_elem}
\alias{reg_elem}
\alias{irreg_elem}
\alias{has_elem}
\alias{get_elem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract List Elements by Type, Function or Regular Expression}
\description{
A suite of functions to subset or extract from (potentially complex) lists and list-like structures. Subsetting may occur according to certain data types, or using identifier functions or regular expresssions to search the list for certain objects.

\itemize{
\item \code{atomic.elem} and \code{list.elem} are non-recursive functions to extract and replace the atomic and sub-list elements at the top-level of the list tree.
\item \code{reg.elem} is the recursive equivalent of \code{atomic.elem} and returns the 'regular' part of the list - with atomic elements in the final nodes. See \code{\link{is.regular}} and \code{\link{is.unlistable}}. \code{irreg.elem} returns all the non-regular elements (i.e. call and terms objects, formulas, etc...). See Examples.
\item \code{get.elem} returns the part of the list responding to either an identifier function, regular expression or exact element name. \code{has.elem} checks for the existence of the searched element and returns \code{TRUE} if a match is found. See Examples.
}
}
\usage{
## Non-recursive (top-level) subsetting and replacing
atomic.elem(l)
atomic.elem(l) <- y
list.elem(l)
list.elem(l) <- y

## Recursive separation of regular (atomic) and irregular (non-atomic) parts
reg.elem(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE)
irreg.elem(l, recursive = TRUE, keep.tree = FALSE, keep.class = FALSE)

## Extract elements using a function or regular expression
get.elem(l, FoR, recursive = TRUE, DF.as.list = TRUE, keep.tree = FALSE,
         keep.class = FALSE, regex = TRUE, ...)

## Check for the existence of elements
has.elem(l, FoR, recursive = TRUE, DF.as.list = TRUE, regex = TRUE, ...)
#
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a list.}
  \item{y}{a list of the same length as the extracted subet of \code{l}.}
  \item{FoR}{a function returning TRUE or FALSE when applied to elements of \code{l}, or a character vector of regular expressions.}
  \item{recursive}{should the list search be recursive (i.e. go though all the elements), or just at the top-level?}
  \item{DF.as.list}{treat data.frame's like (sub-)lists or like atomic elements?}
  \item{keep.tree}{\code{TRUE} always returns the entire list tree leading up to all matched results, while \code{FALSE} drops the top-level part of the tree if possible.}
  \item{keep.class}{for classed objects: Should the class be ratained?}
  \item{regex}{should regular expression search be used on the list names, or only exact matches?}
  \item{...}{further arguments to \code{grep}.}
}
\details{
A list is made up of regular and irregular elements. I defined regular elements as all elements that are either atomic or a list (see \code{\link{is.regular}}). \code{reg.elem} with \code{recursive = TRUE} therefore extracts the subset of the list tree leading up to atomic elements in the final nodes. This part of the list tree is unlistable - calling \code{is.unlistable(reg.elem(l))} will be \code{TRUE} for all lists \code{l}. Conversely, all elements left behind by \code{reg.elem} will be picked up be \code{irreg.elem} (if available). Thus \code{is.unlistable(reg.elem(l))} is always \code{FALSE} for lists with irregular elements (otherwise \code{irreg.elem} returns an empty list). \cr

If \code{keep.tree = TRUE}, \code{reg.elem}, \code{irreg.elem} and \code{get.elem} always return the entire list tree, but cut off all of the branches not leading to the desired result. If \code{keep.tree = FALSE}, top-level parts of the tree are omitted so far this is possible. For example in a nested list with three levels and one data-matrix in one of the final branches, \code{get.elem(l, is.matrix, keep.tree = TRUE)} will return a list (\code{lres}) of depth 3, from which the matrix can be accessed as \code{lres[[1]][[1]][[1]]}. This however does not make much sense. \code{get.elem(l, is.matrix, keep.tree = FALSE)} will therefore figgure out that it can drop the entire tree and return just the matrix. \code{keep.tree = FALSE} makes additional optimizations if matching elements are at far-apart corners in a nested structure, by only perserving the hierarchy if elements are above each other on the same branch. Thus for a list \code{l <- list(list(2,list("a",1)),list(1,list("b",2)))} calling \code{get.elem(l,is.character)} will just return \code{list("a","b")}.
}
\value{

}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
