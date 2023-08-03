\name{list-processing}
\alias{A8-list-processing}
\alias{list-processing}
\title{List Processing} % \emph{collapse}

\description{
\emph{collapse} provides the following set of functions to efficiently work with lists of R objects:
\itemize{
\item \bold{Search and Identification}\itemize{
\item \code{\link{is_unlistable}} checks whether a (nested) list is composed of atomic objects in all final nodes, and thus unlistable to an atomic vector using \code{\link{unlist}}.
\item \code{\link{ldepth}} determines the level of nesting of the list (i.e. the maximum number of nodes of the list-tree).
\item \code{\link{has_elem}} searches elements in a list using element names, regular expressions applied to element names, or a function applied to the elements, and returns \code{TRUE} if any matches were found.
}

\item \bold{Subsetting} \itemize{
\item \code{\link{atomic_elem}} examines the top-level of a list and returns a sublist with the atomic elements. Conversely \code{\link{list_elem}} returns the sublist of elements which are themselves lists or list-like objects.
 \item \code{\link{reg_elem}} and \code{\link{irreg_elem}} are recursive versions of the former. \code{\link{reg_elem}} extracts the 'regular' part of the list-tree leading to atomic elements in the final nodes, while \code{\link{irreg_elem}} extracts the 'irregular' part of the list tree leading to non-atomic elements in the final nodes. (\emph{Tip}: try calling both on an \code{lm} object). Naturally for all lists \code{l}, \code{is_unlistable(reg_elem(l))} evaluates to \code{TRUE}.
\item \code{\link{get_elem}} extracts elements from a list using element names, regular expressions applied to element names, a function applied to the elements, or element-indices used to subset the lowest-level sub-lists. by default the result is presented as a simplified list containing all matching elements. With the \code{keep.tree} option however \code{\link{get_elem}} can also be used to subset lists i.e. maintain the full tree but cut off non-matching branches.
}

\item \bold{Splitting and Transposition} \itemize{
\item \code{\link{rsplit}} recursively splits a vector or data frame into subsets according to combinations of (multiple) vectors / factors - by default returning a (nested) list. If \code{flatten = TRUE}, the list is flattened yielding the same result as \code{\link{split}}. \code{rsplit} is also faster than \code{\link{split}}, particularly for data frames.

\item \code{\link{t_list}} efficiently transposes nested lists of lists, such as those obtained from splitting a data frame by multiple variables using \code{\link{rsplit}}.
}

\item \bold{Apply Functions} \itemize{
\item \code{\link{rapply2d}} is a recursive version of \code{\link{lapply}} with two key differences to \code{\link{rapply}} to apply a function to nested lists of data frames or other list-based objects.
}

\item \bold{Unlisting / Row-Binding} \itemize{
\item \code{\link{unlist2d}} efficiently unlists unlistable lists in 2-dimensions and creates a data frame (or \emph{data.table}) representation of the list. This is done by recursively flattening and row-binding R objects in the list while creating identifier columns for each level of the list-tree and (optionally) saving the row-names of the objects in a separate column. \code{\link{unlist2d}} can thus also be understood as a recursive generalization of \code{do.call(rbind, l)}, for lists of vectors, data frames, arrays or heterogeneous objects. A simpler version for non-recursive row-binding lists of lists / data.frames, is also available by \code{\link{rowbind}}.
}

}
}
\section{Table of Functions}{
  \tabular{lll}{\emph{ Function }  \tab\tab \emph{ Description }  \cr
  % \code{\link{is.regular}} \tab\tab \code{function(x) is.atomic(x) || is.list(x)} \cr
  \code{\link{is_unlistable}} \tab\tab Checks if list is unlistable \cr
  \code{\link{ldepth}} \tab\tab Level of nesting / maximum depth of list-tree \cr
  \code{\link{has_elem}} \tab\tab Checks if list contains a certain element \cr
  \code{\link{get_elem}} \tab\tab Subset list / extract certain elements \cr
  \code{\link{atomic_elem}} \tab\tab Top-level subset atomic elements \cr
  \code{\link{list_elem}} \tab\tab Top-level subset list/list-like elements \cr
  \code{\link{reg_elem}} \tab\tab Recursive version of \code{atomic_elem}: Subset / extract 'regular' part of list \cr
  \code{\link{irreg_elem}} \tab\tab Subset / extract non-regular part of list \cr
  \code{\link{rsplit}} \tab\tab Recursively split vectors or data frames / lists \cr
  \code{\link{t_list}} \tab\tab Transpose lists of lists \cr
  \code{\link{rapply2d}} \tab\tab Recursively apply functions to lists of data objects \cr
  \code{\link{unlist2d}} \tab\tab Recursively unlist/row-bind lists of data objects in 2D, to data frame or \emph{data.table} \cr
  \code{\link{rowbind}} \tab\tab Non-recursive binding of lists or lists / data.frames. \cr
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}
}
\keyword{list}
\keyword{manip}
\keyword{documentation}


