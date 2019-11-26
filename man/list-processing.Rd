\name{list-processing}
\alias{list-processing}
\title{List Processing}

\description{
\code{collapse} provides the following set of functions to work with lists of R objects:
\itemize{
\item \bold{Identification:} \code{\link{is.regular}} checks whether an R object is either atomic or a list. A (nested) list composed of regular objects at each level of the list-tree is unlistable to an atomic vector, checked by \code{\link{is.unlistable}}. In addition, \code{\link{ldepth}} determines the level of nesting of the list (i.e. the maximum number of nodes of the list-tree). \code{\link{has_elem}}

\item \bold{Subsetting:} \itemize{
\item \code{\link{atomic_elem}} examines the top-level of a list and returns a sublist with the atomic elements. Conversely \code{\link{list_elem}} returns the sublist of elements which are themselves lists or list-like objects.
 \item \code{\link{reg_elem}} and \code{\link{irreg_elem}} are recursive versions of the former. \code{\link{reg_elem}} extracts the regular part of the list-tree (leading to atomic elements in the final nodes), while \code{\link{irreg_elem}} extracts the 'irregular' part of the list tree leading to non-atomic elements in the final nodes. (\emph{Tipp}: try calling both on an \code{lm} object). Naturally for all lists \code{l}, \code{is.unlistable(reg_elem(l))} evaluates to \code{TRUE}...
\item \code{\link{get_elem}},
}

\item \bold{Apply Functions:} \code{\link{rapply2d}} is a recursive version of \code{base::lapply} with two key differences to \code{base::rapply}: (1) Data frames are considered as atomic objects, not as (sub-)lists, and (2) the result is not simplified.

\item \bold{Unlisting / Row-Binding:} \code{\link{unlist2d}} efficiently unlists unlistable lists in 2-dimensions and creates a data.frame representation of the list (unlike \code{base::unlist} which returns an atomic vector). This is done by flattening and row-binding R objects in the list while creating identifier columns for each level of the list-tree and (optionally) saving the row-names of the objects in a separate column. \code{\link{unlist2d}} can thus also be understood as a recursive generalization of \code{do.call(rbind, l)}, for lists of vectors, data.frames, arrays or heterogenous objects.


}
}
\section{Table of Functions}{
  \tabular{lll}{\emph{ Function / S3 Generic } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Methods } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
  \code{\link{dapply}} \tab No methods, works with matrices and data frames \tab apply functions to rows or columns \cr
  \code{\link{BY}} \tab \code{default, matrix, data.frame, grouped_df} \tab Split-Apply-Combine computing \cr
  \code{\link{TRA}} \tab \code{default, matrix, data.frame, grouped_df} \tab replace and sweep out statistics \cr
  \code{\link[=fscale]{fscale/STD}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df} \tab scale / standardize data \cr
  \code{\link{W}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab demean / center data \cr
  \code{\link{B}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab compute means / average data \cr
  \code{\link{HDW}} \tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab high-dimensional centering and lm residuals \cr
  \code{\link{HDB}} \tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab high-dimensional averages and lm fitted values \cr
  \code{\link[=flag]{flag/L/F}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of) lags / leads \cr
  \code{\link[=fdiff]{fdiff/D}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of lagged/leaded and iterated) differences \cr
  \code{\link[=fdiff]{fgrowth/G}} \tab \code{default, matrix, data.frame, pseries, pdata.frame, grouped_df}  \tab (sequences of lagged/leaded and iterated) growth rates or log-differences
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}
}



