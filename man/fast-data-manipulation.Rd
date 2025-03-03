\name{fast-data-manipulation}
\alias{A3-fast-data-manipulation}
\alias{fast-data-manipulation}
\title{Fast Data Manipulation}

\description{
\emph{collapse} provides the following functions for fast manipulation of (mostly) data frames.
\itemize{
\item \code{\link{fselect}} is a much faster alternative to \code{dplyr::select} to select columns using expressions involving column names. \code{\link{get_vars}} is a more versatile and programmer friendly function to efficiently select and replace columns by names, indices, logical vectors, regular expressions, or using functions to identify columns.

\item \code{\link{num_vars}}, \code{\link{cat_vars}}, \code{\link{char_vars}}, \code{\link{fact_vars}}, \code{\link{logi_vars}} and \code{\link{date_vars}} are convenience functions to efficiently select and replace columns by data type.

\item \code{\link{add_vars}} efficiently adds new columns at any position within a data frame (default at the end). This can be done vie replacement (i.e. \code{add_vars(data) <- newdata}) or returning the appended data, e.g., \code{add_vars(data, newdata1, newdata2, \dots)}. It is thus also an efficient alternative to \code{\link{cbind.data.frame}}.

\item \code{\link{rowbind}} efficiently combines data frames / lists row-wise. The implementation is derived from \code{data.table::rbindlist}, it is also a fast alternative to \code{\link{rbind.data.frame}}.

\item \code{\link{join}} provides fast, class-agnostic, and verbose table joins.

\item \code{\link{pivot}} efficiently reshapes data, supporting longer, wider and recast pivoting, as well as multi-column-pivots and pivots taking along variable labels.

\item \code{\link{fsubset}} is a much faster version of \code{\link{subset}} to efficiently subset vectors, matrices and data frames. If the non-standard evaluation offered by \code{\link{fsubset}} is not needed, the function \code{\link{ss}} is a much faster and more secure alternative to \code{[.data.frame}.

\item \code{\link[=fslice]{fslice(v)}} is a much faster alternative to \code{dplyr::slice_[head|tail|min|max]} for filtering/deduplicating matrix-like objects (by groups).

\item \code{\link{fsummarise}} is a much faster version of \code{dplyr::summarise}, especially when used together with the \link[=fast-statistical-functions]{Fast Statistical Functions} and \code{\link{fgroup_by}}.

\item \code{\link{fmutate}} is a much faster version of \code{dplyr::mutate}, especially when used together with the \link[=fast-statistical-functions]{Fast Statistical Functions}, the fast \link[=data-transformations]{Data Transformation Functions}, and \code{\link{fgroup_by}}.

\item \code{\link[=ftransform]{ftransform(v)}} is a much faster version of \code{\link{transform}}, which also supports list input and nested pipelines. \code{\link[=ftransform]{settransform(v)}} does all of that by reference, i.e. it assigns to the calling environment. \code{\link[=fcompute]{fcompute(v)}} is similar to \code{\link[=ftransform]{ftransform(v)}} but only returns modified/computed columns. %As a new feature, it is now possible to bulk-process columns with \code{\link{ftransform}}, i.e. \code{ftransform(data, fscale(data[1:2]))} is the same as \code{ftransform(data, col1 = fscale(col1), col2 = fscale(col2))}, and \code{ftransform(data) <- fscale(data[1:2]))} or \code{settransform(data, fscale(data[1:2]))} are both equivalent to \code{data[1:2] <- fscale(data[1:2]))}. Non-matching columns are added to the data.frame.

\item \code{\link{roworder}} is a fast substitute for \code{dplyr::arrange}, but the syntax is inspired by \code{data.table::setorder}.

\item \code{\link{colorder}} efficiently reorders columns in a data frame, see also \code{data.table::setcolorder}.

\item \code{\link{frename}} is a fast substitute for \code{dplyr::rename}, to efficiently rename various objects. \code{\link{setrename}} renames objects by reference. \code{\link{relabel}} and \code{\link{setrelabel}} do the same thing for variable labels (see also \code{\link{vlabels}}).

}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr
                 \code{\link[=fselect]{fselect(<-)}} \tab\tab No methods, for data frames  \tab\tab Fast select or replace columns (non-standard evaluation) \cr
                 \code{\link[=get_vars]{get_vars(<-)}}, \code{\link[=num_vars]{num_vars(<-)}}, \code{\link[=cat_vars]{cat_vars(<-)}}, \code{\link[=char_vars]{char_vars(<-)}}, \code{\link[=fact_vars]{fact_vars(<-)}}, \code{\link[=logi_vars]{logi_vars(<-)}}, \code{\link[=date_vars]{date_vars(<-)}} \tab\tab No methods, for data frames  \tab\tab Fast select or replace columns \cr
                 \code{\link[=add_vars]{add_vars(<-)}} \tab\tab No methods, for data frames  \tab\tab Fast add columns \cr
                 \code{\link{rowbind}} \tab\tab No methods, for lists of lists/data frames \tab\tab Fast row-binding lists \cr
                 \code{\link{join}} \tab\tab No methods, for data frames \tab\tab Fast table joins \cr
                 \code{\link{pivot}} \tab\tab No methods, for data frames \tab\tab Fast reshaping \cr
                 \code{\link{fsubset}} \tab\tab \code{default, matrix, data.frame, pseries, pdata.frame} \tab\tab Fast subset data (non-standard evaluation) \cr
                 \code{\link{ss}} \tab\tab No methods, for data frames \tab\tab Fast subset data frames \cr
                 \code{\link[=fslice]{fslice(v)}} \tab\tab No methods, for matrices and data frames\tab\tab Fast slicing of rows \cr
                 \code{\link{fsummarise}} \tab\tab No methods, for data frames \tab\tab Fast data aggregation \cr
                 \code{\link{fmutate}}, \code{\link[=ftransform]{(f/set)transform(v)(<-)}} \tab\tab No methods, for data frames \tab\tab Compute, modify or delete columns (non-standard evaluation) \cr
                 %\code{\link{settransform}} \tab\tab No methods, for data frames \tab\tab Compute, modify or delete columns by reference (non-standard evaluation) \cr
                 \code{\link[=fcompute]{fcompute(v)}} \tab\tab No methods, for data frames \tab\tab Compute or modify columns, returned in a new data frame (non-standard evaluation) \cr
                 \code{\link[=roworder]{roworder(v)}} \tab\tab No methods, for data frames incl. pdata.frame \tab\tab Reorder rows and return data frame (standard and non-standard evaluation) \cr
                 \code{\link[=colorder]{colorder(v)}} \tab\tab No methods, for data frames \tab\tab Reorder columns and return data frame (standard and non-standard evaluation) \cr
                 \code{\link[=frename]{(f/set)rename}}, \code{\link[=frename]{(set)relabel}} \tab\tab No methods, for all objects with 'names' attribute \tab\tab Rename and return object / relabel columns in a data frame. \cr
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=quick-conversion]{Quick Data Conversion}, \link[=recode-replace]{Recode and Replace Values}
}
\keyword{manip}
\keyword{documentation}


