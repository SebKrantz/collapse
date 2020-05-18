\name{A3-data-frame-manipulation} % \name{Time-Series and Panel Computations}
\alias{A3-data-frame-manipulation}
% \alias{tscomp}
\title{Fast Data Frame Manipulation} %

\description{
\emph{collapse} provides the following functions for fast manipulation of (mostly) data.frames.
\itemize{
\item \code{\link{fselect}} is a much faster alternative to \code{dplyr::select} to select columns using expressions involving column names. \code{\link{get_vars}} is a more versatile and programmer friendly function to efficiently select and replace columns by names, indices, logical vectors, regular expressions or using functions to identify columns.

\item The functions \code{\link{num_vars}}, \code{\link{cat_vars}}, \code{\link{char_vars}}, \code{\link{fact_vars}}, \code{\link{logi_vars}} and \code{\link{Date_vars}} are convenience functions to efficiently select and replace columns by data type.

\item \code{\link{add_vars}} efficiently adds new columns at any position within a data.frame (default at the end). This can be done vie replacement (i.e. \code{add_vars(data) <- newdata}) or returning the appended data (i.e. \code{add_vars(data, newdata1, newdata2, ...)}). Because of the latter, \code{add_vars} is also a more efficient alternative to \code{cbind.data.frame}.

\item \code{\link{fsubset}} is a much faster version of \code{base::subset} for efficiently subset vectors, matrices and data.frames. If the non-standard evaluation offered by \code{\link{fsubset}} is not needed, the function \code{\link{ss}} is a much faster and also more secure alternative to \code{[.data.frame}.

\item \code{\link{ftransform}} is a much faster version of \code{base::transform}, to modify and delete existing columns or append a data frame with new computed columns. \code{\link{settransform}} does all of that by reference, i.e. it modifies the data frame in the global environment. \code{\link{fcompute}} is similar to \code{\link{ftransform}} but only returns modified and computed columns in a new data frame.

}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr
                 \code{\link{fselect}} \tab\tab No methods, for data.frame's  \tab\tab Fast select columns (non-standard evaluation) \cr
                 \code{\link{get_vars}}, \code{\link{num_vars}}, \code{\link{cat_vars}}, \code{\link{char_vars}}, \code{\link{fact_vars}}, \code{\link{logi_vars}}, \code{\link{Date_vars}} \tab\tab No methods, for data.frame's  \tab\tab Fast select or replace columns \cr
                 \code{\link{add_vars}} \tab\tab No methods, for data.frame's  \tab\tab Fast add columns \cr
                 \code{\link{fsubset}} \tab\tab \code{default, matrix, data.frame} \tab\tab Fast subset data (non-standard evaluation) \cr
                 \code{\link{ss}} \tab\tab No methods, for data.frame's \tab\tab Fast subset data.frames's \cr
                 \code{\link{ftransform}} \tab\tab No methods, for data.frame's \tab\tab Compute, modify or delete columns (non-standard evaluation) \cr
                 \code{\link{settransform}} \tab\tab No methods, for data.frame's \tab\tab Compute, modify or delete columns by reference (non-standard evaluation) \cr
                 \code{\link{fcompute}} \tab\tab No methods, for data.frame's \tab\tab Compute or modify columns, returned in a new data.frame (non-standard evaluation) \cr
}
}
\seealso{
\link[=A4-quick-conversion]{Quick Data Conversion}, \link[=collapse-documentation]{Collapse Overview}
}
\keyword{manip}
\keyword{documentation}


