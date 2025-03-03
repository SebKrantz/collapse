\name{fast-grouping-ordering}
\alias{A2-fast-grouping-ordering}
\alias{fast-grouping-ordering}
\title{Fast Grouping and Ordering} % \emph{collapse}

\description{
\emph{collapse} provides the following functions to efficiently group and order data:
\itemize{
\item \code{\link[=radixorder]{radixorder(v)}}, provides fast radix-ordering through direct access to the method \code{\link[=order]{order(..., method = "radix")}}, as well as the possibility to return some attributes very useful for grouping data and finding unique elements. The function \code{\link[=roworder]{roworder(v)}} efficiently reorders a data frame. %The source code for both \code{\link{radixorder}} and \code{\link{order(\dots, method = "radix")}, comes from \code{data.table:::forder}.
%\code{\link{radixorder}} was modified to optionally return either a vector of group starts, a vector of group sizes, or both as an attribute, and also an attribute providing the size of the largest group and a logical statement on whether the input was already ordered. The function \code{\link{radixorderv}} exists as a programmers alternative.

\item \code{\link[=group]{group(v)}} provides fast grouping in first-appearance order of rows, based on a hashing algorithm in C. Objects have class 'qG', see below.

\item \code{\link{GRP}} creates \emph{collapse} grouping objects of class 'GRP' based on \code{\link{radixorder}} or \code{\link{group}}. 'GRP' objects form the central building block for grouped operations and programming in \emph{collapse} and are very efficient inputs to all \emph{collapse} functions supporting grouped operations. %A 'GRP' object provides information about (1) the number of groups, (2) which rows belong to which group, (3) the group sizes, (4) the unique groups, (5) the variables used for grouping, (6) whether the grouping and initial inputs were ordered and (7) (optionally) the output from \code{\link{radixorder}} containing the ordering vector with group starts and maximum group size attributes.

\item \code{\link{fgroup_by}} provides a fast replacement for \code{dplyr::group_by}, creating a grouped data frame (or data.table / tibble etc.) with a 'GRP' object attached. This grouped frame can be used for grouped operations using \emph{collapse}'s fast functions. % \emph{dplyr} functions will treat this tibble like an ordinary (non-grouped) one.

\item \code{\link{fmatch}} is a fast alternative to \code{\link[base]{match}}, which also supports matching of data frame rows.

\item \code{\link{funique}} is a faster version of \code{\link{unique}}. The data frame method also allows selecting unique rows according to a subset of the columns. \code{\link{fnunique}} efficiently calculates the number of unique values/rows. \code{\link{fduplicated}} is a fast alternative to \code{\link{duplicated}}. \code{\link{any_duplicated}} is a simpler and faster alternative to \code{\link{anyDuplicated}}.

\item \code{\link[=fcount]{fcount(v)}} computes group counts based on a subset of columns in the data, and is a fast replacement for \code{dplyr::count}. % \code{\link{fcountv}} is a programmers version of the function.

\item \code{\link{qF}}, shorthand for 'quick-factor' implements very fast factor generation from atomic vectors using either radix ordering \code{method = "radix"} or hashing \code{method = "hash"}. Factors can also be used for efficient grouped programming with \emph{collapse} functions, especially if they are generated using \code{qF(x, na.exclude = FALSE)} which assigns a level to missing values and attaches a class 'na.included' ensuring that no additional missing value checks are executed by \emph{collapse} functions.

\item \code{\link{qG}}, shorthand for 'quick-group', generates a kind of factor-light without the levels attribute but instead an attribute providing the number of levels. Optionally the levels / groups can be attached, but without converting them to character. Objects have a class 'qG', which is also recognized in the \emph{collapse} ecosystem.

\item \code{\link{fdroplevels}} is a substantially faster replacement for \code{\link{droplevels}}.

\item \code{\link{finteraction}} is a fast alternative to \code{\link{interaction}} implemented as a wrapper around \code{as_factor_GRP(GRP(\dots))}. It can be used to generate a factor from multiple vectors, factors or a list of vectors / factors. Unused factor levels are always dropped.

\item \code{\link{groupid}} is a generalization of \code{data.table::rleid} providing a run-length type group-id from atomic vectors. It is generalization as it also supports passing an ordering vector and skipping missing values. For example \code{\link{qF}} and \code{\link{qG}} with \code{method = "radix"} are essentially implemented using \code{groupid(x, radixorder(x))}.

\item \code{\link{seqid}} is a specialized function which creates a group-id from sequences of integer values. For any regular panel dataset \code{groupid(id, order(id, time))} and \code{seqid(time, order(id, time))} provide the same id variable. \code{\link{seqid}} is especially useful for identifying discontinuities in time-sequences.

\item \code{\link{timeid}} is a specialized function to convert integer or double vectors representing time (such as 'Date', 'POSIXct' etc.) to factor or 'qG' object based on the greatest common divisor of elements (thus preserving gaps in time intervals).

}
}
\section{Table of Functions}{
  \tabular{lllll}{\emph{ Function / S3 Generic }   \tab\tab \emph{ Methods }   \tab\tab \emph{ Description }  \cr
                 \code{\link[=radixorder]{radixorder(v)}}  \tab\tab No methods, for data frames and vectors  \tab\tab Radix-based ordering + grouping information \cr
                 \code{\link[=roworder]{roworder(v)}} \tab\tab No methods, for data frames incl. pdata.frame \tab\tab Row sorting/reordering \cr
                 \code{\link[=group]{group(v)}}  \tab\tab No methods, for data frames and vectors  \tab\tab Hash-based grouping + grouping information \cr
                 \code{\link{GRP}} \tab\tab \code{default, GRP, factor, qG, grouped_df, pseries, pdata.frame}  \tab\tab Fast grouping and a flexible grouping object \cr
                 \code{\link{fgroup_by}} \tab\tab No methods, for data frames  \tab\tab Fast grouped data frame \cr
                 \code{\link{fmatch}} \tab\tab No methods, for vectors and data frames  \tab\tab Fast matching \cr
                 \code{\link{funique}}, \code{\link{fnunique}}, \code{\link{fduplicated}}, \code{\link{any_duplicated}} \tab\tab \code{default, data.frame, sf, pseries, pdata.frame, list}  \tab\tab Fast (number of) unique values/rows \cr
                \code{\link[=fcount]{fcount(v)}} \tab\tab Internal generic, supports vectors, matrices, data.frames, lists, grouped_df and pdata.frame \tab\tab Fast group counts \cr
                 \code{\link{qF}} \tab\tab No methods, for vectors \tab\tab Quick factor generation \cr
                 \code{\link{qG}} \tab\tab No methods, for vectors \tab\tab Quick grouping of vectors and a 'factor-light' class \cr
                 \code{\link{fdroplevels}} \tab\tab \code{factor, data.frame, list} \tab\tab Fast removal of unused factor levels \cr
                 \code{\link{finteraction}} \tab\tab No methods, for data frames and vectors \tab\tab Fast interactions \cr
                 \code{\link{groupid}} \tab\tab No methods, for vectors \tab\tab Run-length type group-id \cr
                 \code{\link{seqid}} \tab\tab No methods, for integer vectors \tab\tab Run-length type integer sequence-id \cr
                 \code{\link{timeid}} \tab\tab No methods, for integer or double vectors \tab\tab Integer-id from time/date sequences \cr
}
}
\seealso{
\link[=collapse-documentation]{Collapse Overview}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=time-series-panel-series]{Time Series and Panel Series}
}
\keyword{manip}
\keyword{documentation}


