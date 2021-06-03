\name{NEWS}
\title{\emph{collapse} News}

\section{Changes in version 1.2.1 (end of May 2020)}{

\itemize{

\item Minor fixes for 1.2.0 issues that prevented correct installation on Mac OS X and a vignette rebuilding error on solaris.

\item \code{fmode.grouped_df} with groups and weights now saves the sum of the weights instead of the max (this makes more sense as the max only applies if all elements are unique).

}

}

\section{Changes in version 1.2.0 (mid May 2020)}{

\subsection{Changes to Functionality}{
\itemize{
\item \emph{grouped_df} methods for fast statistical functions now always attach the grouping variables to the output in aggregations, unless argument \code{keep.group_vars = FALSE}. (formerly grouping variables were only attached if also present in the data. Code hinged on this feature should be adjusted.)

\item \code{qF} \code{ordered} argument default was changed to \code{ordered = FALSE}, and the \code{NA} level is now only added if \code{na.exclude = FALSE}. Thus \code{qF} now behaves exactly like \code{as.factor}.

\item \code{Recode} is depreciated in favor of \code{recode_num} and \code{recode_char}, it will be removed soon. Similarly \code{replace_non_finite} was renamed to \code{replace_Inf}.

\item In \code{mrtl} and \code{mctl} the argument \code{ret} was renamed \code{return} and now takes descriptive character arguments (the previous version was a direct C++ export and unsafe, code written with these functions should be adjusted).

\item \code{GRP} argument \code{order} is depreciated in favor of argument \code{decreasing}. \code{order} can still be used but will be removed at some point.
}
}

\subsection{Bug Fixes}{
\itemize{
\item Fixed a bug in \code{flag} where unused factor levels caused a group size error. %It is still recommended to remove unused factor levels when programming with factors, some functions check for them, others not. For example \code{fmean(data, f)} will simply generate a missing row for each unused factor level. If in doubt, use safer \code{GRP} objects for grouped programming. A general level check for all functions will not be implemented as this requires an additional pass in some cases.
}
}

\subsection{Improvements}{
\itemize{
\item Faster grouping with \code{GRP} and faster factor generation with added radix method + automatic dispatch between hash and radix method. \code{qF} is now ~ 5x faster than \code{as.factor} on character and around 30x faster on numeric data. Also \code{qG} was enhanced.

\item Further slight speed tweaks here and there.

\item \code{collap} now provides more control for weighted aggregations with additional arguments \code{w}, \code{keep.w} and \code{wFUN} to aggregate the weights as well. The defaults are \code{keep.w = TRUE} and \code{wFUN = fsum}. A specialty of \code{collap} remains that \code{keep.by = TRUE} and \code{keep.w = TRUE} also work for external objects passed, so code of the form \code{collap(data, by, FUN, catFUN, w = data$weights)} will now have an aggregated \code{weights} vector in the first column. %In such cases use \code{keep.w = FALSE} to omit the weights or \code{collap(data, by, FUN, catFUN, w = ~ weights)} to keep the column order.

\item \code{qsu} now also allows weights to be passed in formula i.e. \code{qsu(data, by = ~ group, pid = ~ panelid, w = ~ weights)}.

\item \code{fdiff} now supports quasi-differences i.e. \eqn{x_t - \rho x_{t-1}} and quasi-log differences i.e. \eqn{log(x_t) - \rho log(x_{t-1})}. an arbitrary \eqn{\rho} can be supplied.

\item Added a \code{Dlog} operator for faster access to log-differences.

\item \code{fgrowth} now has a \code{scale} argument, the default is \code{scale = 100} which provides growth rates in percentage terms (as before), but this may now be changed.

\item All statistical and transformation functions now have a hidden list method, so they can be applied to unclassed list-objects as well. An error is however provided in grouped operations with unequal-length columns.
}
}

\subsection{Additions}{
\itemize{

\item Added a suite of functions for fast data manipulation:
\itemize{
  \item \code{fselect} selects variables from a data frame and is equivalent but much faster than \code{dplyr::select}.
  \item \code{fsubset} is a much faster version of \code{base::subset} to subset vectors, matrices and data.frames. The function \code{ss} was also added as a faster alternative to \code{[.data.frame}.
  \item \code{ftransform} is a much faster update of \code{base::transform}, to transform data frames by adding, modifying or deleting columns. The function \code{settransform} does all of that by reference.
  \item \code{fcompute} is equivalent to \code{ftransform} but returns a new data frame containing only the columns computed from an existing one.
  \item \code{na_omit} is a much faster and enhanced version of \code{base::na.omit}.
  \item \code{replace_NA} efficiently replaces missing values in multi-type data.
}

\item Added function \code{fgroup_by} as a much faster version of \code{dplyr::group_by} based on \emph{collapse} grouping. It attaches a 'GRP' object to a data frame, but only works with \emph{collapse}'s fast functions. This allows \emph{dplyr} like manipulations that are fully \emph{collapse} based and thus significantly faster, i.e. \code{data \%>\% fgroup_by(g1,g2) \%>\% fselect(cola,colb) \%>\% fmean}. Note that \code{data \%>\% dplyr::group_by(g1,g2) \%>\% dplyr::select(cola,colb) \%>\% fmean} still works, in which case the \emph{dplyr} 'group' object is converted to 'GRP' as before. However \code{data \%>\% fgroup_by(g1,g2) \%>\% dplyr::summarize(...)} does not work.

\item Added function \code{varying} to efficiently check the variation of multi-type data over a dimension or within groups.

\item Added function \code{radixorder}, same as \code{base::order(..., method = "radix")} but more accessible and with built-in grouping features.

\item Added functions \code{seqid} and \code{groupid} for generalized run-length type id variable generation from grouping and time variables. \code{seqid} in particular strongly facilitates lagging / differencing irregularly spaced panels using \code{flag}, \code{fdiff} etc.
}
}
}

\section{Changes in version 1.1.0 (early April 2020)}{

  \itemize{

\item Fixed remaining gcc10, LTO and valgrind issues in C/C++ code, and added some more tests (there are now ~ 5300 tests ensuring that \emph{collapse} statistical functions perform as expected).

\item Fixed the issue that supplying an unnamed list to \code{GRP()}, i.e. \code{GRP(list(v1, v2))} would give an error. Unnamed lists are now automatically named 'Group.1', 'Group.2', etc...

\item Fixed an issue where aggregating by a single id in \code{collap()} (i.e. \code{collap(data, ~ id1)}), the id would be coded as factor in the aggregated data.frame. All variables including id's now retain their class and attributes in the aggregated data.

\item  Added weights (\code{w}) argument to \code{fsum} and \code{fprod}. \emph{Note}: \code{fmedian} will also support weights as soon as I am able to implement a sufficiently fast (i.e. linear time) algorithm. I also hope to introduce (weighted) quantiles. I am happy for any help with these features.

\item Added an argument \code{mean = 0} to \code{fwithin / W}. This allows simple and grouped centering on an arbitrary mean, \code{0} being the default. For grouped centering \code{mean = "overall.mean"} can be specified, which will center data on the overall mean of the data. The logical argument \code{add.global.mean = TRUE} used to toggle this in \emph{collapse} 1.0.0 is therefore depreciated.

\item Added arguments \code{mean = 0} (the default) and \code{sd = 1} (the default) to \code{fscale / STD}. These arguments now allow to (group) scale and center data to an arbitrary mean and standard deviation. Setting \code{mean = FALSE} just scales data while preserving the mean(s). Special options for grouped scaling are \code{mean = "overall.mean"}, and \code{sd = "within.sd"}. The latter will scale the data such that the standard deviation of each group is equal to the within- standard deviation (= the standard deviation computed on the group-centered data). The fast but numerically unstable algorithm for variance calculation toggled with \code{stable.algo = FALSE} was removed from \code{fscale}. It is still available for \code{fvar} and \code{fsd}.

\item Added the modulus (\code{\%\%}) and subtract modulus (\code{-\%\%}) operations to \code{TRA()}.

\item Added the function \code{finteraction}, for fast interactions, and \code{as_character_factor} to coerce a factor, or all factors in a list, to character (analogous to \code{as_numeric_factor}). Also exported the function \code{ckmatch}, for matching with an error message showing non-matched elements.

}

}

