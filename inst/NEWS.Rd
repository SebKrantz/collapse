\name{NEWS}
\title{\emph{collapse} News}
\section{Changes in version 1.1.0}{
  \itemize{
\item Fixed issues with 2 tests failing on old windows release, and added some more tests.

\item Fixed the issue that supplying an unnamed list to \code{GRP()}, i.e. \code{GRP(list(v1, v2))} would give an error. Unnamed lists are now automatically named 'Group.1', 'Group.2', etc...

\item Fixed an issue where aggregating by a single id in \code{collap()} (i.e. \code{collap(data, ~ id1)}), the id would be coded as factor in the aggregated data.frame. All variables including id's now retain their class and attributes in the aggregated data.

\item Added an argument \code{mean = 0} to \code{fwithin / W}. This allows simple and grouped centering on an arbitrary mean, \code{0} being the default. For grouped centering \code{mean = "overall.mean"} can be specified, which will center data on the overall mean of the data. The logical argument \code{add.global.mean = TRUE} used to toggle this in *collapse* 1.0.0 is therefore depreciated.

\item Added arguments \code{mean = 0} (the default) and \code{sd = 1} (the default) to \code{fscale / STD}. These arguments now allow to (group) scale and center data to an arbitrary mean and standard deviation. Setting \code{mean = FALSE} just scales data while preserving the mean(s). Special options for grouped scaling are \code{mean = "overall.mean"}, and \code{sd = "within.sd"}. The latter will scale the data such that the standard deviation of each group is equal to the within- standard deviation (= the standard deviation computed on the group-centered data). The fast but numerically unstable algorithm for variance calculation toggled with \code{stable.algo = FALSE} was removed from \code{fscale}. It is still available for \code{fvar} and \code{fsd}.

\item Added the modulus (\code{\%\%}) and subtract modulus (\code{-\%\%}) operations to \code{TRA()}.

\item Added the function \code{finteraction}, for fast interactions, and \code{as.character_factor} to coerce a factor, or all factors in a list, to character (analogous to \code{as.numeric_factor}). Also exported the function \code{ckmatch}, for matching with an error message showing non-matched elements.
  }
}

