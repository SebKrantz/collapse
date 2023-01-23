\name{collapse-options}
\alias{collapse-options}
\alias{set_collapse}
\alias{get_collapse}
\alias{.op}
\title{\emph{collapse} Package Options}
\description{Global options affecting package operation. There are 2 kinds of options, those set using \code{\link{options}} (retrievable using \code{\link{getOption}}), and those set (to avoid the performance overhead of \code{getOption()}) using \code{set_collapse()} (retrievable using \code{get_collapse()}). The latter are implemented using an \link{environment} called \code{.op} contained in the package namespace.}

\usage{
set_collapse(\dots)
get_collapse(opts = NULL)
}
\arguments{
\item{\dots}{either comma separated options, or a list of options. Currently only options \code{nthreads = integer} or \code{na.rm = TRUE|FALSE} are supported. These are added to \code{.op}.}
\item{opts}{character. A vector of options to receive from \code{.op}, or \code{NULL} for all options.}

}

\value{
\code{set_collapse()} returns the old content of \code{.op} invisibly as a list. \code{get_collapse()}, if called with only one option, returns the value of the option, and otherwise a list.
}

\section{Options Set Using \code{options()}}{
\itemize{
\item \code{option("collapse_unused_arg_action")} regulates how generic functions (such as the \link[=fast-statistical-functions]{Fast Statistical Functions}) in the package react when an unknown argument is passed to a method. The default action is \code{"warning"} which issues a warning. Other options are \code{"error"}, \code{"message"} or \code{"none"}, whereby the latter enables silent swallowing of such arguments.

\item \code{option("collapse_mask")} can be used to create additional functions in the \emph{collapse} namespace when loading the package, which will mask some existing base R and \emph{dplyr} functions. In particular, \emph{collapse} provides a large number of functions that start with 'f' e.g. \code{fsubset}, \code{ftransform}, \code{fdroplevels} etc.. Specifying \code{options(collapse_mask = c("fsubset", "ftransform", "fdroplevels"))} before loading the package will make additional functions \code{subset}, \code{transform}, and \code{droplevels} available to the user, and mask the corresponding base R functions when the package is attached. In general, all functions starting with 'f' can be passed to the option. There are also a couple of keywords that you can specify to add groups of functions, and 2 special functions:
\itemize{
\item \code{"manip"} adds data manipulation functions: \code{fsubset, ftransform, ftransform<-, ftransformv, fcompute, fcomputev, fselect, fselect<-, fgroup_by, fgroup_vars, fungroup, fsummarise, fsummarize, fmutate, frename, findex_by, findex}
\item \code{"helper"} adds the functions: \code{fdroplevels}, \code{finteraction}, \code{funique}, \code{fnunique}, \code{fcount}, \code{fcountv}, \code{fquantile}, \code{frange}, \code{fdist}, \code{fnlevels}, \code{fnrow} and \code{fncol}. % fdim not because of infinite recursion
\item \code{"fast-fun"} adds the functions contained in the macro: \code{.FAST_FUN}.
\item \code{"fast-stat-fun"} adds the functions contained in the macro: \code{.FAST_STAT_FUN}.
\item \code{"fast-trfm-fun"} adds the functions contained in: \code{setdiff(.FAST_FUN, .FAST_STAT_FUN)}.
\item \code{"all"} turns on all of the above, and additionally exports a function \code{n()} for use in \code{summarise} and \code{mutate}, and masks \code{base::table()} by the much faster \code{qtab()} function.
\item Since v1.8.8 it is also possible to pass \code{"n"} and \code{"qtab"} directly to the option.
}
Note that none of these options will impact internal \emph{collapse} code, but they may change the way your programs run. \code{"manip"} is probably the safest option to start with.
Specifying \code{"fast-fun"}, \code{"fast-stat-fun"}, \code{"fast-trfm-fun"} or \code{"all"} are ambitious as they replace basic R functions like \code{sum} and \code{max}, introducing \emph{collapse}'s \code{na.rm = TRUE} default (which can now be changed using \code{set_collapse}) and different behavior for matrices and data frames. These options also change some internal macros so that base R functions like \code{sum} or \code{max} called inside \code{fsummarise}, \code{fsummarize}, \code{fmutate} or \code{collap} will also receive vectorized execution. In other words, if you put \code{options(collapse_mask = "all")} before loading the package, and you have a collapse-compatible line of \emph{dplyr} code like \code{wlddev |> group_by(region, income) |> summarise(across(PCGDP:POP, sum))}, this will now receive fully optimized execution. %Note however that because of \code{collapse}'s \code{na.rm = TRUE} default, the result will be different unless you add \code{na.rm = FALSE}.
In General, this option is for your convenience, if you want to write visually more appealing code or you want to translate existing \emph{dplyr} codes to \emph{collapse}. Use with care! % \bold{Note} that the option takes effect upon loading the package (code is in the \code{.onLoad} function), so it needs to be set before any function from the package is accessed in any way by any code you run.



%Thus for production code I generally recommend not using it, unless you can ensure that the option is always set before any code is run, and that \emph{collapse} is always attached after \emph{dplyr}.

% Note also that I have not yet investigated in much detail the impact on internal codes of replacing basic R functions like \code{sum}, \code{min} and \code{max} in the package. In terms of unit testing it seems to be quite ok (a number of tests break where I compare fast statistical functions to base functions because of the \code{na.rm = TRUE} default, but otherwise no major issues arise). In general, while this option is cool and can be of great help in translating existing \emph{dplyr} codes to \emph{collapse}, it is experimental and not the best idea for new production code.

\item \code{option("collapse_export_F")}, if set to \code{TRUE}, exports the lead operator \code{F} in the package namespace. The operator was exported by default until v1.9.0, but is now hidden inside the package due to too many problems with \code{base::F}. Alternatively, the operator can be accessed using \code{collapse:::F}.


% \item \code{option("collapse_DT_alloccol")} sets how many empty columns \emph{collapse} data manipulation functions like \code{ftransform} allocate when taking a shallow copy of \emph{data.table}'s. The default is \code{100L}. Note that the \emph{data.table} default is \code{getOption("datatable.alloccol") = 1024L}. I chose a lower default because shallow copies are taken by each data manipulation function if you manipulate \emph{data.table}'s with collapse, and the cost increases with the number of overallocated columns. With 100 columns, the cost is 2-5 microseconds per copy.

\item \code{option("collapse_nthreads")} and \code{option("collapse_na.rm")} can be set before loading the package to initialize \emph{collapse} with different defaults than \code{nthreads = 1L} and \code{na.rm = TRUE} (e.g. using an \code{\link{.Rprofile}} file). Once loaded, these options have no effect, and users need to use \code{set_collapse()} to change the argument defaults.
}
}
\note{
Options \code{"collapse_mask"} and \code{"collapse_export_F"} (and \code{"collapse_nthreads"}, \code{"collapse_na.rm"}) need to be set before the package is loaded, which means before any component of the package is accessed in any way, including the loading of packages that depend on \emph{collapse}.

A safe way to set them is by using a \code{\link{.Rprofile}} file in your user or project directory (see also \href{https://www.statmethods.net/interface/customizing.html}{here}, the user-level file is located at \code{file.path(Sys.getenv("HOME"), ".Rprofile")} and can be edited using \code{file.edit(Sys.getenv("HOME"), ".Rprofile")}), or by using a \href{https://fastverse.github.io/fastverse/articles/fastverse_intro.html#custom-fastverse-configurations-for-projects}{\code{.fastverse}} configuration file in the project directory.
}

\seealso{
  \link[=collapse-documentation]{Collapse Overview}, \link{collapse-package}
}

\examples{
# Setting new values
oldopts = set_collapse(nthreads = 2, na.rm = FALSE)

# Getting the values
get_collapse()
get_collapse("nthreads")

# Resetting
set_collapse(oldopts)

rm(oldopts)
}

\keyword{documentation}
