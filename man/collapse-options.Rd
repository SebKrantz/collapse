\name{collapse-options}
\alias{collapse-options}
\alias{AA4-collapse-options}
\alias{set_collapse}
\alias{get_collapse}
\alias{.op}
\title{\emph{collapse} Package Options}
\description{\emph{collapse} is globally configurable to an extent few packages are: the default value of key function arguments governing the behavior of its algorithms, and the exported namespace, can be adjusted interactively through the \code{set_collapse()} function.

These options are saved in an internal environment called \code{.op} (for safety and performance reasons) visible in the documentation of some functions such as \code{\link{fmean}}. The contents of this environment can be accessed using \code{get_collapse()}.

There are also a few options that can be set using \code{\link{options}} (retrievable using \code{\link{getOption}}). These options mainly affect package startup behavior.

%Global options affecting package operation. There are 2 kinds of options, those set using \code{\link{options}} (retrievable using \code{\link{getOption}}), and those set (to avoid the performance overhead of \code{getOption()}) using \code{set_collapse()} (retrievable using \code{get_collapse()}). The latter are implemented using an \link{environment} called \code{.op} contained in the package namespace.
}

\usage{
set_collapse(\dots)
get_collapse(opts = NULL)
}
\arguments{
\item{\dots}{either comma separated options, or a single list of options. The available options are:

  \tabular{lll}{
    \code{na.rm} \tab\tab logical, default \code{TRUE}. Sets the default for statistical algorithms such as the \link[=fast-statistical-functions]{Fast Statistical Functions} to skip missing values. If your data does not have missing values, or only in rare cases, it is recommended to change this to \code{FALSE} for performance gains. \emph{Note} that this does not affect other (non-statistical) uses of \code{na.rm} arguments, such as in \code{\link{pivot}}. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{sort} \tab\tab logical, default \code{TRUE}. Sets the default for grouping operations to be sorted. This also applies to factor generation using \code{\link{qF}} and tabulation with \code{\link{qtab}}, but excludes other uses of \code{sort} arguments where grouping is not the objective (such as in \code{\link{funique}} or \code{\link{pivot}}). In general, sorted grouping (internally using \code{\link{radixorder}}) is slower than hash-based direct grouping (internally using \code{\link{group}}). However, if data is pre-sorted, sorted grouping is slightly faster. In general, if records don't need to be sorted or you want to maintain their first-appearance order, changing this to \code{FALSE} is recommended and often brings substantial performance gains. \emph{Note} that this also affects internal grouping applied when atomic vectors (except for factors) or lists are passed to \code{g} arguments in \link[=fast-statistical-functions]{Fast Statistical Functions}. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{nthreads} \tab\tab integer, default 1. Sets the default for OpenMP multithreading, available in certain statistical and data manipulation functions. Setting values greater than 1 is strongly recommended with larger datasets. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{stable.algo} \tab\tab logical, default \code{TRUE}. Option passed to \code{\link[=fvar]{fvar()/fsd()}} and \code{\link[=qsu]{qsu()}}. \code{FALSE} enables one-pass standard deviation calculation, which is very fast, but might incur catastrophic cancellation if numbers are large and the variance is small. see \code{\link{fvar}} for details. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{stub} \tab\tab logical, default \code{TRUE}. Controls whether \link[=.OPERATOR_FUN]{transformation operators} (\code{.OPERATOR_FUN}) such as \code{\link{W}}, \code{\link{L}}, \code{\link{STD}} etc. add prefixes to transformed columns of matrix and data.frame-like objects. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{verbose} \tab\tab integer, default \code{1}. Print additional (diagnostic) information or messages when executing code. Currently only used in \code{\link{join}} and \code{\link{roworder}}. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{digits} \tab\tab integer, default \code{2}. Number of digits to print, e.g. in \code{\link{descr}} or \code{\link{pwcor}}. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

    \code{mask} \tab\tab character, default \code{NULL}. Allows masking existing base R/dplyr functions with faster \emph{collapse} versions, by creating additional functions in the namespace and instantly exporting them: \cr \tab\tab \cr

    \tab\tab For example \code{set_collapse(mask = "unique")} (or, equivalently, \code{set_collapse(mask = "funique")}) will create \code{unique <- funique} in the \emph{collapse} namespace, export \code{unique()}, and silently detach and attach the namespace again so R can find it - all in millisecond. Thus calling \code{unique()} afterwards uses the \emph{collapse} version - which is many times faster. \code{funique} remains available and you can still call \code{base::unique} explicitly. \cr \tab\tab \cr

  \tab\tab All \emph{collapse} functions starting with 'f' can be passed to the option (with or without the 'f') e.g. \code{set_collapse(mask = c("subset", "transform", "droplevels"))} creates \code{subset <- fsubset}, \code{transform <- ftransform} etc. Special functions are \code{"n"} and \code{"table"/"qtab"}, and \code{"\%in\%"}, which create \code{n <- GRPN} (for use in \code{(f)summarise}/\code{(f)mutate}), \code{table <- qtab}, and replace \code{\%in\%} with a fast version using \code{\link{fmatch}}, respectively. \cr \tab\tab \cr

  \tab\tab There are also a couple of convenience keywords that you can use to mask groups of functions: \cr \tab\tab \cr

\tab\tab - \code{"manip"} adds data manipulation functions: \code{fsubset, fslice, fslicev, ftransform, ftransform<-, ftransformv, fcompute, fcomputev, fselect, fselect<-, fgroup_by, fgroup_vars, fungroup, fsummarise, fsummarize, fmutate, frename, findex_by, findex}. \cr \tab\tab \cr
\tab\tab - \code{"helper"} adds the functions: \code{fdroplevels}, \code{finteraction}, \code{fmatch}, \code{funique}, \code{fnunique}, \code{fduplicated}, \code{fcount}, \code{fcountv}, \code{fquantile}, \code{frange}, \code{fdist}, \code{fnlevels}, \code{fnrow} and \code{fncol}. \cr \tab\tab \cr % fdim not because of infinite recursion
\tab\tab - \code{"special"} exports \code{n()}, \code{table()} and \code{\%in\%}. See above. \cr \tab\tab \cr
\tab\tab - \code{"fast-fun"} adds the functions contained in the macro: \code{.FAST_FUN}. See also Note. \cr \tab\tab \cr
\tab\tab - \code{"fast-stat-fun"} adds the functions contained in the macro: \code{.FAST_STAT_FUN}. See also Note. \cr \tab\tab \cr
\tab\tab - \code{"fast-trfm-fun"} adds the functions contained in: \code{setdiff(.FAST_FUN, .FAST_STAT_FUN)}. See also Note. \cr \tab\tab \cr
\tab\tab - \code{"all"} turns on all of the above.\cr \tab\tab \cr


\tab\tab The re-attaching of the namespace places \emph{collapse} at the top of the search path (after the global environment), implying that all its exported functions will take priority over other libraries. Users can use \code{fastverse::fastverse_conflicts()} to check which functions are masked following \code{set_collapse(mask = ...)}. The option can be changed at any time with immediate effect. Using \code{set_collapse(mask = NULL)} removes all masked functions from the namespace, and can also be called simply to place \emph{collapse} at the top of the search path. \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr

 \code{remove} \tab\tab character, default \code{NULL}. Similar to 'mask': allows removing functions from the exported namespace (they are still in the namespace, just no longer exported). All \emph{collapse} functions can be passed here. This argument is always evaluated after 'mask', thus you can also remove masked functions again i.e. after setting a keyword which masks a bunch of functions. There are also a couple of convenience keywords you can specify to bulk-remove certain functions:
\cr \tab\tab \cr

\tab\tab - \code{"shorthand"} removes function shorthands: \code{gv, gv<-, av, av<-, nv, nv<-, gvr, gvr<-, itn, ix, slt, slt<-, sbt, gby, iby, mtt, smr, tfm, tfmv, tfm<-, settfm, settfmv, rnm}. \cr \tab\tab \cr
\tab\tab - \code{"infix"} removes infix functions: \code{\%!=\%, \%[!]in\%, \%[!]iin\%, \%*=\%, \%+=\%, \%-=\%, \%/=\%, \%=\%, \%==\%, \%c*\%, \%c+\%, \%c-\%, \%c/\%, \%cr\%, \%r*\%, \%r+\%, \%r-\%, \%r/\%, \%rr\%}.\cr \tab\tab \cr
\tab\tab - \code{"operator"} removes functions contained in the macro: \code{.OPERATOR_FUN}.\cr \tab\tab \cr
\tab\tab - \code{"old"} removes depreciated functions contained in the macro: \code{.COLLAPSE_OLD}.\cr \tab\tab \cr


\tab\tab Like 'mask', the option is alterable and reversible. Specifying \code{set_collapse(remove = NULL)} restores the exported namespace. Also like 'mask', this option silently detaches and attaches \emph{collapse} again, ensuring that it is at the top of the search path.

 % \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr \tab\tab \cr
  }
}
\item{opts}{character. A vector of options to receive from \code{.op}, or \code{NULL} for a list of all options.}

}

\value{
\code{set_collapse()} returns the old content of \code{.op} invisibly as a list. \code{get_collapse()}, if called with only one option, returns the value of the option, and otherwise a list.
}
\note{
Setting keywords "fast-fun", "fast-stat-fun", "fast-trfm-fun" or "all" with \code{set_collapse(mask = ...)} will also adjust internal optimization flags, e.g. in \code{\link[=fsummarise]{(f)summarise}} and \code{\link[=fmutate]{(f)mutate}}, so that these functions - and all expressions containing them - receive vectorized execution (see examples of \code{\link[=fsummarise]{(f)summarise}} and \code{\link[=fmutate]{(f)mutate}}). Users should be aware of expressions like \code{fmutate(mu = sum(var) / lenth(var))}: this usually gets executed by groups, but with these keywords set,this will be vectorized (like \code{fmutate(mu = fsum(var) / lenth(var))}) implying grouped sum divided by overall length. In this case \code{fmutate(mu = base::sum(var) / lenth(var))} needs to be specified to retain the original result.

\emph{Note} that passing individual functions like \code{set_collapse(mask = "(f)sum")} will \bold{not} change internal optimization flags for these functions. This is to ensure consistency i.e. you can be either all in (by setting appropriate keywords) or all out when it comes to vectorized stats with basic R names.

\emph{Note} also that masking does not change documentation links, so you need to look up the f- version of a function to get the right documentation.

A safe way to set options affecting startup behavior is by using a \code{\link{.Rprofile}} file in your user or project directory (see also \href{https://www.datacamp.com/doc/r/customizing}{here}, the user-level file is located at \code{file.path(Sys.getenv("HOME"), ".Rprofile")} and can be edited using \code{file.edit(Sys.getenv("HOME"), ".Rprofile")}), or by using a \href{https://fastverse.org/fastverse/articles/fastverse_intro.html#custom-fastverse-configurations-for-projects}{\code{.fastverse}} configuration file in the project directory.

\code{options("collapse_remove")} does in fact remove functions from the namespace and cannot be reversed by \code{set_collapse(remove = NULL)} once the package is loaded. It is only reversed by re-loading \emph{collapse}.
}

\section{Options Set Using \code{options()}}{
\itemize{
\item \code{"collapse_unused_arg_action"} regulates how generic functions (such as the \link[=fast-statistical-functions]{Fast Statistical Functions}) in the package react when an unknown argument is passed to a method. The default action is \code{"warning"} which issues a warning. Other options are \code{"error"}, \code{"message"} or \code{"none"}, whereby the latter enables silent swallowing of such arguments.

% \item \code{option("collapse_mask")}

%none of these options will impact internal \emph{collapse} code, but they may change the way your programs run. \code{"manip"} is probably the safest option to start with.
%Specifying \code{"fast-fun"}, \code{"fast-stat-fun"}, \code{"fast-trfm-fun"} or \code{"all"} are ambitious as they replace basic R functions like \code{sum} and \code{max}, introducing \emph{collapse}'s \code{na.rm = TRUE} default (which can now be changed using \code{set_collapse}) and different behavior for matrices and data frames.

% These options also change some internal macros so that base R functions like \code{sum} or \code{max} called inside \code{fsummarise}, \code{fsummarize}, \code{fmutate} or \code{collap} will also receive vectorized execution. In other words, if you put \code{options(collapse_mask = "all")} before loading the package, and you have a collapse-compatible line of \emph{dplyr} code like \code{wlddev |> group_by(region, income) |> summarise(across(PCGDP:POP, sum))}, this will now receive fully optimized execution. %Note however that because of \code{collapse}'s \code{na.rm = TRUE} default, the result will be different unless you add \code{na.rm = FALSE}.
% In General, this option is for your convenience, if you want to write visually more appealing code or you want to translate existing \emph{dplyr} codes to \emph{collapse}. Use with care!

%Thus for production code I generally recommend not using it, unless you can ensure that the option is always set before any code is run, and that \emph{collapse} is always attached after \emph{dplyr}.


\item \code{"collapse_export_F"}, if set to \code{TRUE}, exports the lead operator \code{F} in the package namespace when loading the package. The operator was exported by default until v1.9.0, but is now hidden inside the package due to too many problems with \code{base::F}. Alternatively, the operator can be accessed using \code{collapse:::F}.


% \item \code{option("collapse_DT_alloccol")} sets how many empty columns \emph{collapse} data manipulation functions like \code{ftransform} allocate when taking a shallow copy of \emph{data.table}'s. The default is \code{100L}. Note that the \emph{data.table} default is \code{getOption("datatable.alloccol") = 1024L}. I chose a lower default because shallow copies are taken by each data manipulation function if you manipulate \emph{data.table}'s with collapse, and the cost increases with the number of overallocated columns. With 100 columns, the cost is 2-5 microseconds per copy.

\item \code{"collapse_nthreads"}, \code{"collapse_na_rm"}, \code{"collapse_sort"}, \code{"collapse_stable_algo"}, \code{"collapse_verbose"}, \code{"collapse_digits"}, \code{"collapse_mask"} and \code{"collapse_remove"} can be set before loading the package to initialize \code{.op} with different defaults (e.g. using an \code{\link{.Rprofile}} file). Once loaded, these options have no effect, and users need to use \code{set_collapse()} to change them. See also the Note.
}
}

\seealso{
  \link[=collapse-documentation]{Collapse Overview}, \link{collapse-package}
}

\examples{
# Setting new values
oldopts <- set_collapse(nthreads = 2, na.rm = FALSE)

# Getting the values
get_collapse()
get_collapse("nthreads")

# Resetting
set_collapse(oldopts)
rm(oldopts)

\dontrun{
## This is a typical working setup I use:
library(fastverse)
# Loading other stats packages with fastverse_extend():
# displays versions, checks conflicts, and installs if unavailable
fastverse_extend(qs, fixest, grf, glmnet, install = TRUE)
# Now setting collapse options with some namespace modification
set_collapse(
  nthreads = 4,
  sort = FALSE,
  mask = c("manip", "helper", "special", "mean", "scale"),
  remove = "old"
)
# Final conflicts check (optional)
fastverse_conflicts()

# For some simpler scripts I also use
set_collapse(
  nthreads = 4,
  sort = FALSE,
  mask = "all",
  remove = c("old", "between") # I use data.table::between > fbetween
)

# This is now collapse code
mtcars |>
  subset(mpg > 12) |>
  group_by(cyl) |>
  sum()
}

## Changing what happens with unused arguments
oldopts <- options(collapse_unused_arg_action = "message") # default: "warning"
fmean(mtcars$mpg, bla = 1)

# Now nothing happens, same as base R
options(collapse_unused_arg_action = "none")
fmean(mtcars$mpg, bla = 1)
mean(mtcars$mpg, bla = 1)

options(oldopts)
rm(oldopts)
}

\keyword{documentation}
