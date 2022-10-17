\name{across}
\alias{across}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply Functions Across Multiple Columns
}
\description{
\code{across()} can be used inside \code{\link{fmutate}} and \code{\link{fsummarise}} to apply one or more functions to a selection of columns. It is overall very similar to \code{dplyr::across}, but does not support some \code{rlang} features, has some additional features (arguments), and is optimized to work with \emph{collapse}'s, \code{\link{.FAST_FUN}}, yielding much faster computations.
}
\usage{
across(.cols = NULL, .fns, ..., .names = NULL,
       .apply = "auto", .transpose = "auto")

# acr(...) can be used to abbreviate across(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.cols}{select columns using column names and expressions (e.g. \code{a:b} or \code{c(a, b, c:f)}), column indices, logical vectors, or functions yielding a logical value e.g. \code{is.numeric}. \code{NULL} applies functions to all columns except for grouping columns.}
  \item{.fns}{A function, character vector of functions or list of functions. Vectors / lists can be named to yield alternative names in the result (see \code{.names}). This argument is evaluated inside \code{substitute()}, and the content (not the names of vectors/lists) is checked against \code{.FAST_FUN} and \code{.OPERATOR_FUN}. Matching functions receive vectorized execution, other functions are applied to the data in a standard way.}
  \item{\dots}{further arguments to \code{.fns}. Arguments are evaluated in the data environment and split by groups as well (for non-vectorized functions, if of the same length as the data).}
  \item{.names}{controls the naming of computed columns. \code{NULL} generates names of the form \code{coli_funj} if multiple functions are used. \code{.names = TRUE} enables this for a single function, \code{.names = FALSE} disables it for multiple functions (sensible for functions such as \code{.OPERATOR_FUN} that rename columns (if \code{.apply = FALSE})). Setting \code{.names = "flip"} generates names of the form \code{funj_coli}. It is also possible to supply a function with two arguments for column and function names e.g. \code{function(c, f) paste0(f, "_", c)}. Finally, you can supply a custom vector of names which must match \code{length(.cols) * length(.fns)}.}
  \item{.apply}{controls whether functions are applied column-by-column (\code{TRUE}) or to multiple columns at once (\code{FALSE}). The default, \code{"auto"}, does the latter for vectorized functions, which have an efficient data frame method. It can also be sensible to use \code{.apply = FALSE} for non-vectorized functions, especially multivariate functions like \code{\link{lm}} or \code{\link{pwcor}}, or functions renaming the data. See Examples.  }
  \item{.transpose}{with multiple \code{.fns}, \code{.transpose} controls whether the result is ordered first by column, then by function (\code{TRUE}), or vice-versa (\code{FALSE}). \code{"auto"} does the former if all functions yield results of the same dimensions (dimensions may differ if \code{.apply = FALSE}). See Examples.}
}
\note{
\code{across} does not support \emph{purr}-style lambdas, and does not support \code{dplyr}-style predicate functions e.g. \code{across(where(is.numeric), sum)}, simply use \code{across(is.numeric, sum)}. In contrast to \code{dplyr}, you can also compute on grouping columns.

%In general, my mission with \code{collapse} is not to create a \code{dplyr}-clone, but to take some of the useful features and make them robust and fast using base R and C/C++, with the aim of having a stable API. So don't ask me to implement the latest \emph{dplyr} feature, unless you firmly believe it is very useful and will be around 10 years from now.
}

\seealso{
\code{\link{fsummarise}}, \code{\link{fmutate}}, \link[=fast-data-manipulation]{Fast Data Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# Basic (Weighted) Summaries
fsummarise(wlddev, across(PCGDP:GINI, fmean, w = POP))
\donttest{ % No code relying on suggested package
library(magrittr) # Note: Used because |> is not available on older R versions
wlddev \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(PCGDP:GINI, fmean, w = POP))

# Note that for these we don't actually need across...
fselect(wlddev, PCGDP:GINI) \%>\% fmean(w = wlddev$POP, drop = FALSE)
wlddev \%>\% fgroup_by(region, income) \%>\%
    fselect(PCGDP:GINI, POP) \%>\% fmean(POP, keep.w = FALSE)
collap(wlddev, PCGDP + LIFEEX + GINI ~ region + income, w = ~ POP, keep.w = FALSE)

# But if we want to use some base R function that reguires argument splitting...
wlddev \%>\% na_omit(cols = "POP") \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(PCGDP:GINI, weighted.mean, w = POP, na.rm = TRUE))

# Or if we want to apply different functions...
wlddev \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(PCGDP:GINI, list(mu = fmean, sd = fsd), w = POP),
               POP_sum = fsum(POP), OECD = fmean(OECD))
# Note that the above still detects fmean as a fast function, the names of the list
# are irrelevant, but the function name must be typed or passed as a character vector,
# Otherwise functions will be executed by groups e.g. function(x) fmean(x) won't vectorize

# Same, naming in a different way
wlddev \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(PCGDP:GINI, list(mu = fmean, sd = fsd), w = POP, .names = "flip"),
               sum_POP = fsum(POP), OECD = fmean(OECD))

# Or we want to do more advanced things..
# Such as nesting data frames..
qTBL(wlddev) \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(c(PCGDP, LIFEEX, ODA),
               function(x) list(Nest = list(x)),
               .apply = FALSE))
# Or linear models..
qTBL(wlddev) \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(c(PCGDP, LIFEEX, ODA),
               function(x) list(Mods = list(lm(PCGDP ~., x))),
               .apply = FALSE))
# Or cumputing grouped correlation matrices
qTBL(wlddev) \%>\% fgroup_by(region, income) \%>\%
    fsummarise(across(c(PCGDP, LIFEEX, ODA),
      function(x) qDF(pwcor(x), "Variable"), .apply = FALSE))

# Here calculating 1- and 10-year lags and growth rates of these variables
qTBL(wlddev) \%>\% fgroup_by(country) \%>\%
    fmutate(across(c(PCGDP, LIFEEX, ODA), list(L, G),
                   n = c(1, 10), t = year, .names = FALSE))

# Same but variables in different order
qTBL(wlddev) \%>\% fgroup_by(country) \%>\%
    fmutate(across(c(PCGDP, LIFEEX, ODA), list(L, G), n = c(1, 10),
                   t = year, .names = FALSE, .transpose = FALSE))

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
