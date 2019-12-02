\name{collap}
\alias{advanced-data-aggregation}
\alias{collap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Advanced Data Aggregation
}
\description{
\code{collap} is a fast and easy to use multi-purpose data aggregation command.

It can aggregate data with multiple data types, aggregate (parallelized) using multiple functions to several output formats, and perform (parallelized) fully customized aggregations where the user decides which variables are aggregated using which functions. \code{collap} is made compatible with \code{collapse}'s \link[=fast-statistical-functions]{Fast Statistical Functions}, allowing for extremely fast conventional and weighted aggregation, but even with \code{base} functions it is significantly faster than \code{stats::aggregate}.
}
\usage{
collap(X, by, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
       keep.by = TRUE, keep.col.order = TRUE, sort.row = TRUE, parallel = FALSE,
       mc.cores = 1L, return = c("wide","list","long","long_dupl"),
       give.names = "auto", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data.frame, or an object coercible to data.frame using \code{\link{qDF}}.}
  \item{by}{a one-or two sided formula, i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}, or alternativaley a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{X}.}
\item{FUN}{a function, list of functions (i.e. \code{list(fsum, fmean, fsd)} or \code{list(myfun1 = function(x).., sd = sd)}), or a character vector of function names, which are automatically applied only to numeric variables.}
\item{catFUN}{same as \code{FUN}, but applied only to categorical (non-numeric) typed columns (\code{\link{is.categorical}}).}
\item{cols}{select columns to aggregate using a function, column names or indices. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
\item{custom}{a named list specifying a fully customized aggregation task. The names of the list are function names and the content column-names or column-indices of associated columns. For example \code{custom = list(fmean = 1:6, fsd = 7:9, fmode = 10:11)} tells \code{collap} to aggregate columns 1-6 of \code{X} using the mean, columns 7-9 using the standard deviation etc. \emph{Note}: \code{custom} lets \code{collap} ignore any inputs passed to \code{FUN}, \code{catFUN} or \code{cols}.}
\item{keep.by}{logical. \code{FALSE} will omit grouping variables from the output.}
\item{keep.col.order}{logical. Retain original column order post-aggregtion.}
\item{sort.row}{logical. Sort rows by the groups.}
\item{parallel}{logical. Use \code{parallel::mclapply} instead of \code{lapply} for multi-function or custom aggregation.}
\item{mc.cores}{integer. Argument to \code{parallel::mclapply} setting the number of cores to use.}
\item{return}{character. Control the output format when aggregating with multiple functions or performing custom aggregation. "wide" (default) returns a wider data frame with added columns for each additional function. "list" returns a list of \code{data.frame}'s - one for each function. "long" adds a column "Function" and row-binds the results from differnt functions using \code{data.table::rbindlist}. "long.dupl" is a special option for aggregating multi-type data using multiple \code{FUN} but only one \code{catFUN} or vice-versa. In that case the format is long and data aggregated using only one function is duplicated. See Examples to understand this!}
\item{give.names}{logical. Create unique names of aggregated columnns by adding a prefix 'FUN.'. 'auto' will automatically create such prefixes whenever multiple functions are applied to a column or \code{custom} is used.}
\item{...}{additional arguments passed to all functions supplied to \code{FUN}, \code{catFUN} or \code{custom}. }
}

\details{
\code{collap} automatically checks each function passed to it whether it is a \link[=fast-statistical-functions]{Fast Statistical Function} (i.e. whether the function name is contained in \code{.FAST_STAT_FUN}). If the function is a fast function, \code{collap} only does the grouping and then calls the function to carry out the grouped computations. If the function is not one of \code{.FAST_STAT_FUN}, \code{\link{BY}} is called internally to perform the computation. The resulting computataions from each function are put into a list and recombined to produce the desired output format as controlled by the \code{return} argument. When multiple functions are used with \code{collap}, setting \code{parallel = TRUE} and the number of cores with \code{mc.cores} will instruct \code{collap} to execute these function calls in parallel using \code{parallel::mclapply}. If only a single function is used which is not a \code{.FAST_STAT_FUN}, the \code{parallel} and \code{mc.cores} arguments are handed down to \code{\link{BY}}. See Examples.
}
\value{
\code{X} aggregated by \code{by}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
 %}
% \note{
%%  ~~further notes~~
 %}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BY}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## World Development Panel Data

# Simple and Multi-Type Aggregation ----------------------------
head(collap(wlddev, ~ country + decade))                        # Aggregate this data by country and decade
head(collap(wlddev, ~ country + decade, cols = is.numeric))     # Aggregate only numeric columns
head(collap(wlddev, ~ country + decade, cols = 9:12))           # Only the 4 series
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade))         # Only GDP and life-expactancy
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade, fsum))   # Using the sum instead
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade, sum,     # Same using base::sum -> slower!!
            na.rm = TRUE))
head(collap(wlddev, wlddev[c("country","decade")], fsum,        # same thing, now exploring different input formats
            cols = 9:10))
head(collap(wlddev[9:10], wlddev[c("country","decade")], fsum))

g <- GRP(wlddev, ~ country + decade)                            # Precomputing the grouping can help if I want to do multiple grouped operations on the same data !
head(collap(wlddev, g, keep.by = FALSE))

head(collap(wlddev, ~ country + decade, fmean, flast))          # Aggregate categorical data using not the mode but the last element
head(collap(wlddev, ~ country + decade, catFUN = flast,         # Aggregate only categorical data
            cols = is.categorical))

# Weighted aggregation -----------------------------------------
weights <- abs(rnorm(nrow(wlddev)))                             # Adding a random weight vector - may also have missing values
head(collap(wlddev, ~ country + decade, w = weights))           # Takes weighted mean for numeric and weighted mode for categorical data.

# Multi-Function Aggregation -----------------------------------
head(collap(wlddev, ~ country + decade, list(fmean, fNobs),     # Saving mean and number of observations
            cols = 9:12))

head(collap(wlddev, ~ country + decade,                         # same using base R -> substantially slower
            list(mean = mean,
                 Nobs = function(x,...) sum(!is.na(x))),
            cols = 9:12, na.rm = TRUE))

head(collap(wlddev, ~ country + decade,                         # list output format
            list(fmean, fNobs), cols = 9:12, return = "list"))

head(collap(wlddev, ~ country + decade,                         # long output format
            list(fmean, fNobs), cols = 9:12, return = "long"))

head(collap(wlddev, ~ country + decade,                         # also aggregating the categorical data, and duplicating it 2 times
            list(fmean, fNobs), return = "long_dupl"))

head(collap(wlddev, ~ country + decade,                         # now also using 2 functions on categorical data
            list(fmean, fNobs), list(fmode, flast),
            keep.col.order = FALSE))

head(collap(wlddev, ~ country + decade,                         # more functions, string input format, parallelized execution (choose more than 1 cores, depending on your machine)
            c("fmean","fsum","fNobs","fsd","fvar"),
            c("fmode","ffirst","flast","fNdistinct"),
            parallel = TRUE, mc.cores = 1L,
            keep.col.order = FALSE))

# Custom Aggregation -------------------------------------------
head(collap(wlddev, ~ country + decade,                         # custom aggregation
            custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)))

head(collap(wlddev, ~ country + decade,                         # using column names
            custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"),
                          flast = "date")))

head(collap(wlddev, ~ country + decade,                         # weighted parallelized custom aggregation
            custom = list(fmean = 9:12, fsd = 9:10,
                          fmode = 7:8), w = weights,
            parallel = TRUE, mc.cores = 1L))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
