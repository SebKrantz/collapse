\name{collap}
\alias{collap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Advanced Data Aggregation
}
\description{
\code{collap} is a fast and easy to use multi-purpose data aggregation command.

It can aggregate data with multiple data types, aggregate (parallelized) using multiple functions to several output formats, and perform (parallelized) fully customized aggregations where the user decides which variables are aggregated using which functions. \code{collap} is made compatible with \code{collapse}'s fast statistical functions, allowing for extremely fast conventional and weighted aggregation (significantly faster than \code{data.table} on smaller data).
}
\usage{
collap(X, by, FUN = fmean, catFUN = fmode, cols = NULL, custom = NULL,
       keep.by = TRUE, sort.col = TRUE, sort.row = TRUE, parallel = FALSE,
       mc.cores = 1L, multi.FUN.out = c("wide","list","long","long_dupl"),
       give.names = "auto", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data.frame, or an object coercible to data.frame using \code{qDF}.}
    \item{by}{a one-or two sided formula, i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}, or alternativaley a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{X}.}
\item{FUN}{A function, list of functions (i.e. \code{list(fsum, fmean, fsd)} or \code{list(myfun1 = function(x).., sd = sd)}), or character vector of function names automatically applied only to numeric variables.}
\item{catFUN}{same as \code{FUN}, but applied only to categorical (non-numeric) typed columns.}
\item{cols}{select columns to aggregate using a function, column names or indices. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
\item{custom}{a named list specifying a fully customized aggregation task. The names of the list are function names and the content column-names or column-indices of associated columns. For example \code{custom = list(fmean = 1:6, fsd = 7:9, fmode = 10:11)} tells \code{collap} to aggregate columns 1-6 of \code{x} using the mean, columns 7-9 using the standard deviation etc. \emph{Note}: \code{custom} lets \code{collap} ignore any inputs passed to \code{FUN}, \code{catFUN} or \code{cols}.}
\item{keep.by}{logical. \code{FALSE} will omit grouping variables from the output.}
\item{sort.col}{logical. Retain original column order post-aggregtion.}
\item{sort.row}{logical. Sort rows by the groups.}
\item{parallel}{logical. Use \code{parallel::mclapply} instead of \code{lapply} for multi-function or custom aggregation.}
\item{mc.cores}{integer. Argument to \code{parallel::mclapply} setting the number of cores to use.}
\item{multi.FUN.out}{character. Control the output format when aggregating with multiple functions or performing custom aggregation. "wide" (default) returns a wider data frame with added columns for each additional function. "list" returns a list of \code{data.frames} - one for each function. "long" adds a column "Function" and row-binds the results from differnt functions using \code{data.table::rbindlist}. "long.dupl" is a special option for aggregating multi-type data using multiple \code{FUN} but only one \code{catFUN} or vice-versa. In that case the format is long and data aggregated using only one function is duplicated. See Examples to understand this!}
\item{...}{additional arguments passed to all functions supplied to \code{FUN}, \code{catFUN} or \code{custom}. }
}

\details{
\code{collap} automatically checks each function passed to it whether it is a \code{collapse} function (i.e. whether the function name is one of "fmean", "fmedian", "fsum", "fprod", "fsd", "fvar", "fmin", "fmax", "fMode", "fNobs" or "fNunique") or other function. If the function is a fast \code{collapse} function, \code{collap} only computes the groups and passes them to the fast functions which carry out the grouped computations. Afterwards the results are recombined to the desired output format. If functions are not fast \code{collapse::BY} is called internally. This is still faster than \code{base::aggregate} but significantly slower than using \code{collapse} functions.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## World Development Panel Data
View(wlddev)                             # See also ?wlddev

# Simple and Multi-Type Aggregation ------------
View(collap(wlddev, ~ country + decade))                    # Aggregate this data by country and decade
View(collap(wlddev, ~ country + decade, cols = is.numeric)) # Aggregate only numeric columns
View(collap(wlddev, ~ country + decade, cols = 9:12))       # Only the 4 series
View(collap(wlddev, PCGDP + LIFEEX ~ country + decade))     # Only GDP and life-expactancy
View(collap(wlddev, PCGDP + LIFEEX ~ country + decade, fsum)) # Using the sum instead
View(collap(wlddev, PCGDP + LIFEEX ~ country + decade, sum, na.rm = TRUE)) # Same using base::sum -> slower!!
View(collap(wlddev, wlddev[c("country","decade")], fsum, cols = 9:10)) # same thing, now exploring different input formats
View(collap(wlddev[9:10], wlddev[c("country","decade")], fsum))

g <- GRP(wlddev, ~ country + decade)  # Precomputing the grouping can help if I want to do multiple grouped operations on the same data !
View(collap(wlddev, g, keep.by = FALSE))

View(collap(wlddev, ~ country + decade, fmean, flast)) # Aggregate categorical data using not the mode but the last element
                                                       # Aggregate only categorical data:
View(collap(wlddev, ~ country + decade, catFUN = flast, cols = is.categorical))

# Weighted aggregation ------------
weights <- abs(rnorm(nrow(wlddev)))                   # Adding a random weight vector - may also have missing values
View(collap(wlddev, ~ country + decade, w = weights)) # takes weighted mean for numeric and weighted mode for categorical data.

# Multi-Function Aggregation ------------
View(collap(wlddev, ~ country + decade, list(fmean, fNobs), cols = 9:12)) # Saving mean and number of observations

View(collap(wlddev, ~ country + decade,  # same using base R -> substantially slower
            list(mean = mean, Nobs = function(x,...) sum(!is.na(x))),
            cols = 9:12, na.rm = TRUE))

View(collap(wlddev, ~ country + decade,  # list output format
            list(fmean, fNobs), cols = 9:12,
            multi.FUN.out = "list"))

View(collap(wlddev, ~ country + decade,  # long output format
            list(fmean, fNobs), cols = 9:12,
            multi.FUN.out = "long"))

View(collap(wlddev, ~ country + decade,  # also aggregating the categorical data, and duplicating it 2 times
            list(fmean, fNobs),
            multi.FUN.out = "long_dupl"))

View(collap(wlddev, ~ country + decade,  # now also using 2 functions on categorical data
            list(fmean, fNobs), list(fmode, flast)))

View(collap(wlddev, ~ country + decade,  # parallelized execution (not very useful on small data like this)
            list(fmean, fNobs), list(fmode, flast), parallel = TRUE))

View(collap(wlddev, ~ country + decade,  # more functions, string input format
            c("fmean","fsum","fNobs","fsd","fvar"),
            c("fmode","ffirst","flast","fNdistinct"), parallel = TRUE))

# Custom Aggregation ------------
View(collap(wlddev, ~ country + decade,  # custom aggregation
            custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8)))

View(collap(wlddev, ~ country + decade,  # using column names
            custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"), flast = "date")))

View(collap(wlddev, ~ country + decade,  # parallelized custom aggregation
            custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8), parallel = TRUE))

View(collap(wlddev, ~ country + decade,  # weighted parallelized custom aggregation
            custom = list(fmean = 9:12, fsd = 9:10, fmode = 7:8), w = weights, parallel = TRUE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
