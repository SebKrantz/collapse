\name{collap}
\alias{advanced-aggregation}
\alias{A5-advanced-aggregation}
\alias{collap}
\alias{collapv}
\alias{collapg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Advanced Data Aggregation
}
\description{
\code{collap} is a fast and versatile multi-purpose data aggregation command.

It performs simple and weighted aggregations, multi-type aggregations automatically applying different functions to numeric and categorical columns, multi-function aggregations applying multiple functions to each column, and fully custom aggregations where the user passes a list mapping functions to columns.

% \code{collap} works with \emph{collapse}'s \link[=fast-statistical-functions]{Fast Statistical Functions}, providing extremely fast conventional and weighted aggregation. It also works with other functions but this does not deliver high speeds on large data and does not support weighted aggregations.

% \code{collap} supports formula and data (i.e. grouping vectors or lists of vectors) input to \code{by}, whereas \code{collapv} allows names and indices of grouping columns to be passed to \code{by}.
}
\usage{
# Main function: allows formula and data input to `by` and `w` arguments
collap(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum,
       custom = NULL, \dots, keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE,
       sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE, return.order = sort,
       method = "auto", parallel = FALSE, mc.cores = 2L,
       return = c("wide","list","long","long_dupl"), give.names = "auto")

# Programmer function: allows column names and indices input to `by` and `w` arguments
collapv(X, by, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum,
        custom = NULL, \dots, keep.by = TRUE, keep.w = TRUE, keep.col.order = TRUE,
        sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE, return.order = sort,
        method = "auto", parallel = FALSE, mc.cores = 2L,
        return = c("wide","list","long","long_dupl"), give.names = "auto")

# Auxiliary function: for grouped data ('grouped_df') input + non-standard evaluation
collapg(X, FUN = fmean, catFUN = fmode, cols = NULL, w = NULL, wFUN = fsum,
        custom = NULL, keep.group_vars = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame, or an object coercible to data frame using \code{\link{qDF}}.}
  \item{by}{for \code{collap}: a one-or two sided formula, i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}, or a atomic vector, list of vectors or \code{\link{GRP}} object used to group \code{X}. For \code{collapv}: names or indices of grouping columns, or a logical vector or selector function such as \code{\link{is_categorical}} selecting grouping columns.}
\item{FUN}{a function, list of functions (i.e. \code{list(fsum, fmean, fsd)} or \code{list(sd = fsd, myfun1 = function(x)..)}), or a character vector of function names, which are automatically applied only to numeric variables.}
\item{catFUN}{same as \code{FUN}, but applied only to categorical (non-numeric) typed columns (\code{\link{is_categorical}}).}
\item{cols}{select columns to aggregate using a function, column names, indices or logical vector. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
\item{w}{weights. Can be passed as numeric vector or alternatively as formula i.e. \code{~ weightvar} in \code{collap} or column name / index etc. i.e. \code{"weightvar"} in \code{collapv}. \code{collapg} supports non-standard evaluations so \code{weightvar} can be indicated without quotes. }
\item{wFUN}{same as \code{FUN}: Function(s) to aggregate weight variable if \code{keep.w = TRUE}. By default the sum of the weights is computed in each group.}

\item{custom}{a named list specifying a fully customized aggregation task. The names of the list are function names and the content columns to aggregate using this function (same input as \code{cols}). For example \code{custom = list(fmean = 1:6, fsd = 7:9, fmode = 10:11)} tells \code{collap} to aggregate columns 1-6 of \code{X} using the mean, columns 7-9 using the standard deviation etc. \emph{Notes}: \code{custom} lets \code{collap} ignore any inputs passed to \code{FUN}, \code{catFUN} or \code{cols}. Since v1.6.0 you can also rename columns e.g. \code{custom = list(fmean = c(newname = "col1", "col2"), fmode = c(newname = 3))}.}
\item{keep.by, keep.group_vars}{logical. \code{FALSE} will omit grouping variables from the output. \code{TRUE} keeps the variables, even if passed externally in a list or vector (unlike other \emph{collapse} functions).}
\item{keep.w}{logical. \code{FALSE} will omit weight variable from the output i.e. no aggregation of the weights. \code{TRUE} aggregates and adds weights, even if passed externally as a vector (unlike other \emph{collapse} functions).}
\item{keep.col.order}{logical. Retain original column order post-aggregation.}
\item{sort, decreasing, na.last, return.order, method}{logical / character. Arguments passed to \code{\link{GRP.default}} and affecting the row-order in the aggregated data frame and the grouping algorithm.}
\item{parallel}{logical. Use \code{\link{mclapply}} instead of \code{lapply} to parallelize the computation at the column level. Not available for Windows.}
\item{mc.cores}{integer. Argument to \code{\link{mclapply}} setting the number of cores to use, default is 2.}
\item{return}{character. Control the output format when aggregating with multiple functions or performing custom aggregation. "wide" (default) returns a wider data frame with added columns for each additional function. "list" returns a list of data frames - one for each function. "long" adds a column "Function" and row-binds the results from different functions using \code{data.table::rbindlist}. "long_dupl" is a special option for aggregating multi-type data using multiple \code{FUN} but only one \code{catFUN} or vice-versa. In that case the format is long and data aggregated using only one function is duplicated. See Examples.}
\item{give.names}{logical. Create unique names of aggregated columns by adding a prefix 'FUN.var'. \code{'auto'} will automatically create such prefixes whenever multiple functions are applied to a column. % By default \code{"."} is used as a separator between 'FUN' and 'var'. It is also possible to choose a different separator by specifying \code{give.names = "_"}, for example.
}
\item{\dots}{additional arguments passed to all functions supplied to \code{FUN}, \code{catFUN}, \code{wFUN} or \code{custom}. Since v1.9.0 these are also split by groups for non-\link[=fast-statistical-functions]{Fast Statistical Functions}. The behavior of \link[=fast-statistical-functions]{Fast Statistical Functions} with unused arguments is regulated by \code{option("collapse_unused_arg_action")} and defaults to \code{"warning"}. \code{collapg} also allows other arguments to \code{collap} except for \code{sort, decreasing, na.last, return.order, method} and \code{keep.by}.}
}

\details{
\code{collap} automatically checks each function passed to it whether it is a \link[=fast-statistical-functions]{Fast Statistical Function} (i.e. whether the function name is contained in \code{.FAST_STAT_FUN}). If the function is a fast statistical function, \code{collap} only does the grouping and then calls the function to carry out the grouped computations (vectorized in C/C++), resulting in high aggregation speeds, even with weights. If the function is not one of \code{.FAST_STAT_FUN}, \code{\link{BY}} is called internally to perform the computation. The resulting computations from each function are put into a list and recombined to produce the desired output format as controlled by the \code{return} argument. This is substantially slower, particularly with many groups.

When setting \code{parallel = TRUE} on a non-windows computer, aggregations will efficiently be parallelized at the column level using \code{\link{mclapply}} utilizing \code{mc.cores} cores. Some \link[=fast-statistical-functions]{Fast Statistical Function} support multithreading i.e. have an \code{nthreads} argument that can be passed to \code{collap}. Using C-level multithreading is much more effective than R-level parallelism, and also works on Windows, but the two should never be combined.

When the \code{w} argument is used, the weights are passed to all functions except for \code{wFUN}. This may be undesirable in settings like \code{collap(data, ~ id, custom = list(fsum = ..., fmean = ...), w = ~ weights)} where we wish to aggregate some columns using the weighted mean, and others using a simple sum or another unweighted statistic. %Since many \link[=fast-statistical-functions]{Fast Statistical Functions} including \code{\link{fsum}} support weights, the above computes a weighted mean and a weighted sum. A couple of workarounds were outlined \href{https://github.com/SebKrantz/collapse/issues/96}{here}, but \emph{collapse} 1.5.0 incorporates an easy solution into \code{collap}:
Therefore it is possible to append \link[=fast-statistical-functions]{Fast Statistical Functions} by \code{_uw} to yield an unweighted computation. So for the above example one can specify: \code{collap(data, ~ id, custom = list(fsum_uw = ..., fmean = ...), w = ~ weights)} to get the weighted mean and the simple sum. \emph{Note} that the \code{_uw} functions are not available for use outside collap. Thus one also needs to quote them when passing to the \code{FUN} or \code{catFUN} arguments, e.g. use \code{collap(data, ~ id, fmean, "fmode_uw", w = ~ weights)}. %\emph{Note} also that it is never necessary for functions passed to \code{wFUN} to be appended like this, as the weights are never used to aggregate themselves.

}
\value{
\code{X} aggregated. If \code{X} is not a data frame it is coerced to one using \code{\link{qDF}} and then aggregated.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
 %}
%\note{
% (1) Additional arguments passed are not split by groups. Weighted aggregations with user defined functions should be done with \code{\link{BY}}, \code{\link{fsummarise}}, or using the \emph{data.table} package.

% (2) Move above...


% (3) The dispatch between using optimized \link[=fast-statistical-functions]{Fast Statistical Functions} performing grouped computations internally or calling \code{BY} to perform split-apply-combine computing is done by matching the function name against \code{.FAST_STAT_FUN}. Thus code like \code{collapse::collap(data, ~ id, collapse::fmedian)} does not yield an optimized computation, as \code{"collapse::fmedian" \%!in\% .FAST_STAT_FUN}. It is sufficient to write \code{collapse::collap(data, ~ id, "fmedian")} to get the desired result when the \emph{collapse} namespace is not attached.

%If you want to perform optimized computations with \code{collap} without loading the pacckage, load the functions beforehand as well, e.g. \code{fmedian <- collapse::fmedian; data, ~ id, fmedian)}. Alternatively it is of course also possible to use \code{collapse::fmedian(collapse::fgroup_by(data, id))}, or something similar...

% \code{collap} by default (\code{keep.by = TRUE, keep.w = TRUE}) preserves all arguments passed to the \code{by} or \code{w} arguments, whether passed in a formula or externally. The names of externally passed vectors and lists are intelligently extracted. So it is possible to write \code{collap(iris, iris$Species)}, and obtain an aggregated data frame with two \code{Species} columns, whereas \code{collap(iris, ~ Species)} only has one \code{Species} column. Similarly for weight vectors passed to \code{w}. In this regard \code{collap} is more sophisticated than other \emph{collapse} functions where preservation of grouping and weight variables is restricted to formula use. For example \code{STD(iris, iris$Species)} does not preserve \code{Species} in the output, whereas \code{STD(iris, ~ Species)} does. This \code{collap} feature is there simply for convenience, for example sometimes a survey is disaggregated into several datasets, and this now allows easy pulling of identifiers or weights from other datasets for aggregations. If all information is available in one dataset, just using formulas is highly recommended.
%}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{fsummarise}}, \code{\link{BY}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## A Simple Introduction --------------------------------------
head(iris)
collap(iris, ~ Species)                                        # Default: FUN = fmean for numeric
collapv(iris, 5)                                               # Same using collapv
collap(iris, ~ Species, fmedian)                               # Using the median
collap(iris, ~ Species, fmedian, keep.col.order = FALSE)       # Groups in-front
collap(iris, Sepal.Width + Petal.Width ~ Species, fmedian)     # Only '.Width' columns
collapv(iris, 5, cols = c(2, 4))                               # Same using collapv
collap(iris, ~ Species, list(fmean, fmedian))                  # Two functions
collap(iris, ~ Species, list(fmean, fmedian), return = "long") # Long format
collapv(iris, 5, custom = list(fmean = 1:2, fmedian = 3:4))    # Custom aggregation
collapv(iris, 5, custom = list(fmean = 1:2, fmedian = 3:4),    # Raw output, no column reordering
        return = "list")
collapv(iris, 5, custom = list(fmean = 1:2, fmedian = 3:4),    # A strange choice..
        return = "long")
collap(iris, ~ Species, w = ~ Sepal.Length)                    # Using Sepal.Length as weights, ..
weights <- abs(rnorm(fnrow(iris)))
collap(iris, ~ Species, w = weights)                           # Some random weights..
collap(iris, iris$Species, w = weights)                        # Note this behavior..
collap(iris, iris$Species, w = weights,
       keep.by = FALSE, keep.w = FALSE)

% \donttest{iris |> fgroup_by(Species) |> collapg()}                        # dplyr style, but faster

## Multi-Type Aggregation --------------------------------------
head(wlddev)                                                    # World Development Panel Data
head(collap(wlddev, ~ country + decade))                        # Aggregate by country and decade
head(collap(wlddev, ~ country + decade, fmedian, ffirst))       # Different functions
head(collap(wlddev, ~ country + decade, cols = is.numeric))     # Aggregate only numeric columns
head(collap(wlddev, ~ country + decade, cols = 9:13))           # Only the 5 series
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade))         # Only GDP and life-expactancy
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade, fsum))   # Using the sum instead
head(collap(wlddev, PCGDP + LIFEEX ~ country + decade, sum,     # Same using base::sum -> slower!
            na.rm = TRUE))
head(collap(wlddev, wlddev[c("country","decade")], fsum,        # Same, exploring different inputs
            cols = 9:10))
head(collap(wlddev[9:10], wlddev[c("country","decade")], fsum))
head(collapv(wlddev, c("country","decade"), fsum))              # ..names/indices with collapv
head(collapv(wlddev, c(1,5), fsum))

g <- GRP(wlddev, ~ country + decade)                            # Precomputing the grouping
head(collap(wlddev, g, keep.by = FALSE))                        # This is slightly faster now
# Aggregate categorical data using not the mode but the last element
head(collap(wlddev, ~ country + decade, fmean, flast))
head(collap(wlddev, ~ country + decade, catFUN = flast,         # Aggregate only categorical data
            cols = is_categorical))


## Weighted Aggregation ----------------------------------------
# We aggregate to region level using population weights
head(collap(wlddev, ~ region + year, w = ~ POP))                # Takes weighted mean for numeric..
# ..and weighted mode for categorical data. The weight vector is aggregated using fsum

head(collap(wlddev, ~ region + year, w = ~ POP,                 # Aggregating weights using sum
            wFUN = list(sum = fsum, max = fmax)))               # and max (corresponding to mode)


## Multi-Function Aggregation ----------------------------------
head(collap(wlddev, ~ country + decade, list(mean = fmean, N = fnobs),  # Saving mean and Nobs
            cols = 9:13))

head(collap(wlddev, ~ country + decade,                         # Same using base R -> slower
            list(mean = mean,
                 N = function(x, \dots) sum(!is.na(x))),
            cols = 9:13, na.rm = TRUE))

lapply(collap(wlddev, ~ country + decade,                       # List output format
       list(mean = fmean, N = fnobs), cols = 9:13, return = "list"), head)

head(collap(wlddev, ~ country + decade,                         # Long output format
     list(mean = fmean, N = fnobs), cols = 9:13, return = "long"))

head(collap(wlddev, ~ country + decade,                         # Also aggregating categorical data,
     list(mean = fmean, N = fnobs), return = "long_dupl"))      # and duplicating it 2 times

head(collap(wlddev, ~ country + decade,                         # Now also using 2 functions on
     list(mean = fmean, N = fnobs), list(mode = fmode, last = flast),   # categorical data
            keep.col.order = FALSE))

head(collap(wlddev, ~ country + decade,                         # More functions, string input,
            c("fmean","fsum","fnobs","fsd","fvar"),             # parallelized execution
            c("fmode","ffirst","flast","fndistinct"),           # (choose more than 1 cores,
            parallel = TRUE, mc.cores = 1L,                     # depending on your machine)
            keep.col.order = FALSE))


## Custom Aggregation ------------------------------------------
head(collap(wlddev, ~ country + decade,                         # Custom aggregation
            custom = list(fmean = 11:13, fsd = 9:10, fmode = 7:8)))

head(collap(wlddev, ~ country + decade,                         # Using column names
            custom = list(fmean = "PCGDP", fsd = c("LIFEEX","GINI"),
                          flast = "date")))

head(collap(wlddev, ~ country + decade,                         # Weighted parallelized custom
            custom = list(fmean = 9:12, fsd = 9:10,             # aggregation
                          fmode = 7:8), w = ~ POP,
            wFUN = list(fsum, fmax),
            parallel = TRUE, mc.cores = 1L))

head(collap(wlddev, ~ country + decade,                         # No column reordering
            custom = list(fmean = 9:12, fsd = 9:10,
                          fmode = 7:8), w = ~ POP,
            wFUN = list(fsum, fmax),
            parallel = TRUE, mc.cores = 1L, keep.col.order = FALSE))

## Piped Use --------------------------------------------------
iris |> fgroup_by(Species) |> collapg()
wlddev |> fgroup_by(country, decade) |> collapg() |> head()
wlddev |> fgroup_by(region, year) |> collapg(w = POP) |> head()
wlddev |> fgroup_by(country, decade) |> collapg(fmedian, flast) |> head()
wlddev |> fgroup_by(country, decade) |>
  collapg(custom = list(fmean = 9:12, fmode = 5:7, flast = 3)) |> head()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line
