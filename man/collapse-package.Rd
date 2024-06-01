\name{collapse-package}
\alias{collapse-package}
\alias{collapse}
\docType{package}
\title{ %   \code{collapse:}
Advanced and Fast Data Transformation
}
\description{
\emph{collapse} is a C/C++ based package for data transformation and statistical computing in R. Its aims are:
\itemize{
\item To facilitate complex data transformation, exploration and computing tasks in R.
\item To help make R code fast, flexible, parsimonious and programmer friendly. % \emph{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, its not a misc package
}
It is made compatible with the \emph{tidyverse}, \emph{data.table}, \emph{sf}, \emph{units}, \emph{xts/zoo}, and the \emph{plm} approach to panel data.

}
\section{Getting Started}{
Read the short \href{https://sebkrantz.github.io/collapse/articles/collapse_documentation.html}{vignette} on documentation resources, and check out the built in \link[=collapse-documentation]{documentation}.  % A careful consideration of the \href{https://raw.githubusercontent.com/SebKrantz/cheatsheets/master/collapse.pdf}{cheat sheet} is recommended for quick starters.
% or read the \href{https://sebkrantz.github.io/collapse/articles/collapse_intro.html}{introductory vignette}. % All vignettes can be accessed on the \href{https://sebkrantz.github.io/collapse/}{package website}. A cheatsheet is available \href{https://raw.githubusercontent.com/SebKrantz/cheatsheets/master/collapse.pdf}{here}. A compact introduction for quick-starters is provided in the examples section below.
}


% \section{Key Features} {
% \cr
% \bold{Key Features:}    % (In more detail in \link[=collapse-documentation]{Collapse Overview}) % Key functionality: % Key areas Key topics addressed by \emph{collapse} are: % where \emph{collapse} offers innovative solutions are:

% \tabular{lll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
% \enumerate{  % (1) \tab\tab \dots \cr

%                 \item  \emph{Advanced data programming}: A full set of fast statistical functions
%         supporting grouped and weighted computations on vectors, matrices and
%         data frames. Fast (ordered) and programmable grouping, factor
%         generation, manipulation of data frames and data object conversions. %  select, subset, transform, replace, add and delete data frame columns.
%                  \item \emph{Advanced aggregation}: Fast and easy multi-data-type, multi-function,
%         weighted, parallelized and fully customized data aggregation.
%                  \item \emph{Advanced transformations}: Fast (grouped, weighted) replacing and
%         sweeping out of statistics, scaling / standardizing, centering (i.e.
%         between and within transformations), higher-dimensional centering
%         (i.e. multiple fixed effects transformations), linear
%         prediction and partialling-out.
%                  \item \emph{Advanced time-computations}: Fast (sequences of) lags / leads, and
%         (lagged / leaded, iterated, quasi-, log-) differences and growth
%         rates on (unordered) time series and panel data. Multivariate auto,
%         partial and cross-correlation functions for panel data.
%         Panel data to (ts-)array conversions.
%                  \item \emph{List processing}: (Recursive) list search / identification, extraction /
%         subsetting, data-apply, and generalized row-binding / unlisting in 2D.
%                  \item \emph{Advanced data exploration}: Fast (grouped, weighted, panel-decomposed)
%         summary statistics for complex multilevel / panel data.
% }
% }
% }

\section{Details}{
% Put this below bullet points ?? % \emph{collapse} provides a carefully conceived % \emph{collapse} provides a compact set of functions % organized roughly into several topics
\emph{collapse} provides an integrated suite of statistical and data manipulation functions that greatly extend and enhance the capabilities of base R. In a nutshell, \emph{collapse} provides:

\itemize{
\item Fast C/C++ based (grouped, weighted) computations embedded in highly optimized R code.
\item More complex statistical, time series / panel data and recursive (list-processing) operations.
\item A flexible and generic approach supporting and preserving many R objects. % standard  % (S3 generic statistical functions, class/attribute preservation). % , compatibility with \emph{dplyr}, \emph{plm} and \emph{data.table}
\item Optimized programming in standard and non-standard evaluation.
}
% To explain this a bit:
The statistical functions in \emph{collapse} are S3 generic with core methods for vectors, matrices and data frames, and internally support grouped and weighted computations carried out in C/C++. %Thus functions need only be called once for column-wise and/or grouped computations, providing a lot of extra speed and full support for sampling weights.
%R code is strongly optimized and inputs are swiftly passed to compiled C/C++ code, %, with further checks run at that level.
%where various efficient algorithms are implemented. %This approach enables flexible and parsimonious programming and data manipulation in R at high speeds.

%when passed to a \emph{collapse} statistical function together with a suitable data object % To provide extra speed and programmability
%To facilitate efficient programming, core S3 methods, grouping and ordering functionality and some C-level helper functions can be accessed by the user.

%For example \code{GRP()} creates grouping objects directly passed to C++ by statistical functions. \code{fgroup_by()} attaches these objects to a data frame, yielding efficient chained calls when combined with \emph{magrittr} pipes, \link[=fast-data-manipulation]{fast manipulation functions} and \link[=fast-statistical-functions]{fast statistical functions}.

%Performance gains are also realized when grouping with factors, or computing on grouped (\emph{dplyr}) or panel data (\emph{plm}) frames.

%Hence \emph{collapse} enables optimized programming and data manipulation in both standard and non-standard evaluation.

 %The function \code{fgroup_by} can be used to efficiently create a grouped tibble inside dplyr-like chained calls (constructed with \emph{magrittr} pipes, fast manipulation functions like \code{fselect}, \code{fsubset}, \code{ftransform} and fast statistical functions). Thus \emph{collapse} enables optimized programming in both standard and non-standard evaluation.


% attributes of atomic objects are preserved if the dimensions are unaltered by the computation, and data frame attributes are always preserved,

 Additional methods and C-level features enable broad based compatibility with \emph{dplyr} (grouped tibble), \emph{data.table}, \emph{sf} and \emph{plm} panel data classes. Functions and core methods seek to preserve object attributes (including column attributes such as variable labels), ensuring flexibility and effective workflows with a very broad range of R objects (including most time-series classes). See also the \href{https://sebkrantz.github.io/collapse/articles/collapse_object_handling.html}{vignette} on \emph{collapse}'s handling of R objects.

 Missing values are efficiently skipped at C/C++ level. The package default is \code{na.rm = TRUE}. This can be changed using \code{\link[=set_collapse]{set_collapse(na.rm = FALSE)}}. Missing weights are generally supported. % Core functionality and all statistical functions / computations are tested with > 20,000 unit tests for Base R equivalence, exempting some improvements (e.g. \code{fsum(NA, na.rm = TRUE)} evaluates to \code{NA}, not 0 (unless \code{fill = TRUE}), similarly for \code{fmin} and \code{fmax}; no \code{NaN} values are generated from computations involving \code{NA} values). %Generic functions provide some \link[=collapse-options]{security} against silent swallowing of arguments.

 %Hence they also handle various date and time series classes etc., and can easily be integreated into most approaches to data transformation.

 %A global \code{option("collapse_unused_arg_action")} can be set to regulate the behavior when unused arguments are passed to a generic function. The default is to issue a warning.


\emph{collapse} installs with a built-in hierarchical \link[=collapse-documentation]{documentation} facilitating the use of the package. % The \href{https://sebkrantz.github.io/collapse/articles/index.html}{vignettes} are complimentary and also follow a more structured approach. % and \code{base/stats}. % extra methods warrant, provide ?

%Apart from general performance considerations, \emph{collapse} excels at applications involving fast panel data transformations and techniques, fast weighted computations (e.g. weighted aggregation), fast programming and aggregation with categorical and mixed-type data, fast programming with (multivariate) time series, and programming with lists of data objects.

%Other broad areas are fast grouped and weighted programming to implement new statistical techniques and packages, and fast data manipulation code (i.e. server-side for \code{shiny} apps).

% The package largely avoids non-standard evaluation and exports core methods for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data frames, \dots) % high computation  %(aggregation and transformations ~10x \emph{data.table} on data <1 Mio obs.).

% Beyond speed, flexibility and parsimony in coding, a central objective of \emph{collapse} is to facilitate advanced / complex operations on data.
The package is coded both in C and C++ and built with \emph{Rcpp}, but also uses C/C++ functions from \emph{data.table}, \emph{kit}, \emph{fixest}, \emph{weights}, \emph{stats} and \emph{RcppArmadillo / RcppEigen}. % For the moment \emph{collapse} does not utilize low-level parallelism (such as OpenMP).
% \emph{collapse} is built with \code{Rcpp} and imports \code{C} functions from \emph{data.table}, \emph{lfe} and \emph{stats}. %, and uses \code{ggplot2} visualizations.


}
\section{Author(s)}{ % \author{
\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}

Other contributors from packages \emph{collapse} utilizes:
\itemize{
\item Matt Dowle, Arun Srinivasan and contributors worldwide (\emph{data.table})
\item Dirk Eddelbuettel and contributors worldwide (\emph{Rcpp}, \emph{RcppArmadillo}, \emph{RcppEigen})
\item Morgan Jacob (\emph{kit})
\item Laurent Berge (\emph{fixest})
\item Josh Pasek (\emph{weights})
\item R Core Team and contributors worldwide (\emph{stats})
}
I thank many people from diverse fields for helpful answers on Stackoverflow, Joris Meys for encouraging me and helping to set up the \href{https://github.com/SebKrantz/collapse}{GitHub repository} for \emph{collapse}, and many other people for feature requests and helpful suggestions.
}
\section{Developing / Bug Reporting}{
\itemize{
\item Please report issues at \url{https://github.com/SebKrantz/collapse/issues}.
\item Please send pull-requests to the 'development' branch of the repository.
}
}

%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
%\seealso{
%  Optional links to other man pages
%}

\examples{
## Note: this set of examples is is certainly non-exhaustive and does not
## showcase many recent features, but remains a very good starting point

## Let's start with some statistical programming
v <- iris$Sepal.Length
d <- num_vars(iris)    # Saving numeric variables
f <- iris$Species      # Factor

# Simple statistics
fmean(v)               # vector
fmean(qM(d))           # matrix (qM is a faster as.matrix)
fmean(d)               # data.frame

# Preserving data structure
fmean(qM(d), drop = FALSE)     # Still a matrix
fmean(d, drop = FALSE)         # Still a data.frame

# Weighted statistics, supported by most functions...
w <- abs(rnorm(fnrow(iris)))
fmean(d, w = w)

# Grouped statistics...
fmean(d, f)

# Groupwise-weighted statistics...
fmean(d, f, w)

# Simple Transformations...
head(fmode(d, TRA = "replace"))    # Replacing values with the mode
head(fmedian(d, TRA = "-"))        # Subtracting the median
head(fsum(d, TRA = "\%"))           # Computing percentages
head(fsd(d, TRA = "/"))            # Dividing by the standard-deviation (scaling), etc...

# Weighted Transformations...
head(fnth(d, 0.75, w = w, TRA = "replace"))  # Replacing by the weighted 3rd quartile

# Grouped Transformations...
head(fvar(d, f, TRA = "replace"))  # Replacing values with the group variance
head(fsd(d, f, TRA = "/"))         # Grouped scaling
head(fmin(d, f, TRA = "-"))        # Setting the minimum value in each species to 0
head(fsum(d, f, TRA = "/"))        # Dividing by the sum (proportions)
head(fmedian(d, f, TRA = "-"))     # Groupwise de-median
head(ffirst(d, f, TRA = "\%\%"))     # Taking modulus of first group-value, etc. ...

# Grouped and weighted transformations...
head(fsd(d, f, w, "/"), 3)         # weighted scaling
head(fmedian(d, f, w, "-"), 3)     # subtracting the weighted group-median
head(fmode(d, f, w, "replace"), 3) # replace with weighted statistical mode

## Some more advanced transformations...
head(fbetween(d))                             # Averaging (faster t.: fmean(d, TRA = "replace"))
head(fwithin(d))                              # Centering (faster than: fmean(d, TRA = "-"))
head(fwithin(d, f, w))                        # Grouped and weighted (same as fmean(d, f, w, "-"))
head(fwithin(d, f, w, mean = 5))              # Setting a custom mean
head(fwithin(d, f, w, theta = 0.76))          # Quasi-centering i.e. d - theta*fbetween(d, f, w)
head(fwithin(d, f, w, mean = "overall.mean")) # Preserving the overall mean of the data
head(fscale(d))                               # Scaling and centering
head(fscale(d, mean = 5, sd = 3))             # Custom scaling and centering
head(fscale(d, mean = FALSE, sd = 3))         # Mean preserving scaling
head(fscale(d, f, w))                         # Grouped and weighted scaling and centering
head(fscale(d, f, w, mean = 5, sd = 3))       # Custom grouped and weighted scaling and centering
head(fscale(d, f, w, mean = FALSE,            # Preserving group means
            sd = "within.sd"))                # and setting group-sd to fsd(fwithin(d, f, w), w = w)
head(fscale(d, f, w, mean = "overall.mean",   # Full harmonization of group means and variances,
            sd = "within.sd"))                # while preserving the level and scale of the data.

head(get_vars(iris, 1:2))                      # Use get_vars for fast selecting, gv is shortcut
head(fhdbetween(gv(iris, 1:2), gv(iris, 3:5))) # Linear prediction with factors and covariates
head(fhdwithin(gv(iris, 1:2), gv(iris, 3:5)))  # Linear partialling out factors and covariates
ss(iris, 1:10, 1:2)                            # Similarly fsubset/ss for fast subsetting rows

# Simple Time-Computations..
head(flag(AirPassengers, -1:3))                # One lead and three lags
head(fdiff(EuStockMarkets,                     # Suitably lagged first and second differences
      c(1, frequency(EuStockMarkets)), diff = 1:2))
head(fdiff(EuStockMarkets, rho = 0.87))        # Quasi-differences (x_t - rho*x_t-1)
head(fdiff(EuStockMarkets, log = TRUE))        # Log-differences
head(fgrowth(EuStockMarkets))                  # Exact growth rates (percentage change)
head(fgrowth(EuStockMarkets, logdiff = TRUE))  # Log-difference growth rates (percentage change)
# Note that it is not necessary to use factors for grouping.
fmean(gv(mtcars, -c(2,8:9)), mtcars$cyl) # Can also use vector (internally converted using qF())
fmean(gv(mtcars, -c(2,8:9)),
      gv(mtcars, c(2,8:9)))              # or a list of vector (internally grouped using GRP())
g <- GRP(mtcars, ~ cyl + vs + am)        # It is also possible to create grouping objects
print(g)                                 # These are instructive to learn about the grouping,
plot(g)                                  # and are directly handed down to C++ code
fmean(gv(mtcars, -c(2,8:9)), g)          # This can speed up multiple computations over same groups
fsd(gv(mtcars, -c(2,8:9)), g)

# Factors can efficiently be created using qF()
f1 <- qF(mtcars$cyl)                     # Unlike GRP objects, factors are checked for NA's
f2 <- qF(mtcars$cyl, na.exclude = FALSE) # This can however be avoided through this option
class(f2)                                # Note the added class

library(microbenchmark)
microbenchmark(fmean(mtcars, f1), fmean(mtcars, f2)) # A minor difference, larger on larger data

with(mtcars, finteraction(cyl, vs, am))  # Efficient interactions of vectors and/or factors
finteraction(gv(mtcars, c(2,8:9)))       # .. or lists of vectors/factors

# Simple row- or column-wise computations on matrices or data frames with dapply()
dapply(mtcars, quantile)                 # column quantiles
dapply(mtcars, quantile, MARGIN = 1)     # Row-quantiles
  # dapply preserves the data structure of any matrices / data frames passed
  # Some fast matrix row/column functions are also provided by the matrixStats package
# Similarly, BY performs grouped comptations
BY(mtcars, f2, quantile)
BY(mtcars, f2, quantile, expand.wide = TRUE)
# For efficient (grouped) replacing and sweeping out computed statistics, use TRA()
sds <- fsd(mtcars)
head(TRA(mtcars, sds, "/"))     # Simple scaling (if sd's not needed, use fsd(mtcars, TRA = "/"))

microbenchmark(TRA(mtcars, sds, "/"), sweep(mtcars, 2, sds, "/")) # A remarkable performance gain..

sds <- fsd(mtcars, f2)
head(TRA(mtcars, sds, "/", f2)) # Groupd scaling (if sd's not needed: fsd(mtcars, f2, TRA = "/"))

# All functions above perserve the structure of matrices / data frames
# If conversions are required, use these efficient functions:
mtcarsM <- qM(mtcars)                      # Matrix from data.frame
head(qDF(mtcarsM))                         # data.frame from matrix columns
head(mrtl(mtcarsM, TRUE, "data.frame"))    # data.frame from matrix rows, etc..
head(qDT(mtcarsM, "cars"))                 # Saving row.names when converting matrix to data.table
head(qDT(mtcars, "cars"))                  # Same use a data.frame
\donttest{ % No code relying on suggested packages and the tidyverse, also need to reduce execution time for CRAN
## Now let's get some real data and see how we can use this power for data manipulation
head(wlddev) # World Bank World Development Data: 216 countries, 61 years, 5 series (columns 9-13)

# Starting with some discriptive tools...
namlab(wlddev, class = TRUE)           # Show variable names, labels and classes
fnobs(wlddev)                          # Observation count
pwnobs(wlddev)                         # Pairwise observation count
head(fnobs(wlddev, wlddev$country))    # Grouped observation count
fndistinct(wlddev)                     # Distinct values
descr(wlddev)                          # Describe data
varying(wlddev, ~ country)             # Show which variables vary within countries
qsu(wlddev, pid = ~ country,           # Panel-summarize columns 9 though 12 of this data
    cols = 9:12, vlabels = TRUE)       # (between and within countries)
qsu(wlddev, ~ region, ~ country,       # Do all of that by region and also compute higher moments
    cols = 9:12, higher = TRUE)        # -> returns a 4D array
qsu(wlddev, ~ region, ~ country, cols = 9:12,
    higher = TRUE, array = FALSE) |>                           # Return as a list of matrices..
unlist2d(c("Variable","Trans"), row.names = "Region") |> head()# and turn into a tidy data.frame
pwcor(num_vars(wlddev), P = TRUE)                           # Pairwise correlations with p-value
pwcor(fmean(num_vars(wlddev), wlddev$country), P = TRUE)    # Correlating country means
pwcor(fwithin(num_vars(wlddev), wlddev$country), P = TRUE)  # Within-country correlations
psacf(wlddev, ~country, ~year, cols = 9:12)                 # Panel-data Autocorrelation function
pspacf(wlddev, ~country, ~year, cols = 9:12)                # Partial panel-autocorrelations
psmat(wlddev, ~iso3c, ~year, cols = 9:12) |> plot()         # Convert panel to 3D array and plot

## collapse offers a few very efficent functions for data manipulation:
# Fast selecting and replacing columns
series <- get_vars(wlddev, 9:12)     # Same as wlddev[9:12] but 2x faster
series <- fselect(wlddev, PCGDP:ODA) # Same thing: > 100x faster than dplyr::select
get_vars(wlddev, 9:12) <- series     # Replace, 8x faster wlddev[9:12] <- series + replaces names
fselect(wlddev, PCGDP:ODA) <- series # Same thing

# Fast subsetting
head(fsubset(wlddev, country == "Ireland", -country, -iso3c))
head(fsubset(wlddev, country == "Ireland" & year > 1990, year, PCGDP:ODA))
ss(wlddev, 1:10, 1:10) # This is an order of magnitude faster than wlddev[1:10, 1:10]

# Fast transforming
head(ftransform(wlddev, ODA_GDP = ODA / PCGDP, ODA_LIFEEX = sqrt(ODA) / LIFEEX))
settransform(wlddev, ODA_GDP = ODA / PCGDP, ODA_LIFEEX = sqrt(ODA) / LIFEEX) # by reference
head(ftransform(wlddev, PCGDP = NULL, ODA = NULL, GINI_sum = fsum(GINI)))
head(ftransformv(wlddev, 9:12, log))                   # Can also transform with lists of columns
head(ftransformv(wlddev, 9:12, fscale, apply = FALSE)) # apply = FALSE invokes fscale.data.frame
settransformv(wlddev, 9:12, fscale, apply = FALSE)     # Changing the data by reference
ftransform(wlddev) <- fscale(gv(wlddev, 9:12))         # Same thing (using replacement method)

library(magrittr) # Same thing, using magrittr
wlddev \%<>\% ftransformv(9:12, fscale, apply = FALSE)
wlddev \%>\% ftransform(gv(., 9:12) |>              # With compound pipes: Scaling and lagging
                        fscale() |> flag(0:2, iso3c, year)) |> head()

# Fast reordering
head(roworder(wlddev, -country, year))
head(colorder(wlddev, country, year))

# Fast renaming
head(frename(wlddev, country = Ctry, year = Yr))
setrename(wlddev, country = Ctry, year = Yr)     # By reference
head(frename(wlddev, tolower, cols = 9:12))

# Fast grouping
fgroup_by(wlddev, Ctry, decade) |> fgroup_vars() |> head()
rm(wlddev)                                       # .. but only works with collapse functions

## Now lets start putting things together
wlddev |> fsubset(year > 1990, region, income, PCGDP:ODA) |>
  fgroup_by(region, income) |> fmean()         # Fast aggregation using the mean

# Same thing using dplyr manipulation verbs
library(dplyr)
wlddev |> filter(year > 1990) |> select(region, income, PCGDP:ODA) |>
  group_by(region,income) |> fmean()       # This is already a lot faster than summarize_all(mean)

wlddev |> fsubset(year > 1990, region, income, PCGDP:POP) |>
  fgroup_by(region, income) |> fmean(POP)     # Weighted group means

wlddev |> fsubset(year > 1990, region, income, PCGDP:POP) |>
  fgroup_by(region, income) |> fsd(POP)       # Weighted group standard deviations

wlddev |> na_omit(cols = "POP") |> fgroup_by(region, income) |>
  fselect(PCGDP:POP) |> fnth(0.75, POP)       # Weighted group third quartile

wlddev |> fgroup_by(country) |> fselect(PCGDP:ODA) |>
  fwithin() |> head()                         # Within transformation
wlddev |> fgroup_by(country) |> fselect(PCGDP:ODA) |>
  fmedian(TRA = "-") |> head()                # Grouped centering using the median
# Replacing data points by the weighted first quartile:
wlddev |> na_omit(cols = "POP") |> fgroup_by(country) |>
  fselect(country, year, PCGDP:POP) \%>\%
  ftransform(fselect(., -country, -year) |>
             fnth(0.25, POP, "fill")) |> head()

wlddev |> fgroup_by(country) |> fselect(PCGDP:ODA) |> fscale() |> head() # Standardizing
wlddev |> fgroup_by(country) |> fselect(PCGDP:POP) |>
   fscale(POP) |> head()  # Weighted..

wlddev |> fselect(country, year, PCGDP:ODA) |>  # Adding 1 lead and 2 lags of each variable
  fgroup_by(country) |> flag(-1:2, year) |> head()
wlddev |> fselect(country, year, PCGDP:ODA) |>  # Adding 1 lead and 10-year growth rates
  fgroup_by(country) |> fgrowth(c(0:1,10), 1, year) |> head()

# etc...

# Aggregation with multiple functions
wlddev |> fsubset(year > 1990, region, income, PCGDP:ODA) |>
  fgroup_by(region, income) \%>\% {
    add_vars(fgroup_vars(., "unique"),
             fmedian(., keep.group_vars = FALSE) |> add_stub("median_"),
             fmean(., keep.group_vars = FALSE) |> add_stub("mean_"),
             fsd(., keep.group_vars = FALSE) |> add_stub("sd_"))
  } |> head()

# Transformation with multiple functions
wlddev |> fselect(country, year, PCGDP:ODA) |>
  fgroup_by(country) \%>\% {
    add_vars(fdiff(., c(1,10), 1, year) |> flag(0:2, year),  # Sequence of lagged differences
             ftransform(., fselect(., PCGDP:ODA) |> fwithin() |> add_stub("W.")) |>
               flag(0:2, year, keep.ids = FALSE))            # Sequence of lagged demeaned vars
  } |> head()

# With ftransform, can also easily do one or more grouped mutations on the fly..
settransform(wlddev, median_ODA = fmedian(ODA, list(region, income), TRA = "fill"))

settransform(wlddev, sd_ODA = fsd(ODA, list(region, income), TRA = "fill"),
                     mean_GDP = fmean(PCGDP, country, TRA = "fill"))

wlddev \%<>\% ftransform(fmedian(list(median_ODA = ODA, median_GDP = PCGDP),
                               list(region, income), TRA = "fill"))

# On a groped data frame it is also possible to grouped transform certain columns
# but perform aggregate operatins on others:
wlddev |> fgroup_by(region, income) \%>\%
    ftransform(gmedian_GDP = fmedian(PCGDP, GRP(.), TRA = "replace"),
               omedian_GDP = fmedian(PCGDP, TRA = "replace"),  # "replace" preserves NA's
               omedian_GDP_fill = fmedian(PCGDP)) |> tail()

rm(wlddev)

## For multi-type data aggregation, the function collap() offers ease and flexibility
# Aggregate this data by country and decade: Numeric columns with mean, categorical with mode
head(collap(wlddev, ~ country + decade, fmean, fmode))

# taking weighted mean and weighted mode:
head(collap(wlddev, ~ country + decade, fmean, fmode, w = ~ POP, wFUN = fsum))

# Multi-function aggregation of certain columns
head(collap(wlddev, ~ country + decade,
            list(fmean, fmedian, fsd),
            list(ffirst, flast), cols = c(3,9:12)))

# Customized Aggregation: Assign columns to functions
head(collap(wlddev, ~ country + decade,
            custom = list(fmean = 9:10, fsd = 9:12, flast = 3, ffirst = 6:8)))

# For grouped data frames use collapg
wlddev |> fsubset(year > 1990, country, region, income, PCGDP:ODA) |>
  fgroup_by(country) |> collapg(fmean, ffirst) |>
  ftransform(AMGDP = PCGDP > fmedian(PCGDP, list(region, income), TRA = "fill"),
             AMODA = ODA > fmedian(ODA, income, TRA = "replace_fill")) |> head()

## Additional flexibility for data transformation tasks is offerend by tidy transformation operators
# Within-transformation (centering on overall mean)
head(W(wlddev, ~ country, cols = 9:12, mean = "overall.mean"))
# Partialling out country and year fixed effects
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year)))
# Same, adding ODA as continuous regressor
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year) + ODA))
# Standardizing (scaling and centering) by country
head(STD(wlddev, ~ country, cols = 9:12))
# Computing 1 lead and 3 lags of the 4 series
head(L(wlddev, -1:3, ~ country, ~year, cols = 9:12))
# Computing the 1- and 10-year first differences
head(D(wlddev, c(1,10), 1, ~ country, ~year, cols = 9:12))
head(D(wlddev, c(1,10), 1:2, ~ country, ~year, cols = 9:12))     # ..first and second differences
# Computing the 1- and 10-year growth rates
head(G(wlddev, c(1,10), 1, ~ country, ~year, cols = 9:12))
# Adding growth rate variables to dataset
add_vars(wlddev) <- G(wlddev, c(1, 10), 1, ~ country, ~year, cols = 9:12, keep.ids = FALSE)
get_vars(wlddev, "G1.", regex = TRUE) <- NULL # Deleting again

# These operators can conveniently be used in regression formulas:
# Using a Mundlak (1978) procedure to estimate the effect of OECD on LIFEEX, controlling for PCGDP
lm(LIFEEX ~ log(PCGDP) + OECD + B(log(PCGDP), country),
   wlddev |> fselect(country, OECD, PCGDP, LIFEEX) |> na_omit())

# Adding 10-year lagged life-expectancy to allow for some convergence effects (dynamic panel model)
lm(LIFEEX ~ L(LIFEEX, 10, country) + log(PCGDP) + OECD + B(log(PCGDP), country),
   wlddev |> fselect(country, OECD, PCGDP, LIFEEX) |> na_omit())

# Tranformation functions and operators also support indexed data classes:
wldi <- findex_by(wlddev, country, year)
head(W(wldi$PCGDP))                      # Country-demeaning
head(W(wldi, cols = 9:12))
head(W(wldi$PCGDP, effect = 2))          # Time-demeaning
head(W(wldi, effect = 2, cols = 9:12))
head(HDW(wldi$PCGDP))                    # Country- and time-demeaning
head(HDW(wldi, cols = 9:12))
head(STD(wldi$PCGDP))                    # Standardizing by country
head(STD(wldi, cols = 9:12))
head(L(wldi$PCGDP, -1:3))                # Panel-lags
head(L(wldi, -1:3, 9:12))
head(G(wldi$PCGDP))                      # Panel-Growth rates
head(G(wldi, 1, 1, 9:12))

lm(Dlog(PCGDP) ~ L(Dlog(LIFEEX), 0:3), wldi)   # Panel data regression
rm(wldi)
}
# Remove all objects used in this example section
rm(v, d, w, f, f1, f2, g, mtcarsM, sds, series, wlddev)
}


\keyword{package}
\keyword{manip}
