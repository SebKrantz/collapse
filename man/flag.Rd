\name{flag}
\alias{flag}
\alias{flag.default}
\alias{flag.matrix}
\alias{flag.data.frame}
\alias{flag.pseries}
\alias{flag.pdata.frame}
\alias{flag.grouped_df}
\alias{L}
\alias{L.default}
\alias{L.matrix}
\alias{L.data.frame}
\alias{L.pseries}
\alias{L.pdata.frame}
\alias{L.grouped_df}
\alias{F}
\alias{F.default}
\alias{F.matrix}
\alias{F.data.frame}
\alias{F.pseries}
\alias{F.pdata.frame}
\alias{F.grouped_df}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Lags and Leads for Time Series and Panel Data
}
\description{
\code{flag} is an S3 generic to compute (sequences of) lags and leads. \code{L} and \code{F} are wrappers around \code{flag} representing the lag- and lead-operators, such that \code{L(x,-1) = F(x,1) = F(x)} and \code{L(x,-3:3) = F(x,3:-3)}. \code{L} and \code{F} provide more flexibility than \code{flag} when applied to data frames (i.e.  column subsetting, formula input and id-variable-preservation capabilities\dots), but are otherwise identical.

\emph{Note:} Since v1.9.0, \code{F} is no longer exported, but can be accessed using \code{collapse:::F}, or through setting \code{options(collapse_export_F = TRUE)} before loading the package. The syntax is the same as \code{L}.

% (\code{flag} is more of a programmers function in style of the \link[=fast-statistical-functions]{Fast Statistical Functions} while \code{L} and \code{F} are more practical to use in regression formulas or for computations on data frames.)
}
\usage{
flag(x, n = 1, \dots)
   L(x, n = 1, \dots)

\method{flag}{default}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{L}{default}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)

\method{flag}{matrix}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, \dots)
\method{L}{matrix}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)

\method{flag}{data.frame}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = length(n) > 1L, \dots)
\method{L}{data.frame}(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{flag}{pseries}(x, n = 1, fill = NA, stubs = length(n) > 1L, shift = "time", \dots)
\method{L}{pseries}(x, n = 1, fill = NA, stubs = TRUE, shift = "time", \dots)

\method{flag}{pdata.frame}(x, n = 1, fill = NA, stubs = length(n) > 1L, shift = "time", \dots)
\method{L}{pdata.frame}(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE,
  shift = "time", keep.ids = TRUE, \dots)

# Methods for grouped data frame / compatibility with dplyr:

\method{flag}{grouped_df}(x, n = 1, t = NULL, fill = NA, stubs = length(n) > 1L, keep.ids = TRUE, \dots)
\method{L}{grouped_df}(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector / time series, (time series) matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df'). Data must not be numeric.}
  \item{n}{integer. A vector indicating the lags / leads to compute (passing negative integers to \code{flag} or \code{L} computes leads, passing negative integers to \code{F} computes lags).}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{t}{a time vector or list of vectors. Data frame methods also allows one-sided formula i.e. \code{~time}. grouped_df method supports lazy-evaluation i.e. \code{time} (no quotes). Either support wrapping a transformation function e.g. \code{~timeid(time)}, \code{qG(time)} etc.. See also Details on how \code{t} is processed.}
  \item{cols}{\emph{data.frame method}: Select columns to lag using a function, column names, indices or a logical vector. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{stubs}{logical. \code{TRUE} will rename all lagged / leaded columns by adding a stub or prefix "L\code{n}." / "F\code{n}.".}
  \item{shift}{\emph{pseries / pdata.frame methods}: character. \code{"time"} performs a fully identified time-lag (if the index contains a time variable), whereas \code{"row"} performs a simple (group) lag, where observations are shifted based on the present order of rows (in each group). The latter is significantly faster, but requires time series / panels to be regularly spaced and sorted by time within each group.}
  \item{keep.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: Logical. Drop all identifiers from the output (which includes all variables passed to \code{by} or \code{t}). \emph{Note}: For 'grouped_df' / 'pdata.frame' identifiers are dropped, but the \code{"groups"} / \code{"index"} attributes are kept.}
    \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
If a single integer is passed to \code{n}, and \code{g/by} and \code{t} are left empty, \code{flag/L/F} just returns \code{x} with all columns lagged / leaded by \code{n}. If \code{length(n)>1}, and \code{x} is an atomic vector (time series), \code{flag/L/F} returns a (time series) matrix with lags / leads computed in the same order as passed to \code{n}. If instead \code{x} is a matrix / data frame, a matrix / data frame with \code{ncol(x)*length(n)} columns is returned where columns are sorted first by variable and then by lag (so all lags computed on a variable are grouped together). \code{x} can be of any standard data type.

With groups/panel-identifiers supplied to \code{g/by}, \code{flag/L/F} efficiently computes a panel-lag/lead by shifting the entire vector(s) but inserting \code{fill} elements in the right places. If \code{t} is left empty, the data needs to be ordered such that all values belonging to a group are consecutive and in the right order. It is not necessary that the groups themselves are alphabetically ordered. If a time-variable is supplied to \code{t} (or a list of time-variables uniquely identifying the time-dimension), the series / panel is fully identified and lags / leads can be securely computed even if the data is unordered / irregular.

%It is also possible to lag unordered or irregular time series utilizing only the \code{t} argument to identify the temporal dimension of the data.

% Since v1.5.0 \code{flag/L/F} provide full built-in support for irregular time series and unbalanced panels. The suggested workaround using the \code{\link{seqid}} function is therefore no longer necessary.
%\code{flag/L/F} supports balanced panels and unbalanced panels where various individuals are observed for different time-sequences (both start, end and duration of observation can differ for each individual). \code{flag/L/F} does not natively support irregularly spaced time series and panels, that is situations where there are either gaps in time and/or repeated observations in the same time-period for some individual (see also computational details below). For such cases the function \code{\link{seqid}} can be used to generate an appropriate panel-identifier (i.e. splitting individuals with an irregular time-sequence into multiple individuals with regular time-sequences before applying \code{flag/L/F}).

%(in that case data is shifted around and \code{fill} values are inserted in such a way that if the data were sorted afterwards the result would be identical to computing lags / leads on sorted data). Internally this works by using the grouping- and time-variable(s) to create an ordering and then accessing the panel-vector(s) through this ordering. If the data is just a bit unordered, such computations are nearly as fast as computations on ordered data (without \code{t}), however, if the data is very unordered, it can take significantly longer. Since most panel data come perfectly or pretty ordered, I recommend always supplying \code{t} to be on the safe-side.

% It is also possible to compute lags / leads on unordered time series (thus utilizing \code{t} but leaving \code{g/by} empty), although this is probably more rare to encounter than unordered panels. Irregularly spaced time series can also be lagged using a panel- identifier generated with \code{\link{seqid}}.

\bold{Note} that the \code{t} argument is processed as follows: If \code{is.factor(t) || (is.numeric(t) && !is.object(t))} (i.e. \code{t} is a factor or plain numeric vector), it is assumed to represent unit timesteps (e.g. a 'year' variable in a typical dataset), and thus coerced to integer using \code{as.integer(t)} and directly passed to C++ without further checks or transformations at the R-level. Otherwise, if \code{is.object(t) && is.numeric(unclass(t))} (i.e. \code{t} is a numeric time object, most likely 'Date' or 'POSIXct'), this object is passed through \code{\link{timeid}} before going to C++. Else (e.g. \code{t} is character), it is passed through \code{\link{qG}} which performs ordered grouping. If \code{t} is a list of multiple variables, it is passed through \code{\link{finteraction}}. You can customize this behavior by calling any of these functions (including \code{unclass/as.integer}) on your time variable beforehand.

At the C++ level, if both \code{g/by} and \code{t} are supplied, \code{flag} works as follows: Use two initial passes to create an ordering through which the data are accessed. First-pass: Calculate minimum and maximum time-value for each individual. Second-pass: Generate an internal ordering vector (\code{o}) by placing the current element index into the vector slot obtained by adding the cumulative group size and the current time-value subtracted its individual-minimum together. This method of computation is faster than any sort-based method and delivers optimal performance if the panel-id supplied to \code{g/by} is already a factor variable, and if \code{t} is an integer/factor variable. For irregular time/panel series, \code{length(o) > length(x)}, and \code{o} represents the unobserved 'complete series'. If \code{length(o) > 1e7 && length(o) > 3*length(x)}, a warning is issued to make you aware of potential performance implications of the oversized ordering vector.

%If \code{t} is not factor or integer but instead \code{is.double(t) && !is.object(t)}, it is assumed to be integer represented by double and converted using \code{as.integer(t)}. For other objects such as dates, \code{t} is grouped using \code{\link{qG}} or \code{\link{GRP}} (for multiple time identifiers). Similarly, if \code{g/by} is not factor or 'GRP' object, \code{\link{qG}} or \code{\link{GRP}} will be called to group the respective identifier. Since grouping is more expensive than computing lags, prepare the data for optimal performance (or use \emph{plm} classes). See also the Note.

%A caveat of not using sort-based methods is that gaps or repeated values in time are only recognized towards the end of the second pass where they cannot be rectified anymore, and thus \code{flag/L/F} does not natively support irregular panels but throws an error.


The 'indexed_series' ('pseries') and 'indexed_frame' ('pdata.frame')  methods automatically utilize the identifiers attached to these objects, which are already factors, thus lagging is quite efficient. However, the internal ordering vector still needs to be computed, thus if data are known to be ordered and regularly spaced, using \code{shift = "row"} to toggle a simple group-lag (same as utilizing \code{g} but not \code{t} in other methods) can yield a significant performance gain.

%and thus securely and efficiently compute fully identified panel-lags. If these objects have > 2 panel-identifiers attached to them, the last identifier is assumed to be the time-variable, and the others are taken as grouping-variables and interacted. Note that \code{flag/L/F} is significantly faster than \code{plm::lag/plm::lead} since the latter is written in R and based on a Split-Apply-Combine logic.
}
\value{
\code{x} lagged / leaded \code{n}-times, grouped by \code{g/by}, ordered by \code{t}. See Details and Examples.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%Since v1.7.0, if \code{is.double(t) && !is.object(t)}, it is coerced to integer using \code{as.integer(t)}. This is to avoid the inefficiency of ordered grouping, and owes to the fact that in most data imported into R, the time (year) variables are coded as double although they should be integer.
% }

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{fdiff}}, \code{\link{fgrowth}}, \link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple Time Series: AirPassengers
L(AirPassengers)                      # 1 lag
flag(AirPassengers)                   # Same
L(AirPassengers, -1)                  # 1 lead

head(L(AirPassengers, -1:3))          # 1 lead and 3 lags - output as matrix

## Time Series Matrix of 4 EU Stock Market Indicators, 1991-1998
tsp(EuStockMarkets)                                     # Data is recorded on 260 days per year
freq <- frequency(EuStockMarkets)
plot(stl(EuStockMarkets[,"DAX"], freq))                 # There is some obvious seasonality
head(L(EuStockMarkets, -1:3 * freq))                    # 1 annual lead and 3 annual lags
summary(lm(DAX ~., data = L(EuStockMarkets,-1:3*freq))) # DAX regressed on its own annual lead,
                                                        # lags and the lead/lags of the other series
## World Development Panel Data
head(flag(wlddev, 1, wlddev$iso3c, wlddev$year))        # This lags all variables,
head(L(wlddev, 1, ~iso3c, ~year))                       # This lags all numeric variables
head(L(wlddev, 1, ~iso3c))                              # Without t: Works because data is ordered
head(L(wlddev, 1, PCGDP + LIFEEX ~ iso3c, ~year))       # This lags GDP per Capita & Life Expectancy
head(L(wlddev, 0:2, ~ iso3c, ~year, cols = 9:10))       # Same, also retaining original series
head(L(wlddev, 1:2, PCGDP + LIFEEX ~ iso3c, ~year,      # Two lags, dropping id columns
       keep.ids = FALSE))

# Regressing GDP on its's lags and life-Expectancy and its lags
summary(lm(PCGDP ~ ., L(wlddev, 0:2, ~iso3c, ~year, 9:10, keep.ids = FALSE)))

## Indexing the data: facilitates time-based computations
wldi <- findex_by(wlddev, iso3c, year)
head(L(wldi, 0:2, cols = 9:10))                              # Again 2 lags of GDP and LIFEEX
head(L(wldi$PCGDP))                                          # Lagging an indexed series
summary(lm(PCGDP ~ L(PCGDP,1:2) + L(LIFEEX,0:2), wldi))      # Running the lm again
summary(lm(PCGDP ~ ., L(wldi, 0:2, 9:10, keep.ids = FALSE))) # Same thing
\donttest{ % No code relying on suggested package
## Using grouped data:
library(magrittr)
wlddev \%>\% fgroup_by(iso3c) \%>\% fselect(PCGDP,LIFEEX) \%>\% flag(0:2)
wlddev \%>\% fgroup_by(iso3c) \%>\% fselect(year,PCGDP,LIFEEX) \%>\% flag(0:2,year) # Also using t (safer)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{ts} % __ONLY ONE__ keyword per line
