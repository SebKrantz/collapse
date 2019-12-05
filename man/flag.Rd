\name{flag-L-F}
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
Fast Lags and Leads for Time-Series and Panel Data
}
\description{
\code{flag} is an S3 generic to compute (sequences of) lags and leads. \code{L} and \code{F} are wrappers around \code{flag} representing the lag- and lead-operators, such that \code{L(x,-1) = F(x,1) = F(x)} and \code{L(x,-3:3) = F(x,3:-3)}. \code{L} & \code{F} provide more flexibility than \code{flag} when applied to data frames (i.e.  column subsetting, formula input and id-variable-preservation capabilities...), but are otherwise identical.

(\code{flag} is more of a programmers function in style of the \link[=fast-statistical-functions]{Fast Statistical Functions} while \code{L} & \code{F} are more practical to use in regression formulas or for computations on data frames.)
}
\usage{
flag(x, n = 1, \dots)
   L(x, n = 1, \dots)
   F(x, n = 1, \dots)
%## Default S3 methods:
\method{flag}{default}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{L}{default}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{F}{default}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
%## S3 methods for class 'matrix'
\method{flag}{matrix}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{L}{matrix}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{F}{matrix}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
%## S3 methods for class 'data.frame'
\method{flag}{data.frame}(x, n = 1, g = NULL, t = NULL, fill = NA, stubs = TRUE, \dots)
\method{L}{data.frame}(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
\method{F}{data.frame}(x, n = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)

# Methods for compatibility with plm:
%## S3 methods for class 'pseries'
\method{flag}{pseries}(x, n = 1, fill = NA, stubs = TRUE, \dots)
\method{L}{pseries}(x, n = 1, fill = NA, stubs = TRUE, \dots)
\method{F}{pseries}(x, n = 1, fill = NA, stubs = TRUE, \dots)
%## S3 methods for class 'pdata.frame'
\method{flag}{pdata.frame}(x, n = 1, fill = NA, stubs = TRUE, \dots)
\method{L}{pdata.frame}(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE,
  keep.ids = TRUE, \dots)
\method{F}{pdata.frame}(x, n = 1, cols = is.numeric, fill = NA, stubs = TRUE,
  keep.ids = TRUE, \dots)

# Methods for compatibility with dplyr:
%## S3 methods for class 'grouped_df'
\method{flag}{grouped_df}(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
\method{L}{grouped_df}(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
\method{F}{grouped_df}(x, n = 1, t = NULL, fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, matrix, data.frame, panel-series (\code{plm::pseries}), panel-data.frame (\code{plm::pdata.frame}) or grouped tibble (\code{dplyr::grouped_df}). \emph{Note}: Data must not be numeric.}
  \item{n}{an integer vector indicating the lags/leads to compute.}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{t}{same input as \code{g}, to indicate the time-variable. For safe computation of lags/leads on unordered time-series and panels. \emph{Note}: Data frame method also allows one-sided formula i.e. \code{~time}, and \code{grouped_df} method also allows lazy-evaluation i.e. \code{time} (no quotes).}
  \item{cols}{\emph{data.frame method}: Select columns to lag/lead using a function, column names or indices. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{stubs}{logical. \code{TRUE} will rename all lagged / leaded columns by adding a stub or prefix "L\code{n}." / "F\code{n}.".}
  \item{keep.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: Logical. Drop all panel-identifiers from the output (which includes all variables passed to \code{by} or \code{t}). \emph{Note}: For panel-data.frame's and grouped tibbles identifiers are dropped, but the 'index' / 'groups' attributes are kept.}
    \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
If a single integer is passed to \code{n}, and \code{g/by} and \code{t} are left empty, \code{flag/L/F} just returs \code{x} with all columns lagged / leaded by \code{n}. If \code{length(n)>1}, and \code{x} is an atomic vector, \code{flag/L/F} returns a matrix with lags / leads computed in the same order as passed to \code{n}. If instead \code{x} is a matrix / data.frame, a matrix / data.frame with \code{ncol(x)*length(n)} columns is returned where columns are sorted first by variable and then by lag (so all lags computed on a variable are grouped together). \code{x} can be of any standard data type.

With groups/panel-identifiers supplied to \code{g/by}, \code{flag/L/F} efficiently computes a panel-lag by shifting the entire vector(s) but inserting \code{fill} elements in the right places. If \code{t} is left empty, the data needs to be ordered such that all values belonging to a group are consecutive and in the right order. It is not necessary that the groups themselves occur in the right order. If a time-variable is supplied to \code{t} (or a list of time-variables uniquely identifying the time-dimension), the panel is fully identified and lags / leads can be securely computed even if the data is completely unordered (in that case data is shifted around and \code{fill} values are inserted in such a way that if the data were sorted afterwards the result would be identical to computing lags / leads on sorted data). Internally this works by using the grouping- and time-variable(s) to create an ordering and then accessing the panel-vector(s) through this ordering. If the data is just a bit unordered, such computations are nearly as fast as computations on ordered data (without \code{t}), however, if the data is very unordered, it can take significantly longer. Since most panel-data come perfectly or pretty ordered, I recommend always supplying \code{t} to be on the safe-side.

It is also possible to compute lags / leads on unordered time-series (thus utilizing \code{t} but leaving \code{g/by} empty), although this is probably more rare to encounter than unordered panels.

The methods applying to \code{plm} objects (panel-series and panel-data.frames) automatically utilize the panel-identifiers attached to these objects and thus securely compute fully identified panel-lags. If these objects have > 2 panel-identifiers attached to them, the last identifier is assumed to be the time-variable, and the others are taken as grouping-variables and interacted. I note that \code{flag/L/F} is significantly faster than \code{plm::lag/plm::lead} since the latter is written in R and based on a Split-Apply-Combine logic.
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
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=fdiff]{fdiff/D}}, \code{\link[=fgrowth]{fgrowth/G}}, \link[=time-series-panel-series]{Time-Series and Panel-Series}
}
\examples{
## Simple Time-Series: Airpassengers
L(AirPassengers)                      # 1 lag
F(AirPassengers)                      # 1 lead

all.identical(L(AirPassengers),       # 3 identical ways of computing 1 lag
              flag(AirPassengers),
              F(AirPassengers,-1))

L(AirPassengers,-1:3)                 # 1 lead and 3 lags - output as matrix

## Time-Series Matrix of 4 EU Stock Market Indicators, 1991-1998
tsp(EuStockMarkets)                                     # Data is recorded on 260 days per year
freq <- frequency(EuStockMarkets)
plot(stl(EuStockMarkets[,"DAX"], freq))                 # There is some obvious seasonality
L(EuStockMarkets,-1:3*freq)                             # 1 annual lead and 3 annual lags
summary(lm(DAX ~., data = L(EuStockMarkets,-1:3*freq))) # DAX regressed on it's own annual lead,
                                                        # lags and the lead/lags of the other series

## World Development Panel Data
head(flag(wlddev, 1, wlddev$iso3c, wlddev$year))        # This lags all variables,
head(L(wlddev, 1, ~iso3c, ~year))                       # This lags all numeric variables
head(L(wlddev, 1, ~iso3c))                              # Without t: Works because data is ordered
head(L(wlddev, 1, PCGDP + LIFEEX ~ iso3c, ~year))       # This lags GDP per Capita & Life Expectancy
head(L(wlddev, 0:2, ~ iso3c, ~year, cols = 9:10))       # Same, also retaining original series
head(L(wlddev, 1:2, PCGDP + LIFEEX ~ iso3c, ~year,      # Two lags, dropping id columns
       keep.ids = FALSE))

# Different ways of regressing GDP on its's lags and life-Expectancy and it's lags
summary(lm(PCGDP ~ ., L(wlddev, 0:2, ~iso3c, ~year, 9:10, keep.ids = FALSE)))     # 1 - Precomputing
summary(lm(PCGDP ~ L(PCGDP,1:2,iso3c,year) + L(LIFEEX,0:2,iso3c,year), wlddev))   # 2 - Ad-hoc
summary(lm(PCGDP ~ L(PCGDP,1:2,iso3c) + L(LIFEEX,0:2,iso3c), wlddev))             # 3 - same no year
g = qF(wlddev$iso3c); t = qF(wlddev$year)                                         # 4- Precomputing
summary(lm(PCGDP ~ L(PCGDP,1:2,g,t) + L(LIFEEX,0:2,g,t), wlddev))                 # panel-id's

## Using plm:
pwlddev <- plm::pdata.frame(wlddev, index = c("iso3c","year"))
head(L(pwlddev, 0:2, 9:10))                                     # Again 2 lags of GDP and LIFEEX
PCGDP <- pwlddev$PCGDP                                          # A panel-Series of GDP per Capita
L(PCGDP)                                                        # Lagging the panel series
summary(lm(PCGDP ~ ., L(pwlddev, 0:2, 9:10, keep.ids = FALSE))) # Running the lm again: WORKS!
# THIS DOES NOT WORK: Unfortunately lm drops the attributes of the columns,
# so L.default is used here and ordinary lags are computed. (with and attach don't retain attr.)
summary(lm(PCGDP ~ L(PCGDP,1:2) + L(LIFEEX,0:2), pwlddev))
LIFEEX <- pwlddev$LIFEEX                                        # To make it work, create pseries
summary(lm(PCGDP ~ L(PCGDP,1:2) + L(LIFEEX,0:2)))               # THIS WORKS !!

## Using dplyr:
library(dplyr)
wlddev \%>\% group_by(iso3c) \%>\% select(PCGDP,LIFEEX) \%>\% L(0:2)
wlddev \%>\% group_by(iso3c) \%>\% select(year,PCGDP,LIFEEX) \%>\% L(0:2,year) # Also using t (safer)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{ts}% __ONLY ONE__ keyword per line
