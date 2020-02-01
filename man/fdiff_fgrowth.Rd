\name{fdiff-fgrowth-D-G}
\alias{fdiff}
\alias{fdiff.default}
\alias{fdiff.matrix}
\alias{fdiff.data.frame}
\alias{fdiff.pseries}
\alias{fdiff.pdata.frame}
\alias{fdiff.grouped_df}
\alias{fgrowth}
\alias{fgrowth.default}
\alias{fgrowth.matrix}
\alias{fgrowth.data.frame}
\alias{fgrowth.pseries}
\alias{fgrowth.pdata.frame}
\alias{fgrowth.grouped_df}
\alias{D}
\alias{D.default}
\alias{D.matrix}
\alias{D.data.frame}
\alias{D.pseries}
\alias{D.pdata.frame}
\alias{D.grouped_df}
\alias{G}
\alias{G.default}
\alias{G.matrix}
\alias{G.data.frame}
\alias{G.pseries}
\alias{G.pdata.frame}
\alias{G.grouped_df}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ % Lagged and Iterated
Fast Differences and Growth Rates for Time-Series and Panel Data
}
\description{
\code{fdiff} and \code{fgrowth} are S3 generics to compute (sequences of) suitably lagged / leaded and iterated differences and growth rates / log-differences, respectively. \code{D} and \code{G} are wrappers around \code{fdiff} and \code{fgrowth} representing the 'difference-operator' and the 'growth-operator'. \code{D} / \code{G} provide more flexibility than \code{fdiff} / \code{fgrowth} when applied to data frames, but are otherwise identical.

(\code{fdiff} and \code{fgrowth} are programmers functions in style of the \link[=A1-fast-statistical-functions]{Fast Statistical Functions} while \code{D} & \code{G} are more practical to use in regression formulas or for computations on data frames.)
}
\usage{
  fdiff(x, n = 1, diff = 1, \dots)
fgrowth(x, n = 1, diff = 1, \dots)
      D(x, n = 1, diff = 1, \dots)
      G(x, n = 1, diff = 1, \dots)
% ## Default S3 methods:
\method{fdiff}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
      stubs = TRUE, \dots)
\method{fgrowth}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, stubs = TRUE, \dots)
\method{D}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
  stubs = TRUE, \dots)
\method{G}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
  logdiff = FALSE, stubs = TRUE, \dots)
%## S3 methods for class 'matrix'
\method{fdiff}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
      stubs = TRUE, \dots)
\method{fgrowth}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, stubs = TRUE, \dots)
\method{D}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
  stubs = TRUE, \dots)
\method{G}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
  logdiff = FALSE, stubs = TRUE, \dots)
%## S3 methods for class 'data.frame'
\method{fdiff}{data.frame}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
      stubs = TRUE, \dots)
\method{fgrowth}{data.frame}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, stubs = TRUE, \dots)
\method{D}{data.frame}(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, stubs = TRUE, keep.ids = TRUE, \dots)
\method{G}{data.frame}(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, \dots)

# Methods for compatibility with plm:
%## S3 methods for class 'pseries'
\method{fdiff}{pseries}(x, n = 1, diff = 1, fill = NA, stubs = TRUE, \dots)
\method{fgrowth}{pseries}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, \dots)
\method{D}{pseries}(x, n = 1, diff = 1, fill = NA, stubs = TRUE, \dots)
\method{G}{pseries}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, \dots)
%## S3 methods for class 'pdata.frame'
\method{fdiff}{pdata.frame}(x, n = 1, diff = 1, fill = NA, stubs = TRUE, \dots)
\method{fgrowth}{pdata.frame}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, stubs = TRUE, \dots)
\method{D}{pdata.frame}(x, n = 1, diff = 1, cols = is.numeric, fill = NA, stubs = TRUE,
  keep.ids = TRUE, \dots)
\method{G}{pdata.frame}(x, n = 1, diff = 1, cols = is.numeric, fill = NA,
  logdiff = FALSE, stubs = TRUE, keep.ids = TRUE, \dots)

# Methods for compatibility with dplyr:
%## S3 methods for class 'grouped_df'
\method{fdiff}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, stubs = TRUE,
      keep.ids = TRUE, \dots)
\method{fgrowth}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE,
        stubs = TRUE, keep.ids = TRUE, \dots)
\method{D}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, stubs = TRUE,
  keep.ids = TRUE, \dots)
\method{G}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE,
  stubs = TRUE, keep.ids = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (\code{plm::pseries}), panel-data.frame (\code{plm::pdata.frame}) or grouped tibble (\code{dplyr::grouped_df}).}
  \item{n}{a integer vector indicating the number of lags or leads.}
  \item{diff}{a vector of integers > 1 indicating the order of differencing / taking growth rates or log-differences.}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{t}{same input as \code{g}, to indicate the time-variable. For safe computation of differences/growth rates on unordered time-series and panels. \emph{Notes}: data.frame method also allows name, index or one-sided formula i.e. \code{~time}. grouped_df method also allows lazy-evaluation i.e. \code{time} (no quotes).}
  \item{cols}{\emph{data.frame method}: Select columns to difference/compute growth rates using a function, column names or indices. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{logdiff}{logical. compute log-differences instead of exact growth rates. See Details.}
  \item{stubs}{logical. \code{TRUE} will rename all differenced columns by adding a prefix "L\code{n}D\code{diff}." / "F\code{n}D\code{diff}." and a prefix "L\code{n}G\code{diff}." / "F\code{n}G\code{diff}." for growth rates.}
  \item{keep.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: Logical. Drop all panel-identifiers from the output (which includes all variables passed to \code{by} or \code{t}). \emph{Note}: For panel-data.frame's and grouped tibbles identifiers are dropped, but the 'index' / 'groups' attributes are kept.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
By default, \code{fdiff/D|fgrowth/G} return \code{x} with all columns differenced | converted to growth rates. Differences are computed as \code{repeat(diff){x[i] - x[i-n]}}, growth rates as \code{repeat(diff){(x[i] - x[i-n])/x[i-n]*100}} and log-differences as \code{repeat(diff){(log(x[i]) - log(x[i-n]))*100}}. Setting \code{diff = 2} thus returns differences of differences | growth rates of growth rates etc... and setting \code{n = 2} returns simple differences computed by subtracting twice-lagged \code{x} from \code{x}. It is also possible to compute forward differences | growth rates by passing negative \code{n} values. \code{n} also supports sequences of integers (lags), and \code{diff} supports positive sequences of integers (differences):

If more than one value is passed to \code{n} and/or \code{diff}, the data is expanded-wide as follows: If \code{x} is an atomic vector or time-series, a (time-series) matrix is returned with columns ordered first by lag, then by difference. If \code{x} is a matrix or data.frame, each column is expanded in like manor such that the output has \code{ncol(x)*length(n)*length(diff)} columns ordered first by column name, then by lag, then by difference.

With groups/panel-identifiers supplied to \code{g/by}, \code{fdiff/D|fgrowth/G} efficiently compute panel-differences | growth rates by inserting \code{fill} elements in the right places. If \code{t} is left empty, the data needs to be ordered such that all values belonging to a group are consecutive and in the right order. It is not necessary that the groups themselves occur in the right order. If time-variable(s) are supplied to \code{t}, the panel is fully identified and differences | growth rates can be securely computed even if the data is completely unordered (in that case data is shifted around and \code{fill} values are inserted in such a way that if the data were sorted afterwards the result would be identical to computing on sorted data). Internally this works by using the grouping- and time-variables to create an ordering and then accessing the panel-vector(s) through this ordering. If the data is just a bit unordered, such computations are nearly as fast as computations on ordered data (without \code{t}), however, if the data is very unordered, it can take significantly longer. Since most panel-data come perfectly or pretty ordered, I recommend always supplying \code{t} to be on the safe-side.

It is also possible to compute differences | growth rates on unordered vectors / time-series (thus utilizing \code{t} but leaving \code{g/by} empty).

The methods applying to \code{plm} objects (panel-series and panel-data.frames) automatically utilize the panel-identifiers attached to these objects and thus securely compute fully identified panel-differences. If these objects have > 2 panel-identifiers attached to them, the last identifier is assumed to be the time-variable, and the others are taken as grouping-variables and interacted.
}
\value{
\code{fdiff/D} returns \code{x} differenced \code{diff} times using lags \code{n} of itself.
\code{fgrowth/G} returns \code{x} where the growth rate or log-difference was taken \code{diff} times using lags \code{n} of itself. Computations can be grouped by \code{g/by} and/or ordered by \code{t}. See Details and Examples.
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
\code{\link[=flag]{flag/L/F}}, \link[=A7-time-series-panel-series]{Time-Series and Panel-Series}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple Time-Series: Airpassengers
D(AirPassengers)                      # 1st difference, same as fdiff(AirPassengers)
D(AirPassengers,-1)                   # forward difference
G(AirPassengers)                      # growth rate, same as fgrowth(AirPassengers)
G(AirPassengers, logdiff = TRUE)      # log-difference
D(AirPassengers,1,2)                  # second difference
G(AirPassengers,1,2)                  # growth rate of growth rate
D(AirPassengers,12)                   # seasonal difference (data is monthly)
G(AirPassengers,12)                   # seasonal growth rate (data is monthly)

D(AirPassengers,-2:2,1:3)             # sequence of leaded/lagged and iterated differences

# let's do some visual analysis
plot(AirPassengers)                   # plot the series - seasonal pattern is evident
plot(stl(AirPassengers, "periodic"))  # Seasonal decomposition
plot(D(AirPassengers,c(1,12),1:2))    # plotting ordinary and seasonal first and second differences
plot(G(AirPassengers,c(1,12),1:2))    # same using growth rates
plot(stl(window(G(AirPassengers,12),  # Taking seasonal growth rate removes most seasonal variation
                1950), "periodic"))


## Time-Series Matrix of 4 EU Stock Market Indicators, recorded 260 days per year
plot(G(EuStockMarkets,c(0,260)))                    # Plot series and annual growth rates
summary(lm(L260G1.DAX ~., G(EuStockMarkets,260)))   # Annual growth rate of DAX regressed on the
                                                    # growth rates of the other indicators

## World Development Panel Data
head(fgrowth(num_vars(wlddev), 1, 1,                    # Computes growth rates of numeric variables
             wlddev$country, wlddev$year))              # fgrowth/fdiff require externally inputs...
head(G(wlddev, 1, 1, ~country, ~year))                  # Growth of numeric variables, id's attached
head(G(wlddev, 1, 1, ~country))                         # Without t: Works because data is ordered
head(G(wlddev, 1, 1, PCGDP + LIFEEX ~ country, ~year))  # Growth of GDP per Capita & Life Expectancy
head(G(wlddev, 0:1, 1, ~ country, ~year, cols = 9:10))  # Same, also retaining original series
head(G(wlddev, 0:1, 1, ~ country, ~year, 9:10,          # Dropping id columns
       keep.ids = FALSE))

# Dynamic Panel-Data Models:
summary(lm(G(PCGDP,1,1,iso3c,year) ~                    # GDP growth regressed on it's lagged level
             L(PCGDP,1,iso3c,year) +                    # and the growth rate of Life Expanctancy
             G(LIFEEX,1,1,iso3c,year), data = wlddev))

g = qF(wlddev$country)                                  # Omitting t and precomputing g allows for a
summary(lm(G(PCGDP,1,1,g) ~ L(PCGDP,1,g) +              # bit more parsimonious specification
                            G(LIFEEX,1,1,g), wlddev))

summary(lm(G1.PCGDP ~.,                                 # Now adding level and lagged level of
L(G(wlddev,0:1,1, ~ country, ~year,9:10),0:1,           # LIFEEX and lagged growth rates
  ~ country, ~year, keep.ids = FALSE)[-1]))

## Using plm can make things easier, but avoid attaching or 'with' calls:
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))
head(G(pwlddev, 0:1, 1, 9:10))                          # Again growth rates of LIFEEX and PCGDP
PCGDP <- pwlddev$PCGDP                                  # A panel-Series of GDP per Capita
D(PCGDP)                                                # Differencing the panel series.
summary(lm(G1.PCGDP ~.,                                 # Running the dynamic model again ->
           data = L(G(pwlddev,0:1,1,9:10),0:1,          # code becomes a bit simpler
                    keep.ids = FALSE)[-1]))

# One could be tempted to also do something like this, but THIS DOES NOT WORK!!!:
# lm drops the attributes (-> with(pwlddev, PCGDP) drops attr. so G.default and L.matrix are used)
summary(lm(G(PCGDP) ~ L(G(PCGDP,0:1)) + L(G(LIFEEX,0:1),0:1), pwlddev))

# To make it work, one needs to create pseries (note: attach(pwlddev) also won't work)
LIFEEX <- pwlddev$LIFEEX
summary(lm(G(PCGDP) ~ L(G(PCGDP,0:1)) + L(G(LIFEEX,0:1),0:1))) # THIS WORKS !!

## Using dplyr:
library(dplyr)
wlddev \%>\% group_by(country) \%>\%
             select(PCGDP,LIFEEX) \%>\% D(0:1,1:2)           # Adding a first and second difference
wlddev \%>\% group_by(country) \%>\%
             select(year,PCGDP,LIFEEX) \%>\% D(0:1,1:2,year) # Also using t (safer)
wlddev \%>\% group_by(country) \%>\%                         # Growth rates, dropping id's
             select(year,PCGDP,LIFEEX) \%>\% G(0:1,1:2,year, keep.ids = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{ts}
