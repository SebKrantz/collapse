\name{lags & leads}
\alias{lags&leads}
\alias{L}
\alias{F}
\alias{flag}
\alias{flead}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Lags and Leads for Time-Series and Panel Data
}
\description{
\code{flag} is an S3 generic to compute (sequences of) lags. \code{flead} computes leads and is a wrapper around \code{flag} i.e. \code{flag(x,-1) = flead(x,1)} and \code{flag(x,-3:3) = flead(x,3:-3)}. \code{L} is also a wrapper around \code{flag} and represents the 'lag-operator'. \code{L} provides more options than \code{flag} when applied to data.frame's, but is otherwise identical. To round things off, \code{F}, the 'lead-operator', is a wrapper around \code{L} such that \code{L(x,-1) = F(x,1)}, etc... 
}
\usage{
## Default S3 methods:
flag/L|flead/F(x, n = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)

## S3 methods for class 'matrix'
flag/L|flead/F(x, n = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)

## S3 methods for class 'data.frame'
flag|flead(x, n = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)
L|F(x, n = 1, by = NULL, t = NULL, cols = is.numeric, 
    fill = NA, give.names = TRUE, drop.ids = FALSE, ...)

# Methods for compatibility with plm:

## S3 methods for class 'pseries'
flag/L|flead/F(x, n = 1, fill = NA, give.names = TRUE, ...)

## S3 methods for class 'pdata.frame'
flag|flead(x, n = 1, fill = NA, give.names = TRUE, ...)
L|F(x, n = 1, cols = is.numeric, fill = NA, give.names = TRUE, drop.ids = FALSE, ...)

# Methods for compatibility with dplyr:

## S3 methods for class 'grouped_df'
flag/L|flead/F(x, n = 1, t = NULL, fill = NA, give.names = TRUE, drop.ids = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, matrix, data.frame, panel-series (\code{\link{plm::pseries}}), panel-data.frame (\code{\link{plm::pdata.frame}}) or grouped tibble (\code{\link{dplyr::grouped_df}}).}
  \item{n}{an integer vector indicating the lags/leads to compute.}
  \item{g}{a factor, \code{GRP()} object, atomic vector (> internally converted to factor) or a list of vectors / factors (> internally converted to a \code{GRP()} object) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: same as g, but also allows for column-names, column-indices and one- or two-sided formulas. See Examples.}
  \item{t}{same input as g, to indicate the time-variable. For safe computation of lags/leads on unordered time-series and panels. \emph{Notes}: data.frame method also allows name, index or one-sided formula i.e. \code{~time}. grouped_df method also allows lazy-evaluation i.e. \code{time} (no quotes).}
  \item{cols}{\emph{data.frame method}: Select columns to lag/lead using a function, column names or indices. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{give.names}{\code{TRUE} will rename all lagged / leaded columns by adding a prefix "L\code{n}." / "F\code{n}.".}
  \item{drop.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: drop all panel-identifiers from the output (which includes all variables passed to 'by' or 't'). For \code{pdata.frame}'s identifiers are dropped, but the 'index' attribute is always kept.}
}
\details{
If a single integer is passed to \code{n}, and \code{g/by} and \code{t} are left empty, \code{flag/L|flead/F} just return \code{x} with all columns lagged | leaded by \code{n}. If \code{length(n)>1}, and \code{x} is an atomic vector, \code{flag/L|flead/F} return a matrix with lags | leads computed in the same order as passed to \code{n}. If instead \code{x} is a matrix | data.frame, a matrix | data.frame with \code{ncol(x)*length(n)} columns is returned where columns are sorted first by variable and then by lag (so all lags computed on a variable are grouped together). \code{x} can be of any standard data type.

With groups/panel-identifiers supplied to \code{g/by}, \code{flag/L|flead/F} efficiently compute a panel-lag by shifting the entire vector(s) but inserting \code{fill} elements in the right places. If \code{t} is left empty, the data needs to be ordered such that all values belonging to a group are consecutive and in the right order. It is not necessary that the groups themselves occur in the right order. If time-variable(s) are supplied to \code{t}, the panel is fully identified and lags | leads can be securely computed even if the data is completely unordered (in that case data is shifted around and \code{fill} values are inserted in such a way that if the data were sorted afterwards the result would be identical to computing lags | leads on sorted data). Internally this works by using the grouping- and time-variables to create an ordering and then accessing the panel-vector(s) through this ordering. If the data is just a bit unordered, such computations are nearly as fast as computations on ordered data (without \code{t}), however, if the data is very unordered, it can take significantly longer. Since most panel-data come perfectly or pretty ordered, I recommend always supplying \code{t} to be on the safe-side. 

It is also possible to compute lags | leads on unordered time-series (thus utilizing \code{t} but leaving \code{g/by} empty), although this is probably more rare to encounter than unordered panels. 

The methods applying to \code{plm} objects (panel-series and panel-data.frames) automatically utilize the panel-identifiers attached to these objects and thus securely compute fully identified panel-lags. If these objects have > 2 panel-identifiers attached to them, the last identifier is assumed to be the time-variable, and the others are taken as grouping-variables and interacted. I note that \code{flag/L|flead/F} are significantly faster than \code{plm::lag|plm::lead} since the latter are written in R (vs C++) and based on a split-apply-combine logic which can become very slow on large panels with many groups. 
}
\value{
\code{x} lagged | leaded by \code{n}. See Details and Examples. 
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
## Simple Time-Series: Airpassengers
L(AirPassengers)                      # 1 lag
F(AirPassengers)                      # 1 lead

all.identical(L(AirPassengers),       # 4 identical ways of computing 1 lag
              flag(AirPassengers),    
              F(AirPassengers,-1), 
              flead(AirPassengers,-1)) 

L(AirPassengers,-1:3)                 # 1 lead and 3 lags - output as matrix

## Time-Series Matrix of 4 EU Stock Market Indicators
head(L(EuStockMarkets,-1:3))                       # 1 lead and 3 lags
summary(lm(DAX ~., data = L(EuStockMarkets,-1:3))) # DAX regressed on it's own lead, lags and the lead and lags of the other indicators

## World Development Panel Data
head(flag(wlddev, 1, wlddev$country, wlddev$year))   # This lags all variables. With flag all inputs need to be passed externally, so L is more convenient: 
head(L(wlddev, 1, ~country, ~year))                  # This lags all numeric variables
head(L(wlddev, 1, ~country))                         # Without t: Also works because data is ordered 
head(L(wlddev, 1, PCGDP + LIFEEX ~ country, ~year))  # This lags GDP per Capita and Life Expectancy
head(L(wlddev, 0:2, ~ country, ~year, cols = 9:10))  # Same, also retaining original series
head(L(wlddev, 0:2, 1, 4, 9:10))                     # Same, even more parsimonious
head(L(wlddev, 1:2, PCGDP + LIFEEX ~ country, ~year, drop.ids = TRUE)) # Two lags, dropping id columns

# Different ways of regressing GDP on its's lags and life-Expectancy and it's lags
summary(lm(PCGDP ~ ., L(wlddev, 0:2, PCGDP + LIFEEX ~ country, ~year, drop.ids = TRUE)))# 1 - Precomputing lags
summary(lm(PCGDP ~ L(PCGDP,1:2,country,year) + L(LIFEEX,0:2,country,year), wlddev))     # 2 - Ad-hoc computation in lm formula
summary(lm(PCGDP ~ L(PCGDP,1:2,country) + L(LIFEEX,0:2,country), wlddev))               # 3 - same but without year
g = qF(wlddev$country); t = qF(wlddev$year)
summary(lm(PCGDP ~ L(PCGDP,1:2,g,t) + L(LIFEEX,0:2,g,t), wlddev))                       # 4- Precomputing panel-identifiers

## Using plm:
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))
head(L(pwlddev, 0:2, 9:10))                                   # Again computing 2 lags of GDP and LIFEEX
PCGDP <- pwlddev$PCGDP                                        # A panel-Series of GDP per Capita
L(PCGDP)                                                      # Lagging the panel series. About 25x faster than plm::lag(PCGDP) on my computer
summary(lm(PCGDP ~ ., L(pwlddev, 0:2, 9:10, drop.ids = TRUE)))# Running the lm again: WORKS!
summary(lm(PCGDP ~ L(PCGDP,1:2) + L(LIFEEX,0:2), pwlddev))    # THIS DOES NOT WORK: Unfortunately lm drops the attributes of the columns, so L.default is used here and ordinary lags are computed. (generally with and attach don't retain column attributes)
LIFEEX <- pwlddev$LIFEEXA                                     # To make it work, one could create pseries
summary(lm(PCGDP ~ L(PCGDP,1:2) + L(LIFEEX,0:2)))             # THIS WORKS !!

## Using dplyr:
library(dplyr)
wlddev \%>\% group_by(country) \%>\% select(PCGDP,LIFEEX) \%>\% L(0:2)      
wlddev \%>\% group_by(country) \%>\% select(year,PCGDP,LIFEEX) \%>\% L(0:2,year) # Also using t (safer)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
