\name{differences & growth rates}
\alias{diff&growth}
\alias{D}
\alias{G}
\alias{fdiff}
\alias{fgrowth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Lagged and Iterated Differences and Growth Rates for Time-Series and Panel Data
}
\description{
\code{fdiff} and \code{fgrowth} are S3 generics to compute (sequences of) suitably lagged or leaded and iterated differences and growth rates / log-differences, respectively. \code{D} and \code{G} are parsimonious wrappers around \code{fdiff} and \code{fgrowth} representing the 'difference-operator' and the 'growth-operator'. \code{D}/\code{G} provide more options than \code{fdiff}/\code{fgrowth} when applied to data.frame's, but are otherwise identical.
}
\usage{
## Default S3 methods:
  fdiff/D(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)
fgrowth/G(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, give.names = TRUE, ...)

## S3 methods for class 'matrix'
  fdiff/D(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)
fgrowth/G(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, give.names = TRUE, 

## S3 methods for class 'data.frame'
  fdiff(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, give.names = TRUE, ...)
fgrowth(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE, give.names = TRUE, 
      D(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric, 
        fill = NA, give.names = TRUE, drop.ids = FALSE, ...) 
      G(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric, 
        fill = NA, logdiff = FALSE, give.names = TRUE, drop.ids = FALSE, ...) 

# Methods for compatibility with plm:

## S3 methods for class 'pseries'
  fdiff/D(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...)
fgrowth/G(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, give.names = TRUE, ...)

## S3 methods for class 'pdata.frame'
  fdiff(x, n = 1, diff = 1, fill = NA, give.names = TRUE, ...)
fgrowth(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, give.names = TRUE, ...)
      D(x, n = 1, diff = 1, cols = is.numeric, fill = NA, give.names = TRUE, drop.ids = FALSE, ...)
      G(x, n = 1, diff = 1, cols = is.numeric, fill = NA, logdiff = FALSE, give.names = TRUE, drop.ids = FALSE, ...)
      
# Methods for compatibility with dplyr:

## S3 methods for class 'grouped_df'
  fdiff/D(x, n = 1, diff = 1, t = NULL, fill = NA, give.names = TRUE, drop.ids = FALSE, ...)
fgrowth/G(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE, give.names = TRUE, drop.ids = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (\code{\link{plm::pseries}}), panel-data.frame (\code{\link{plm::pdata.frame}}) or grouped tibble (\code{\link{dplyr::grouped_df}}).}
  \item{n}{a integer vector indicating the number of lags or leads.}
  \item{diff}{a vector of integers > 1 indicating the order of differencing / growth rates.}
  \item{g}{a factor, \code{GRP()} object, atomic vector (> internally converted to factor) or a list of vectors / factors (> internally converted to a \code{GRP()} object) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: same as g, but also allows for column-names, column-indices and one- or two-sided formulas. See Examples.}
  \item{t}{same input as g, to indicate the time-variable. For safe computation of differnces/growth rates on unordered time-series and panels. \emph{Notes}: data.frame method also allows name, index or one-sided formula i.e. \code{~time}. grouped_df method also allows lazy-evaluation i.e. \code{time} (no quotes).}
  \item{cols}{\emph{data.frame method}: Select columns to difference/compute growth rates using a function, column names or indices. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{logdiff}{compute log-differences instead of exact growth rates. See Details.}
  \item{give.names}{\code{TRUE} will rename all differenced columns by adding a prefix "L\code{n}D\code{diff}." / "F\code{n}D\code{diff}." and a prefix "L\code{n}G\code{diff}." / "F\code{n}G\code{diff}." for growth rates.}
  \item{drop.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: drop all panel-identifiers from the output (which includes all variables passed to 'by' or 't'). For \code{pdata.frame}'s identifiers are dropped, but the 'index' attribute is always kept.}
}
\details{
By default, \code{fdiff/D|fgrowth/G} return \code{x} with all columns differenced | converted to growth rates. Differences are computed as \code{repeat(diff){x[i] - x[i-n]}}, growth rates as \code{repeat(diff){(x[i] - x[i-n])/x[i-n]*100}} and log-differences as \code{repeat(diff){(log(x[i]) - log(x[i-n]))*100}}. Setting \code{diff = 2} thus returns differences of differences | growth rates of growth rates etc... and setting \code{n = 2} returns simple differences computed by subtracting twice-lagged \code{x} from \code{x}. It is also possible to compute forward differences | growth rates by passing negative \code{n} values. \code{n} also supports sequences of integers (lags), and \code{diff} supports positive and continuous sequences of integers (differences):

If more than one value is passed to \code{n} and/or \code{diff}, the data is expanded-wide as follows: If \code{x} is an atomic vector or time-series, a (time-series) matrix is returned with columns ordered first by lag, then by difference. If \code{x} is a matrix or data.frame, each column is expanded in like manor such that the output has \code{length(n)*length(diff)*ncol(x)} with columns ordered first by column name, then by lag, then by difference.

With groups/panel-identifiers supplied to \code{g/by}, \code{fdiff/D|fgrowth/G} efficiently compute panel-differences | growth rates by inserting \code{fill} elements in the right places. If \code{t} is left empty, the data needs to be ordered such that all values belonging to a group are consecutive and in the right order. It is not necessary that the groups themselves occur in the right order. If time-variable(s) are supplied to \code{t}, the panel is fully identified and differences | growth rates can be securely computed even if the data is completely unordered (in that case data is shifted around and \code{fill} values are inserted in such a way that if the data were sorted afterwards the result would be identical to computing on sorted data). Internally this works by using the grouping- and time-variables to create an ordering and then accessing the panel-vector(s) through this ordering. If the data is just a bit unordered, such computations are nearly as fast as computations on ordered data (without \code{t}), however, if the data is very unordered, it can take significantly longer. Since most panel-data come perfectly or pretty ordered, I recommend always supplying \code{t} to be on the safe-side. 

It is also possible to compute differences | growth rates on unordered vectors / time-series (thus utilizing \code{t} but leaving \code{g/by} empty). 

The methods applying to \code{plm} objects (panel-series and panel-data.frames) automatically utilize the panel-identifiers attached to these objects and thus securely compute fully identified panel-differences. If these objects have > 2 panel-identifiers attached to them, the last identifier is assumed to be the time-variable, and the others are taken as grouping-variables and interacted.  
}
\value{
\code{fdiff/D} returns \code{x} differenced \code{diff} times using lags \code{n} of itself. 
\code{fgrowth/G} returns \code{x} where the growth rate or log-difference was taken \code{diff} times using lags \code{n} of itself. See Details and Examples. 
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
D(AirPassengers)                      # 1st difference, same as fdiff(AirPassengers)
D(AirPassengers,-1)                   # forward difference
G(AirPassengers)                      # growth rate, same as fgrowth(AirPassengers)
G(AirPassengers, logdiff = TRUE)      # log-difference
D(AirPassengers,1,2)                  # second difference
G(AirPassengers,1,2)                  # growth rate of growth rate
D(AirPassengers,12)                   # seasonal difference (data is monthly)
G(AirPassengers,12)                   # seasonal growth rate (data is monthly)

D(AirPassengers,-2:2,1:3)             # sequence of leaded/lagged and iterated differences - output as matrix

# let's do some visual analysis
plot(AirPassengers)                   # plot the series - seasonal pattern is evident
plot(D(AirPassengers,c(1,12),1:2))    # plotting the ordinary and seasonal first and second differences (second differences could be useful to spot periods of rapid expansion of the industry)
plot(G(AirPassengers,c(1,12),1:2))    # same using growth rates

## Time-Series Matrix of 4 EU Stock Market Indicators
D(EuStockMarkets,1:2,1:2)                          # 2 lags and 2 differences for each indicator
summary(lm(G1.DAX ~., data = G(EuStockMarkets)))   # growth rate of DAX regressed on the growth rates of the other indicators

## World Development Panel Data
head(fgrowth(num.vars(wlddev), 1, 1, wlddev$country, wlddev$year))   # This computes the growth rate of all variables. With fgrowth/fdiff all inputs need to be passed externally, so G/D is more convenient: 
head(G(wlddev, 1, 1, ~country, ~year))                  # Growth rate of all Numeric Variables, identifiers attached
head(G(wlddev, 1, 1, ~country))                         # Without t: Also works because data is ordered 
head(G(wlddev, 1, 1, PCGDP + LIFEEX ~ country, ~year))  # Growth rate of GDP per Capita and Life Expectancy
head(G(wlddev, 0:1, 1, ~ country, ~year, cols = 9:10))  # Same, also retaining original series
head(G(wlddev, 0:1, 1, 1, 4, 9:10))                     # Same, even more parsimonious
head(G(wlddev, 0:1, 1, 1, 4, 9:10, drop.ids = TRUE))    # Adding growth rate of growth rate and dropping id columns

# Dynamic Panel-Data Models: 
summary(lm(G(PCGDP,1,1,iso3c,year) ~ # GDP growth regressed on it's lagged level and 
           L(PCGDP,1,iso3c,year) +   # the growth rate of Life Expanctancy
           G(LIFEEX,1,1,iso3c,year), data = wlddev)) 

g = qF(wlddev$country)               # omitting t and precomputing g allows for a bit more parsimonious specification
summary(lm(G(PCGDP,1,1,g) ~ L(PCGDP,1,g) + G(LIFEEX,1,1,g), wlddev)) 
           
summary(lm(G1.PCGDP ~.,              # now adding level and lagged level of LIFEEX and lagged growth rates
        data = L(G(wlddev,0:1,1,1,4,9:10),0:1,1,2, drop.ids = TRUE)[-1]))           
           
## Using plm can make things easier, but be careful what you are doing:
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))
head(G(pwlddev, 0:1, 1, 9:10))       # Again computing growth rates of LIFEEX and PCGDP
PCGDP <- pwlddev$PCGDP               # A panel-Series of GDP per Capita
D(PCGDP)                             # Differencing the panel series.
summary(lm(G1.PCGDP ~.,              # running the dynamic model again -> code becomes a bit simpler
        data = L(G(pwlddev,0:1,1,9:10),0:1,drop.ids = TRUE)[-1]))           
        
# One could be tempted to also do something like this, but THIS DOES NOT WORK!!!: lm drops the attributes (-> with(pwlddev, PCGDP) returns a vector and not a panel-series so G.default and L.matrix are used) 
summary(lm(G(PCGDP) ~ L(G(PCGDP,0:1)) + L(G(LIFEEX,0:1),0:1), pwlddev)) 

# To make it work, one needs to create pseries (note: attach(pwlddev) also won't work as this also removes attributes)
LIFEEX <- pwlddev$LIFEEXA                                     
summary(lm(G(PCGDP) ~ L(G(PCGDP,0:1)) + L(G(LIFEEX,0:1),0:1)))# THIS WORKS !!

## Using dplyr:
library(dplyr)
wlddev \%>\% group_by(country) \%>\% select(PCGDP,LIFEEX) \%>\% D(0:1,1:2)           # Adding a first and second difference
wlddev \%>\% group_by(country) \%>\% select(year,PCGDP,LIFEEX) \%>\% D(0:1,1:2,year) # Also using t (safer)
wlddev \%>\% group_by(country) \%>\% select(year,PCGDP,LIFEEX) \%>\% G(0:1,1:2,year, drop.ids = TRUE) # Growth rates, dropping id's

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
