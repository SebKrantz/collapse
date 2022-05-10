\name{fgrowth}
\alias{fgrowth}
\alias{fgrowth.default}
\alias{fgrowth.matrix}
\alias{fgrowth.data.frame}
\alias{fgrowth.list}
\alias{fgrowth.pseries}
\alias{fgrowth.pdata.frame}
\alias{fgrowth.grouped_df}
\alias{G}
\alias{G.default}
\alias{G.matrix}
\alias{G.data.frame}
\alias{G.list}
\alias{G.pseries}
\alias{G.pdata.frame}
\alias{G.grouped_df}
\title{ % Lagged and Iterated
Fast Growth Rates for Time Series and Panel Data
}
\description{
\code{fgrowth} is a S3 generic to compute (sequences of) suitably lagged / leaded, iterated and compounded growth rates, obtained with via the exact method of computation or through log differencing. By default growth rates are provided in percentage terms, but any scale factor can be applied. The growth operator \code{G} is a parsimonious wrapper around \code{fgrowth}, and also provides more flexibility when applied to data frames.

}
\usage{
fgrowth(x, n = 1, diff = 1, \dots)
      G(x, n = 1, diff = 1, \dots)

\method{fgrowth}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, scale = 100, power = 1, stubs = TRUE, \dots)
\method{G}{default}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE,
  scale = 100, power = 1, stubs = TRUE, \dots)

\method{fgrowth}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, scale = 100, power = 1,
        stubs = length(n) + length(diff) > 2L, \dots)
\method{G}{matrix}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA, logdiff = FALSE,
  scale = 100, power = 1, stubs = TRUE, \dots)

\method{fgrowth}{data.frame}(x, n = 1, diff = 1, g = NULL, t = NULL, fill = NA,
        logdiff = FALSE, scale = 100, power = 1,
        stubs = length(n) + length(diff) > 2L, \dots)
\method{G}{data.frame}(x, n = 1, diff = 1, by = NULL, t = NULL, cols = is.numeric,
  fill = NA, logdiff = FALSE, scale = 100, power = 1, stubs = TRUE,
  keep.ids = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{fgrowth}{pseries}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100,
        power = 1, stubs = length(n) + length(diff) > 2L, shift = "time", \dots)
\method{G}{pseries}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100,
  power = 1, stubs = TRUE, shift = "time", \dots)

\method{fgrowth}{pdata.frame}(x, n = 1, diff = 1, fill = NA, logdiff = FALSE, scale = 100,
        power = 1, stubs = length(n) + length(diff) > 2L, shift = "time", \dots)
\method{G}{pdata.frame}(x, n = 1, diff = 1, cols = is.numeric, fill = NA, logdiff = FALSE,
  scale = 100, power = 1, stubs = TRUE, shift = "time", keep.ids = TRUE, \dots)

# Methods for grouped data frame / compatibility with dplyr:

\method{fgrowth}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE,
        scale = 100, power = 1, stubs = length(n) + length(diff) > 2L,
        keep.ids = TRUE, \dots)
\method{G}{grouped_df}(x, n = 1, diff = 1, t = NULL, fill = NA, logdiff = FALSE,
  scale = 100, power = 1, stubs = TRUE, keep.ids = TRUE, \dots)
}
\arguments{
  \item{x}{a numeric vector / time series, (time series) matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df').}
  \item{n}{integer. A vector indicating the number of lags or leads.}
  \item{diff}{integer. A vector of integers > 1 indicating the order of taking growth rates, e.g. \code{diff = 2} means computing the growth rate of the growth rate.}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{t}{a time vector or list of vectors. See \code{\link{flag}}.}
  \item{cols}{\emph{data.frame method}: Select columns to compute growth rates using a function, column names, indices or a logical vector. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{fill}{value to insert when vectors are shifted. Default is \code{NA}. }
  \item{logdiff}{logical. Compute log-difference growth rates instead of exact growth rates. See Details.}
  \item{scale}{logical. Scale factor post-applied to growth rates, default is 100 which gives growth rates in percentage terms. See Details.}
  \item{power}{numeric. Apply a power to annualize or compound growth rates e.g. \code{fgrowth(AirPassengers, 12, power = 1/12)} is equivalent to \code{((AirPassengers/flag(AirPassengers, 12))^(1/12)-1)*100}.}
  \item{stubs}{logical. \code{TRUE} will rename all computed columns by adding a prefix "L\code{n}G\code{diff}." / "F\code{n}G\code{diff}.", or "L\code{n}Dlog\code{diff}." / "F\code{n}Dlog\code{diff}." if \code{logdiff = TRUE}.}
    \item{shift}{\emph{pseries / pdata.frame methods}: character. \code{"time"} or \code{"row"}. See \code{\link{flag}} for details.}
  \item{keep.ids}{\emph{data.frame / pdata.frame / grouped_df methods}: Logical. Drop all identifiers from the output (which includes all variables passed to \code{by} or \code{t}). \emph{Note}: For 'grouped_df' / 'pdata.frame' identifiers are dropped, but the \code{"groups"} / \code{"index"} attributes are kept.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
\code{fgrowth/G} by default computes exact growth rates using \code{repeat(diff) ((x[i]/x[i-n])^power - 1)*scale}, so for \code{diff > 1} it computes growth rate of growth rates. If \code{logdiff = TRUE}, approximate growth rates are computed using \code{log(x[i]/x[i-n])*scale} for \code{diff = 1} and \code{repeat(diff-1) x[i] - x[i-n]} thereafter (usually \code{diff = 1} for log-differencing). For further details see the help pages of \code{\link{fdiff}} and \code{\link{flag}}.

}
\value{
\code{x} where the growth rate was taken \code{diff} times using lags \code{n} of itself, scaled by \code{scale}. Computations can be grouped by \code{g/by} and/or ordered by \code{t}. See Details and Examples.
}
\seealso{
\code{\link[=flag]{flag/L/F}}, \code{\link[=fdiff]{fdiff/D/Dlog}}, \link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple Time Series: AirPassengers
G(AirPassengers)                      # Growth rate, same as fgrowth(AirPassengers)
G(AirPassengers, logdiff = TRUE)      # Log-difference
G(AirPassengers, 1, 2)                # Growth rate of growth rate
G(AirPassengers, 12)                  # Seasonal growth rate (data is monthly)

head(G(AirPassengers, -2:2, 1:3))     # Sequence of leaded/lagged and iterated growth rates

# let's do some visual analysis
plot(G(AirPassengers, c(0, 1, 12)))
plot(stl(window(G(AirPassengers, 12), # Taking seasonal growth rate removes most seasonal variation
                1950), "periodic"))


## Time Series Matrix of 4 EU Stock Market Indicators, recorded 260 days per year
plot(G(EuStockMarkets,c(0,260)))                    # Plot series and annual growth rates
summary(lm(L260G1.DAX ~., G(EuStockMarkets,260)))   # Annual growth rate of DAX regressed on the
                                                    # growth rates of the other indicators

## World Development Panel Data
head(fgrowth(num_vars(wlddev), 1, 1,                    # Computes growth rates of numeric variables
             wlddev$country, wlddev$year))              # fgrowth requires external inputs..
head(G(wlddev, 1, 1, ~country, ~year))                  # Growth of numeric variables, id's attached
head(G(wlddev, 1, 1, ~country))                         # Without t: Works because data is ordered
head(G(wlddev, 1, 1, PCGDP + LIFEEX ~ country, ~year))  # Growth of GDP per Capita & Life Expectancy
head(G(wlddev, 0:1, 1, ~ country, ~year, cols = 9:10))  # Same, also retaining original series
head(G(wlddev, 0:1, 1, ~ country, ~year, 9:10,          # Dropping id columns
       keep.ids = FALSE))

}
\keyword{manip}
\keyword{ts}
