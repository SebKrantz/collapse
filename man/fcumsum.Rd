\name{fcumsum}
\alias{fcumsum}
\alias{fcumsum.default}
\alias{fcumsum.matrix}
\alias{fcumsum.data.frame}
\alias{fcumsum.pseries}
\alias{fcumsum.pdata.frame}
\alias{fcumsum.grouped_df}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Grouped, Ordered) Cumulative Sum for Matrix-Like Objects
}
\description{
\code{fcumsum} is a generic function that computes the (column-wise) cumulative sum of \code{x}, (optionally) grouped by \code{g} and/or ordered by \code{o}. Several options to deal with missing values are provided.
}
\usage{
fcumsum(x, \dots)

\method{fcumsum}{default}(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, \dots)

\method{fcumsum}{matrix}(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, \dots)

\method{fcumsum}{data.frame}(x, g = NULL, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{fcumsum}{pseries}(x, na.rm = .op[["na.rm"]], fill = FALSE, shift = "time", \dots)

\method{fcumsum}{pdata.frame}(x, na.rm = .op[["na.rm"]], fill = FALSE, shift = "time", \dots)

# Methods for grouped data frame / compatibility with dplyr:

\method{fcumsum}{grouped_df}(x, o = NULL, na.rm = .op[["na.rm"]], fill = FALSE, check.o = TRUE,
        keep.ids = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector / time series, (time series) matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df').}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{o}{a vector or list of vectors providing the order in which the elements of \code{x} are cumulatively summed. Will be passed to \code{\link{radixorderv}} unless \code{check.o = FALSE}.}
\item{na.rm}{logical. Skip missing values in \code{x}. Defaults to \code{TRUE} and implemented at very little computational cost.}
\item{fill}{if \code{na.rm = TRUE}, setting \code{fill = TRUE} will overwrite missing values with the previous value of the cumulative sum, starting from 0.}
\item{check.o}{logical. Programmers option: \code{FALSE} prevents passing \code{o} to \code{\link{radixorderv}}, requiring \code{o} to be a valid ordering vector that is integer typed with each element in the range \code{[1, length(x)]}. This gives some extra speed, but will terminate R if any element of \code{o} is too large or too small. }
    \item{shift}{\emph{pseries / pdata.frame methods}: character. \code{"time"} or \code{"row"}. See \code{\link{flag}} for details. The argument here does not control 'shifting' of data but rather the order in which elements are summed.}
  \item{keep.ids}{\emph{pdata.frame / grouped_df methods}: Logical. Drop all identifiers from the output (which includes all grouping variables and variables passed to \code{o}). \emph{Note}: For grouped / panel data frames identifiers are dropped, but the \code{"groups"} / \code{"index"} attributes are kept.}
\item{\dots}{arguments to be passed to or from other methods.}
}
\details{
If \code{na.rm = FALSE}, \code{fcumsum} works like \code{\link{cumsum}} and propagates missing values. The default \code{na.rm = TRUE} skips missing values and computes the cumulative sum on the non-missing values. Missing values are kept. If \code{fill = TRUE}, missing values are replaced with the previous value of the cumulative sum (starting from 0), computed on the non-missing values.

By default the cumulative sum is computed in the order in which elements appear in \code{x}. If \code{o} is provided, the cumulative sum is computed in the order given by \code{radixorderv(o)}, without the need to first sort \code{x}. This applies as well if groups are used (\code{g}), in which case the cumulative sum is computed separately in each group.

The \emph{pseries} and \emph{pdata.frame} methods assume that the last factor in the \link[=findex]{index} is the time-variable and the rest are grouping variables. The time-variable is passed to \code{radixorderv} and used for ordered computation, so that cumulative sums are accurately computed regardless of whether the panel-data is ordered or balanced.

\code{fcumsum} explicitly supports integers. Integers in R are bounded at bounded at +-2,147,483,647, and an integer overflow error will be provided if the cumulative sum (within any group) exceeds +-2,147,483,647. In that case data should be converted to double beforehand.

}
\value{
the cumulative sum of values in \code{x}, (optionally) grouped by \code{g} and/or ordered by \code{o}. See Details and Examples.
}

\seealso{
\code{\link{fdiff}}, \code{\link{fgrowth}}, \link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Non-grouped
fcumsum(AirPassengers)
head(fcumsum(EuStockMarkets))
fcumsum(mtcars)

# Non-grouped but ordered
o <- order(rnorm(nrow(EuStockMarkets)))
all.equal(copyAttrib(fcumsum(EuStockMarkets[o, ], o = o)[order(o), ], EuStockMarkets),
          fcumsum(EuStockMarkets))

## Grouped
head(with(wlddev, fcumsum(PCGDP, iso3c)))

## Grouped and ordered
head(with(wlddev, fcumsum(PCGDP, iso3c, year)))
head(with(wlddev, fcumsum(PCGDP, iso3c, year, fill = TRUE)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{ts} % __ONLY ONE__ keyword per line
