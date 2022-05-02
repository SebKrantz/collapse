\name{psacf}
\alias{psacf}
\alias{psacf.default}
\alias{psacf.pseries}
\alias{psacf.data.frame}
\alias{psacf.pdata.frame}
\alias{pspacf}
\alias{pspacf.default}
\alias{pspacf.pseries}
\alias{pspacf.data.frame}
\alias{pspacf.pdata.frame}
\alias{psccf}
\alias{psccf.default}
\alias{psccf.pseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auto- and Cross- Covariance and Correlation Function Estimation for Panel Series
}
\description{
\code{psacf}, \code{pspacf} and \code{psccf} compute (and by default plot) estimates of the auto-, partial auto- and cross- correlation or covariance functions for panel series. They are analogues to \code{\link{acf}}, \code{\link{pacf}} and \code{\link{ccf}}.
}
\usage{
psacf(x, \dots)
pspacf(x, \dots)
psccf(x, y, \dots)

\method{psacf}{default}(x, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance","partial"),
      plot = TRUE, gscale = TRUE, \dots)
\method{pspacf}{default}(x, g, t = NULL, lag.max = NULL, plot = TRUE, gscale = TRUE, \dots)
\method{psccf}{default}(x, y, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance"),
      plot = TRUE, gscale = TRUE, \dots)

\method{psacf}{data.frame}(x, by, t = NULL, cols = is.numeric, lag.max = NULL,
      type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, \dots)
\method{pspacf}{data.frame}(x, by, t = NULL, cols = is.numeric, lag.max = NULL,
       plot = TRUE, gscale = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{psacf}{pseries}(x, lag.max = NULL, type = c("correlation", "covariance","partial"),
      plot = TRUE, gscale = TRUE, \dots)
\method{pspacf}{pseries}(x, lag.max = NULL, plot = TRUE, gscale = TRUE, \dots)
\method{psccf}{pseries}(x, y, lag.max = NULL, type = c("correlation", "covariance"),
      plot = TRUE, gscale = TRUE, \dots)

 \method{psacf}{pdata.frame}(x, cols = is.numeric, lag.max = NULL,
      type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, \dots)
\method{pspacf}{pdata.frame}(x, cols = is.numeric, lag.max = NULL, plot = TRUE, gscale = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{a numeric vector, 'indexed_series' ('pseries'), data frame or 'indexed_frame' ('pdata.frame').}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{by}{\emph{data.frame method}: Same input as \code{g}, but also allows one- or two-sided formulas using the variables in \code{x}, i.e. \code{~ idvar} or \code{var1 + var2 ~ idvar1 + idvar2}.}
  \item{t}{a time vector or list of vectors. See \code{\link{flag}}.}
    \item{cols}{\emph{data.frame method}: Select columns using a function, column names, indices or a logical vector. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
\item{lag.max}{integer. Maximum lag at which to calculate the acf. Default is \code{2*sqrt(length(x)/ng)} where \code{ng} is the number of groups in the panel series / supplied to \code{g}.}
\item{type}{character. String giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance" or "partial".}
\item{plot}{logical. If \code{TRUE} (default) the acf is plotted.}
\item{gscale}{logical. Do a groupwise scaling / standardization of \code{x, y} (using \code{\link{fscale}} and the groups supplied to \code{g}) before computing panel-autocovariances / correlations. See Details.}
\item{\dots}{further arguments to be passed to \code{\link{plot.acf}}.}
}
\details{
If \code{gscale = TRUE} data are standardized within each group (using \code{\link{fscale}}) such that the group-mean is 0 and the group-standard deviation is 1. This is strongly recommended for most panels to get rid of individual-specific heterogeneity which would corrupt the ACF computations.

After scaling, \code{psacf}, \code{pspacf} and \code{psccf} compute the ACF/CCF by creating a matrix of panel-lags of the series using \code{\link{flag}} and then correlating this matrix with the series (\code{x, y}) using \code{\link{cor}} and pairwise-complete observations. This may require a lot of memory on large data, but is done because passing a sequence of lags to \code{\link{flag}} and thus calling \code{\link{flag}} and \code{\link{cor}} one time is much faster than calling them \code{lag.max} times. The partial ACF is computed from the ACF using a Yule-Walker decomposition, in the same way as in \code{\link{pacf}}.
}
\value{
An object of class 'acf', see \code{\link{acf}}. The result is returned invisibly if \code{plot = TRUE}.}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
%\note{
% \code{psacf} does not compute the true ACF defined as dividing the autocorvariance function of \code{x} by the variance of \code{x}. Rather, for reasons having to do with computational efficiency and assuming use on larger panel-vectors, \code{psacf} simply uses \code{\link{cor} to correlate \code{x} with it's lags (thus dividing the autocovariance by the product of the standard deviations of \code{x} and it's lag).

% For \code{plm::pseries} and \code{plm::pdata.frame}, the first index variable is assumed to be the group-id and the second the time variable. If more than 2 index variables are attached to \code{plm::pseries}, the last one is taken as the time variable and the others are taken as group-id's and interacted.

%The \code{pdata.frame} method only works for properly subsetted objects of class 'pdata.frame'. A list of 'pseries' will not work.
%}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## World Development Panel Data
head(wlddev)                                                    # See also help(wlddev)
psacf(wlddev$PCGDP, wlddev$country, wlddev$year)                # ACF of GDP per Capita
psacf(wlddev, PCGDP ~ country, ~year)                           # Same using data.frame method
psacf(wlddev$PCGDP, wlddev$country)                             # The Data is sorted, can omit t
pspacf(wlddev$PCGDP, wlddev$country)                            # Partial ACF
psccf(wlddev$PCGDP, wlddev$LIFEEX, wlddev$country)              # CCF with Life-Expectancy at Birth

psacf(wlddev, PCGDP + LIFEEX + ODA ~ country, ~year)            # ACF and CCF of GDP, LIFEEX and ODA
psacf(wlddev, ~ country, ~year, c(9:10,12))                     # Same, using cols argument
pspacf(wlddev, ~ country, ~year, c(9:10,12))                    # Partial ACF

## Using indexed data:
wldi <- findex_by(wlddev, iso3c, year)  # Creating a indexed frame
PCGDP <- wldi$PCGDP                     # Indexed Series of GDP per Capita
LIFEEX <- wldi$LIFEEX                   # Indexed Series of Life Expectancy
psacf(PCGDP)                            # Same as above, more parsimonious
pspacf(PCGDP)
psccf(PCGDP, LIFEEX)
psacf(wldi[c(9:10,12)])
pspacf(wldi[c(9:10,12)])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{ts}
