\name{psacf}
\alias{psacf}
% \alias{psacf.default}
% \alias{psacf.pseries}
% \alias{psacf.data.frame}
% \alias{psacf.pdata.frame}
\alias{pspacf}
% \alias{pspacf.default}
% \alias{pspacf.pseries}
% \alias{pspacf.data.frame}
% \alias{pspacf.pdata.frame}
\alias{psccf}
% \alias{psccf.default}
% \alias{psccf.pseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auto- and Cross- Covariance and -Correlation Function Estimation for Panel-Series
}
\description{
\code{psacf}, \code{pspacf} and \code{psccf} compute (and by default plot) estimates of the auto-, partial auto- and cross- correlation or covariance functions for panel-vectors and \code{plm::pseries}. They are faithful analogues to \code{stats::acf}, \code{stats::pacf} and \code{stats::ccf}.
}
\usage{
psacf(x, ...)
pspacf(x, ...)
psccf(x, y, ...)

## Default S3 methods:
 psacf(x, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance","partial"),
       plot = TRUE, gscale = TRUE, ...)
pspacf(x, g, t = NULL, lag.max = NULL, plot = TRUE, gscale = TRUE, ...)
 psccf(x, y, g, t = NULL, lag.max = NULL, type = c("correlation", "covariance"),
       plot = TRUE, gscale = TRUE, ...)

## S3 methods for class 'pseries'
 psacf(x, lag.max = NULL, type = c("correlation", "covariance","partial"),
       plot = TRUE, gscale = TRUE, ...)
pspacf(x, lag.max = NULL, plot = TRUE, gscale = TRUE, ...)
 psccf(x, y, lag.max = NULL, type = c("correlation", "covariance"),
       plot = TRUE, gscale = TRUE, ...)

## S3 methods for class 'data.frame'
 psacf(x, by, t = NULL, cols = is.numeric, lag.max = NULL,
       type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...)
pspacf(x, by, t = NULL, cols = is.numeric, lag.max = NULL,
       type = c("correlation", "covariance","partial"), plot = TRUE, gscale = TRUE, ...)

## S3 methods for class 'pdata.frame'
 psacf(x, lag.max = NULL, type = c("correlation", "covariance","partial"),
       plot = TRUE, gscale = TRUE, ...)
pspacf(x, lag.max = NULL, plot = TRUE, gscale = TRUE, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x, y}{a numeric vector, panel-series (\code{plm::pseries}), data.frame or
panel-data-frame (\code{plm::pdata.frame}).}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x, y}. }
  \item{by}{\emph{data frame method}: Same input as \code{g}, but also allows one- or two-sided formulas using the variables in \code{x}, i.e. \code{~ idvar} or \code{var1 + var2 ~ idvar1 + idvar2}.}

    \item{t}{same input as g, to indicate the time-variable. For secure computations on unordered panel-vectors. Data frame method also takes one-sided formula i.e. \code{~time}.}
    \item{cols}{\emph{data.frame method}: Select columns using a function, column names or indices. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
\item{lag.max}{maximum lag at which to calculate the acf. Default is \code{2*sqrt(length(x)/ng)} where \code{ng} is the number of groups in the panel-series / supplied to \code{g}.}
\item{type}{character string giving the type of acf to be computed. Allowed values are "correlation" (the default), "covariance" or "partial".}
\item{plot}{logical. If \code{TRUE} (the default) the acf is plotted.}
\item{gscale}{logical. Do a groupwise scaling / standardization of \code{x, y} (using \code{collapse::fscale} and the groups supplied to \code{g}) before computing panel-autocorrelations. Recommended and \code{TRUE} by default.}
\item{...}{further arguments to be passed to \code{stats:::plot.acf}.}
}
\details{
If \code{gscale = TRUE} data are standardized within each group (using \code{collapse::fscale}) such that the group-mean is 0 and the group-standard deviation is 1. This is strongly recommended for most panels to get rid of individual-specific heterogeneity which would corrupt the ACF computations.

After scaling, \code{psacf}, \code{pspacf} and \code{psccf} compute the ACF/CCF by creating a matrix of panel-lags of the series using \code{collapse::flag} and then correlating this matrix with the series (\code{x, y}) using \code{stats::cor} and pairwise-complete observations. This may require a lot of memory on large data, but is done because passing a sequence of lags to \code{collapse::flag} and thus calling \code{collapse::flag} and \code{stats::cor} one time is much faster than calling them \code{lag.max} times. The partial ACF is computed from the ACF in the same way as in \code{stats::pacf}.
}
\value{
An object of class "acf", see ?stats::acf. The result is returned invisibly if plot is TRUE.}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
\note{
For \code{plm::pseries} and \code{plm::pdata.frame}, the first index variable is taken to be the group-id and the second the time variable. If more than 2 index variables are attached to \code{plm::pseries}, the last one is taken as the time variable and the others are taken as group-id's and interacted.

The \code{pdata.frame} method only works for properly subsetted objects of class 'pdata.frame'. A list of 'pseries' won't work.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[=time-series-panel-series]{Time-Series and Panel-Series} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## World Development Panel Data
head(wlddev)                                                    # see also help(wlddev)
psacf(wlddev$PCGDP, wlddev$country, wlddev$year)                # ACF of GDP per Capita
psacf(wlddev, PCGDP ~ country, ~year)                           # Same using data.frame method
psacf(wlddev$PCGDP, wlddev$country)                             # The Data is sorted, so we could leave out t
pspacf(wlddev$PCGDP, wlddev$country)                            # Partial ACF
psccf(wlddev$PCGDP, wlddev$LIFEEX, wlddev$country)              # CCF with Life-Expectancy at Birth

psacf(wlddev, PCGDP + LIFEEX + ODA ~ country, ~year)            # ACF and CCF of GDP, Life Expactancy and ODA received
psacf(wlddev, ~ country, ~year, c(9:10,12))                     # Same, using cols argument
pspacf(wlddev, ~ country, ~year, c(9:10,12))                    # Partial ACF

## Using plm:
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))# Creating a Panel-Data Frame of this data
PCGDP <- pwlddev$PCGDP                                          # Panel-Series of GDP per Capita
LIFEEX <- pwlddev$LIFEEX                                        # Panel-Series of Life Expectancy
psacf(PCGDP)                                                    # Same as above, but more parsimonious
pspacf(PCGDP)
psccf(PCGDP, LIFEEX)
psacf(pwlddev[c(9:10,12)])
pspacf(pwlddev[c(9:10,12)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
