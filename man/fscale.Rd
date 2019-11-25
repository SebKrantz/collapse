\name{fscale-STD}
\alias{fscale}
\alias{fscale.default}
\alias{fscale.matrix}
\alias{fscale.data.frame}
\alias{fscale.pseries}
\alias{fscale.pdata.frame}
\alias{fscale.grouped_df}
\alias{standardize}
\alias{STD}
\alias{STD.default}
\alias{STD.matrix}
\alias{STD.data.frame}
\alias{STD.pseries}
\alias{STD.pdata.frame}
\alias{STD.grouped_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Grouped, Weighted) Scaling and Centering of Matrix-like Objects
}
\description{
\code{fscale} is a generic function to efficiently standardize (scale and center) data. \code{STD} is a wrapper around \code{fscale} representing the 'standardization operator', with more options than \code{fscale} when applied to matrices and data frames.  Standardization can be simple or groupwise, ordinary or weighted.

\emph{Note}: For centering without scaling see \code{\link{W}}, for scaling without centering use \code{\link[=fsd]{fsd(..., TRA = "/")}}.
}
\usage{
\method{fscale}{default}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)

\method{fscale}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE,
       stub = "STD.", ...)

\method{fscale}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
       keep.by = TRUE, keep.w = TRUE, stable.algo = TRUE, stub = "STD.", ...)

# Methods for compatibility with plm:

\method{fscale}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)

\method{fscale}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
       keep.ids = TRUE, keep.w = TRUE, stable.algo = TRUE, stub = "STD.", ...)

# Methods for compatibility with dplyr:

\method{fscale}{grouped_df}(x, w = NULL, na.rm = TRUE, keep.group_vars = TRUE,
       keep.w = TRUE, stable.algo = TRUE, ...)
   STD(x, w = NULL, na.rm = TRUE, keep.group_vars = TRUE,
       keep.w = TRUE, stable.algo = TRUE, stub = "STD.", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (\code{\link{plm::pseries}}), panel-data.frame (\code{\link{plm::pdata.frame}}) or grouped tibble (\code{\link{dplyr::grouped_df}}).}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{by}{\emph{STD data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{w}{a numeric vector of (non-negative) weights. \code{STD} \code{data.frame} and \code{pdata.frame} methods also allow a one-sided formula i.e. \code{~ weightcol}. The \code{grouped_df} (\code{dplyr}) method supports lazy-evaluation. See Examples.}
  \item{na.rm}{logical. skip missing values in \code{x} when computing means and sd's.}
    \item{effect}{\code{plm} methods: Select which panel identifier should be used as grouping variable. 1L means first variable in the \code{\link{index}}, 2L the second etc. if more than one integer is supplied, the corresponding index-variables are interacted. }

  \item{stub}{a prefix or stub to rename all transformed columns. \code{FALSE} will not rename columns.}
  \item{stable.algo}{logical. TRUE uses a faster but numerically unstable algorithm to compute standard deviations. The default is Welfords numerically stable online algorithm. See Details.}
  \item{keep.by\cr keep.ids \cr keep.group_vars}{\emph{data.frame, pdata.frame and grouped_df methods}: Logical. Retain grouping / panel-identifier columns in the output. For \code{STD.data.frame} this only works if grouping variables were passed in a formula.}
  \item{keep.w}{\emph{data.frame, pdata.frame and grouped_df methods}: Logical. Retain column containing the weights in the output. Only works if \code{w} is passed as formula / lazy-expression.}
  \item{...}{arguments to be passed to or from other methods.}
}
\details{
If \code{g = NULL}, \code{fscale} subracts from all data points in \code{x} the mean or weighted mean (if \code{w} is supplied), and then divides this difference by the standard deviation or frequency-weighted standard deviation (if \code{w} is supplied). The result is that all vectors in \code{x} will have mean 0 and standard deviation 1. \cr

With groups supplied to \code{g}, this standardizing becomes groupwise, so that in each group the data points will have mean 0 and standard deviation 1.

If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the mean and sd for that group will be \code{NA}, and all data points belonging to that group will also be \code{NA} in the output.

If \code{na.rm = TRUE}, means and sd's are computed on the available data points, and also the weight vector can have missing values.

By default means and standard deviations are computed using Wedford's numerically stable online algorithm. If \code{stable.algo = FALSE}, a faster but numerically unstable default algorithm is used. See \code{\link{fsd}} for more details regarding the algorithms.

}
\value{
\code{x} standardized (mean = 0, sd = 1), grouped by \code{g/by}. See Details.
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
\code{\link[=B]{B/W}}, \link{Fast Statistical Functions}, \code{\link{TRA}}, \link{Data Transformations}
}
\examples{
## Simple Scaling & Centering / Standardizing
fscale(mtcars)             # Doesn't rename columns
STD(mtcars)                # By default adds a prefix
qsu(STD(mtcars))           # See that is works

## Panel-Data
head(fscale(get_vars(wlddev,9:12), wlddev$iso3c))   # Standardizing 4 series within each country (Note: get_vars(wlddev,9:12) is 2x faster than wlddev[9:12])
head(STD(wlddev, ~iso3c, cols = 9:12))              # Same thing using STD, with identifier column added unless keep.by = FALSE
pwcor(fscale(get_vars(wlddev,9:12), wlddev$iso3c))  # Correlaing these panel-series after standardizing values within each country

## Using plm
pwlddev <- plm::pdata.frame(wlddev, index = c("iso3c","year"))
head(STD(pwlddev))                                  # Standardizing all numeric variables by country
head(STD(pwlddev, effect = 2L))                     # Standardizing all numeric variables by year

## Weighted Standardizing
weights = abs(rnorm(nrow(wlddev)))
head(fscale(get_vars(wlddev,9:12), wlddev$iso3c, weights))
head(STD(wlddev, ~iso3c, weights, 9:12))

# Using dplyr
library(dplyr)
wlddev \%>\% group_by(iso3c) \%>\% select(PCGDP,LIFEEX) \%>\% STD
wlddev \%>\% group_by(iso3c) \%>\% select(PCGDP,LIFEEX) \%>\% STD(weights) # weighted standardizing
wlddev \%>\% group_by(iso3c) \%>\% select(PCGDP,LIFEEX,ODA) \%>\% STD(ODA) # weighting by ODA -> keeps the weight column unless keep.w = FALSE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
