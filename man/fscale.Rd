\name{fscale}
\alias{fscale}
\alias{fscale.default}
\alias{fscale.matrix}
\alias{fscale.data.frame}
\alias{fscale.pseries}
\alias{fscale.pdata.frame}
\alias{fscale.grouped_df}
% \alias{standardize}
\alias{STD}
\alias{STD.default}
\alias{STD.matrix}
\alias{STD.data.frame}
\alias{STD.pseries}
\alias{STD.pdata.frame}
\alias{STD.grouped_df}
% - Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Grouped, Weighted) Scaling and Centering of Matrix-like Objects
}
\description{
\code{fscale} is a generic function to efficiently standardize (scale and center) data. \code{STD} is a wrapper around \code{fscale} representing the 'standardization operator', with more options than \code{fscale} when applied to matrices and data frames.  Standardization can be simple or groupwise, ordinary or weighted. Arbitrary target means and standard deviations can be set, with special options for grouped scaling and centering. It is also possible to scale data without centering i.e. perform mean-preserving scaling.
}
\usage{
fscale(x, \dots)
   STD(x, \dots)

\method{fscale}{default}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)
\method{STD}{default}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)

\method{fscale}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)
\method{STD}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1,
    stub = "STD.", \dots)

\method{fscale}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)
\method{STD}{data.frame}(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
    mean = 0, sd = 1, stub = "STD.", keep.by = TRUE, keep.w = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{fscale}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)
\method{STD}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)

\method{fscale}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, sd = 1, \dots)
\method{STD}{pdata.frame}(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
    mean = 0, sd = 1, stub = "STD.", keep.ids = TRUE, keep.w = TRUE, \dots)

# Methods for grouped data frame / compatibility with dplyr:

\method{fscale}{grouped_df}(x, w = NULL, na.rm = TRUE, mean = 0, sd = 1,
       keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{STD}{grouped_df}(x, w = NULL, na.rm = TRUE, mean = 0, sd = 1,
    stub = "STD.", keep.group_vars = TRUE, keep.w = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{a numeric vector, matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df').}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{by}{\emph{STD data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
    \item{cols}{\emph{STD (p)data.frame method}: Select columns to scale using a function, column names, indices or a logical vector. Default: All numeric columns. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{w}{a numeric vector of (non-negative) weights. \code{STD} data frame and \code{pdata.frame} methods also allow a one-sided formula i.e. \code{~ weightcol}. The \code{grouped_df} (\emph{dplyr}) method supports lazy-evaluation. See Examples.}
  \item{na.rm}{logical. Skip missing values in \code{x} or \code{w} when computing means and sd's.}
  \item{effect}{\emph{plm} methods: Select which panel identifier should be used as group-id. 1L takes the first variable in the \link[=indexing]{index}, 2L the second etc.. Index variables can also be called by name using a character string. More than one variable can be supplied. }
  \item{stub}{a prefix or stub to rename all transformed columns. \code{FALSE} will not rename columns.}
  \item{mean}{the mean to center on (default is 0). If \code{mean = FALSE}, no centering will be performed. In that case the scaling is mean-preserving. A numeric value different from 0 (i.e. \code{mean = 5}) will be added to the data after subtracting out the mean(s), such that the data will have a mean of 5. A special option when performing grouped scaling and centering is \code{mean = "overall.mean"}. In that case the overall mean of the data will be added after subtracting out group means.}
 \item{sd}{the standard deviation to scale the data to (default is 1). A numeric value different from 0 (i.e. \code{sd = 3}) will scale the data to have a standard deviation  of 3. A special option when performing grouped scaling is \code{sd = "within.sd"}. In that case the within standard deviation (= the standard deviation of the group-centered series) will be calculated and applied to each group. The results is that the variance of the data within each group is harmonized without forcing a certain variance (such as 1).}
  \item{keep.by, keep.ids, keep.group_vars}{\emph{data.frame, pdata.frame and grouped_df methods}: Logical. Retain grouping / panel-identifier columns in the output. For \code{STD.data.frame} this only works if grouping variables were passed in a formula.}
  \item{keep.w}{\emph{data.frame, pdata.frame and grouped_df methods}: Logical. Retain column containing the weights in the output. Only works if \code{w} is passed as formula / lazy-expression.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
If \code{g = NULL}, \code{fscale} by default (column-wise) subtracts the mean or weighted mean (if \code{w} is supplied) from all data points in \code{x}, and then divides this difference by the standard deviation or frequency-weighted standard deviation. The result is that all columns in \code{x} will have a (weighted) mean 0 and (weighted) standard deviation 1. Alternatively, data can be scaled to have a mean of \code{mean} and a standard deviation of \code{sd}. If \code{mean = FALSE} the data is only scaled (not centered) such that the mean of the data is preserved. \cr

Means and standard deviations are computed using Welford's numerically stable online algorithm.

With groups supplied to \code{g}, this standardizing becomes groupwise, so that in each group (in each column) the data points will have mean \code{mean} and standard deviation \code{sd}. Naturally if \code{mean = FALSE} then each group is just scaled and the mean is preserved. For centering without scaling see \code{\link{fwithin}}.

If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the mean and sd for that group will be \code{NA}, and all data points belonging to that group will also be \code{NA} in the output.

If \code{na.rm = TRUE}, means and sd's are computed (column-wise) on the available data points, and also the weight vector can have missing values. In that case, the weighted mean an sd are computed on (column-wise) \code{complete.cases(x, w)}, and \code{x} is scaled using these statistics. \emph{Note} that \code{fscale} will not insert a missing value in \code{x} if the weight for that value is missing, rather, that value will be scaled using a weighted mean and standard-deviated computed without itself! (The intention here is that a few (randomly) missing weights shouldn't break the computation when \code{na.rm = TRUE}, but it is not meant for weight vectors with many missing values. If you don't like this behavior, you should prepare your data using \code{x[is.na(w), ] <- NA}, or impute your weight vector for non-missing \code{x}).

Special options for grouped scaling are \code{mean = "overall.mean"} and \code{sd = "within.sd"}. The former group-centers vectors on the overall mean of the data (see \code{\link{fwithin}} for more details) and the latter scales the data in each group to have the within-group standard deviation (= the standard deviation of the group-centered data). Thus scaling a grouped vector with options \code{mean = "overall.mean"} and \code{sd = "within.sd"} amounts to removing all differences in the mean and standard deviations between these groups. In weighted computations, \code{mean = "overall.mean"} will subtract weighted group-means from the data and add the overall weighted mean of the data, whereas \code{sd = "within.sd"} will compute the weighted within- standard deviation and apply it to each group.

}
\value{
\code{x} standardized (mean = mean, standard deviation = sd), grouped by \code{g/by}, weighted with \code{w}. See Details.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
\note{
For centering without scaling see \code{\link[=fwithin]{fwithin/W}}. For simple not mean-preserving scaling use \code{\link[=fsd]{fsd(..., TRA = "/")}}. To sweep pre-computed means and scale-factors out of data see \code{\link{TRA}}.
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{fwithin}}, \code{\link{fsd}}, \code{\link{TRA}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple Scaling & Centering / Standardizing
head(fscale(mtcars))               # Doesn't rename columns
head(STD(mtcars))                  # By default adds a prefix
qsu(STD(mtcars))                   # See that is works
qsu(STD(mtcars, mean = 5, sd = 3)) # Assigning a mean of 5 and a standard deviation of 3
qsu(STD(mtcars, mean = FALSE))     # No centering: Scaling is mean-preserving

## Panel Data
head(fscale(get_vars(wlddev,9:12), wlddev$iso3c))   # Standardizing 4 series within each country
head(STD(wlddev, ~iso3c, cols = 9:12))              # Same thing using STD, id's added
pwcor(fscale(get_vars(wlddev,9:12), wlddev$iso3c))  # Correlaing panel series after standardizing

fmean(get_vars(wlddev, 9:12))                       # This calculates the overall means
fsd(fwithin(get_vars(wlddev, 9:12), wlddev$iso3c))  # This calculates the within standard deviations
head(qsu(fscale(get_vars(wlddev, 9:12),             # This group-centers on the overall mean and
    wlddev$iso3c,                                   # group-scales to the within standard deviation
    mean = "overall.mean", sd = "within.sd"),       # -> data harmonized in the first 2 moments
    by = wlddev$iso3c))

## Indexed data
wldi <- findex_by(wlddev, iso3c, year)
head(STD(wldi))                                  # Standardizing all numeric variables by country
head(STD(wldi, effect = 2L))                     # Standardizing all numeric variables by year

## Weighted Standardizing
weights = abs(rnorm(nrow(wlddev)))
head(fscale(get_vars(wlddev,9:12), wlddev$iso3c, weights))
head(STD(wlddev, ~iso3c, weights, 9:12))
\donttest{ % No code relying on suggested package
# Grouped data
library(magrittr)
wlddev \%>\% fgroup_by(iso3c) \%>\% fselect(PCGDP,LIFEEX) \%>\% STD()
wlddev \%>\% fgroup_by(iso3c) \%>\% fselect(PCGDP,LIFEEX) \%>\% STD(weights) # weighted standardizing
wlddev \%>\% fgroup_by(iso3c) \%>\% fselect(PCGDP,LIFEEX,POP) \%>\% STD(POP) # weighting by POP ->
# ..keeps the weight column unless keep.w = FALSE
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")

