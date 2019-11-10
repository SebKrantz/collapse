\name{fscale}
\alias{fscale}
\alias{standardize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Grouped, Weighted) Scaling and Centering of Matrix-like Objects
}
\description{
\code{fscale} is a generic function to efficiently standardize (scale and center) data. \code{STD} is a wrapper around \code{fscale} representing the 'standardization operator', with more options than \code{fscale} when applied to data frames.  Standardization can be simple or groupwise, ordinary or weighted. 

Note: For centering without scaling see \code{\link{W}}, for scaling without centering see \code{\link{fsd}}.
}
\usage{
\method{fscale}{default}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)

\method{fscale}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, 
       give.names = TRUE, ...)
       
\method{fscale}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE, 
       stable.algo = TRUE, give.names = TRUE, keep.by = TRUE, 
       keep.w = TRUE, ...)
       
# Methods for compatibility with plm:

\method{fscale}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)

\method{fscale}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, stable.algo = TRUE, ...)
   STD(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE, 
       stable.algo = TRUE, give.names = TRUE, keep.pids = TRUE, 
       keep.w = TRUE, ...)

# Methods for compatibility with dplyr:

\method{fscale}{grouped_df}(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, 
       keep.groups = TRUE, keep.w = TRUE, ...)
   STD(x, w = NULL, na.rm = TRUE, stable.algo = TRUE, 
       give.names = TRUE, keep.groups = TRUE, keep.w = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (\code{\link{plm::pseries}}), panel-data.frame (\code{\link{plm::pdata.frame}}) or grouped tibble (\code{\link{dplyr::grouped_df}}).}
  \item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}
  \item{by}{\emph{STD data.frame method}: Same as g, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{w}{a numeric vector of (non-negative) weights. \code{STD} \code{data.frame} and \code{pdata.frame} methods also allow a one-sided formula i.e. \code{~ weightcol}. The \code{grouped_df} (\code{dplyr}) method supports lazy-evaluation. See Examples.}
  \item{na.rm}{skip missing values in \code{x} when computing means and sd's.}
    \item{effect}{\code{plm} methods: Select which panel identifier should be used as grouping variable. 1L means first variable in the \code{\link{index}}, 2L the second etc. if more than one integer is supplied, the corresponding index-variables are interacted. }

  \item{give.names}{\code{TRUE} will rename all transformed columns by adding a prefix "STD.".}
  \item{stable.algo}{TRUE uses a faster but numerically unstable algorithm to compute standard deviations. The default is Welfords numerically stable online algorithm. See Details.}
  \item{keep.by\cr keep.pids \cr keep.groups}{\emph{data.frame, pdata.frame and grouped_df methods}: Retain grouping / panel-identifier columns in the output. For data frames this only works if groups were passed in a formula.}
  \item{keep.w}{\emph{data.frame, pdata.frame and grouped_df methods}: Retain column containing the weights in the output. Only works if \code{w} is passed as formula / lazy-expression.}
}
\details{
If \code{g = NULL}, \code{fscale} will subract from all data points in \code{x} the mean or weighted mean (if \code{w} is supplied), and then divides this difference by the standard deviation or frequency-weighted standard deviation (if \code{w} is supplied). The result is that all vectors in \code{x} will have mean 0 and standard deviation 1. \cr

With groups supplied to \code{g}, this standardizing becomes groupwise, so that in each group the data points will have mean 0 and standard deviation 1. 

If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the mean and sd for that group will be \code{NA}, and all data points belonging to that group will also be \code{NA} in the output.

If \code{na.rm = TRUE}, the weight vector can have missing values.

By default means and standard deviations are computed using Wedford's numerically stable online algorithm. If \code{stable.algo = FALSE}, a faster but numerically unstable default algorithm is used. See \code{\link{fsum}} for more details regarding the algorithms. 

}
\value{
\code{x} standardized (mean = 0, sd = 1). See Details. 
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
## Simple centering and replacement
B(mtcars)                     
W(mtcars)
B(mtcars) + W(mtcars) == mtcars # This should be true for all elements apart from rounding errors

## Groupwise centering and replacement
B(mtcars, mtcars$cyl)                     
W(mtcars, mtcars$cyl)
B(mtcars, mtcars$cyl) + W(mtcars, mtcars$cyl) == mtcars 

## Using B() and W() in regressions
  # Several ways of running the same regression with cyl-fixed effects
  lm(W(mpg,cyl) ~ W(carb,cyl), data = mtcars)                     # Centering each variable individualls
  lm(mpg ~ carb, data = W(mtcars, ~ cyl))                         # Centering the entire data
  lm(mpg ~ carb, data = W(mtcars, ~ cyl, add.global.mean = TRUE)) # Here only the intercept changes
  lm(mpg ~ carb + B(carb,cyl), data = mtcars)                     # Procedure suggested by Mundlack (1978) - partialling out group averages amounts to the same as demeaning the data
  
  # Now with cyl, vs and am fixed effects
  lm(W(mpg,list(cyl,vs,am)) ~ W(carb,list(cyl,vs,am)), data = mtcars)
  lm(mpg ~ carb, data = W(mtcars, ~ cyl + vs + am))
  lm(mpg ~ carb + B(carb,list(cyl,vs,am)), data = mtcars) 
  
  # Now with cyl, vs and am fixed effects weighted by hp:
  lm(W(mpg,list(cyl,vs,am),hp) ~ W(carb,list(cyl,vs,am),hp), data = mtcars)
  lm(mpg ~ carb, data = W(mtcars, ~ cyl + vs + am, ~ hp))
  lm(mpg ~ carb + B(carb,list(cyl,vs,am),hp), data = mtcars) # This now gives a slightly different coefficient

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
