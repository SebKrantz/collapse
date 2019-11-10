\name{B&W}
\alias{B&W}
\alias{B}
\alias{W}
\alias{between}
\alias{within}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Between (Averaging) and Within (Centering) Transformations
}
\description{
\code{B} is a generic function to efficiently obtain averaged data. \code{W} is also generic and centers data on those averages. Averages can be simple or groupwise, ordinary or weighted. 
}
\usage{
\method{B}{default}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, ...)
W(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...)

\method{B}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, give.names = FALSE, ...)
W(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, give.names = FALSE, ...)

\method{B}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, give.names = FALSE, ...)
W(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, give.names = FALSE, ...)

# Methods for compatibility with plm:

\method{B}{pseries}(x, w = NULL, na.rm = TRUE, fill = FALSE, ...)
W(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE, ...)

\method{B}{pdata.frame}(x, w = NULL, na.rm = TRUE, fill = FALSE, give.names = FALSE, ...)
W(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE, give.names = FALSE, ...)

# Methods for compatibility with dplyr:

\method{B}{grouped_df}(x, w = NULL, na.rm = TRUE, fill = FALSE, give.names = FALSE, drop.groups = FALSE, drop.w = TRUE, ...)
W(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE, give.names = FALSE, drop.groups = FALSE, drop.w = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (pseries, see \code{\link{plm}}), panel-data.frame (pdata.frame, see \code{\link{plm}}) or grouped tibble (grouped_df, see \code{\link{dplyr}}).}
  \item{g}{a factor, \code{GRP()} object, atomic vector (> internally converted to factor) or a list of vectors / factors (> internally converted to a \code{GRP()} object) used to group \code{x}.}
  \item{w}{a numeric vector of (non-negative) weights. if \code{na.rm = TRUE}, the weight vector can have missing values, and those observations will be skipped in the computation of averages.}
  \item{na.rm}{skip missing values in \code{x} when computing averages. If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the average for that group will be      \code{NA}, and all data points belonging to that group will also be \code{NA}.}
  \item{give.names}{\code{TRUE} will rename all between (\code{B()}) transformed columns by adding a prefix "B.", and all within (\code{W()}) transformed columns with a prefix "W.".}
  \item{fill}{option to \code{B}: \code{TRUE} will overrwite missing values in \code{x} with the respective average. By default missing values in \code{x} are preserved.}
  \item{add.global.mean}{option to \code{W}: \code{TRUE} will add back the global mean to all data values after subtracting out group-means.}
  \item{drop.groups}{for grouped tibbles: should grouping columns be dropped from the output?}
  \item{drop.w}{for grouped tibbles: should a column in \code{x} used for weighting be dropped from the output?}
}
\details{
Without groups, \code{B} will replace all data points in \code{x} with their mean or weighted mean (if \code{w} is supplied). Similarly \code{W} will subtract the mean from all data points i.e. center the data on the mean. \cr

With groups supplied to \code{g}, the replacment / centering performed by \code{B} / \code{W} becomes groupwise. I like to think of this in terms of panel data: If \code{x} is a vector in such a dataset, \code{xit} denotes a single data-point belonging to group \code{i} in time-period \code{t} (\code{t} need not be a time-period). Then \code{xi.} denotes \code{x}, averaged over \code{t}. \code{B} now returns \code{xi.} and \code{W} returns \code{x - xi.}. Thus for any data \code{x} and any grouping vector \code{g}: \code{B(x,g) + W(x,g) = xi. + x - xi. = x}. In terms of variance, \code{B} only retains the variance between group averages, while \code{W}, by subtracting out group means, only retains the variance within those groups. \cr

The data replacement performed by \code{B} can keep (default) or overwrite missing values (option \code{fill}) in \code{x}. \code{W} can center data simply (default), or add back the global / overall mean in groupwise computations (option \code{add.global.mean}). Let \code{x..} denote the global mean of \code{x}, then \code{W} with \code{add.global.mean = TRUE} returns \code{x - xi. + x..} instead of \code{x - xi.}. This is useful to get rid of group-differences but preserve the overall level of the data (as simple groupwise centering will set the overall mean of the data to 0). In regression analysis, centering with \code{add.global.mean = TRUE} will only change the constant term. See Examples.
}
\value{
\code{B} returns \code{x} with every element replaced by its (groupwise) mean (\code{xi.}). \code{W} returns \code{x} where every element was subtracted its (groupwise) mean (\code{x - xi.}). See Details. 
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
