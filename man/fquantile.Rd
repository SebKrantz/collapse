\name{fquantile}
\alias{fquantile}
\alias{.quantile}
\alias{frange}
\alias{.range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Weighted) Sample Quantiles and Range}
\description{
A faster alternative to \code{\link{quantile}} (written fully in C), that supports sampling weights, and can also quickly compute quantiles from an ordering vector (e.g. \code{order(x)}). \code{frange} provides a fast alternative to \code{\link{range}}.
}
\usage{
fquantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), w = NULL,
          o = if(length(x) > 1e5L && length(probs) > log(length(x)))
              radixorder(x) else NULL,
          na.rm = .op[["na.rm"]], type = 7L, names = TRUE,
          check.o = is.null(attr(o, "sorted")))

# Programmers version: no names, intelligent defaults, or checks
.quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), w = NULL, o = NULL,
          na.rm = TRUE, type = 7L, names = FALSE, check.o = FALSE)

# Fast range (min and max)
frange(x, na.rm = .op[["na.rm"]], finite = FALSE)
.range(x, na.rm = TRUE, finite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric or integer vector.}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
  \item{w}{a numeric vector of strictly positive sampling weights. Missing weights are only supported if \code{x} is also missing.}
  \item{o}{integer. An vector giving the ordering of the elements in \code{x}, such that \code{identical(x[o], sort(x))}. If available this considerably speeds up the estimation.}
  \item{na.rm}{logical. Remove missing values, default \code{TRUE}. }
  \item{finite}{logical. Omit all non-finite values.}
  \item{type}{integer. Quantile types 4-9. See \code{\link{quantile}}. Further details are provided in \href{https://www.tandfonline.com/doi/abs/10.1080/00031305.1996.10473566}{Hyndman and Fan (1996)} who recommended type 8. The default method is type 7.}
  \item{names}{logical. Generates names of the form \code{paste0(round(probs * 100, 1), "\%")} (in C). Set to \code{FALSE} for speedup. }
  \item{check.o}{logical. If \code{o} is supplied, \code{TRUE} runs through \code{o} once and checks that it is valid, i.e. that each element is in \code{[1, length(x)]}. Set to \code{FALSE} for significant speedup if \code{o} is known to be valid. }
}
\details{
\code{fquantile} is implemented using a quickselect algorithm in C, inspired by \emph{data.table}'s \code{gmedian}. The algorithm is applied incrementally to different sections of the array to find individual quantiles. If many quantile probabilities are requested, sorting the whole array with the fast \code{\link{radixorder}} algorithm is more efficient. The default threshold for this (\code{length(x) > 1e5L && length(probs) > log(length(x))}) is conservative, given that quickselect is generally more efficient on longitudinal data with similar values repeated by groups. With random data, my investigations yield that a threshold of \code{length(probs) > log10(length(x))} would be more appropriate.

\code{frange} is considerably more efficient than \code{\link{range}}, requiring only one pass through the data instead of two. For probabilities 0 and 1, \code{fquantile} internally calls \code{frange}.

Following \href{https://www.tandfonline.com/doi/abs/10.1080/00031305.1996.10473566}{Hyndman and Fan (1996)}, the quantile type-\eqn{i} quantile function of the sample \eqn{X} can be written as a weighted average of two order statistics:

\deqn{\hat{Q}_{X,i}(p) = (1 - \gamma) X_{(j)} + \gamma X_{(j + 1)}}

where \eqn{j = \lfloor pn + m \rfloor,\ m \in \mathbb{R}} and \eqn{\gamma = pn + m - j,\ 0 \le \gamma \le 1}, with \eqn{m} differing by quantile type (\eqn{i}). For example, the default type 7 quantile estimator uses \eqn{m = 1 - p}, see \code{\link{quantile}}.

For weighted data with normalized weights \eqn{w = \{w_1, ..., w_n\}}, where \eqn{w_k > 0} and \eqn{\sum_k w_k = 1}, let \eqn{\{w_{(1)}, ..., w_{(n)}\}} be the weights for each order statistic and \eqn{W_{(k)} = \operatorname{Weight}[X_j \le X_{(k)}] = \sum_{j=1}^k w_{(j)}} the cumulative weight for each order statistic.

We can then first find the largest value \eqn{l} such that the cumulative normalized weight \eqn{W_{(l)} \leq p}, and replace \eqn{pn} with \eqn{l + (p - W_{(l)})/w_{(l+1)}}, where \eqn{w_{(l+1)}} is the weight of the next observation. This gives:

\deqn{j = \lfloor l + \frac{p - W_{(l)}}{w_{(l+1)}} + m \rfloor}
\deqn{\gamma = l + \frac{p - W_{(l)}}{w_{(l+1)}} + m - j}

For a more detailed exposition \href{https://htmlpreview.github.io/?https://github.com/mjskay/uncertainty-examples/blob/master/weighted-quantiles.html}{see these excellent notes} by Matthew Kay. See also the R implementation of weighted quantiles type 7 in the Examples below.
}
\note{
The new weighted quantile algorithm from v2.1.0 does not skip zero weights anymore as this is technically very difficult (it is not clear if \eqn{j} hits a zero weight element whether one should move forward or backward to find an alternative). Thus, all non-missing elements are considered and weights should be strictly positive.
}
\value{
A vector of quantiles. If \code{names = TRUE}, \code{fquantile} generates names as \code{paste0(round(probs * 100, 1), "\%")} (in C).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\author{
Sebastian Krantz based on \href{https://htmlpreview.github.io/?https://github.com/mjskay/uncertainty-examples/blob/master/weighted-quantiles.html}{notes} by Matthew Kay.
}
\references{
Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, \emph{American Statistician} 50, 361–365. doi:10.2307/2684934.

Wicklin, R. (2017) Sample quantiles: A comparison of 9 definitions; SAS Blog. https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.html

Wikipedia: https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample

Weighted Quantiles by Matthew Kay: https://htmlpreview.github.io/?https://github.com/mjskay/uncertainty-examples/blob/master/weighted-quantiles.html
}
\seealso{
\code{\link{fnth}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Basic range and quantiles
frange(mtcars$mpg)
fquantile(mtcars$mpg)

## Checking computational equivalence to stats::quantile()
w = alloc(abs(rnorm(1)), 32)
o = radixorder(mtcars$mpg)
for (i in 5:9) print(all_obj_equal(fquantile(mtcars$mpg, type = i),
                                   fquantile(mtcars$mpg, type = i, w = w),
                                   fquantile(mtcars$mpg, type = i, o = o),
                                   fquantile(mtcars$mpg, type = i, w = w, o = o),
                                    quantile(mtcars$mpg, type = i)))

## Demonstaration: weighted quantiles type 7 in R
wquantile7R <- function(x, w, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = TRUE) {
  if(na.rm && anyNA(x)) {             # Removing missing values (only in x)
    cc = whichNA(x, invert = TRUE)    # The C code first calls radixorder(x), which places
    x = x[cc]; w = w[cc]              # missing values last, so removing = early termination
  }
  o = radixorder(x)                   # Ordering
  wo = proportions(w[o])
  Wo = cumsum(wo)                     # Cumulative sum
  res = sapply(probs, function(p) {
    l = which.max(Wo > p) - 1L        # Lower order statistic
    s = l + (p - Wo[l])/wo[l+1L] + 1 - p
    j = floor(s)
    gamma = s - j
    (1 - gamma) * x[o[j]] + gamma * x[o[j+1L]]  # Weighted quantile
  })
  if(names) names(res) = paste0(as.integer(probs * 100), "\%")
  res
} # Note: doesn't work for min and max.

wquantile7R(mtcars$mpg, mtcars$wt)

all.equal(wquantile7R(mtcars$mpg, mtcars$wt),
          fquantile(mtcars$mpg, c(0.25, 0.5, 0.75), mtcars$wt))

## Efficient grouped quantile estimation: use .quantile for less call overhead
BY(mtcars$mpg, mtcars$cyl, .quantile, names = TRUE, expand.wide = TRUE)
BY(mtcars, mtcars$cyl, .quantile, names = TRUE)
mtcars |> fgroup_by(cyl) |> BY(.quantile)

## With weights
BY(mtcars$mpg, mtcars$cyl, .quantile, w = mtcars$wt, names = TRUE, expand.wide = TRUE)
BY(mtcars, mtcars$cyl, .quantile, w = mtcars$wt, names = TRUE)
mtcars |> fgroup_by(cyl) |> fselect(-wt) |> BY(.quantile, w = mtcars$wt)
mtcars |> fgroup_by(cyl) |> fsummarise(across(-wt, .quantile, w = wt))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{univar}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
