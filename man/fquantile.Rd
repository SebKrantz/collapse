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
  \item{w}{a numeric vector of sampling weights. Missing weights are only supported if \code{x} is also missing.}
  \item{o}{integer. An vector giving the ordering of the elements in \code{x}, such that \code{identical(x[o], sort(x))}. If available this considerably speeds up the estimation.}
  \item{na.rm}{logical. Remove missing values, default \code{TRUE}. }
  \item{finite}{logical. Omit all non-finite values.}
  \item{type}{integer. Quantile types 5-9. See \code{\link{quantile}}. Further details are provided in Hyndman and Fan (1996) who recommended type 8. The default method is type 7.}
  \item{names}{logical. Generates names of the form \code{paste0(round(probs * 100, 1), "\%")} (in C). Set to \code{FALSE} for speedup. }
  \item{check.o}{logical. If \code{o} is supplied, \code{TRUE} runs through \code{o} once and checks that it is valid, i.e. that each element is in \code{[1, length(x)]}. Set to \code{FALSE} for significant speedup if \code{o} is known to be valid. }
}
\details{
\code{fquantile} is implemented using a quickselect algorithm in C, inspired by \emph{data.table}'s \code{gmedian}. The algorithm is applied incrementally to different sections of the array to find individual quantiles. If many quantile probabilities are requested, sorting the whole array with the fast \code{\link{radixorder}} algorithm is more efficient. The default threshold for this (\code{length(x) > 1e5L && length(probs) > log(length(x))}) is conservative, given that quickselect is generally more efficient on longitudinal data with similar values repeated by groups. With random data, my investigations yield that a threshold of \code{length(probs) > log10(length(x))} would be more appropriate.

Weighted quantile estimation, in a nutshell, is done by internally calling \code{radixorder(x)} (unless \code{o} is supplied), and summing the weights in order until the lowest required order statistic \code{j} is found, which corresponds to exceeding a target sum of weights that is a function of the probability \code{p}, the quantile method (see \code{\link{quantile}}), the total sum of weights, and the smallest (non-zero) weight. For quantile type 7 the target sum is \code{sumwp = (sum(w) - min(w)) * p} (resembling \code{(n - 1) * p} in the unweighted case). Then, a continuous index \code{h} in [0, 1] is determined as one minus the difference between the sum of weights associated with \code{j} and the target sum, divided by the weight of element \code{j}, that is \code{h = 1 - (sumwj - sumwp) / w[j]}. A weighted quantile can then be computed as a weighted average of 2 order statistics, exactly as in the unweighted case: \code{WQ[i](p) = (1 - h) x[j] + h x[j+1]}. If the order statistic \code{j+1} has a zero weight, \code{j+2} is taken (or \code{j+3} if \code{j+2} also has zero weight etc..). The Examples section provides a demonstration in R that is roughly equivalent to the algorithm just outlined.

\code{frange} is considerably more efficient than \code{\link{range}}, which calls both \code{\link{min}} and \code{\link{max}}, and thus requires 2 full passes instead of 1 required by \code{frange}. If only probabilities \code{0} and \code{1} are requested, \code{fquantile} internally calls \code{frange}.

}
\value{
A vector of quantiles. If \code{names = TRUE}, \code{fquantile} generates names as \code{paste0(round(probs * 100, 1), "\%")} (in C).
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fnth}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
frange(mtcars$mpg)

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
  if(anyv(w, 0)) {                    # Removing zero weights
    nzw = whichv(w, 0, invert = TRUE) # In C, skipping zero weight order statistics is built
    x = x[nzw]; w = w[nzw]            # into the quantile algorithm, as outlined above
  }
  o = radixorder(x)                   # Ordering
  wo = w[o]
  w_cs = cumsum(wo)                   # Cumulative sum
  sumwp = sum(w)                      # Computing sum(w) - min(w)
  sumwp = sumwp - min(w)
  sumwp = sumwp * probs               # Target sums of weights for quantile type 7
  res = sapply(sumwp, function(tsump) {
    j = which.max(w_cs > tsump)           # Lower order statistic
    hl = (w_cs[j] - tsump) / wo[j]        # Index weight of x[j]  (h = 1 - hl)
    hl * x[o[j]] + (1 - hl) * x[o[j+1L]]  # Weighted quantile
  })
  if(names) names(res) = paste0(as.integer(probs * 100), "\%")
  res
} # Note: doesn't work for min and max. Overall the C code is significantly more rigorous.

wquantile7R(mtcars$mpg, mtcars$wt)

all.equal(wquantile7R(mtcars$mpg, mtcars$wt),
          fquantile(mtcars$mpg, c(0.25, 0.5, 0.75), mtcars$wt))

## Efficient grouped quantile estimation: use .quantile for less call overhead
BY(mtcars$mpg, mtcars$cyl, .quantile, names = TRUE, expand.wide = TRUE)
BY(mtcars, mtcars$cyl, .quantile, names = TRUE)
library(magrittr)
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
