\name{fquantile}
\alias{fquantile}
\alias{frange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Weighted) Sample Quantiles and Range of Vector}
\description{
A faster alternative to \code{\link{quantile}} (written fully in C), that supports sampling weights, and can also quickly compute quantiles from an ordering vector (e.g. \code{order(x)}). \code{frange} provides a fast alternative to \code{\link{range}}.
}
\usage{
fquantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), w = NULL, o = NULL,
          na.rm = TRUE, type = 7L, names = TRUE, check.o = TRUE)

frange(x, na.rm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric or integer vector.}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
  \item{w}{a numeric vector of sampling weights. Missing weights are only supported if \code{x} is also missing.}
  \item{o}{integer. An vector giving the ordering of the elements in \code{x}, such that \code{identical(x[o], sort(x))}. If available this considerably speeds up the estimation.}
  \item{na.rm}{logical. Remove missing values, default \code{TRUE}. }
  \item{type}{integer. Quantile types 5-9. See \code{\link{quantile}}. }
  \item{names}{logical. Generates names of the form \code{paste0(as.integer(probs * 100), "\%")} (in C). Set to \code{FALSE} for speedup, especially during repeated executions. }
  \item{check.o}{logical. If \code{o} is supplied, \code{TRUE} runs through \code{o} once and checks that it is valid, i.e. that each element is in \code{[1, length(x)]}. Set to \code{FALSE} for significant speedup if \code{o} is known to be valid. }
}
\details{
\code{fquantile} is implemented using a quickselect algorithm in C, inspired by \emph{data.table}'s \code{gmedian}.

Weighted quantile estimation, in a nutshell, is done by internally calling \code{radixorder(x)} (unless \code{o} is supplied), summing the weights in order until the highest required order statistic \code{j+1} is found, which corresponds to reaching a target sum of weights that is a function of the probability \code{p}, the quantile method (see \code{\link{quantile}}), the total sum of weights, and the average weight. Then, a continuous index \code{h} in [0, 1] is determined from the difference between the sum of weights associated with \code{j+1} and the target sum, divided by the weight of element \code{j+1}. A quantile can then be computed as a weighted average of 2 order statistics: \code{Q[i](p) = (1 - h) x[j] + h x[j+1]} (as in the unweighted case). If the lower order statistic \code{x[j]} has a zero weight, \code{x[j-1]} is taken (or \code{x[j-2]} if \code{x[j-1]} also has zero weight etc..). Taking into account the average weight in the determination of the target sum (and therefore \code{h}) ensures that weighted quantile estimation gives the same results as unweighted estimation if the weights are all equal.

\code{frange} is considerably more efficient than \code{\link{range}}, which calls both \code{\link{min}} and \code{\link{max}}, and thus requires 2 full passes instead of 1 required by \code{frange}.
}
\value{
A vector of quantiles. If \code{names = TRUE}, \code{fquantile} generates names as \code{paste0(as.integer(probs * 100), "\%")} (in C).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
frange(mtcars$mpg)

w = rep(1, 32)
o = radixorder(mtcars$mpg)
for (i in 5:9) print(all_obj_equal(fquantile(mtcars$mpg, type = i),
                                   fquantile(mtcars$mpg, type = i, w = w),
                                   fquantile(mtcars$mpg, type = i, o = o),
                                   fquantile(mtcars$mpg, type = i, w = w, o = o),
                                    quantile(mtcars$mpg, type = i)))

# Efficient grouped quantile estimation
BY(mtcars$mpg, mtcars$cyl, fquantile, expand.wide = TRUE)
BY(mtcars, mtcars$cyl, fquantile)

# With weights
BY(mtcars$mpg, mtcars$cyl, fquantile, w = mtcars$wt, expand.wide = TRUE)
BY(mtcars, mtcars$cyl, fquantile, w = mtcars$wt)

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
