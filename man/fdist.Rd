\name{fdist}
\alias{fdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast and Flexible Distance Computations
}
\description{
A fast and flexible replacement for \code{\link{dist}}, to compute euclidean distances.
}
\usage{
fdist(x, v = NULL, ..., method = "euclidean", nthreads = .op[["nthreads"]])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector or matrix. Data frames/lists can be passed but will be converted to matrix using \code{\link{qM}}. Non-numeric (double) inputs will be coerced. }
  \item{v}{an (optional) numeric (double) vector such that \code{length(v) == NCOL(x)}, to compute distances with (the rows of) \code{x}. Other vector types will be coerced.}
  \item{\dots}{not used. A placeholder for possible future arguments.}
  \item{method}{an integer or character string indicating the method of computing distances.
  \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab \code{"euclidean"}   \tab\tab euclidean distance \cr
                 2 \tab\tab \code{"euclidean_squared"} \tab\tab squared euclidean distance (more efficient) \cr
    }
    %\emph{Note:} The mahalanobis distance can be computed using: \code{x_mahal = t(forwardsolve(t(chol(cov(x))), t(x)))}. See Examples.
  }
  \item{nthreads}{integer. The number of threads to use. If \code{v = NULL} (full distance matrix), multithreading is along the distance matrix columns (decreasing thread loads as matrix is lower triangular). If \code{v} is supplied, multithreading is at the sub-column level (across elements).}
}
\value{
If \code{v = NULL}, a full lower-triangular distance matrix between the rows of \code{x} is computed and returned as a 'dist' object (all methods apply, see \code{\link{dist}}). Otherwise, a numeric vector of distances of each row of \code{x} with \code{v} is returned. See Examples.
}
\note{
\code{fdist} does not check for missing values, so \code{NA}'s will result in \code{NA} distances.

\code{kit::topn} is a suitable complimentary function to find nearest neighbors. It is very efficient and skips missing values by default.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{flm}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# Distance matrix
m = as.matrix(mtcars)
str(fdist(m)) # Same as dist(m)

# Distance with vector
d = fdist(m, fmean(m))
kit::topn(d, 5)  # Index of 5 nearest neighbours

# Mahalanobis distance
m_mahal = t(forwardsolve(t(chol(cov(m))), t(m)))
fdist(m_mahal, fmean(m_mahal))
sqrt(unattrib(mahalanobis(m, fmean(m), cov(m))))
\donttest{
# Distance of two vectors
x <- rnorm(1e6)
y <- rnorm(1e6)
microbenchmark::microbenchmark(
  fdist(x, y),
  fdist(x, y, nthreads = 2),
  sqrt(sum((x-y)^2))
)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{multivariate}
\keyword{nonparametric}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
