\name{pwcor, pwcov, pwNobs}
\alias{pwcor}
\alias{pwcov}
\alias{pwNobs}
\alias{print.pwcov}
\alias{print.pwcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Correlations, Covariances and Observation Count
}
\description{
Computes pairwise Pearsons correlations, covariances and observation counts. Pairwise correlations and covariances can be computed together with observation counts and p-values, and output as 3D array (default) or list of matrices. For an equivalent and faster implementation of \code{pwcor} see \code{Hmisc::rcorr} (written in Fortran). A major feature \code{pwcor} and \code{pwcov} is their sophisticated print method.
}
\usage{
pwcor(X, ..., N = FALSE, P = FALSE, array = TRUE)

pwcov(X, ..., N = FALSE, P = FALSE, array = TRUE)

pwNobs(X)

\method{print}{pwcor}(x, digits = 2L, sig.level = 0.05, show = c("all","lower.tri","upper.tri"),
      spacing = 1L, ...)

\method{print}{pwcov}(x, digits = 2L, sig.level = 0.05, show = c("all","lower.tri","upper.tri"),
      spacing = 1L, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or data.frame, for \code{pwcor} and \code{pwcov} all columns must be numeric.}
  \item{x}{an object of class 'pwcor' / 'pwcov'. }
  \item{N}{logical. \code{TRUE} also computes pairwise observation counts.}
  \item{P}{logical. \code{TRUE} also computes pairwise p-values (same as \code{\link{cor.test}}).}
  \item{array}{logical. If \code{N = TRUE} or \code{P = TRUE}, \code{TRUE} (default) returns output as 3D array whereas \code{FALSE} returns a list of matrices.}
  \item{digits}{integer. The number of digits to round to in print. }
  \item{sig.level}{numeric. P-value threshold below which a \code{'*'} is displayed above significant coefficients if \code{P = TRUE}. }
  \item{show}{character. The part of the correlation / covariance matrix to display. }
  \item{spacing}{integer. Controls the spacing between different reported quantities in the printout of the matrix: 0 - compressed, 1 - single space, 2 - double space.}
  \item{...}{other arguments passed to \code{\link{cor}} or \code{\link{cov}}. Only sensible if \code{P = FALSE}. }
}

\value{
a numeric matrix, 3D array or list of matrices of the computed statistics. For \code{pwcor} and \code{pwcov} the object has a class 'pwcor' and 'pwcov', respectively.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{qsu}}, \link[=collapse-documentation]{Collapse Overview} %% ~~objects to See Also as
}
\examples{
mna <- na_insert(mtcars)
pwcor(mna)
pwcov(mna)
pwNobs(mna)
pwcor(mna, N = TRUE)
pwcor(mna, P = TRUE)
pwcor(mna, N = TRUE, P = TRUE)
aperm(pwcor(mna, N = TRUE, P = TRUE))
print(pwcor(mna, N = TRUE, P = TRUE), digits = 3, sig.level = 0.01, show = "lower.tri")
pwcor(mna, N = TRUE, P = TRUE, array = FALSE)
print(pwcor(mna, N = TRUE, P = TRUE, array = FALSE), show = "lower.tri")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate} % use one of  RShowDoc("KEYWORDS")

