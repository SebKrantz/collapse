\name{pwcor-pwcov-pwnobs}
\alias{pwcor}
\alias{pwcov}
\alias{pwnobs}
\alias{print.pwcov}
\alias{print.pwcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
(Pairwise, Weighted) Correlations, Covariances and Observation Counts
}
\description{
Computes (pairwise, weighted) Pearson's correlations, covariances and observation counts. Pairwise correlations and covariances can be computed together with observation counts and p-values, and output as 3D array (default) or list of matrices. \code{pwcor} and \code{pwcov} offer an elaborate print method.
}
\usage{
pwcor(X, \dots, w = NULL, N = FALSE, P = FALSE, array = TRUE, use = "pairwise.complete.obs")

pwcov(X, \dots, w = NULL, N = FALSE, P = FALSE, array = TRUE, use = "pairwise.complete.obs")

pwnobs(X)

\method{print}{pwcor}(x, digits = .op[["digits"]], sig.level = 0.05,
      show = c("all","lower.tri","upper.tri"), spacing = 1L, return = FALSE, \dots)

\method{print}{pwcov}(x, digits = .op[["digits"]], sig.level = 0.05,
      show = c("all","lower.tri","upper.tri"), spacing = 1L, return = FALSE, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a matrix or data.frame, for \code{pwcor} and \code{pwcov} all columns must be numeric. All functions are faster on matrices, so converting is advised for large data (see \code{\link{qM}}).}
  \item{x}{an object of class 'pwcor' / 'pwcov'. }
  \item{w}{numeric. A vector of (frequency) weights. }
  \item{N}{logical. \code{TRUE} also computes pairwise observation counts.}
  \item{P}{logical. \code{TRUE} also computes pairwise p-values (same as \code{\link{cor.test}} and \code{Hmisc::rcorr}).}
  \item{array}{logical. If \code{N = TRUE} or \code{P = TRUE}, \code{TRUE} (default) returns output as 3D array whereas \code{FALSE} returns a list of matrices.}
  \item{use}{argument passed to \code{\link{cor}} / \code{\link{cov}}. If \code{use != "pairwise.complete.obs"}, \code{sum(complete.cases(X))} is used for \code{N}, and p-values are computed accordingly. }
  \item{digits}{integer. The number of digits to round to in print. }
  \item{sig.level}{numeric. P-value threshold below which a \code{'*'} is displayed above significant coefficients if \code{P = TRUE}. }
  \item{show}{character. The part of the correlation / covariance matrix to display. }
  \item{spacing}{integer. Controls the spacing between different reported quantities in the printout of the matrix: 0 - compressed, 1 - single space, 2 - double space.}
  \item{return}{logical. \code{TRUE} returns the formatted object from the print method for exporting. The default is to return \code{x} invisibly.}
  \item{\dots}{other arguments passed to \code{\link{cor}} or \code{\link{cov}}. Only sensible if \code{P = FALSE}. }
}

\value{
a numeric matrix, 3D array or list of matrices with the computed statistics. For \code{pwcor} and \code{pwcov} the object has a class 'pwcor' and 'pwcov', respectively.
}

\note{
\code{weights::wtd.cors} is imported for weighted pairwise correlations (written in C for speed). For weighted correlations with bootstrap SE's see \code{weights::wtd.cor} (bootstrap can be slow). Weighted correlations for complex surveys are implemented in \code{jtools::svycor}. An equivalent and faster implementation of \code{pwcor} (without weights) is provided in \code{Hmisc::rcorr} (written in Fortran).
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{qsu}}, \link[=summary-statistics]{Summary Statistics}, \link[=collapse-documentation]{Collapse Overview} %% ~~objects to See Also as
}
\examples{
mna <- na_insert(mtcars)
pwcor(mna)
pwcov(mna)
pwnobs(mna)
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

