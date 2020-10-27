\name{flm}
\alias{flm}

\title{
Fast (Sparse, Weighted) Linear Model Fitting
}
\description{
\code{flm} is a fast linear model command that takes matrices as input and (by default) only returns a coefficient matrix. 6 different efficient fitting methods are implemented: 4 using base R linear algebra, and 2 utilizing the \emph{RcppArmadillo} and \emph{RcppEigen} packages. The function itself only has an overhead of 5-10 microseconds, and is thus well suited as a bootstrap workhorse.
}
\usage{
flm(y, x, w = NULL, add.icpt = FALSE,
    sparse = FALSE, return.raw = FALSE,
    method = c("lm", "solve", "qr", "arma", "chol", "eigen"),
    eigen.method = 3L, tol = 1e-07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a response vector or matrix. Multiple dependent variables are only supported by methods "lm", "solve", "qr" and "chol".}
  \item{x}{a matrix or regressors.}
  \item{w}{a weight vector.}
  \item{add.icpt}{logical. \code{TRUE} adds an intercept column named '(Intercept)' to \code{x}.}
  \item{sparse}{logical. \code{TRUE} coerces \code{x} to a sparse matrix using \code{as(x, "dgCMatrix")}.}
  \item{return.raw}{logical. \code{TRUE} returns the original output from the different methods. For 'lm', 'arma' and 'eigen', this includes additional statistics such as residuals, fitted values or standard errors. The other methods just return coefficients but in different formats. }
  \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "lm"   \tab\tab uses \code{\link{.lm.fit}}. \cr
                 2 \tab\tab "solve" \tab\tab \code{solve(crossprod(x), crossprod(x, y))}. \cr
                 3 \tab\tab "qr"   \tab\tab \code{qr.coef(qr(x, LAPACK = TRUE), y)}. \cr
                 4 \tab\tab "arma"   \tab\tab uses \code{RcppArmadillo::fastLmPure} (directly called from C). \cr
                 5 \tab\tab "chol"   \tab\tab \code{chol2inv(chol(crossprod(x))) \%*\% crossprod(x, y)} (quite fast but requires \code{crossprod(x)} to be positive definite).  \cr
                 6 \tab\tab "eigen"   \tab\tab uses \code{RcppEigen::fastLmPure} (directly called from C, very fast but potentially unstable, depending on the method). \cr
  }
}
\item{eigen.method}{integer. Select the method of computated used by \code{RcppEigen::fastLmPure}:
        \tabular{lll}{\emph{ Int. } \tab\tab \emph{ Description }  \cr
                 0 \tab\tab column-pivoted QR decomposition. \cr
                 1 \tab\tab unpivoted QR decomposition. \cr
                 2 \tab\tab LLT Cholesky. \cr
                 3 \tab\tab LDLT Cholesky. \cr
                 4 \tab\tab Jacobi singular value decomposition (SVD).  \cr
                 5 \tab\tab method based on the eigenvalue-eigenvector decomposition of X'X. \cr
  }
  See \code{vignette("RcppEigen-Introduction", package = "RcppEigen")} for details on these methods and benchmark results. Run \code{source(system.file("examples", "lmBenchmark.R", package = "RcppEigen"))} to re-run the benchmark on your machine.
}
\item{tol}{numeric. A numerical tolerance for the solution - only applicable with methods "lm", "solve" and "qr".}
}

% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
\value{
If \code{return.raw = FALSE}, a matrix of coefficients with the rows corresponding to the columns of \code{x}, otherwise the raw results from the various methods are returned.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=HDW]{fHDwithin/HDW}}, \code{\link{fFtest}}, \link[=A6-data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
coef <- flm(mtcars$mpg, qM(mtcars[c("hp","carb")]),
            mtcars$wt, add.icpt = TRUE)
coef

lmcoef <- coef(lm(mpg ~ hp + carb, weights = wt, mtcars))
lmcoef

all.equal(drop(coef), lmcoef)

all_obj_equal(lapply(1:6, function(i)
 flm(mtcars$mpg, qM(mtcars[c("hp","carb")]),
     mtcars$wt, add.icpt = TRUE, method = i)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
