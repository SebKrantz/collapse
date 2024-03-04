\name{arithmetic}
\alias{arithmetic}
\alias{\%rr\%}
\alias{\%r+\%}
\alias{\%r-\%}
\alias{\%r*\%}
\alias{\%r/\%}
\alias{\%cr\%}
\alias{\%c+\%}
\alias{\%c-\%}
\alias{\%c*\%}
\alias{\%c/\%}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Row/Column Arithmetic for Matrix-Like Objects
}
\description{
Fast operators to perform row- or column-wise replacing and sweeping operations of vectors on matrices, data frames, lists. See also \code{\link{setop}} for math by reference and \code{\link{setTRA}} for sweeping by reference.
}
\usage{
## Perform the operation with v and each row of X

X \%rr\% v    # Replace rows of X with v
X \%r+\% v    # Add v to each row of X
X \%r-\% v    # Subtract v from each row of X
X \%r*\% v    # Multiply each row of X with v
X \%r/\% v    # Divide each row of X by v

## Perform a column-wise operation between V and X

X \%cr\% V    # Replace columns of X with V
X \%c+\% V    # Add V to columns of X
X \%c-\% V    # Subtract V from columns of X
X \%c*\% V    # Multiply columns of X with V
X \%c/\% V    # Divide columns of X by V
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix, data frame or list like object (with rows (r) columns (c) matching \code{v} / \code{V}).}
  \item{v}{for row operations: an atomic vector of matching \code{NCOL(X)}. If \code{X} is a data frame, \code{v} can also be a list of scalar atomic elements. It is also possible to sweep lists of vectors \code{v} out of lists of matrices or data frames \code{X}.}

  \item{V}{for column operations: a suitable scalar, vector, or matrix / data frame matching \code{NROW(X)}. \code{X} can also be a list of vectors / matrices in which case \code{V} can be a scalar / vector / matrix or matching list of scalars / vectors / matrices.}

}
\details{
With a matrix or data frame \code{X}, the default behavior of R when calling \code{X op v} (such as multiplication \code{X * v}) is to perform the operation of \code{v} with each column of \code{X}. The equivalent operation is performed by \code{X \%cop\% V}, with the difference that it computes significantly faster if \code{X}/\code{V} is a data frame / list. A more complex but frequently required task is to perform an operation with \code{v} on each row of \code{X}. This is provided based on efficient C++ code by the \code{\%rop\%} set of functions, e.g. \code{X \%r*\% v} efficiently multiplies \code{v} to each row of \code{X}.
}
\value{
\code{X} where the operation with \code{v} / \code{V} was performed on each row or column. All attributes of \code{X} are preserved.
}
\note{
\emph{Computations and Output:} These functions are all quite simple, they only work with \code{X} on the LHS i.e. \code{v \%op\% X} will likely fail. The row operations are simple wrappers around \code{\link{TRA}} which provides more operations including grouped replacing and sweeping (where \code{v} would be a matrix or data frame with less rows than \code{X} being mapped to the rows of \code{X} by grouping vectors). One consequence is that just like \code{\link{TRA}}, row-wise mathematical operations (+, -, *, /) always yield numeric output, even if both \code{X} and \code{v} may be integer. This is different for column- operations which depend on base R and may also preserve integer data.

\emph{Rules of Arithmetic:} Since these operators are defined as simple infix functions, the normal rules of arithmetic are not respected. So \code{a \%c+\% b \%c*\% c} evaluates as \code{(a \%c+\% b) \%c*\% c}. As with all chained infix operations, they are just evaluated sequentially from left to right.

\emph{Performance Notes:} The function \code{\link{setop}} and a related set of \code{\%op=\%} operators as well as the \code{\link{setTRA}} function can be used to perform these operations by reference, and are faster if copies of the output are not required!! Furthermore, for Fast Statistical Functions, using \code{fmedian(X, TRA = "-")} will be a tiny bit faster than \code{X \%r-\% fmedian(X)}. Also use \code{fwithin(X)} for fast centering using the mean, and \code{fscale(X)} for fast scaling and centering or mean-preserving scaling.



}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{setop}}, \code{\link{TRA}}, \code{\link{dapply}}, \link[=efficient-programming]{Efficient Programming}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}

\examples{
## Using data frame's / lists
v <- mtcars$cyl
mtcars \%cr\% v
mtcars \%c-\% v
mtcars \%r-\% seq_col(mtcars)
mtcars \%r-\% lapply(mtcars, quantile, 0.28)

mtcars \%c*\% 5       # Significantly faster than mtcars * 5
mtcars \%c*\% mtcars  # Significantly faster than mtcars * mtcars

## Using matrices
X <- qM(mtcars)
X \%cr\% v
X \%c-\% v
X \%r-\% dapply(X, quantile, 0.28)

## Chained Operations
library(magrittr) # Needed here to evaluate infix operators in sequence
mtcars \%>\% fwithin() \%r-\% rnorm(11) \%c*\% 5 \%>\%
    tfm(mpg = fsum(mpg)) \%>\% qsu()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
\keyword{math}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
