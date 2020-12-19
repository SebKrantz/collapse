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
Fast operators to perform row- or column-wise replacing and sweeping operations of vectors on matrices, data frames, lists.
}
\usage{
## Perform the operation with v and each row of X

X \%rr\% v    # Replace rows of X with v
X \%r+\% v    # Add v to each row of X
X \%r-\% v    # Subtract v from each row of X
X \%r*\% v    # Multiply each row of X with v
X \%r/\% v    # Divide each row of X by v

## Perform the operation with v and each column of X

X \%cr\% v    # Replace columns of X with v
X \%c+\% v    # Add v to each column of X
X \%c-\% v    # Subtract v from each column of X
X \%c*\% v    # Multiply each column of X with v
X \%c/\% v    # Divide each column of X by v
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix, data frame or list like object with rows (r) ro columns (c) matching \code{v}.}
  \item{v}{a suitable atomic vector. If \code{X} is a list, \code{v} can also be a scalar (for column-operations) or a matching list of atomic elements (for row-operations on data frame's / lists).}
}
\details{
With a matrix or data frame \code{X}, the default behavior of R when calling \code{X op v} (such as multiplication \code{X * v}) is to perform the operation of the vector \code{v} with each column of \code{X}. The equivalent operation is performed by \code{X \%cop\% v}, with the difference that it computes significantly faster if \code{X} is a data frame / list. A more complex but frequently required task is to perform an operation with \code{v} on each row of \code{X}. This is provided based on efficient C++ code by the \code{\%rop\%} set of functions, e.g. \code{X \%r*\% v} efficiently multiplies \code{v} to each row of \code{X}.
}
\value{
\code{X} where the operation with \code{v} was performed on each row or column. All attributes of \code{X} are preserved.
}
\note{
These functions are all quite simple, they only work with \code{X} on the LHS i.e. \code{v \%op\% X} will likely fail. The row operations are simple wrappers around \code{\link{TRA}} which provides more operations including grouped replacing and sweeping (where \code{V} would be a matrix or data frame with less rows than \code{X} being mapped to the rows of \code{X} by grouping vectors). One consequence is that just like \code{\link{TRA}}, row-wise mathematical operations (+, -, *, /) always yield numeric output, even if both \code{X} and \code{v} may be integer. This is different for column- operations which depend on base R and may also preserve integer data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{TRA}}, \code{\link{dapply}}, \link[=A6-data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}

\examples{
# Using data frame's / lists
v <- mtcars$cyl
mtcars \%cr\% v
mtcars \%c-\% v
mtcars \%r-\% lapply(mtcars, quantile, 0.28)

# Using matrices
X <- qM(mtcars)
X \%cr\% v
X \%c-\% v
X \%r-\% dapply(X, quantile, 0.28)

# General Notes: For Fast Statistical Functions, using
# fmedian(X, TRA = "-") will be faster than X \%r-\% fmedian(X).
# Also use fwithin(X) for fast centering using the mean, and
# fscale(X) for fast scaling and centering or mean-preserving scaling

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
