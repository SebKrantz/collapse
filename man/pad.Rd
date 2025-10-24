\name{pad}
\alias{pad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pad Matrix-Like Objects with a Value
}
\description{
The \code{pad} function inserts elements / rows filled with \code{value} into a vector matrix or data frame \code{X} at positions given by \code{i}. It is particularly useful to expand objects returned by statistical procedures which remove missing values to the original data dimensions.
}
\usage{
pad(X, i, value = NA, method = c("auto", "xpos", "vpos"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  a vector, matrix, data frame or list of equal-length columns.
}
  \item{i}{
  either an integer (positive or negative) or logical vector giving positions / rows of \code{X} into which \code{value}'s should be inserted, or, alternatively, a positive integer vector with \code{length(i) == NROW(X)}, but with some gaps in the indices into which \code{value}'s can be inserted, or a logical vector with \code{sum(i) == NROW(X)} such that \code{value}'s can be inserted for \code{FALSE} values in the logical vector. See also \code{method} and Examples.
}
  \item{value}{
 a scalar value to be replicated and inserted into \code{X} at positions / rows given by \code{i}. Default is \code{NA}.
}
  \item{method}{
 an integer or string specifying the use of \code{i}. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto" \tab\tab automatic method selection: If \code{i} is positive integer and \code{length(i) == NROW(X)} or if \code{i} is logical and \code{sum(i) == NROW(X)}, choose method "xpos", else choose "vpos". \cr \tab\tab\tab\tab \cr \tab\tab\tab\tab \cr
                 1 \tab\tab "xpos"   \tab\tab \code{i} is a vector of positive integers or a logical vector giving the positions of the the elements / rows of \code{X}. \code{values}'s are inserted where there are gaps / \code{FALSE} values in \code{i}. \cr \tab\tab\tab\tab \cr \tab\tab\tab\tab \cr
                 2 \tab\tab "vpos" \tab\tab \code{i} is a vector of positive / negative integers or a logical vector giving the positions at which \code{values}'s / rows should be inserted into \code{X}.
  }
}
}
\value{
\code{X} with elements / rows filled with \code{value} inserted at positions given by \code{i}.
}
\seealso{
\code{\link{append}}, \link[=recode-replace]{Recode and Replace Values}, \link[=small-helpers]{Small (Helper) Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
v <- 1:3

pad(v, 1:2)       # Automatic selection of method "vpos"
pad(v, -(1:2))    # Same thing
pad(v, c(TRUE, TRUE, FALSE, FALSE, FALSE)) # Same thing

pad(v, c(1, 3:4)) # Automatic selection of method "xpos"
pad(v, c(TRUE, FALSE, TRUE, TRUE, FALSE))  # Same thing

head(pad(wlddev, 1:3)) # Insert 3 missing rows at the beginning of the data
head(pad(wlddev, 2:4)) # ... at rows positions 2-4

# pad() is mostly useful for statistical models which only use the complete cases:
mod <- lm(LIFEEX ~ PCGDP, wlddev)
# Generating a residual column in the original data (automatic selection of method "vpos")
settfm(wlddev, resid = pad(resid(mod), mod$na.action))
# Another way to do it:
r <- resid(mod)
i <- as.integer(names(r))
resid2 <- pad(r, i)        # automatic selection of method "xpos"
# here we need to add some elements as flast(i) < nrow(wlddev)
resid2 <- c(resid2, rep(NA, nrow(wlddev)-length(resid2)))
# See that these are identical:
identical(unattrib(wlddev$resid), resid2)

# Can also easily get a model matrix at the dimensions of the original data
mm <- pad(model.matrix(mod), mod$na.action)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
