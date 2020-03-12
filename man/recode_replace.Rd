\name{AA1-recode-replace}
\alias{AA1-recode-replace}
\alias{Recode}
\alias{replace_non_finite}
\alias{replace_outliers}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recode and Replace Values in Matrix-Like Objects
}
\description{
\itemize{
\item \code{Recode} can be used to replace multiple values in vectors, matrices or data.frames, using either exact (\code{==}) or regular expression matching.
\item \code{replace_non_finite} replaces \code{NaN/Inf/-Inf} (or optionally only \code{Inf/-Inf}) with a value (default is \code{NA}).
\item \code{replace_outliers} replaces values falling outside a 1- or 2-sided numeric threshold or outside a certain number of column- standard deviations with a value (default is \code{NA}).
}
}
\usage{
Recode(X, ..., copy = FALSE, reserve.na.nan = TRUE, regex = FALSE)

replace_non_finite(X, value = NA, replace.nan = TRUE)

replace_outliers(X, limits, value = NA, single.limit = c("SDs","min","max"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix or data.frame.}
  \item{...}{comma-separated recode arguments of the form: \code{name = newname, `2` = 0, `NaN` = 0, `NA` = 0, `Inf` = NA, `-Inf` = NA}, etc...}
  \item{limits}{either a vector of two-numeric values \code{c(minval, maxval)} constituting a two-sided outlier threshold, or a single numeric value constituting either factor of standard deviations (default), or the minimum or maximum of a one-sided outlier threshold. See also \code{single.limit}.}
  \item{value}{a single (scalar) value to replace matching elements with. Default is \code{NA}.}
  \item{copy}{logical. For reciprocal or sequential replacements of the form \code{a = b, b = c} make a copy of \code{X} to prevent \code{a} being replaced with \code{b} and then all \code{b}-values being replaced with \code{c} again. In general \code{Recode} does the replacements one-after the other, starting with the first. }
  \item{reserve.na.nan}{logical. \code{TRUE} identifies \code{NA} and \code{NaN} as special numeric values and does the correct replacement. \code{FALSE} will treat \code{NA/NaN} as strings, and thus not match numeric \code{NA/NaN}. \emph{Note}: This is not an issue for \code{Inf/-Inf}, which are matched in both numeric and character variables. }
  \item{regex}{logical. If \code{TRUE}, all recode-argument names are (sequentially) passed to \code{\link{grepl}} as a pattern to search \code{X}. All matches are replaced.}
  \item{replace.nan}{logical. \code{TRUE} (default) replaces \code{NaN/Inf/-Inf}. \code{FALSE} replaces only \code{Inf/-Inf}.}
  \item{single.limit}{a character or integer (only applies if \code{length(limits) == 1}):
  \itemize{
  \item \code{1 - "SDs"} specifies that \code{limits} will be interpreted as a (two-sided) threshold in column standard-deviations. The underlying code is equivalent to \code{X[abs(fscale(X)) > limits] <- value}.
  \item \code{2 - "min"} specifies that \code{limits} will be interpreted as a (one-sided) minimum threshold. The underlying code is equivalent to \code{X[X < limits] <- value}.
  \item \code{3 - "max"} specifies that \code{limits} will be interpreted as a (one-sided) maximum threshold. The underlying code is equivalent to \code{X[X > limits] <- value}.
  }
}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
% }
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
\note{
\code{Recode} is not suitable for recoding factors or other classed objects / columns, it simply does \code{X[X == value] <- replacement} in a more efficient way. For classed objects, see for example \code{dplyr::recode}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[=AA2-small-helpers]{Small (Helper) Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
Recode(c("a","b","c"), a = "b", b = "c")
Recode(c("a","b","c"), a = "b", b = "c", copy = TRUE)
Recode(c("a","b","c"), a = "b", b = "a", copy = TRUE)
Recode(month.name, ber = NA, regex = TRUE)
mtcr <- Recode(mtcars, `0` = 2, `4` = Inf, `1` = NaN)
replace_non_finite(mtcr)
replace_non_finite(mtcr, replace.nan = FALSE)
replace_outliers(mtcars, c(2, 100))                 # replace all values below 2 and above 100 w. NA
replace_outliers(mtcars, 2, single.limit = "min")   # replace all value smaller than 2 with NA
replace_outliers(mtcars, 100, single.limit = "max") # replace all value larger than 100 with NA
replace_outliers(mtcars, 2)                         # replace all values above or below 2 column-
                                                    # standard-deviations from the column-mean w. NA
replace_outliers(                                   # Passing a grouped_df, pseries or pdata.frame
num_vars(dplyr::group_by(iris, Species)), 2)        # allows to remove outliers according to
                                                    # in-group standard-deviation. see ?fscale

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{documentation}
