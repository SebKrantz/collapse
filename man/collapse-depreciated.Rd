\name{collapse-depreciated}
\alias{collapse-depreciated}
\alias{Recode}
\alias{replace_non_finite}
\alias{is.regular}
\title{
 Depreciated \emph{collapse} Functions
}
\description{
The functions \code{Recode} and \code{replace_non_finite} available until \emph{collapse} v1.1.0 will be removed soon. Since v1.2.0, \code{Recode} is replaced by \code{\link{recode_num}} and \code{\link{recode_char}} and \code{replace_non_finite} is replaced by \code{\link{replace_Inf}}. Since version 1.5.1, \code{is.regular} is depreciated - the function is not very useful and clashes with a more important one in the \emph{zoo} package.

%The function \code{as.factor.GRP} was renamed to \code{\link{as.factor_GRP}} to make it clear that this is not a method, and \code{group_names.GRP} was renamed to \code{\link{GRPnames}} for the same reason and to increase parsimony.
}
\usage{
Recode(X, \dots, copy = FALSE, reserve.na.nan = TRUE, regex = FALSE)

replace_non_finite(X, value = NA, replace.nan = TRUE)

is.regular(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix or data frame.}
  \item{x}{an R object. }
  \item{\dots}{comma-separated recode arguments of the form: \code{name = newname, `2` = 0, `NaN` = 0, `NA` = 0, `Inf` = NA, `-Inf` = NA}, etc\dots}
  \item{value}{a single (scalar) value to replace matching elements with. Default is \code{NA}.}
  \item{copy}{logical. For reciprocal or sequential replacements of the form \code{a = b, b = c} make a copy of \code{X} to prevent \code{a} being replaced with \code{b} and then all \code{b}-values being replaced with \code{c} again. In general \code{Recode} does the replacements one-after the other, starting with the first. }
  \item{reserve.na.nan}{logical. \code{TRUE} identifies \code{NA} and \code{NaN} as special numeric values and does the correct replacement. \code{FALSE} will treat \code{NA/NaN} as strings, and thus not match numeric \code{NA/NaN}. \emph{Note}: This is not an issue for \code{Inf/-Inf}, which are matched in both numeric and character variables. }
  \item{regex}{logical. If \code{TRUE}, all recode-argument names are (sequentially) passed to \code{\link{grepl}} as a pattern to search \code{X}. All matches are replaced.}
  \item{replace.nan}{logical. \code{TRUE} (default) replaces \code{NaN/Inf/-Inf}. \code{FALSE} replaces only \code{Inf/-Inf}.}
}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }
% \value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% \dots
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

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\link[=AA1-recode-replace]{Recode and Replace Values}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
\dontrun{
Recode(c("a","b","c"), a = "b", b = "c")
Recode(c("a","b","c"), a = "b", b = "c", copy = TRUE)
Recode(c("a","b","c"), a = "b", b = "a", copy = TRUE)
Recode(month.name, ber = NA, regex = TRUE)
mtcr <- Recode(mtcars, `0` = 2, `4` = Inf, `1` = NaN)
replace_non_finite(mtcr)
replace_non_finite(mtcr, replace.nan = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{documentation}
