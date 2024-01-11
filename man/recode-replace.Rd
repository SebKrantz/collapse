\name{recode-replace}
\alias{AA1-recode-replace}
\alias{recode-replace}
\alias{recode_num}
\alias{recode_char}
\alias{replace_na}
\alias{replace_inf}
\alias{replace_outliers}
\title{
Recode and Replace Values in Matrix-Like Objects
}
\description{
A small suite of functions to efficiently perform common recoding and replacing tasks in matrix-like objects.
}
\usage{
recode_num(X, \dots, default = NULL, missing = NULL, set = FALSE)

recode_char(X, \dots, default = NULL, missing = NULL, regex = FALSE,
            ignore.case = FALSE, fixed = FALSE, set = FALSE)

replace_na(X, value = 0, cols = NULL, set = FALSE, type = "const")

replace_inf(X, value = NA, replace.nan = FALSE, set = FALSE)

replace_outliers(X, limits, value = NA,
                 single.limit = c("sd", "mad", "min", "max"),
                 ignore.groups = FALSE, set = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a vector, matrix, array, data frame or list of atomic objects. \code{replace_outliers} has internal methods for \link[=fgroup_by]{grouped} and \link[=findex_by]{indexed} data.}
  \item{\dots}{comma-separated recode arguments of the form: \code{value = replacement, `2` = 0, Secondary = "SEC"} etc. \code{recode_char} with \code{regex = TRUE} also supports regular expressions i.e. \code{`^S|D$` = "STD"} etc.}
\item{default}{optional argument to specify a scalar value to replace non-matched elements with.}
\item{missing}{optional argument to specify a scalar value to replace missing elements with. \emph{Note} that to increase efficiency this is done before the rest of the recoding i.e. the recoding is performed on data where missing values are filled!}
\item{set}{logical. \code{TRUE} does replacements by reference (i.e. in-place modification of the data) and returns the result invisibly.}
\item{type}{character. One of \code{"const"}, \code{"locf"} (last non-missing observation carried forward) or \code{"focb"} (first non-missing observation carried back). The latter two ignore \code{value}.}
\item{regex}{logical. If \code{TRUE}, all recode-argument names are (sequentially) passed to \code{\link{grepl}} as a pattern to search \code{X}. All matches are replaced. \emph{Note} that \code{NA}'s are also matched as strings by \code{grepl}. }
\item{value}{a single (scalar) value to replace matching elements with. In \code{replace_outliers} setting \code{value = "clip"} will replace outliers with the corresponding threshold values. See Examples.}
\item{cols}{select columns to replace missing values in using a function, column names, indices or a logical vector.}
  \item{replace.nan}{logical. \code{TRUE} replaces \code{NaN/Inf/-Inf}. \code{FALSE} (default) replaces only \code{Inf/-Inf}.}
  \item{limits}{either a vector of two-numeric values \code{c(minval, maxval)} constituting a two-sided outlier threshold, or a single numeric value:}

  \item{single.limit}{character, controls the behavior if \code{length(limits) == 1}:
  \itemize{
  \item \code{"sd"/"mad":} \code{limits} will be interpreted as a (two-sided) outlier threshold in terms of (column) standard deviations/median absolute deviations. For the standard deviation this is equivalent to \code{X[abs(fscale(X)) > limits] <- value}. Since \code{fscale} is S3 generic with methods for 'grouped_df', 'pseries' and 'pdata.frame', the standardizing will be grouped if such objects are passed (i.e. the outlier threshold is then measured in within-group standard deviations) unless \code{ignore.groups = TRUE}. The same holds for median absolute deviations.
  \item \code{"min"/"max":} \code{limits} will be interpreted as a (one-sided) minimum/maximum threshold. The underlying code is equivalent to \code{X[X </> limits] <- value}.
}
}
\item{ignore.groups}{logical. If \code{length(limits) == 1} and \code{single.limit \%in\% c("sd", "mad")} and \code{X} is a 'grouped_df', 'pseries' or 'pdata.frame', \code{TRUE} will ignore the grouped nature of the data and calculate outlier thresholds on the entire dataset rather than within each group.}

\item{ignore.case, fixed}{logical. Passed to \code{\link{grepl}} and only applicable if \code{regex = TRUE}.}
}
\details{
\itemize{
\item \code{recode_num} and \code{recode_char} can be used to efficiently recode multiple numeric or character values, respectively. The syntax is inspired by \code{dplyr::recode}, but the functionality is enhanced in the following respects: (1) when passed a data frame / list, all appropriately typed columns will be recoded. (2) They preserve the attributes of the data object and of columns in a data frame / list, and (3) \code{recode_char} also supports regular expression matching using \code{\link{grepl}}.
\item \code{replace_na} efficiently replaces \code{NA/NaN} with a value (default is \code{0}). data can be multi-typed, in which case appropriate columns can be selected through the \code{cols} argument. For numeric data a more versatile alternative is provided by \code{data.table::nafill} and \code{data.table::setnafill}.
\item \code{replace_inf} replaces \code{Inf/-Inf} (or optionally \code{NaN/Inf/-Inf}) with a value (default is \code{NA}). It skips non-numeric columns in a data frame.
\item \code{replace_outliers} replaces values falling outside a 1- or 2-sided numeric threshold or outside a certain number of standard deviations or median absolute deviation with a value (default is \code{NA}). It skips non-numeric columns in a data frame.
}
}
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
These functions are not generic and do not offer support for factors or date(-time) objects. see \code{dplyr::recode_factor}, \emph{forcats} and other appropriate packages for dealing with these classes.

Simple replacing tasks on a vector can also effectively be handled by, \code{\link{setv}} / \code{\link{copyv}}. Fast vectorized switches are offered by package \emph{kit} (functions \code{iif}, \code{nif}, \code{vswitch}, \code{nswitch}) as well as \code{data.table::fcase} and \code{data.table::fifelse}. Using switches is more efficient than \code{recode_*}, as \code{recode_*} creates an internal copy of the object to enable cross-replacing.

Function \code{\link{TRA}}, and the associated \code{TRA} ('transform') argument to \link[=fast-statistical-functions]{Fast Statistical Functions} also has option \code{"replace_na"}, to replace missing values with a statistic computed on the non-missing observations, e.g. \code{fmedian(airquality, TRA = "replace_na")} does median imputation.
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{pad}}, \link[=efficient-programming]{Efficient Programming}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
recode_char(c("a","b","c"), a = "b", b = "c")
recode_char(month.name, ber = NA, regex = TRUE)
mtcr <- recode_num(mtcars, `0` = 2, `4` = Inf, `1` = NaN)
replace_inf(mtcr)
replace_inf(mtcr, replace.nan = TRUE)
replace_outliers(mtcars, c(2, 100))                 # Replace all values below 2 and above 100 w. NA
replace_outliers(mtcars, c(2, 100), value = "clip") # Clipping outliers to the thresholds
replace_outliers(mtcars, 2, single.limit = "min")   # Replace all value smaller than 2 with NA
replace_outliers(mtcars, 100, single.limit = "max") # Replace all value larger than 100 with NA
replace_outliers(mtcars, 2)                         # Replace all values above or below 2 column-
                                                    # standard-deviations from the column-mean w. NA
replace_outliers(fgroup_by(iris, Species), 2)       # Passing a grouped_df, pseries or pdata.frame
                                                    # allows to remove outliers according to
                                                    # in-group standard-deviation. see ?fscale
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{documentation}
