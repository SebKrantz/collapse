\name{fbetween-fwithin-B-W}
\alias{B}
\alias{B.default}
\alias{B.matrix}
\alias{B.data.frame}
\alias{B.pseries}
\alias{B.pdata.frame}
\alias{B.grouped_df}
\alias{W}
\alias{W.default}
\alias{W.matrix}
\alias{W.data.frame}
\alias{W.pseries}
\alias{W.pdata.frame}
\alias{W.grouped_df}
\alias{fbetween}
\alias{fbetween.default}
\alias{fbetween.matrix}
\alias{fbetween.data.frame}
\alias{fbetween.pseries}
\alias{fbetween.pdata.frame}
\alias{fbetween.grouped_df}
\alias{fwithin}
\alias{fwithin.default}
\alias{fwithin.matrix}
\alias{fwithin.data.frame}
\alias{fwithin.pseries}
\alias{fwithin.pdata.frame}
\alias{fwithin.grouped_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Between (Averaging) and Within (Centering) Transformations
}
\description{
\code{fbetween} and \code{fwithin} are S3 generics to efficiently obtain between-transformed (averaged) or within-transformed (demeaned) data. These operations can be performed groupwise and/or weighted. \code{B} and \code{W} are wrappers around \code{fbetween} and \code{fwithin} representing the 'between-operator' and the 'within-operator'. \code{B} / \code{W} provide more flexibility than \code{fbetween} / \code{fwithin} when applied to data frames (i.e.  column subsetting, formula input, auto-renaming and id-variable-preservation capabilities...), but are otherwise identical.

(\code{fbetween} and \code{fwithin} are simple programmers functions in style of the \link[=A1-fast-statistical-functions]{Fast Statistical Functions} while \code{B} and \code{W} are more practical to use in regression formulas or for ad-hoc computations on data frames.)
}
\usage{
fbetween(x, \dots)
 fwithin(x, \dots)
       B(x, \dots)
       W(x, \dots)
% ## Default S3 methods:
\method{fbetween}{default}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{default}(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
\method{B}{default}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{W}{default}(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
% ## S3 methods for class 'matrix'
\method{fbetween}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
\method{B}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, stub = "B.", \dots)
\method{W}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, stub = "W.", \dots)
%## S3 methods for class 'data.frame'
\method{fbetween}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
\method{B}{data.frame}(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
  fill = FALSE, stub = "B.", keep.by = TRUE, keep.w = TRUE, \dots)
\method{W}{data.frame}(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
  add.global.mean = FALSE, stub = "W.", keep.by = TRUE, keep.w = TRUE, \dots)

# Methods for compatibility with plm:
%## S3 methods for class 'pseries'
\method{fbetween}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
\method{B}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{W}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
%## S3 methods for class 'pdata.frame'
\method{fbetween}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, add.global.mean = FALSE, \dots)
\method{B}{pdata.frame}(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
  fill = FALSE, stub = "B.", keep.ids = TRUE, keep.w = TRUE, \dots)
\method{W}{pdata.frame}(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
  add.global.mean = FALSE, stub = "W.", keep.ids = TRUE, keep.w = TRUE, \dots)

# Methods for compatibility with dplyr:
%## S3 methods for class 'grouped_df'
\method{fbetween}{grouped_df}(x, w = NULL, na.rm = TRUE, fill = FALSE,
         keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{fwithin}{grouped_df}(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
        keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{B}{grouped_df}(x, w = NULL, na.rm = TRUE, fill = FALSE,
  stub = "B.", keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{W}{grouped_df}(x, w = NULL, na.rm = TRUE, add.global.mean = FALSE,
  stub = "W.", keep.group_vars = TRUE, keep.w = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector, matrix, data.frame, panel-series (\code{plm::pseries}), panel-data.frame (\code{plm::pdata.frame}) or grouped tibble (\code{dplyr::grouped_df}).}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{by}{\emph{B and W data.frame method}: Same as g, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{w}{a numeric vector of (non-negative) weights. \code{B / W} \code{data.frame} and \code{pdata.frame} methods also allow a one-sided formula i.e. \code{~ weightcol}. The \code{grouped_df} (\code{dplyr}) method supports lazy-evaluation. See Examples.}
  %\item{w}{a numeric vector of (non-negative) weights. if \code{na.rm = TRUE}, the weight vector can have missing values, and those observations will be skipped in the computation of averages.}
    \item{cols}{\emph{data.frame method}: Select columns to center/average using a function, column names or indices. Default: All numeric variables. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
  \item{na.rm}{logical. skip missing values in \code{x} when computing averages. If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the average for that group will be \code{NA}, and all data points belonging to that group will also be \code{NA}.}
      \item{effect}{\code{plm} methods: Select which panel identifier should be used as grouping variable. 1L means first variable in the \code{plm::index}, 2L the second etc. if more than one integer is supplied, the corresponding index-variables are interacted. }
  \item{stub}{a prefix or stub to rename all transformed columns. \code{FALSE} will not rename columns.}
  \item{fill}{\emph{option to \code{fbetween/B}}: Logical. \code{TRUE} will overrwite missing values in \code{x} with the respective average. By default missing values in \code{x} are preserved.}
  \item{add.global.mean}{\emph{option to \code{fwithin/W}}: Logical. \code{TRUE} will add back the global mean to all data values after subtracting out group-means.}
  \item{keep.by, keep.ids, keep.group_vars}{\emph{B and W data.frame, pdata.frame and grouped_df methods}: Logical. Retain grouping / panel-identifier columns in the output. For data frames this only works if grouping variables were passed in a formula.}
  \item{keep.w}{\emph{B and W data.frame, pdata.frame and grouped_df methods}: Logical. Retain column containing the weights in the output. Only works if \code{w} is passed as formula / lazy-expression.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
Without groups, \code{fbetween/B} replaces all data points in \code{x} with their mean or weighted mean (if \code{w} is supplied). Similarly \code{fwithin/W} subtracts the mean from all data points i.e. centers the data on the mean. \cr

With groups supplied to \code{g}, the replacment / centering performed by \code{fbetween/B} | \code{fwithin/W} becomes groupwise. I like to think of this in terms of panel data: If \code{x} is a vector in such a dataset, \code{xit} denotes a single data-point belonging to group \code{i} in time-period \code{t} (\code{t} need not be a time-period). Then \code{xi.} denotes \code{x}, averaged over \code{t}. \code{fbetween/B} now returns \code{xi.} and \code{fwithin/W} returns \code{x - xi.}. Thus for any data \code{x} and any grouping vector \code{g}: \code{B(x,g) + W(x,g) = xi. + x - xi. = x}. In terms of variance, \code{fbetween/B} only retains the variance between group averages, while \code{fwithin/W}, by subtracting out group means, only retains the variance within those groups. \cr

The data replacement performed by \code{fbetween/B} can keep (default) or overwrite missing values (option \code{fill}) in \code{x}. \code{fwithin/W} can center data simply (default), or add back the global / overall mean in groupwise computations (option \code{add.global.mean}). Let \code{x..} denote the global mean of \code{x}, then \code{fwithin/W} with \code{add.global.mean = TRUE} returns \code{x - xi. + x..} instead of \code{x - xi.}. This is useful to get rid of group-differences but preserve the overall level of the data (as simple groupwise centering will set the overall mean of the data to 0). In regression analysis, centering with \code{add.global.mean = TRUE} will only change the constant term. See Examples.
}
\value{
\code{fbetween/B} returns \code{x} with every element replaced by its (groupwise) mean (\code{xi.}). \code{fwithin/W} returns \code{x} where every element was subtracted its (groupwise) mean (\code{x - xi.} or \code{x - xi. + x..}). See Details.
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
\code{\link[=HDW]{fHDbetween/HDB and fHDwithin/HDW}}, \code{\link[=fscale]{fscale/STD}}, \code{\link{TRA}}, \link[=A6-data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple centering and averaging
fbetween(mtcars)
B(mtcars)
fwithin(mtcars)
W(mtcars)
fbetween(mtcars) + fwithin(mtcars) == mtcars # This should be true apart from rounding errors

## Groupwise centering and averaging
fbetween(mtcars, mtcars$cyl)
 fwithin(mtcars, mtcars$cyl)
fbetween(mtcars, mtcars$cyl) + fwithin(mtcars, mtcars$cyl) == mtcars

W(wlddev, ~ iso3c, cols = 9:12)    # Center the 4 series in this dataset by country
cbind(get_vars(wlddev,"iso3c"),    # Same thing done manually using fwithin...
      add_stub(fwithin(get_vars(wlddev,9:12), wlddev$iso3c), "W."))

## Using B() and W() in regressions:

# Several ways of running the same regression with cyl-fixed effects
lm(W(mpg,cyl) ~ W(carb,cyl), data = mtcars)                     # Centering each individually
lm(mpg ~ carb, data = W(mtcars, ~ cyl, stub = FALSE))           # Centering the entire data
lm(mpg ~ carb, data = W(mtcars, ~ cyl, stub = FALSE,            # Here only the intercept changes
                        add.global.mean = TRUE))
lm(mpg ~ carb + B(carb,cyl), data = mtcars)                     # Procedure suggested by
# ...Mundlack (1978) - partialling out group averages amounts to the same as demeaning the data

# Now with cyl, vs and am fixed effects
lm(W(mpg,list(cyl,vs,am)) ~ W(carb,list(cyl,vs,am)), data = mtcars)
lm(mpg ~ carb, data = W(mtcars, ~ cyl + vs + am, stub = FALSE))
lm(mpg ~ carb + B(carb,list(cyl,vs,am)), data = mtcars)

# Now with cyl, vs and am fixed effects weighted by hp:
lm(W(mpg,list(cyl,vs,am),hp) ~ W(carb,list(cyl,vs,am),hp), data = mtcars)
lm(mpg ~ carb, data = W(mtcars, ~ cyl + vs + am, ~ hp, stub = FALSE))
lm(mpg ~ carb + B(carb,list(cyl,vs,am),hp), data = mtcars)       # Gives a different coefficient!!

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line
