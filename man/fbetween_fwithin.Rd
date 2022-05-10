\name{fbetween-fwithin}
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
Fast Between (Averaging) and (Quasi-)Within (Centering) Transformations
}
\description{
\code{fbetween} and \code{fwithin} are S3 generics to efficiently obtain between-transformed (averaged) or (quasi-)within-transformed (demeaned) data. These operations can be performed groupwise and/or weighted. \code{B} and \code{W} are wrappers around \code{fbetween} and \code{fwithin} representing the 'between-operator' and the 'within-operator'.

(\code{B} / \code{W} provide more flexibility than \code{fbetween} / \code{fwithin} when applied to data frames (i.e.  column subsetting, formula input, auto-renaming and id-variable-preservation capabilities\dots), but are otherwise identical.)

%(\code{fbetween} and \code{fwithin} are simple programmers functions in style of the \link[=fast-statistical-functions]{Fast Statistical Functions} while \code{B} and \code{W} are more practical to use in regression formulas or for ad-hoc computations on data frames.)
}
\usage{
fbetween(x, \dots)
 fwithin(x, \dots)
       B(x, \dots)
       W(x, \dots)

\method{fbetween}{default}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{default}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)
\method{B}{default}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{W}{default}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)

\method{fbetween}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)
\method{B}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, stub = "B.", \dots)
\method{W}{matrix}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, stub = "W.", \dots)

\method{fbetween}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{data.frame}(x, g = NULL, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)
\method{B}{data.frame}(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
  fill = FALSE, stub = "B.", keep.by = TRUE, keep.w = TRUE, \dots)
\method{W}{data.frame}(x, by = NULL, w = NULL, cols = is.numeric, na.rm = TRUE,
  mean = 0, theta = 1, stub = "W.", keep.by = TRUE, keep.w = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{fbetween}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)
\method{B}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{W}{pseries}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)

\method{fbetween}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, fill = FALSE, \dots)
\method{fwithin}{pdata.frame}(x, effect = 1L, w = NULL, na.rm = TRUE, mean = 0, theta = 1, \dots)
\method{B}{pdata.frame}(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
  fill = FALSE, stub = "B.", keep.ids = TRUE, keep.w = TRUE, \dots)
\method{W}{pdata.frame}(x, effect = 1L, w = NULL, cols = is.numeric, na.rm = TRUE,
  mean = 0, theta = 1, stub = "W.", keep.ids = TRUE, keep.w = TRUE, \dots)

# Methods for grouped data frame / compatibility with dplyr:

\method{fbetween}{grouped_df}(x, w = NULL, na.rm = TRUE, fill = FALSE,
         keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{fwithin}{grouped_df}(x, w = NULL, na.rm = TRUE, mean = 0, theta = 1,
        keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{B}{grouped_df}(x, w = NULL, na.rm = TRUE, fill = FALSE,
  stub = "B.", keep.group_vars = TRUE, keep.w = TRUE, \dots)
\method{W}{grouped_df}(x, w = NULL, na.rm = TRUE, mean = 0, theta = 1,
  stub = "W.", keep.group_vars = TRUE, keep.w = TRUE, \dots)
}

\arguments{
    \item{x}{a numeric vector, matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df').}
  \item{g}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{group}}) used to group \code{x}.}
  \item{by}{\emph{B and W data.frame method}: Same as g, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{w}{a numeric vector of (non-negative) weights. \code{B}/\code{W} data frame and \code{pdata.frame} methods also allow a one-sided formula i.e. \code{~ weightcol}. The \code{grouped_df} (\emph{dplyr}) method supports lazy-evaluation. See Examples.}
    \item{cols}{\emph{B/W (p)data.frame methods}: Select columns to scale using a function, column names, indices or a logical vector. Default: All numeric columns. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{na.rm}{logical. Skip missing values in \code{x} and \code{w} when computing averages. If \code{na.rm = FALSE} and a \code{NA} or \code{NaN} is encountered, the average for that group will be \code{NA}, and all data points belonging to that group in the output vector will also be \code{NA}.}
      \item{effect}{\emph{plm} methods: Select which panel identifier should be used as grouping variable. 1L takes the first variable in the \link[=indexing]{index}, 2L the second etc. Index variables can also be called by name using a character string. If more than one variable is supplied, the corresponding index-factors are interacted. }
  \item{stub}{a prefix or stub to rename all transformed columns. \code{FALSE} will not rename columns.}
  \item{fill}{\emph{option to \code{fbetween}/\code{B}}: Logical. \code{TRUE} will overwrite missing values in \code{x} with the respective average. By default missing values in \code{x} are preserved.}
  \item{mean}{\emph{option to \code{fwithin}/\code{W}}: The mean to center on, default is 0, but a different mean can be supplied and will be added to the data after the centering is performed. A special option when performing grouped centering is \code{mean = "overall.mean"}. In that case the overall mean of the data will be added after subtracting out group means.}
  \item{theta}{\emph{option to \code{fwithin}/\code{W}}: Double. An optional scalar parameter for quasi-demeaning i.e. \code{x - theta * xi.}. This is useful for variance components ('random-effects') estimators. see Details.}
  \item{keep.by, keep.ids, keep.group_vars}{\emph{B and W data.frame, pdata.frame and grouped_df methods}: Logical. Retain grouping / panel-identifier columns in the output. For data frames this only works if grouping variables were passed in a formula.}
  \item{keep.w}{\emph{B and W data.frame, pdata.frame and grouped_df methods}: Logical. Retain column containing the weights in the output. Only works if \code{w} is passed as formula / lazy-expression.}
  \item{\dots}{arguments to be passed to or from other methods.}
}
\details{
Without groups, \code{fbetween}/\code{B} replaces all data points in \code{x} with their mean or weighted mean (if \code{w} is supplied). Similarly \code{fwithin/W} subtracts the (weighted) mean from all data points i.e. centers the data on the mean. \cr

With groups supplied to \code{g}, the replacement / centering performed by \code{fbetween/B} | \code{fwithin/W} becomes groupwise. In terms of panel data notation: If \code{x} is a vector in such a panel dataset, \code{xit} denotes a single data-point belonging to group \code{i} in time-period \code{t} (\code{t} need not be a time-period). Then \code{xi.} denotes \code{x}, averaged over \code{t}. \code{fbetween}/\code{B} now returns \code{xi.} and \code{fwithin}/\code{W} returns \code{x - xi.}. Thus for any data \code{x} and any grouping vector \code{g}: \code{B(x,g) + W(x,g) = xi. + x - xi. = x}. In terms of variance, \code{fbetween/B} only retains the variance between group averages, while \code{fwithin}/\code{W}, by subtracting out group means, only retains the variance within those groups. \cr

The data replacement performed by \code{fbetween}/\code{B} can keep (default) or overwrite missing values (option \code{fill = TRUE}) in \code{x}. \code{fwithin/W} can center data simply (default), or add back a mean after centering (option \code{mean = value}), or add the overall mean in groupwise computations (option \code{mean = "overall.mean"}). Let \code{x..} denote the overall mean of \code{x}, then \code{fwithin}/\code{W} with \code{mean = "overall.mean"} returns \code{x - xi. + x..} instead of \code{x - xi.}. This is useful to get rid of group-differences but preserve the overall level of the data. In regression analysis, centering with \code{mean = "overall.mean"} will only change the constant term. See Examples.

If \code{theta != 1}, \code{fwithin}/\code{W} performs quasi-demeaning \code{x - theta * xi.}. If \code{mean = "overall.mean"}, \code{x - theta * xi. + theta * x..} is returned, so that the mean of the partially demeaned data is still equal to the overall data mean \code{x..}. A numeric value passed to \code{mean} will simply be added back to the quasi-demeaned data i.e. \code{x - theta * xi. + mean}.

Now in the case of a linear panel model \eqn{y_{it} = \beta_0 + \beta_1 X_{it} + u_{it}} with \eqn{u_{it} = \alpha_i + \epsilon_{it}}. If \eqn{\alpha_i \neq \alpha = const.} (there exists individual heterogeneity), then pooled OLS is at least inefficient and inference on \eqn{\beta_1} is invalid. If \eqn{E[\alpha_i|X_{it}] = 0} (mean independence of individual heterogeneity \eqn{\alpha_i}), the variance components or 'random-effects' estimator provides an asymptotically efficient FGLS solution by estimating a transformed model \eqn{y_{it}-\theta y_{i.}  = \beta_0 + \beta_1 (X_{it} - \theta X_{i.}) + (u_{it} - \theta u_{i.}}), where \eqn{\theta = 1 - \frac{\sigma_\alpha}{\sqrt(\sigma^2_\alpha + T \sigma^2_\epsilon)}}. An estimate of \eqn{\theta} can be obtained from the an estimate of \eqn{\hat{u}_{it}} (the residuals from the pooled model). If \eqn{E[\alpha_i|X_{it}] \neq 0}, pooled OLS is biased and inconsistent, and taking \eqn{\theta = 1} gives an unbiased and consistent fixed-effects estimator of \eqn{\beta_1}. See Examples.
}
\value{
\code{fbetween}/\code{B} returns \code{x} with every element replaced by its (groupwise) mean (\code{xi.}). Missing values are preserved if \code{fill = FALSE} (the default). \code{fwithin/W} returns \code{x} where every element was subtracted its (groupwise) mean (\code{x - theta * xi. + mean} or, if \code{mean = "overall.mean"}, \code{x - theta * xi. + theta * x..}). See Details.
}
\references{
Mundlak, Yair. 1978. On the Pooling of Time Series and Cross Section Data. \emph{Econometrica} 46 (1): 69-85.
}
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link[=HDW]{fhdbetween/HDB and fhdwithin/HDW}}, \code{\link[=fscale]{fscale/STD}}, \code{\link{TRA}}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple centering and averaging
head(fbetween(mtcars))
head(B(mtcars))
head(fwithin(mtcars))
head(W(mtcars))
all.equal(fbetween(mtcars) + fwithin(mtcars), mtcars)

## Groupwise centering and averaging
head(fbetween(mtcars, mtcars$cyl))
head(fwithin(mtcars, mtcars$cyl))
all.equal(fbetween(mtcars, mtcars$cyl) + fwithin(mtcars, mtcars$cyl), mtcars)

head(W(wlddev, ~ iso3c, cols = 9:13))    # Center the 5 series in this dataset by country
head(cbind(get_vars(wlddev,"iso3c"),     # Same thing done manually using fwithin..
      add_stub(fwithin(get_vars(wlddev,9:13), wlddev$iso3c), "W.")))

## Using B() and W() for fixed-effects regressions:

# Several ways of running the same regression with cyl-fixed effects
lm(W(mpg,cyl) ~ W(carb,cyl), data = mtcars)                     # Centering each individually
lm(mpg ~ carb, data = W(mtcars, ~ cyl, stub = FALSE))           # Centering the entire data
lm(mpg ~ carb, data = W(mtcars, ~ cyl, stub = FALSE,            # Here only the intercept changes
                        mean = "overall.mean"))
lm(mpg ~ carb + B(carb,cyl), data = mtcars)                     # Procedure suggested by
# ..Mundlak (1978) - partialling out group averages amounts to the same as demeaning the data
\donttest{ % No code relying on suggested package
plm::plm(mpg ~ carb, mtcars, index = "cyl", model = "within")   # "Proof"..
}
# This takes the interaction of cyl, vs and am as fixed effects
lm(W(mpg) ~ W(carb), data = iby(mtcars, id = finteraction(cyl, vs, am)))
lm(mpg ~ carb, data = W(mtcars, ~ cyl + vs + am, stub = FALSE))
lm(mpg ~ carb + B(carb,list(cyl,vs,am)), data = mtcars)

# Now with cyl fixed effects weighted by hp:
lm(W(mpg,cyl,hp) ~ W(carb,cyl,hp), data = mtcars)
lm(mpg ~ carb, data = W(mtcars, ~ cyl, ~ hp, stub = FALSE))
lm(mpg ~ carb + B(carb,cyl,hp), data = mtcars)       # WRONG ! Gives a different coefficient!!

## Manual variance components (random-effects) estimation
res <- HDW(mtcars, mpg ~ carb)[[1]]  # Get residuals from pooled OLS
sig2_u <- fvar(res)
sig2_e <- fvar(fwithin(res, mtcars$cyl))
T <- length(res) / fndistinct(mtcars$cyl)
sig2_alpha <- sig2_u - sig2_e
theta <- 1 - sqrt(sig2_alpha) / sqrt(sig2_alpha + T * sig2_e)
lm(mpg ~ carb, data = W(mtcars, ~ cyl, theta = theta, mean = "overall.mean", stub = FALSE))
\donttest{ % No code relying on suggested package
# A slightly different method to obtain theta...
plm::plm(mpg ~ carb, mtcars, index = "cyl", model = "random")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line
