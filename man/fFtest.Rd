\name{fFtest}
\alias{fFtest}
\alias{fFtest.default}
\alias{fFtest.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Weighted) F-test for Linear Models (with Factors)
}
\description{
\code{fFtest} computes an R-squared based F-test for the exclusion of the variables in \code{exc}, where the full (unrestricted) model is defined by variables supplied to both \code{exc} and \code{X}. The test is efficient and designed for cases where both \code{exc} and \code{X} may contain multiple factors and continuous variables. There is also an efficient 2-part formula method.
}
\usage{
fFtest(...) # Internal method dispatch: formula if is.call(..1) || is.call(..2)

\method{fFtest}{default}(y, exc, X = NULL, w = NULL, full.df = TRUE, \dots)

\method{fFtest}{formula}(formula, data = NULL, weights = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a numeric vector: the dependent variable.}
  \item{exc}{a numeric vector, factor, numeric matrix or list / data frame of numeric vectors and/or factors: variables to test / exclude.}
  \item{X}{a numeric vector, factor, numeric matrix or list / data frame of numeric vectors and/or factors: covariates to include in both the restricted (without \code{exc}) and unrestricted model. If left empty (\code{X = NULL}), the test amounts to the F-test of the regression of \code{y} on \code{exc}.}
  \item{w}{numeric. A vector of (frequency) weights.}
  \item{formula}{a 2-part formula: \code{y ~ exc | X}, where both \code{exc} and \code{X} are expressions connected with \code{+}, and \code{X} can be omitted. \emph{Note} that other operators (\code{:}, \code{*}, \code{^}, \code{-}, etc.) are not supported, you can interact variables using standard functions like \code{\link[=itn]{finteraction/itn}} or \code{magrittr::multiply_by} inside the formula e.g. \code{log(y) ~ x1 + itn(x2, x3) | x4} or \code{log(y) ~ x1 + multiply_by(x2, x3) | x4}.}
 \item{data}{a named list or data frame.}
\item{weights}{a weights vector or expression that results in a vector when evaluated in the \code{data} environment.}
  \item{full.df}{logical. If \code{TRUE} (default), the degrees of freedom are calculated as if both restricted and unrestricted models were estimated using \code{lm()} (i.e. as if factors were expanded to matrices of dummies). \code{FALSE} only uses one degree of freedom per factor.  }
\item{\dots}{other arguments passed to \code{fFtest.default} or to \code{fhdwithin}. Sensible options might be the \code{lm.method} argument or further control parameters to \code{fixest::demean}, the workhorse function underlying \code{fhdwithin} for higher-order centering tasks. }

}
\details{
Factors and continuous regressors are efficiently projected out using \code{\link{fhdwithin}}, and the option \code{full.df} regulates whether a degree of freedom is subtracted for each used factor level (equivalent to dummy-variable estimator / expanding factors), or only one degree of freedom per factor (treating factors as variables). The test automatically removes missing values and considers only the complete cases of \code{y, exc} and \code{X}. Unused factor levels in \code{exc} and \code{X} are dropped.

\emph{Note} that an intercept is always added by \code{\link{fhdwithin}}, so it is not necessary to include an intercept in data supplied to \code{exc} / \code{X}.
}
\value{
A 5 x 3 numeric matrix of statistics. The columns contain statistics:
\enumerate{
\item the R-squared of the model
\item the numerator degrees of freedom i.e. the number of variables (k) and used factor levels if \code{full.df = TRUE}
\item the denominator degrees of freedom: N - k - 1.
\item the F-statistic
\item the corresponding P-value
}
The rows show these statistics for:
\enumerate{
\item the Full (unrestricted) Model (\code{y ~ exc + X})
\item the Restricted Model (\code{y ~ X})
\item the Exclusion Restriction of \code{exc}. The R-squared shown is simply the difference of the full and restricted R-Squared's, not the R-Squared of the model \code{y ~ exc}.
}
If \code{X = NULL}, only a vector of the same 5 statistics testing the model (\code{y ~ exc}) is shown.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% \dots
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{flm}}, \code{\link{fhdwithin}}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## We could use fFtest as a simple seasonality test:
fFtest(AirPassengers, qF(cycle(AirPassengers)))         # Testing for level-seasonality
fFtest(AirPassengers, qF(cycle(AirPassengers)),         # Seasonality test around a cubic trend
        poly(seq_along(AirPassengers), 3))
fFtest(fdiff(AirPassengers), qF(cycle(AirPassengers)))  # Seasonality in first-difference

## A more classical example with only continuous variables
fFtest(mpg ~ cyl + vs | hp + carb, mtcars)
fFtest(mtcars$mpg, mtcars[c("cyl","vs")], mtcars[c("hp","carb")])
\donttest{ % requires fixest package
## Now encoding cyl and vs as factors
fFtest(mpg ~ qF(cyl) + qF(vs) | hp + carb, mtcars)
fFtest(mtcars$mpg, lapply(mtcars[c("cyl","vs")], qF), mtcars[c("hp","carb")])
}
## Using iris data: A factor and a continuous variable excluded
fFtest(Sepal.Length ~ Petal.Width + Species | Sepal.Width + Petal.Length, iris)
fFtest(iris$Sepal.Length, iris[4:5], iris[2:3])

## Testing the significance of country-FE in regression of GDP on life expectancy
fFtest(log(PCGDP) ~ iso3c | LIFEEX, wlddev)
fFtest(log(wlddev$PCGDP), wlddev$iso3c, wlddev$LIFEEX)
\donttest{ % requires fixest package
## Ok, country-FE are significant, what about adding time-FE
fFtest(log(PCGDP) ~ qF(year) | iso3c + LIFEEX, wlddev)
fFtest(log(wlddev$PCGDP), qF(wlddev$year), wlddev[c("iso3c","LIFEEX")])
}
# Same test done using lm:
data <- na_omit(get_vars(wlddev, c("iso3c","year","PCGDP","LIFEEX")))
full <- lm(PCGDP ~ LIFEEX + iso3c + qF(year), data)
rest <- lm(PCGDP ~ LIFEEX + iso3c, data)
anova(rest, full)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest} % use one of  RShowDoc("KEYWORDS")

