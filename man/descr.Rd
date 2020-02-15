\name{descr}
\alias{descr}
\alias{print.descr}
\alias{as.data.frame.descr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detailed Statistical Description of Data Frame
}
\description{
\code{descr} offers concise description of each variable in a \code{data.frame}. It is built as a wrapper around \code{qsu}, but by default also computes frequency tables with percentages for categorical variables, and quantiles and the number of distinct values for numeric variables (next to the mean, sd, min, max, skewness and kurtosis computed by \code{qsu}).

}
\usage{
descr(X, Ndistinct = TRUE, higher = TRUE, table = TRUE,
      Qprobs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
      label.attr = "label", ...)

\method{print}{descr}(x, n = 6, perc = TRUE, summary = TRUE, ...)

\method{as.data.frame}{descr}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a data.frame or list of atomic vectors.
}
  \item{Ndistinct}{
logical. \code{TRUE} (default) computes the number of distinct values on all variables using \code{\link{fNdistinct}}.
}
  \item{higher}{
logical. Argument is passed down to \code{\link{qsu}}: \code{TRUE} (default) computes the skewness and the kurtosis.
}
  \item{table}{
logical. \code{TRUE} (default) calls \code{\link{table}} on all categorical variables (excluding \link[=is.Date]{Date} variables).
}
  \item{Qprobs}{
probabilities for quantiles to compute on numeric variables, passed down to \code{\link{quantile}}. If something non-numeric is passed (i.e. \code{NULL}, \code{FALSE}, \code{NA}, \code{""} etc.), no quantiles are computed.
}
  \item{label.attr}{
character. The name of a label attribute to display for each variable (if variables are labeled).
}
  \item{\dots}{
other arguments passed to \code{\link{qsu.default}}.
}
\item{x}{an object of class 'descr'.}
\item{n}{integer. The number of first and last entries to display of the table computed for categorical variables.}
\item{perc}{logical. \code{TRUE} (default) adds percentages in brackets behind the frequencies in the table for categorical variables.}
\item{summary}{logical. \code{TRUE} (default) computes and displays a summary of the frequencies if the size of the table for a categorical variables exceeds \code{2*n}.}
}
\details{
\code{descr} was heavily inspired by \code{Hmisc::describe}, but computes about 10x faster. The performance is comparable to \code{base::summary}. \code{descr} was built as a wrapper around \code{\link{qsu}}, to enrich the set of statistics computed by \code{\link{qsu}} for both numeric and categorical variables.

\code{\link{qsu}} itself is yet about 10x faster than \code{descr}, and is optimized for grouped, panel-data and weighted statistics. It is possible to also compute grouped, panel-data and/or weighted statistics with \code{desc} by passing group-ids to \code{g}, panel-ids to \code{pid} or a weight vector to \code{w}. These arguments are handed down to \code{\link{qsu.default}} and only affect the statistics natively computed by \code{qsu}, i.e. passing a weight vector produces a weighted mean, sd, skewness and kurtosis but not weighted quantiles.

The list-object returned from \code{descr} can be converted to a tidy data.frame using \code{as.data.frame}. This representation will not include frequency tables computed for categorical variables, and the method cannot handle arrays of statistics (applicable when \code{g} or \code{pid} arguments are passed to \code{descr}, in that case \code{as.data.frame.descr} will throw an appropriate error).
}
\value{
A 2-level nested list, the top-level containing the statistics computed for each variable, which are themselves stored in a list containing the class, the label, the basic statistics and quantiles / tables computed for the variable. The object is given a class 'descr' and also has the number of observations in the dataset attached as an 'N' attribute, as well as an attribute 'arstat' indicating whether the object contains arrays of statistics.
}

\seealso{
\code{\link{qsu}}, \code{\link{pwcor}}, \link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Standard Use
descr(iris)
descr(wlddev)
descr(GGDC10S)

as.data.frame(descr(wlddev))

## Passing Arguments down to qsu: For Panel-Data Statistics
descr(iris, pid = iris$species)
descr(wlddev, pid = wlddev$iso3c)

## Grouped Statistics
descr(iris, g = iris$species)
descr(GGDC10S, g = GGDC10S$Region)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }% use one of  RShowDoc("KEYWORDS")
