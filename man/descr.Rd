\name{descr}
\alias{descr}
\alias{[.descr}
\alias{print.descr}
\alias{as.data.frame.descr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detailed Statistical Description of Data Frame
}
\description{
\code{descr} offers a concise description of each variable in a data frame. It is built as a wrapper around \code{\link{qsu}}, but also computes frequency tables for categorical variables, and quantiles and the number of distinct values for numeric variables.

}
\usage{
descr(X, Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
      Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
      cols = NULL, label.attr = "label", stepwise = FALSE, \dots)

\method{[}{descr}(x, ...)

\method{as.data.frame}{descr}(x, \dots)

\method{print}{descr}(x, n = 14, perc = TRUE, digits = 2, t.table = TRUE,
      summary = TRUE, reverse = FALSE, stepwise = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a data frame or list of atomic vectors. Atomic vectors, matrices or arrays can be passed but will first be coerced to data frame using \code{\link{qDF}}.
}
  \item{Ndistinct}{
logical. \code{TRUE} (default) computes the number of distinct values on all variables using \code{\link{fndistinct}}.
}
  \item{higher}{
logical. Argument is passed down to \code{\link{qsu}}: \code{TRUE} (default) computes the skewness and the kurtosis.
}
  \item{table}{
logical. \code{TRUE} (default) computes a (sorted) frequency table for all categorical variables (excluding \link[=is_date]{Date} variables).
}

\item{sort.table}{an integer or character string specifying how the frequency table should be presented:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "value"   \tab\tab sort table by values. \cr
                 2 \tab\tab "freq" \tab\tab sort table by frequencies. \cr
                 3 \tab\tab "none"   \tab\tab return table in first-appearance order of values, or levels for factors (most efficient). \cr
  }
  }

  \item{Qprobs}{
double. Probabilities for quantiles to compute on numeric variables, passed down to \code{\link{quantile}}. If something non-numeric is passed (i.e. \code{NULL}, \code{FALSE}, \code{NA}, \code{""} etc.), no quantiles are computed.
}
  \item{cols}{
select columns to describe using column names, indices, a logical vector or a function (e.g. \code{is.numeric}).
}
  \item{label.attr}{
character. The name of a label attribute to display for each variable (if variables are labeled).
}
  \item{\dots}{for \code{descr}: other arguments passed to \code{\link{qsu.default}}. For \code{[.descr}: variable names or indices passed to \code{[.list}. The argument is unused in the \code{print} and \code{as.data.frame} methods.}

\item{x}{an object of class 'descr'.}
\item{n}{integer. The maximum number of table elements to print for categorical variables. If the number of distinct elements is \code{<= n}, the whole table is printed. Otherwise the remaining items are grouped into an '... \%s Others' category.}
\item{perc}{logical. \code{TRUE} (default) adds percentages to the frequencies in the table for categorical variables.}
\item{digits}{integer. The number of decimals to print in statistics and percentage tables.}
\item{t.table}{logical. \code{TRUE} (default) prints a transposed table.}
\item{summary}{logical. \code{TRUE} (default) computes and displays a summary of the frequencies, if the size of the table for a categorical variable exceeds \code{n}.}
\item{reverse}{logical. \code{TRUE} prints contents in reverse order, starting with the last column, so that the dataset can be analyzed by scrolling up the console after calling \code{descr}.}
\item{stepwise}{logical. \code{TRUE} prints one variable at a time. The user needs to press [enter] to see the printout for the next variable. If called from \code{descr}, the computation is also done one variable at a time, and the finished 'descr' object is returned invisibly. This is recommended for larger datasets, where precomputing the statistics for all variables can be time consuming. }
}
\details{
\code{descr} was heavily inspired by \code{Hmisc::describe}, but computes about 10x faster. The performance is comparable to \code{\link{summary}}. \code{descr} was built as a wrapper around \code{\link{qsu}}, to enrich the set of statistics computed by \code{\link{qsu}} for both numeric and categorical variables.

\code{\link{qsu}} itself is yet about 10x faster than \code{descr}, and is optimized for grouped, panel data and weighted statistics. It is possible to also compute grouped, panel data and/or weighted statistics with \code{descr} by passing group-ids to \code{g}, panel-ids to \code{pid} or a weight vector to \code{w}. These arguments are handed down to \code{\link{qsu.default}} and only affect the statistics natively computed by \code{qsu}, i.e. passing a weight vector produces a weighted mean, sd, skewness and kurtosis but not weighted quantiles.

The list-object returned from \code{descr} can be converted to a tidy data frame using \code{as.data.frame}. This representation will not include frequency tables computed for categorical variables, and the method cannot handle arrays of statistics (applicable when \code{g} or \code{pid} arguments are passed to \code{descr}, in that case \code{as.data.frame.descr} will throw an appropriate error).
}
\value{
A 2-level nested list, the top-level containing the statistics computed for each variable, which are themselves stored in a list containing the class, the label, the basic statistics and quantiles / tables computed for the variable. The object is given a class 'descr' and also has the number of observations in the dataset attached as an 'N' attribute, as well as an attribute 'arstat' indicating whether the object contains arrays of statistics, and an attribute 'table' indicating whether \code{table = TRUE} (i.e. the object could contain tables for categorical variables).
}

\seealso{
\code{\link{qsu}}, \code{\link{pwcor}}, \link[=summary-statistics]{Summary Statistics}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Standard Use
descr(iris)
descr(wlddev)
descr(GGDC10S)

# Some useful print options (also try stepwise argument)
print(descr(GGDC10S), reverse = TRUE, t.table = FALSE)
# For bigger data consider: descr(big_data, stepwise = TRUE)

# Generating a data frame
as.data.frame(descr(wlddev, table = FALSE))

## Passing Arguments down to qsu.default: For Panel Data Statistics
descr(iris, pid = iris$Species)
descr(wlddev, pid = wlddev$iso3c)

## Grouped Statistics
descr(iris, g = iris$Species)
descr(GGDC10S, g = GGDC10S$Region)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }% use one of  RShowDoc("KEYWORDS")
