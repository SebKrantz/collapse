\name{descr}
\alias{descr}
\alias{descr.default}
\alias{descr.grouped_df}
\alias{[.descr}
\alias{print.descr}
\alias{as.data.frame.descr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detailed Statistical Description of Data Frame
}
\description{
\code{descr} offers a fast and detailed description of each variable in a data frame. Since v1.9.0 it fully supports grouped and weighted computations.

}
\usage{
descr(X, \dots)

\method{descr}{default}(X, by = NULL, w = NULL, cols = NULL,
      Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
      Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), Qtype = 7L,
      label.attr = "label", stepwise = FALSE, \dots)

\method{descr}{grouped_df}(X, w = NULL,
      Ndistinct = TRUE, higher = TRUE, table = TRUE, sort.table = "freq",
      Qprobs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), Qtype = 7L,
      label.attr = "label", stepwise = FALSE, \dots)

\method{as.data.frame}{descr}(x, \dots, gid = "Group")

\method{print}{descr}(x, n = 14, perc = TRUE, digits = 2, t.table = TRUE, total = TRUE,
      compact = FALSE, summary = !compact, reverse = FALSE, stepwise = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a (grouped) data frame or list of atomic vectors. Atomic vectors, matrices or arrays can be passed but will first be coerced to data frame using \code{\link{qDF}}.
}
  \item{by}{a factor, \code{\link{GRP}} object, or atomic vector / list of vectors (internally grouped with \code{\link{GRP}}), or a one- or two-sided formula e.g. \code{~ group1} or \code{var1 + var2 ~ group1 + group2} to group \code{X}. See Examples.}
    \item{w}{a numeric vector of (non-negative) weights. the default method also supports a one-sided formulas i.e. \code{~ weightcol} or \code{~ log(weightcol)}. The \code{grouped_df} method supports lazy-expressions (same without \code{~}). See Examples.}
    \item{cols}{select columns to describe using column names, indices a logical vector or selector function (e.g. \code{is.numeric}). \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
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
double. Probabilities for quantiles to compute on numeric variables, passed down to \code{\link{.quantile}}. If something non-numeric is passed (i.e. \code{NULL}, \code{FALSE}, \code{NA}, \code{""} etc.), no quantiles are computed.
}
  \item{Qtype}{
integer. Quantile types 5-9 following Hyndman and Fan (1996) who recommended type 8, default 7 as in \code{\link{quantile}}.
}
  \item{label.attr}{
character. The name of a label attribute to display for each variable (if variables are labeled).
}
  \item{\dots}{for \code{descr}: other arguments passed to \code{\link{qsu.default}}. For \code{[.descr}: variable names or indices passed to \code{[.list}. The argument is unused in the \code{print} and \code{as.data.frame} methods.}

\item{x}{an object of class 'descr'.}
\item{n}{integer. The maximum number of table elements to print for categorical variables. If the number of distinct elements is \code{<= n}, the whole table is printed. Otherwise the remaining items are summed into an '... \%s Others' category.}
\item{perc}{logical. \code{TRUE} (default) adds percentages to the frequencies in the table for categorical variables, and, if \code{!is.null(by)}, the percentage of observations in each group.}
\item{digits}{integer. The number of decimals to print in statistics, quantiles and percentage tables.}
\item{t.table}{logical. \code{TRUE} (default) prints a transposed table.}
\item{total}{logical. \code{TRUE} (default) adds a 'Total' column for grouped tables (when using \code{by} argument).}
\item{compact}{logical. \code{TRUE} combines statistics and quantiles to generate a more compact printout. Especially useful with groups (\code{by}).}
\item{summary}{logical. \code{TRUE} (default) computes and displays a summary of the frequencies, if the size of the table for a categorical variable exceeds \code{n}.}
\item{reverse}{logical. \code{TRUE} prints contents in reverse order, starting with the last column, so that the dataset can be analyzed by scrolling up the console after calling \code{descr}.}
\item{stepwise}{logical. \code{TRUE} prints one variable at a time. The user needs to press [enter] to see the printout for the next variable. If called from \code{descr}, the computation is also done one variable at a time, and the finished 'descr' object is returned invisibly. }
\item{gid}{character. Name assigned to the group-id column, when describing data by groups.}
%\item{wsum}{\code{TRUE} sums the weights by groups and adds them as a 'WeightSum' column. Alternatively a name for the column can be supplied. }
%\item{stringsAsFactors}{logical. Make factors from dimension names of 'qsu' array. Same as option to \code{\link{as.data.frame.table}}.}
}
\details{
\code{descr} was heavily inspired by \code{Hmisc::describe}, but is much faster and has more advanced statistical capabilities. It is principally a wrapper around \code{\link{qsu}}, \code{\link{fquantile}} (\code{.quantile}), and \code{\link{fndistinct}} for numeric variables, and computes frequency tables for categorical variables using \code{\link{qtab}}. Date variables are summarized with \code{\link{fnobs}}, \code{\link{fndistinct}} and \code{\link{frange}}.

Since v1.9.0 grouped and weighted computations are fully supported. The use of sampling weights will produce a weighted mean, sd, skewness and kurtosis, and weighted quantiles for numeric data. For categorical data, tables will display the sum of weights instead of the frequencies, and percentage tables as well as the percentage of missing values indicated next to 'Statistics' in print, be relative to the total sum of weights. All this can be done by groups. Grouped (weighted) quantiles are computed using \code{\link{BY}}.

For larger datasets, calling the \code{stepwise} option directly from \code{descr()} is recommended, as precomputing the statistics for all variables before digesting the results can be time consuming.
%\code{\link{qsu}} itself is yet about 10x faster than \code{descr}, and is optimized for grouped, panel data and weighted statistics. It is possible to also compute grouped, panel data and/or weighted statistics with \code{descr} by passing group-ids to \code{g}, panel-ids to \code{pid} or a weight vector to \code{w}. These arguments are handed down to \code{\link{qsu.default}} and only affect the statistics natively computed by \code{qsu}, i.e. passing a weight vector produces a weighted mean, sd, skewness and kurtosis but not weighted quantiles.

The list-object returned from \code{descr} can efficiently be converted to a tidy data frame using the \code{as.data.frame} method. This representation will not include frequency tables computed for categorical variables. %, and the method cannot handle arrays of statistics (applicable when \code{g} or \code{pid} arguments are passed to \code{descr}, in that case \code{as.data.frame.descr} will throw an appropriate error).
}
\value{
A 2-level nested list-based object of class 'descr'. The list has the same size as the dataset, and contains the statistics computed for each variable, which are themselves stored in a list containing the class, the label, the basic statistics and quantiles / tables computed for the variable (in matrix form).

The object has attributes attached providing the 'name' of the dataset, the number of rows in the dataset ('N'), an attribute 'arstat' indicating whether arrays of statistics where generated by passing arguments (e.g. \code{pid}) down to \code{qsu.default}, an attribute 'table' indicating whether \code{table = TRUE} (i.e. the object could contain tables for categorical variables), and attributes 'groups' and/or 'weights' providing a \code{\link{GRP}} object and/or weight vector for grouped and/or weighted data descriptions.
}

\seealso{
\code{\link{qsu}}, \code{\link{qtab}}, \code{\link{fquantile}}, \code{\link{pwcor}}, \link[=summary-statistics]{Summary Statistics}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Simple Use
descr(iris)
descr(wlddev)
descr(GGDC10S)

# Some useful print options (also try stepwise argument)
print(descr(GGDC10S), reverse = TRUE, t.table = FALSE)
# For bigger data consider: descr(big_data, stepwise = TRUE)

# Generating a data frame
as.data.frame(descr(wlddev, table = FALSE))

## Weighted Desciptions
descr(wlddev, w = ~ replace_NA(POP)) # replacing NA's with 0's for fquantile()

## Grouped Desciptions
descr(wlddev, ~ income)
print(descr(wlddev, ~ income), compact = TRUE)

## Grouped & Weighted Desciptions
descr(wlddev, ~ income, w = ~ replace_NA(POP))

## Passing Arguments down to qsu.default: for Panel Data Statistics
descr(iris, pid = iris$Species)
descr(wlddev, pid = wlddev$iso3c)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }% use one of  RShowDoc("KEYWORDS")
