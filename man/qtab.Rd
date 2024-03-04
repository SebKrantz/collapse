\name{qtab}
\alias{qtab}
\alias{qtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Weighted) Cross Tabulation
}
\description{
A versatile and computationally more efficient replacement for \code{\link{table}}. Notably, it also supports tabulations with frequency weights, and computation of a statistic over combinations of variables.
}
\usage{
qtab(..., w = NULL, wFUN = NULL, wFUN.args = NULL,
     dnn = "auto", sort = .op[["sort"]], na.exclude = TRUE,
     drop = FALSE, method = "auto")

qtable(...) # Long-form. Use set_collapse(mask = "table") to replace table()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
atomic vectors or factors spanning the table dimensions, (optionally) with tags for the dimension names, or a data frame / list of these. See Examples.
}
  \item{w}{
  a single vector to aggregate over the table dimensions e.g. a vector of frequency weights.
}
  \item{wFUN}{
a function used to aggregate \code{w} over the table dimensions. The default \code{NULL} computes the sum of the non-missing weights via an optimized internal algorithm. \link[=fast-statistical-functions]{Fast Statistical Functions} also receive vectorized execution.
}
\item{wFUN.args}{
a list of (optional) further arguments passed to \code{wFUN}. See Examples.
}
\item{dnn}{
the names of the table dimensions. Either passed directly as a character vector or list (internally \code{\link{unlist}}'ed), a function applied to the \code{\dots} list (e.g. \code{\link{names}}, or \code{\link{vlabels}}), or one of the following options:
\itemize{
\item \code{"auto"} constructs names based on the \code{\dots} arguments, or calls \code{\link{names}} if a single list is passed as input.
\item \code{"namlab"} does the same as \code{"auto"}, but also calls \code{\link{vlabels}} on the list and appends the names by the variable labels.
}
\code{dnn = NULL} will return a table without dimension names.
}
  \item{sort, na.exclude, drop, method}{
arguments passed down to \code{\link{qF}}:
\itemize{
\item \code{sort = FALSE} orders table dimensions in first-appearance order of items in the data (can be more efficient if vectors are not factors already). Note that for factors this option will both recast levels in first-appearance order and drop unused levels.
\item \code{na.exclude = FALSE} includes \code{NA}'s in the table (equivalent to \code{\link{table}}'s \code{useNA = "ifany"}).
\item \code{drop = TRUE} removes any unused factor levels (= zero frequency rows or columns).
\item \code{method \%in\% c("radix", "hash")} provides additional control over the algorithm used to convert atomic vectors to factors.
}
}

}

\value{
An array of class 'qtab' that inherits from 'table'. Thus all 'table' methods apply to it.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{descr}}, \link[=summary-statistics]{Summary Statistics}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Basic use
qtab(iris$Species)
with(mtcars, qtab(vs, am))
qtab(mtcars[.c(vs, am)])

library(magrittr)
iris \%$\% qtab(Sepal.Length > mean(Sepal.Length), Species)
iris \%$\% qtab(AMSL = Sepal.Length > mean(Sepal.Length), Species)

## World after 2015
wlda15 <- wlddev |> fsubset(year >= 2015) |> collap(~ iso3c)

# Regions and income levels (country frequency)
wlda15 \%$\% qtab(region, income)
wlda15 \%$\% qtab(region, income, dnn = vlabels)
wlda15 \%$\% qtab(region, income, dnn = "namlab")

# Population (millions)
wlda15 \%$\% qtab(region, income, w = POP) |> divide_by(1e6)

# Life expectancy (years)
wlda15 \%$\% qtab(region, income, w = LIFEEX, wFUN = fmean)

# Life expectancy (years), weighted by population
wlda15 \%$\% qtab(region, income, w = LIFEEX, wFUN = fmean,
                  wFUN.args = list(w = POP))

# GDP per capita (constant 2010 US$): median
wlda15 \%$\% qtab(region, income, w = PCGDP, wFUN = fmedian,
                  wFUN.args = list(na.rm = TRUE))

# GDP per capita (constant 2010 US$): median, weighted by population
wlda15 \%$\% qtab(region, income, w = PCGDP, wFUN = fmedian,
                  wFUN.args = list(w = POP))

# Including OECD membership
tab <- wlda15 \%$\% qtab(region, income, OECD)
tab

# Various 'table' methods
tab |> addmargins()
tab |> marginSums(margin = c("region", "income"))
tab |> proportions()
tab |> proportions(margin = "income")
as.data.frame(tab) |> head(10)
ftable(tab, row.vars = c("region", "OECD"))

# Other options
tab |> fsum(TRA = "\%")    # Percentage table (on a matrix use fsum.default)
tab \%/=\% (sum(tab)/100)    # Another way (division by reference, preserves integers)
tab

rm(tab, wlda15)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{multivariate}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
