\name{varying}
\alias{varying}
\alias{varying.default}
\alias{varying.matrix}
\alias{varying.data.frame}
\alias{varying.pseries}
\alias{varying.pdata.frame}
\alias{varying.grouped_df}
\alias{varying.sf}
\title{Fast Check of Variation in Data}  % Vectors, Matrix and Data Frame Columns}
\description{
\code{varying} is a generic function that (column-wise) checks for variation in the values of \code{x}, (optionally) within the groups \code{g} (e.g. a panel-identifier).
}
\usage{
varying(x, ...)

\method{varying}{default}(x, g = NULL, any_group = TRUE, use.g.names = TRUE, ...)

\method{varying}{matrix}(x, g = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...)

\method{varying}{data.frame}(x, by = NULL, cols = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...)

# Methods for indexed data / compatibility with plm:

\method{varying}{pseries}(x, effect = 1L, any_group = TRUE, use.g.names = TRUE, ...)

\method{varying}{pdata.frame}(x, effect = 1L, cols = NULL, any_group = TRUE, use.g.names = TRUE,
        drop = TRUE, ...)

# Methods for grouped data frame / compatibility with dplyr:

\method{varying}{grouped_df}(x, any_group = TRUE, use.g.names = FALSE, drop = TRUE,
        keep.group_vars = TRUE, ...)

# Methods for grouped data frame / compatibility with sf:

\method{varying}{sf}(x, by = NULL, cols = NULL, any_group = TRUE, use.g.names = TRUE, drop = TRUE, ...)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{a vector, matrix, data frame, 'indexed_series' ('pseries'), 'indexed_frame' ('pdata.frame') or grouped data frame ('grouped_df'). Data must not be numeric.}

\item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}.}

\item{by}{same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1 + group2} or \code{var1 + var2 ~ group1 + group2}. See Examples}

\item{any_group}{logical. If \code{!is.null(g)}, \code{FALSE} will check and report variation in all groups, whereas the default \code{TRUE} only checks if there is variation within any group. See Examples.}

  \item{cols}{select columns using column names, indices or a function (e.g. \code{is.numeric}). Two-sided formulas passed to \code{by} overwrite \code{cols}.}

\item{use.g.names}{logical. Make group-names and add to the result as names (default method) or row-names (matrix and data frame methods). No row-names are generated for \emph{data.table}'s.}

\item{drop}{\emph{matrix and data.frame methods:} Logical. \code{TRUE} drops dimensions and returns an atomic vector if the result is 1-dimensional.}

\item{effect}{\emph{plm} methods: Select the panel identifier by which variation in the data should be examined. 1L takes the first variable in the \link[=indexing]{index}, 2L the second etc.. Index variables can also be called by name. More than one index variable can be supplied, which will be interacted.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation.}

\item{...}{arguments to be passed to or from other methods.}

}
\details{
Without groups passed to \code{g}, \code{varying} simply checks if there is any variation in the columns of \code{x} and returns \code{TRUE} for each column where this is the case and \code{FALSE} otherwise. A set of data points is defined as varying if it contains at least 2 distinct non-missing values (such that a non-0 standard deviation can be computed on numeric data). \code{varying} checks for variation in both numeric and non-numeric data.

If groups are supplied to \code{g} (or alternatively a \emph{grouped_df} to \code{x}), \code{varying} can operate in one of 2 modes:

\itemize{
\item If \code{any_group = TRUE} (the default), \code{varying} checks each column for variation in any of the groups defined by \code{g}, and returns \code{TRUE} if such within-variation was detected and \code{FALSE} otherwise. Thus only one logical value is returned for each column and the computation on each column is terminated as soon as any variation within any group was found.
\item If \code{any_group = FALSE}, \code{varying} runs through the entire data checking each group for variation and returns, for each column in \code{x}, a logical vector reporting the variation check for all groups. If a group contains only missing values, a \code{NA} is returned for that group.
}

The \emph{sf} method simply ignores the geometry column.
}
\value{
A logical vector or (if \code{!is.null(g)} and \code{any_group = FALSE}), a matrix or data frame of logical vectors indicating whether the data vary (over the dimension supplied by \code{g}).
}

\seealso{
\link[=summary-statistics]{Summary Statistics}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Checks overall variation in all columns
varying(wlddev)

## Checks whether data are time-variant i.e. vary within country
varying(wlddev, ~ country)

## Same as above but done for each country individually, countries without data are coded NA
head(varying(wlddev, ~ country, any_group = FALSE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % use one of  RShowDoc("KEYWORDS")
