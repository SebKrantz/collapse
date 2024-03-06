\name{fcount}
\alias{fcount}
\alias{fcountv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Efficiently Count Observations by Group
}
\description{
A much faster replacement for \code{dplyr::count}.
}
\usage{
fcount(x, ..., w = NULL, name = "N", add = FALSE,
      sort = FALSE, decreasing = FALSE)

fcountv(x, cols = NULL, w = NULL, name = "N", add = FALSE,
        sort = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a data frame or list-like object, including 'grouped_df' or 'indexed_frame'. Atomic vectors or matrices can also be passed, but will be sent through \code{\link{qDF}}. }
  \item{\dots}{for \code{fcount}: names or sequences of columns to count cases by - passed to \code{\link{fselect}}. For \code{fcountv}: further arguments passed to \code{\link{GRP}} (such as \code{decreasing}, \code{na.last}, \code{method}, \code{effect} etc.). Leaving this empty will count on all columns. }
  \item{cols}{select columns to count cases by, using column names, indices, a logical vector or a selector function (e.g. \code{is_categorical}).}
  \item{w}{a numeric vector of weights, may contain missing values. In \code{fcount} this can also be the (unquoted) name of a column in the data frame. \code{fcountv} also supports a single character name. \emph{Note} that the corresponding argument in \code{dplyr::count} is called \code{wt}, but \emph{collapse} has a global default for weights arguments to be called \code{w}.}
  \item{name}{character. The name of the column containing the count or sum of weights. \code{dplyr::count} it is called \code{"n"}, but \code{"N"} is more consistent with the rest of \emph{collapse} and \emph{data.table}.}
  \item{add}{\code{TRUE} adds the count column to \code{x}. Alternatively \code{add = "group_vars"} (or \code{add = "gv"} for parsimony) can be used to retain only the variables selected for counting in \code{x} and the count.}
  \item{sort, decreasing}{arguments passed to \code{\link{GRP}} affecting the order of rows in the output (if \code{add = FALSE}), and the algorithm used for counting. In general, \code{sort = FALSE} is faster unless data is already sorted by the columns used for counting.
}
}
\value{
If \code{x} is a list, an object of the same type as \code{x} with a column (\code{name}) added at the end giving the count. Otherwise, if \code{x} is atomic, a data frame returned from \code{\link[=qDF]{qDF(x)}} with the count column added. By default (\code{add = FALSE}) only the unique rows of \code{x} of the columns used for counting are returned.
}
\seealso{
\code{\link{GRPN}}, \code{\link{fnobs}}, \code{\link{fndistinct}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
fcount(mtcars, cyl, vs, am)
fcountv(mtcars, cols = .c(cyl, vs, am))
fcount(mtcars, cyl, vs, am, sort = TRUE)
fcount(mtcars, cyl, vs, am, add = TRUE)
fcount(mtcars, cyl, vs, am, add = "group_vars")

## With grouped data
mtcars |> fgroup_by(cyl, vs, am) |> fcount()
mtcars |> fgroup_by(cyl, vs, am) |> fcount(add = TRUE)
mtcars |> fgroup_by(cyl, vs, am) |> fcount(add = "group_vars")

## With indexed data: by default counting on the first index variable
wlddev |> findex_by(country, year) |> fcount()
wlddev |> findex_by(country, year) |> fcount(add = TRUE)
# Use fcountv to pass additional arguments to GRP.pdata.frame,
# here using the effect argument to choose a different index variable
wlddev |> findex_by(country, year) |> fcountv(effect = "year")
wlddev |> findex_by(country, year) |> fcountv(add = "group_vars", effect = "year")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
