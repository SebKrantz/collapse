\name{A2-quick-grouping}
\alias{A2-quick-grouping}
\alias{GRP}
\alias{GRP.default}
\alias{GRP.factor}
\alias{GRP.qG}
\alias{GRP.pseries}
\alias{GRP.pdata.frame}
\alias{GRP.grouped_df}
\alias{is.GRP}
\alias{print.GRP}
\alias{plot.GRP}
\alias{group_names.GRP}
\alias{as.factor.GRP}
\title{Fast Grouping / \code{collapse} Grouping Objects}
\description{
  \code{GRP} performs fast, ordered and unordered, groupings of vectors and data.frames (or lists of vectors) using \code{data.table}'s fast grouping and ordering \code{C} routine (\code{forder}). The output is a list-like object of class 'GRP' which can be printed, plotted and used as an efficient input to all of \code{collapse}'s fast functions, operators, as well as \code{\link{collap}}, \code{\link{BY}} and \code{\link{TRA}}.
}

\usage{
GRP(X, ...)

\method{GRP}{default}(X, by = NULL, sort = TRUE, order = 1L, na.last = TRUE,
    return.groups = TRUE, return.order = FALSE, ...)

\method{GRP}{factor}(X, ...)
\method{GRP}{qG}(X, ...)
\method{GRP}{pseries}(X, ...)
\method{GRP}{pdata.frame}(X, ...)
\method{GRP}{grouped_df}(X, ...)

is.GRP(x)
group_names.GRP(x, force.char = TRUE)
as.factor.GRP(x)

\method{print}{GRP}(x, n = 6, ...)

\method{plot}{GRP}(x, breaks = "auto", type = "s", horizontal = FALSE, ...)
}

\arguments{
  \item{X}{a vector, list of columns or data.frame (default method), or a classed object (conversion/extractor methods).}

  \item{x}{a GRP object.}

  \item{by}{if \code{X} is a data.frame or list, \code{by} can indicate columns to use for the grouping (by default all columns are used). Columns must be passed using a vector of column names, indices, or using a one-sided formula i.e. \code{~ col1 + col2}.}

  \item{sort}{logical. sort the groups (argument passed to \code{data.table:::forderv}, \code{TRUE} is like using \code{keyby} with \code{data.table}, vs. \code{by}).}

  \item{order}{integer. sort the groups in ascending (1L, default) or descending (-1L) order (argument passed to \code{data.table:::forderv}).}

  \item{na.last}{logical. if missing values are encountered in grouping vector/columns, assign them to the last group (argument passed to \code{data.table:::forderv}).}

  \item{return.groups}{logical. include the unique groups in the created 'GRP' object.}

  \item{return.order}{logical. include the output from \code{data.table:::forderv} in the created 'GRP' object.}

  \item{force.char}{logical. Always output group names as character vector, even if a single numeric vector was passed to \code{GRP.default}.}

  \item{n}{integer. Number of groups to print out.}

  \item{breaks}{integer. Number of breaks in the histogram of group-sizes.}

  \item{type}{linetype for plot.}

  \item{horizontal}{logical. \code{TRUE} arranges plots next to each other, instead of above each other.}

  \item{...}{arguments to be passed to or from other methods.}
}
\details{
\code{GRP} is a central function in the \code{collapse} package because it provides the key inputs to facilitate easy and efficient groupwise-programming at the \code{C/C++} level: Information about (1) the number of groups (2) an integer group-id indicating which values / rows belong to which group and (3) information about the size of each group. Provided with these informations, \code{collapse}'s \link[=A1-fast-statistical-functions]{Fast Statistical Functions} pre-allocate intermediate and result vectors of the right sizes and (in most cases) perform grouped statistical computations in a single pass through the data.

The sorting and ordering functionality for \code{GRP} only affects (2), that is groups receive different integer-id's depending on whether the groups are sorted \code{sort = TRUE}, and in which order (\code{order = 1} ascending or \code{order = -1} descending). This in-turn changes the order of values/rows in the output of \code{collapse} functions (the row/value corresponding to group 1 always comes out on top). The default setting with \code{sort = TRUE} and \code{order = 1} results in groups being sorted in ascending order. This is equivalent to performing grouped operations in \code{data.table} using \code{keyby}, whereas \code{sort = FALSE} is equivalent to \code{data.table} grouping with \code{by}.

Evidently \code{GRP} is an S3 generic function with one default method supporting vector and list input and several conversion methods. The most important of these is the conversion of factors to 'GRP' objects and vice-versa. To obtain a 'GRP' object from a factor, one simply gets the number of groups calling \code{ng <- length(levels(f))} (1) and then computes the count of each level (3) using \code{\link[=tabulate]{tabulate(f, ng)}}. The integer group-id (2) is already given by the factor itself after removing the levels and class attributes. The levels are put in a list and moved to position (4) in the 'GRP' object, which is reserved for the unique groups. Going from factor to 'GRP' object thus only requires a tabulation of the levels, whereas creating a factor from a 'GRP' object using \code{as.factor.GRP} does not involve any computations, but may involve interactions if multiple grouping columns were used (which are then interacted to produce unique factor levels) or \code{\link{as.character}} conversions if the grouping column(s) were numeric (which are potentially expensive).

\emph{Note}: For faster factor generation and a factor-light class 'qG' which avoids the coercion of factor levels to character also see \code{\link{qF}} and \code{\link{qG}}.
}
\value{
  A list-like object of class `GRP' containing information about the number of groups, the observations (rows) belonging to each group, the size of each group, the unique group names / definitions, whether the groups are ordered or not and (optionally) the ordering vector used to perform the ordering. The object is structured as follows:
  \tabular{lllllll}{\emph{ List-index }  \tab\tab \emph{ Element-name }   \tab\tab \emph{ Content type } \tab\tab \emph{ Content description} \cr

                 [[1]] \tab\tab N.groups   \tab\tab \code{integer(1)} \tab\tab Number of Groups \cr

                 [[2]] \tab\tab group.id \tab\tab \code{integer(NROW(X))} \tab\tab An integer group-identifier \cr

                 [[3]] \tab\tab group.sizes    \tab\tab \code{integer(N.groups)} \tab\tab Vector of group sizes \cr

                 [[4]] \tab\tab groups    \tab\tab \code{unique(X)} or \code{NULL} \tab\tab Unique groups (same format as input, sorted if \code{sort = TRUE}), or \code{NULL} if \code{return.groups = FALSE} \cr

                 [[5]] \tab\tab group.vars   \tab\tab \code{character} \tab\tab The names of the grouping variables \cr
                 [[6]] \tab\tab ordered   \tab\tab \code{logical(2)} \tab\tab \code{[1]- TRUE} if \code{sort = TRUE}, \code{[2]- TRUE} if \code{X} already sorted \cr

                 [[7]] \tab\tab order     \tab\tab \code{integer(NROW(X))} or \code{NULL} \tab\tab Ordering vector from \code{data.table:::forderv} or \code{NULL} if \code{return.order = FALSE} (the default) \cr

                 [[8]] \tab\tab call \tab\tab \code{call} \tab\tab The \code{GRP()} call, obtained from \code{match.call()}
                 }
}
\seealso{
\code{\link{qF}}, \code{\link{qG}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default method
GRP(mtcars$cyl)
GRP(mtcars, ~ cyl + vs + am)      # or GRP(mtcars, c("cyl","vs","am")) or GRP(mtcars, c(2,8:9))
g <- GRP(mtcars, ~ cyl + vs + am) # saving the object
plot(g)                           # plotting it
group_names.GRP(g)                # retain group names
fsum(mtcars, g)                   # compute the sum of mtcars, grouped by variables cyl, vs and am.

## convert factor to GRP object
GRP(iris$Species)

## get GRP object from a dplyr grouped tibble
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% GRP

}
\keyword{manip}
\keyword{documentation}
