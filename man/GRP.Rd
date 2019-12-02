\name{GRP}
\alias{quick-grouping}
\alias{GRP}
\alias{GRP.default}
\alias{GRP.factor}
\alias{GRP.qG}
\alias{GRP.pseries}
\alias{GRP.pdata.frame}
\alias{GRP.grouped_df}
\alias{is.GRP}
% \alias{print.GRP}
% \alias{plot.GRP}
\alias{group_names.GRP}
\alias{as.factor.GRP}
\title{Fast Grouping / \code{collapse} Grouping Objects}
\description{
  \code{GRP} performs fast, ordered and unordered, groupings of vectors and data.frames (or lists of vectors) using \code{data.table}'s fast grouping and ordering \code{C} routine (\code{forder}). The output is a list-like object of class 'GRP' which can be printed, plotted and used as an efficient input to all of \code{collapse}'s statistical functions, operators, and \code{collap} and \code{BY}.
}

\usage{
GRP(X, ...)

\method{GRP}{default}(X, by = NULL, sort = TRUE, order = 1L, na.last = FALSE,
    return.groups = TRUE, return.order = FALSE, ...)

\method{GRP}{factor}(X, ...)
\method{GRP}{qG}(X, ...)
\method{GRP}{pseries}(X, ...)
\method{GRP}{pdata.frame}(X, ...)
\method{GRP}{grouped_df}(X, ...)

group_names.GRP(x, force.char = TRUE)
as.factor.GRP(x)

\method{print}{GRP}(x, n = 6, ...)

\method{plot}{GRP}(x, breaks = "auto", type = "s", horizontal = FALSE, ...)
}

\arguments{
  \item{X}{a vector, list of columns or data.frame (default method), or a classed object (conversion/extractor methods).}

  \item{x}{a GRP object.}

  \item{by}{if \code{X} is a data.frame or list, \code{by} can indicate columns to use for the grouping (by default all columns are used). columns must be passed using a vector of column names, indices, or using a one-sided formula i.e. \code{~ col1 + col2}.}

  \item{sort}{logical. sort the groups (argument passed to \code{data.table:::forderv}, \code{TRUE} is like using \code{keyby} with \code{data.table}, vs. \code{by}).}

  \item{order}{integer. sort the groups in ascending (1L, default) or descending (-1L) order (argument passed to \code{data.table:::forderv}).}

  \item{na.last}{logical. if missing values are encountered in grouping vector/columns, assign them to the last group (argument passed to \code{data.table:::forderv}).}

  \item{return.groups}{logical. include the unique groups in the created 'GRP' object.}

  \item{return.order}{logical. include the output from \code{data.table:::forderv} in the created 'GRP' object.}

  \item{force.char}{logical. Always output group names as character vector, even if a single numeric vector was passed to \code{GRP.default}.}

  \item{...}{arguments to be passed to or from other methods.}
}
\value{
  A list-like object of class 'GRP' containing information about the number of groups, the observations (rows) belonging to each group, the size of each group, the unique group names / definitions, whether the groups are ordered or not and (optionally) the ordering vector used to perform the ordering. The object is simple and structured as follows:
  \tabular{llll}{\emph{ List-index } \Sexpr{"\u200B"}  \tab \emph{ Element-name } \Sexpr{"\u200B"}  \tab \emph{ Content type } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Content description} \cr

                 [[1]] \tab N.groups   \tab \code{integer(1)} \tab Number of Groups \cr

                 [[2]] \tab group.id \tab \code{integer(NROW(X))} \tab An integer group-identifier \cr

                 [[3]] \tab group.sizes    \tab \code{integer(N.groups)} \tab Vector of group sizes \cr

                 [[4]] \tab groups    \tab \code{unique(X)} or \code{NULL} \tab Unique groups (same format as input, sorted if \code{sort = TRUE}), or \code{NULL} if \code{return.groups = FALSE} \cr

                 [[5]] \tab group.vars   \tab \code{character} \tab The names of the grouping variables \cr
                 [[6]] \tab ordered   \tab \code{logical(2)} \tab \code{[1]- TRUE} if \code{sort = TRUE}, \code{[2]- TRUE} if \code{X} already sorted \cr

                 [[7]] \tab order     \tab \code{integer(NROW(X))} or \code{NULL} \tab Ordering vector from \code{data.table:::forderv} or \code{NULL} if \code{return.order = FALSE} (the default) \cr

                 [[8]] \tab call \tab \code{call} \tab The \code{GRP()} call, obtained from \code{match.call()}
                 }
}
\seealso{
\code{\link{qF}}, \code{\link{qG}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default method
GRP(mtcars$cyl)
GRP(mtcars, ~ cyl + vs + am)      # alternatively use GRP(mtcars, c("cyl","vs","am")) or GRP(mtcars, c(2,8:9))
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
