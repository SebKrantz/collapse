\name{GRP}
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
\title{Fast Grouping / \code{collapse} Grouping Objects}
\usage{
\method{GRP}{default}(X, by = if(is.atomic(X)) NULL else seq_along(X), sort = TRUE, order = 1L,
    na.last = FALSE, return.groups = TRUE, return.order = FALSE, ...)

\method{GRP}{factor}(X, ...)
\method{GRP}{qG}(X, ...)
\method{GRP}{pseries}(X, ...)
\method{GRP}{pdata.frame}(X, ...)
\method{GRP}{grouped_df}(X, ...)

\method{group_names}{GRP}(g)
\method{print}{GRP}(g, n = 6, ...)
\method{plot}{GRP}(g, breaks = "auto", type = "s", horizontal = FALSE, ...)
}
\arguments{
  \item{X}{a vector or data.frame (default method), or classed object (conversion methods).}

  \item{g}{a GRP object.}

  \item{by}{if X is a data.frame or list, \code{by} can indicate columns to use for the grouping (by default all columns are used). columns must be passed using a vector of column names or column indices, or using a one-sided formula i.e. \code{~ col1 + col2}.}

  \item{sort}{sort the groups (argument passed to \code{data.table:::forderv})}

  \item{order}{sort the groups in ascending (1L, default) or descending (-1L) order (argument passed to \code{data.table:::forderv})}

  \item{na.last}{if missing values are encountered in grouping vector/columns, put them last (argument passed to \code{data.table:::forderv})}

  \item{return.groups}{include the unique groups in the created GRP object}

  \item{return.order}{include the output from \code{data.table:::forderv} in the created GRP object}
}
\value{
  A list-like object of class 'GRP' containing information about the number of groups, the observations (rows) belonging to each group, the size of each group, the unique group names / definitions, whether the groups are ordered or not and (optionally) the ordering vector used to perform the ordering. The object is simple and structured as follows:
  \tabular{rllll}{\tab \emph{ List-index } \Sexpr{"\u200B"}  \tab \emph{ Element-name } \Sexpr{"\u200B"}  \tab \emph{ Content type } \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \Sexpr{"\u200B"} \tab \emph{ Content description} \cr
                 \tab [[1]] \tab N.groups   \tab integer(1) \tab Number of Groups \cr
                 \tab [[2]] \tab group.id \tab integer(NROW(X)) \tab An integer group-identifier \cr
                 \tab [[3]] \tab group.sizes    \tab integer(N.groups) \tab Vector of group sizes \cr
                 (default) \Sexpr{"\u200B"}  \tab [[4]] \tab groups    \tab unique(X) \tab Unique groups (same format as input, sorted if \code{sort = TRUE}) \cr
                 \tab [[5]] \tab ordered   \tab logical(2) \tab \code{[1]- TRUE} if \code{sort = TRUE}, \code{[2]- TRUE} if \code{X} already sorted \cr
                 (optional) \Sexpr{"\u200B"}  \tab [[6]] \tab order     \tab integer(NROW(X)) \tab Ordering vector from \code{data.table:::forderv} \cr
                 \tab [[7]] \tab call \tab call \tab The \code{GRP()} call, obtained from \code{match.call()}
                 }
}
\description{
  \code{GRP} performs fast, ordered and unordered, groupings of vectors and data.frames (or lists of vectors) using \code{data.table}'s fast grouping and ordering \code{C} routine (\code{forder}). The output is an object of class 'GRP' which can be printed, plotted and used as an efficient input to all of \code{collapse}'s statistical and analysis functions.
}
\examples{
## default method
GRP(mtcars$cyl)
GRP(mtcars, ~ cyl + vs + am)      # alternatively use GRP(mtcars, c("cyl","vs","am")) or GRP(mtcars, c(2,8:9))
g <- GRP(mtcars, ~ cyl + vs + am) # saving the object
plot(g)                           # plotting it
names(g)                          # retain group names
fsum(mtcars, g)                   # compute the sum of mtcars, grouped by variables cyl, vs and am.

## convert factor to GRP object
GRP(iris$Species)

## get GRP object from a grouped tibble
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% GRP

}
