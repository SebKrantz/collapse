\name{GRP}
\alias{GRP}
\alias{GRP.default}
\alias{GRP.factor}
\alias{GRP.qG}
\alias{GRP.pseries}
\alias{GRP.pdata.frame}
\alias{GRP.grouped_df}
\alias{fgroup_by}
\alias{fgroup_vars}
\alias{is.GRP}
\alias{print.GRP}
\alias{plot.GRP}
\alias{group_names.GRP}
\alias{as.factor.GRP}
\title{Fast Grouping / \emph{collapse} Grouping Objects}
\description{
  \code{GRP} performs fast, ordered and unordered, groupings of vectors and data.frames (or lists of vectors) using \code{\link{radixorderv}}. The output is a list-like object of class 'GRP' which can be printed, plotted and used as an efficient input to all of \emph{collapse}'s fast functions, operators, as well as \code{\link{collap}}, \code{\link{BY}} and \code{\link{TRA}}.

\code{fgroup_by} is similar to \code{dplyr::group_by} but faster. It creates a 'grouped_df', but with a 'GRP' object attached - for faster dplyr-like programming with \emph{collapse}'s fast functions.

There are also several conversion methods to convert to and from 'GRP' objects. The most important of these is \code{GRP.grouped_df}, which returns a 'GRP' object from a 'grouped_df' created with \code{fgroup_by} or \code{dplyr::group_by}.
}

\usage{
GRP(X, ...)

\method{GRP}{default}(X, by = NULL, sort = TRUE, decreasing = FALSE, na.last = TRUE,
    return.groups = TRUE, return.order = FALSE, ...)

\method{GRP}{factor}(X, ...)
\method{GRP}{qG}(X, ...)
\method{GRP}{pseries}(X, effect = 1L, ...)
\method{GRP}{pdata.frame}(X, effect = 1L, ...)
\method{GRP}{grouped_df}(X, ...)

# Identify, get group names, and convert GRP object to factor
is.GRP(x)
group_names.GRP(x, force.char = TRUE)
as.factor.GRP(x, ordered = FALSE)

# Fast version of dplyr::group_by for use with fast functions, see details
fgroup_by(X, ..., sort = TRUE, decreasing = FALSE, na.last = TRUE, return.order = FALSE)

# This gets grouping columns from a grouped_df created with dplyr::group_by or fgroup_by
fgroup_vars(X, return = "data")

\method{print}{GRP}(x, n = 6, ...)

\method{plot}{GRP}(x, breaks = "auto", type = "s", horizontal = FALSE, ...)
}

\arguments{
  \item{X}{a vector, list of columns or data.frame (default method), or a classed object (conversion/extractor methods).}

  \item{x}{a GRP object.}

  \item{by}{if \code{X} is a data.frame or list, \code{by} can indicate columns to use for the grouping (by default all columns are used). Columns must be passed using a vector of column names, indices, or using a one-sided formula i.e. \code{~ col1 + col2}.}

  \item{sort}{logical. This argument only affects character vectors / columns passed. If \code{FALSE}, these are not ordered but simply grouped in the order of first appearance of unique elements. This provides a slight performance gain if only grouping but not alphabetic ordering is required (argument passed to \code{\link{radixorderv}}).}

  \item{ordered}{logical. \code{TRUE} adds a class 'ordered' i.e. generates an ordered factor.}

  \item{decreasing}{logical. Should the sort order be increasing or decreasing? Can be a vector of length equal to the number of arguments in \code{X} / \code{by} (argument passed to \code{\link{radixorderv}}).}

  \item{na.last}{logical. if missing values are encountered in grouping vector/columns, assign them to the last group (argument passed to \code{\link{radixorderv}}).}

  \item{return.groups}{logical. include the unique groups in the created 'GRP' object.}

  \item{return.order}{logical. include the output from \code{\link{radixorderv}} in the created 'GRP' object.}

  \item{force.char}{logical. Always output group names as character vector, even if a single numeric vector was passed to \code{GRP.default}.}

 \item{effect}{\emph{plm} methods: Select which panel identifier should be used as grouping variable. 1L means first variable in the \code{plm::index}, 2L the second etc.. More than one variable can be supplied. }

   \item{return}{an integer or string specifying what \code{fgroup_vars} should return. The options are:
      \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "data"   \tab\tab full grouping columns (default) \cr
                 2 \tab\tab "unique"   \tab\tab unique rows of grouping columns \cr
                 3 \tab\tab "names" \tab\tab names of grouping columns \cr
                 4 \tab\tab "indices"   \tab\tab integer indices of grouping columns \cr
                 5 \tab\tab "named_indices"  \tab\tab named integer indices of grouping columns \cr
                 6 \tab\tab "logical"   \tab\tab logical selection vector of grouping columns \cr
                 7 \tab\tab "named_logical"     \tab\tab named logical selection vector of grouping columns \cr
      }
  }


  \item{n}{integer. Number of groups to print out.}

  \item{breaks}{integer. Number of breaks in the histogram of group-sizes.}

  \item{type}{linetype for plot.}

  \item{horizontal}{logical. \code{TRUE} arranges plots next to each other, instead of above each other.}

  \item{...}{for \code{fgroup_by}: unquoted comma-separated column names of grouping columns. Otherwise: arguments to be passed to or from other methods.}
}
\details{
\code{GRP} is a central function in the \emph{collapse} package because it provides the key inputs to facilitate easy and efficient groupwise-programming at the \code{C/C++} level: Information about (1) the number of groups (2) an integer group-id indicating which values / rows belong to which group and (3) information about the size of each group. Provided with these informations, \emph{collapse}'s \link[=A1-fast-statistical-functions]{Fast Statistical Functions} pre-allocate intermediate and result vectors of the right sizes and (in most cases) perform grouped statistical computations in a single pass through the data.

The sorting and ordering functionality for \code{GRP} only affects (2), that is groups receive different integer-id's depending on whether the groups are sorted \code{sort = TRUE}, and in which order (argument \code{decreasing}). This in-turn changes the order of values/rows in the output of \emph{collapse} functions. \emph{Note} that \code{sort = FALSE} is only effective on character vectors. Numeric grouping vectors will always produce ordered groupings.

%This in-turn changes the order of values/rows in the output of \emph{collapse} functions (the row/value corresponding to group 1 always comes out on top). The default setting with \code{sort = TRUE} and \code{decreasing = FALSE} results in groups being sorted in ascending order. This is equivalent to performing grouped operations in \emph{data.table} using \code{keyby}, whereas \code{sort = FALSE} is equivalent to \emph{data.table} grouping with \code{by}, however this only works if the \code{by} columns are character, numeric grouping columns are always ordered.

Next to group, there is the function \code{fgroup_by} as a significantly faster alternative to \code{dplyr::group_by}. It creates a grouped tibble by attaching a 'GRP' object to a data frame. \emph{collapse} functions with a grouped_df method applied to that data frame will yield grouped computations. Note that \code{fgroup_by} can only be used in combination with \emph{collapse} functions, not with \emph{dplyr} verbs such as \code{summarize} or \code{mutate}.

\code{GRP} is an S3 generic function with one default method supporting vector and list input and several conversion methods:

The conversion of factors to 'GRP' objects by \code{GRP.factor} involves obtaining the number of groups calling \code{ng <- flevels(f)} and then computing the count of each level using \code{\link[=tabulate]{tabulate(f, ng)}}. The integer group-id (2) is already given by the factor itself after removing the levels and class attributes and replacing any missing values with \code{ng + 1L}. The levels are put in a list and moved to position (4) in the 'GRP' object, which is reserved for the unique groups. Going from factor to 'GRP' object thus only requires a tabulation of the levels, whereas creating a factor from a 'GRP' object using \code{as.factor.GRP} does not involve any computations, but may involve interactions if multiple grouping columns were used (which are then interacted to produce unique factor levels) or \code{\link{as.character}} conversions if the grouping column(s) were numeric (which are potentially expensive).

The method \code{GRP.grouped_df} takes the 'groups' attribute from a grouped tibble and converts it to a 'GRP' object. If the grouped tibble was generated using \code{fgroup_by}, all work is done already. If it was created using \code{dplyr::group_by}, a C++ routine is called to efficiently convert the grouping object.

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
\code{\link{qF}}, \code{\link{qG}}, \code{\link{finteraction}}, \link[=collapse-documentation]{Collapse Overview}
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

## dplyr integration
library(dplyr)
mtcars \%>\% group_by(cyl,vs,am) \%>\% GRP    # get GRP object from a dplyr grouped tibble
mtcars \%>\% group_by(cyl,vs,am) \%>\% fmean  # grouped mean using dplyr grouping
mtcars \%>\% fgroup_by(cyl,vs,am) \%>\% fmean # faster alternative with collapse grouping

}
\keyword{manip}
\keyword{documentation}
