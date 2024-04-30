\name{GRP}
\alias{GRP}
\alias{GRP.GRP}
\alias{GRP.default}
\alias{GRP.factor}
\alias{GRP.qG}
\alias{GRP.pseries}
\alias{GRP.pdata.frame}
\alias{GRP.grouped_df}
\alias{fgroup_by}
\alias{gby}
\alias{group_by_vars}
\alias{fgroup_vars}
\alias{fungroup}
\alias{gsplit}
\alias{greorder}
\alias{is_GRP}
\alias{length.GRP}
\alias{print.GRP}
\alias{plot.GRP}
\alias{GRPnames}
\alias{GRPid}
\alias{GRPN}
\alias{as_factor_GRP}
\title{Fast Grouping / \emph{collapse} Grouping Objects}
\description{
  \code{GRP} performs fast, ordered and unordered, groupings of vectors and data frames (or lists of vectors) using \code{\link{radixorderv}} or \code{\link{group}}. The output is a list-like object of class 'GRP' which can be printed, plotted and used as an efficient input to all of \emph{collapse}'s fast statistical and transformation functions and operators (see macros \code{.FAST_FUN} and \code{.OPERATOR_FUN}), as well as to \code{\link{collap}}, \code{\link{BY}} and \code{\link{TRA}}.

\code{fgroup_by} is similar to \code{dplyr::group_by} but faster and class-agnostic. It creates a grouped data frame with a 'GRP' object attached - for fast dplyr-like programming with \emph{collapse}'s fast functions.

There are also several conversion methods to and from 'GRP' objects. Notable among these is \code{GRP.grouped_df}, which returns a 'GRP' object from a grouped data frame created with \code{dplyr::group_by} or \code{fgroup_by}, and the duo \code{GRP.factor} and \code{as_factor_GRP}.

\code{gsplit} efficiently splits a vector based on a 'GRP' object, and \code{greorder} helps to recombine the results. These are the workhorses behind functions like \code{\link{BY}}, and \code{\link{collap}}, \code{\link{fsummarise}} and \code{\link{fmutate}} when evaluated with base R and user-defined functions.
}
\usage{
GRP(X, \dots)

\method{GRP}{default}(X, by = NULL, sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE,
    return.groups = TRUE, return.order = sort, method = "auto",
    call = TRUE, \dots)

\method{GRP}{factor}(X, \dots, group.sizes = TRUE, drop = FALSE, return.groups = TRUE,
    call = TRUE)

\method{GRP}{qG}(X, \dots, group.sizes = TRUE, return.groups = TRUE, call = TRUE)

\method{GRP}{pseries}(X, effect = 1L, \dots, group.sizes = TRUE, return.groups = TRUE,
    call = TRUE)

\method{GRP}{pdata.frame}(X, effect = 1L, \dots, group.sizes = TRUE, return.groups = TRUE,
    call = TRUE)

\method{GRP}{grouped_df}(X, \dots, return.groups = TRUE, call = TRUE)

# Identify 'GRP' objects
is_GRP(x)

\method{length}{GRP}(x)                          # Length of data being grouped
GRPN(x, expand = TRUE, \dots)        # Group sizes (default: expanded to match data length)
GRPid(x, sort = FALSE, \dots)        # Group id (data length, same as GRP(.)$group.id)
GRPnames(x, force.char = TRUE, sep = ".")  # Group names

as_factor_GRP(x, ordered = FALSE, sep = ".") # 'GRP'-object to (ordered) factor conversion

# Efficiently split a vector using a 'GRP' object
gsplit(x, g, use.g.names = FALSE, \dots)

# Efficiently reorder y = unlist(gsplit(x, g)) such that identical(greorder(y, g), x)
greorder(x, g, \dots)

# Fast, class-agnostic pendant to dplyr::group_by for use with fast functions, see details
fgroup_by(.X, \dots, sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE,
          return.groups = TRUE, return.order = sort, method = "auto")
# Standard-evaluation analogue (slim wrapper around GRP.default(), for programming)
group_by_vars(X, by = NULL, ...)
# Shorthand for fgroup_by
gby(.X, \dots, sort = .op[["sort"]], decreasing = FALSE, na.last = TRUE,
    return.groups = TRUE, return.order = sort, method = "auto")

# Get grouping columns from a grouped data frame created with dplyr::group_by or fgroup_by
fgroup_vars(X, return = "data")

# Ungroup grouped data frame created with dplyr::group_by or fgroup_by
fungroup(X, \dots)

\method{print}{GRP}(x, n = 6, \dots)

\method{plot}{GRP}(x, breaks = "auto", type = "l", horizontal = FALSE, \dots)
}
\arguments{
  \item{X}{a vector, list of columns or data frame (default method), or a suitable object (conversion / extractor methods).}

  \item{.X}{a data frame or list.}

  \item{x, g}{a 'GRP' object. For \code{gsplit/greorder}, \code{x} can be a vector of any type, or \code{NULL} to return the integer indices of the groups. \code{gsplit/greorder/GRPN/GRPid} also support vectors or data frames to be passed to \code{g/x}.}

  \item{by}{if \code{X} is a data frame or list, \code{by} can indicate columns to use for the grouping (by default all columns are used). Columns must be passed using a vector of column names, indices, a one-sided formula i.e. \code{~ col1 + col2}, a logical vector (converted to indices) or a selector function e.g. \code{is_categorical}.}

  \item{sort}{logical. If \code{FALSE}, groups are not ordered but simply grouped in the order of first appearance of unique elements / rows. This often provides a performance gain if the data was not sorted beforehand. See also \code{method}.}
  \item{ordered}{logical. \code{TRUE} adds a class 'ordered' i.e. generates an ordered factor.}

  \item{decreasing}{logical. Should the sort order be increasing or decreasing? Can be a vector of length equal to the number of arguments in \code{X} / \code{by} (argument passed to \code{\link{radixorderv}}).}

  \item{na.last}{logical. If missing values are encountered in grouping vector/columns, assign them to the last group (argument passed to \code{\link{radixorderv}}).}

  \item{return.groups}{logical. Include the unique groups in the created GRP object.}

  \item{return.order}{logical. If \code{sort = TRUE}, include the output from \code{\link{radixorderv}} in the created GRP object. This brings performance improvements in \code{gsplit} (and thus also benefits grouped execution of base R functions). }

  \item{method}{character. The algorithm to use for grouping: either \code{"radix"}, \code{"hash"} or \code{"auto"}. \code{"auto"} will chose \code{"radix"} when \code{sort = TRUE}, yielding ordered grouping via \code{\link{radixorderv}}, and \code{"hash"}-based grouping in first-appearance order via \code{\link{group}} otherwise. It is possibly to put \code{method = "radix"} and \code{sort = FALSE}, which will group character data in first appearance order but sort numeric data (a good hybrid option). \code{method = "hash"} currently does not support any sorting, thus putting \code{sort = TRUE} will simply be ignored.}

  \item{group.sizes}{logical. \code{TRUE} tabulates factor levels using \code{\link{tabulate}} to create a vector of group sizes; \code{FALSE} leaves that slot empty when converting from factors.}

  \item{drop}{logical. \code{TRUE} efficiently drops unused factor levels beforehand using \code{\link{fdroplevels}}.}

  \item{call}{logical. \code{TRUE} calls \code{\link{match.call}} and saves it in the final slot of the GRP object.}

  \item{expand}{logical. \code{TRUE} returns a vector the same length as the data. \code{FALSE} returns the group sizes (computed in first-appearance-order of groups if \code{x} is not already a 'GRP' object). }

  \item{force.char}{logical. Always output group names as character vector, even if a single numeric vector was passed to \code{GRP.default}.}

  \item{sep}{character. The separator passed to \code{\link{paste}} when creating group names from multiple grouping variables by pasting them together.}

 \item{effect}{\emph{plm} / indexed data methods: Select which panel identifier should be used as grouping variable. 1L takes the first variable in the \link[=findex]{index}, 2L the second etc., identifiers can also be passed as a character string. More than one variable can be supplied. }

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

  \item{use.g.names}{logical. \code{TRUE} returns a named list, like \code{\link{split}}. \code{FALSE} is slightly more efficient. }

 %  \item{untibble}{logical. \code{TRUE} also removes classes \code{c("tbl_df", "tbl")} from \code{X}. \code{fgroup_by} attaches an attribute \code{"was.tibble"} indicating if \code{X} was a tibble prior to grouping. The argument thus defaults to \code{TRUE} if this attribute is attached and equal to \code{FALSE}, so that \code{identical(X, X |> fgroup_by(...) |> fungroup())}, regardless of the nature of \code{X}. }

  \item{n}{integer. Number of groups to print out.}

  \item{breaks}{integer. Number of breaks in the histogram of group-sizes.}

  \item{type}{linetype for plot.}

  \item{horizontal}{logical. \code{TRUE} arranges plots next to each other, instead of above each other. \emph{Note} that the size of each group is only plotted for objects with less than 10,000 groups.}

  \item{\dots}{for \code{fgroup_by}: unquoted comma-separated column names, sequences of columns, expressions involving columns, and column names, indices, logical vectors or selector functions. See Examples. For \code{group_by_vars}, \code{gsplit}, \code{greorder}, \code{GRPN} and \code{GRPid}: further arguments passed to \code{GRP} (if \code{g/x} is not already a 'GRP' object). For example the \code{by} argument could be used if a data frame is passed.}
}
\details{
\code{GRP} is a central function in the \emph{collapse} package because it provides, in the form of integer vectors, some key pieces of information to efficiently perform grouped operations at the \code{C/C++} level.

Most statistical function require information about (1) the number of groups (2) an integer group-id indicating which values / rows belong to which group and (3) information about the size of each group. Provided with these, \emph{collapse}'s \link[=fast-statistical-functions]{Fast Statistical Functions} pre-allocate intermediate and result vectors of the right sizes and (in most cases) perform grouped statistical computations in a single pass through the data.

The sorting functionality of \code{GRP.default} lets groups receive different integer-id's depending on whether the groups are sorted \code{sort = TRUE} (\code{FALSE} gives first-appearance order), and in which order (argument \code{decreasing}). This affects the order of values/rows in the output whenever an aggregation is performed. % \emph{Note} that \code{sort = FALSE} is only effective on character vectors, numeric grouping vectors will always produce ordered groupings.

%This in-turn changes the order of values/rows in the output of \emph{collapse} functions (the row/value corresponding to group 1 always comes out on top). The default setting with \code{sort = TRUE} and \code{decreasing = FALSE} results in groups being sorted in ascending order. This is equivalent to performing grouped operations in \emph{data.table} using \code{keyby}, whereas \code{sort = FALSE} is equivalent to \emph{data.table} grouping with \code{by}, however this only works if the \code{by} columns are character, numeric grouping columns are always ordered.

Other elements in the object provide information about whether the data was sorted by the variables defining the grouping (6) and the ordering vector (7). These also feed into optimizations in \code{gsplit/greorder} that benefit the execution of base R functions across groups.

Complimentary to \code{GRP}, the function \code{fgroup_by} is a significantly faster and class-agnostic alternative to \code{dplyr::group_by} for programming with \emph{collapse}. It creates a grouped data frame with a 'GRP' object attached in a \code{"groups"} attribute. This data frame has classes 'GRP_df', \dots, 'grouped_df' and 'data.frame', where \dots stands for any other classes the input frame inherits such as 'data.table', 'sf', 'tbl_df', 'indexed_frame' etc.. \emph{collapse} functions with a 'grouped_df' method respond to 'grouped_df' objects created with either \code{fgroup_by} or \code{dplyr::group_by}. The method \code{GRP.grouped_df} takes the \code{"groups"} attribute from a 'grouped_df' and converts it to a 'GRP' object if created with \code{dplyr::group_by}. %If the grouped data frame was generated using \code{fgroup_by}, all work is done already. If it was created using \code{dplyr::group_by}, a C routine is called to efficiently convert the grouping object.

The 'GRP_df' class in front responds to \code{print.GRP_df} which first calls \code{print(fungroup(x), ...)} and prints one line below the object indicating the grouping variables, followed, in square brackets, by some statistics on the group sizes: \code{[N | Mean (SD) Min-Max]}. The mean is rounded to a full number and the standard deviation (SD) to one digit. Minimum and maximum are only displayed if the SD is non-zero. There also exist a method \code{[.GRP_df} which calls \code{\link{NextMethod}} but makes sure that the grouping information is preserved or dropped depending on the dimensions of the result (subsetting rows or aggregation with \emph{data.table} drops the grouping object).

%Note that \code{fgroup_by} can only be used in combination with \emph{collapse} functions, not with \code{dplyr::summarize} or \code{dplyr::mutate} (the grouping object and method of computing results is different). The converse is not true, you can group data with \code{dplyr::group_by} and then apply \emph{collapse} functions. \code{fgroup_by} is class-agnostic, i.e. the classes of the data frame or list passed are preserved, and all standard methods (like subsetting with \code{`[`} or \code{print} methods) apply to the grouped object.

% Apart from the class 'grouped_df' which is added behind any classes the object might inherit (apart from 'data.frame'), a class 'GRP_df' is added in front. This class responds to a \code{print} method. Both first call the corresponding method for the object and then print / attach the grouping information.

\code{GRP.default} supports vector and list input and will also return 'GRP' objects if passed. There is also a hidden method \code{GRP.GRP} which simply returns grouping objects (no re-grouping functionality is offered).

Apart from \code{GRP.grouped_df} there are several further conversion methods:

The conversion of factors to 'GRP' objects by \code{GRP.factor} involves obtaining the number of groups calling \code{ng <- fnlevels(f)} and then computing the count of each level using \code{\link[=tabulate]{tabulate(f, ng)}}. The integer group-id (2) is already given by the factor itself after removing the levels and class attributes and replacing any missing values with \code{ng + 1L}. The levels are put in a list and moved to position (4) in the 'GRP' object, which is reserved for the unique groups. Finally, a sortedness check \code{!is.unsorted(id)} is run on the group-id to check if the data represented by the factor was sorted (6). \code{GRP.qG} works similarly (see also \code{\link{qG}}), and the 'pseries' and 'pdata.frame' methods simply group one or more factors in the \link[=indexing]{index} (selected using the \code{effect} argument) .

Creating a factor from a 'GRP' object using \code{as_factor_GRP} does not involve any computations, but may involve interacting multiple grouping columns using the \code{paste} function to produce unique factor levels. %  or \code{\link{as.character}} conversions if the grouping column(s) were numeric (which are potentially expensive).



%\emph{Note}: For faster factor generation and a factor-light class 'qG' which avoids the coercion of factor levels to character also see \code{\link{qF}} and \code{\link{qG}}.
}
\value{
  A list-like object of class `GRP' containing information about the number of groups, the observations (rows) belonging to each group, the size of each group, the unique group names / definitions, whether the groups are ordered and data grouped is sorted or not, the ordering vector used to perform the ordering and the group start positions. The object is structured as follows:
  \tabular{lllllll}{\emph{ List-index }  \tab\tab \emph{ Element-name }   \tab\tab \emph{ Content type } \tab\tab \emph{ Content description} \cr\cr

                 [[1]] \tab\tab N.groups   \tab\tab \code{integer(1)} \tab\tab Number of Groups \cr\cr

                 [[2]] \tab\tab group.id \tab\tab \code{integer(NROW(X))} \tab\tab An integer group-identifier \cr\cr

                 [[3]] \tab\tab group.sizes    \tab\tab \code{integer(N.groups)} \tab\tab Vector of group sizes \cr\cr

                 [[4]] \tab\tab groups    \tab\tab \code{unique(X)} or \code{NULL} \tab\tab Unique groups (same format as input, except for \code{fgroup_by} which uses a plain list, sorted if \code{sort = TRUE}), or \code{NULL} if \code{return.groups = FALSE} \cr\cr

                 [[5]] \tab\tab group.vars   \tab\tab \code{character} \tab\tab The names of the grouping variables \cr\cr
                 [[6]] \tab\tab ordered   \tab\tab \code{logical(2)} \tab\tab \code{[1]} Whether the groups are ordered: equal to the \code{sort} argument in the default method, or \code{TRUE} if converted objects inherit a class \code{"ordered"} and \code{NA} otherwise, \code{[2]} Whether the data (\code{X}) is already sorted: the result of \code{!is.unsorted(group.id)}. If \code{sort = FALSE} (default method) the second entry is \code{NA}. \cr\cr

                 [[7]] \tab\tab order     \tab\tab \code{integer(NROW(X))} or \code{NULL} \tab\tab Ordering vector from \code{radixorderv} (with \code{"starts"} attribute), or \code{NULL} if \code{return.order = FALSE} \cr\cr
                 [[8]] \tab\tab group.starts \tab\tab \code{integer(N.groups)} or \code{NULL} \tab\tab The first-occurrence positions/rows of the groups. Useful e.g. with \code{ffirst(x, g, na.rm = FALSE)}. \code{NULL} if \code{return.groups = FALSE}. \cr\cr

                 [[9]] \tab\tab call \tab\tab \code{match.call()} or \code{NULL} \tab\tab The \code{GRP()} call, obtained from \code{match.call()}, or \code{NULL} if \code{call = FALSE}
                 }
}
\seealso{
\code{\link{radixorder}}, \code{\link{group}}, \code{\link{qF}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## default method
GRP(mtcars$cyl)
GRP(mtcars, ~ cyl + vs + am)       # Or GRP(mtcars, c("cyl","vs","am")) or GRP(mtcars, c(2,8:9))
g <- GRP(mtcars, ~ cyl + vs + am)  # Saving the object
print(g)                           # Printing it
plot(g)                            # Plotting it
GRPnames(g)                        # Retain group names
GRPid(g)                           # Retain group id (same as g$group.id), useful inside fmutate()
fsum(mtcars, g)                    # Compute the sum of mtcars, grouped by variables cyl, vs and am
gsplit(mtcars$mpg, g)              # Use the object to split a vector
gsplit(NULL, g)                    # The indices of the groups
identical(mtcars$mpg,              # greorder and unlist undo the effect of gsplit
          greorder(unlist(gsplit(mtcars$mpg, g)), g))

## Convert factor to GRP object and vice-versa
GRP(iris$Species)
as_factor_GRP(g)
\donttest{ % The tidyverse regularly causes havoc to CRAN tests in other packages, therefore this is not tested
## dplyr integration
library(dplyr)
mtcars |> group_by(cyl,vs,am) |> GRP()    # Get GRP object from a dplyr grouped tibble
mtcars |> group_by(cyl,vs,am) |> fmean()  # Grouped mean using dplyr grouping
mtcars |> fgroup_by(cyl,vs,am) |> fmean() # Faster alternative with collapse grouping

mtcars |> fgroup_by(cyl,vs,am)            # Print method for grouped data frame

## Adding a column of group sizes.
mtcars |> fgroup_by(cyl,vs,am) |> fsummarise(Sizes = GRPN())
# Note: can also set_collapse(mask = "n") to use n() instead, see help("collapse-options")
# Other usage modes:
mtcars |> fgroup_by(cyl,vs,am) |> fmutate(Sizes = GRPN())
mtcars |> fmutate(Sizes = GRPN(list(cyl,vs,am)))  # Same thing, slightly more efficient

## Various options for programming and interactive use
fgroup_by(GGDC10S, Variable, Decade = floor(Year / 10) * 10) |> head(3)
fgroup_by(GGDC10S, 1:3, 5) |> head(3)
fgroup_by(GGDC10S, c("Variable", "Country")) |> head(3)
fgroup_by(GGDC10S, is.character) |> head(3)
fgroup_by(GGDC10S, Country:Variable, Year) |> head(3)
fgroup_by(GGDC10S, Country:Region, Var = Variable, Year) |> head(3)

## Note that you can create a grouped data frame without materializing the unique grouping columns
fgroup_by(GGDC10S, Variable, Country, return.groups = FALSE) |> fmutate(across(AGR:SUM, fscale))
fgroup_by(GGDC10S, Variable, Country, return.groups = FALSE) |> fselect(AGR:SUM) |> fmean()

## Note also that setting sort = FALSE on unsorted data can be much faster... if not required...
library(microbenchmark)
microbenchmark(gby(GGDC10S, Variable, Country), gby(GGDC10S, Variable, Country, sort = FALSE))
}
}
\keyword{manip}
