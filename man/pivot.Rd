\name{pivot}
\alias{pivot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast and Easy Data Frame Reshaping
}
\description{
\code{pivot()} will become \emph{collapse}'s data reshaping command, to be formally released with v2.0.0. It combines longer-, wider-, and recast-pivoting functionality in a single parsimonious API. Notably, it can also accommodate variable labels. Please give me feedback on the API.
}
\usage{
pivot(data,               # Summary of Documentation:
      ids = NULL,         # identifier cols to preserve
      values = NULL,      # cols containing the data
      names = NULL,       # name(s) of new col(s) | col(s) containing names
      labels = NULL,      # name of new labels col | col(s) containing labels
      how = "longer",     # method: "longer"/"l", "wider"/"w" or "recast"/"r"
      na.rm = FALSE,      # remove rows missing 'values' in reshaped data
      factor = c("names", "labels"), # create new id col(s) as factor variable(s)?
      check.dups = FALSE, # detect duplicate 'ids'+'names' combinations
      # Options below only apply if how = "wider" or "recast"
      fill = NULL,        # value to insert (default NA/NULL)
      drop = TRUE,        # drop unused levels (=columns) if 'names' is factor
      sort = "",          # "ids": sort 'ids' and/or "names": alphabetic casting
      # This only applies when casting multiple long columns ('values')
      transpose = ""      # "cols": applies t_list() before flattening, and/or
)                         # "names": sets names nami_colj. default: colj_nami
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data frame-like object (list of equal-length columns).}
  \item{ids}{identifier columns to keep. Specified using column names, indices, a logical vector or an identifier function e.g. \code{\link{is_categorical}}.}
  \item{values}{columns containing the data to be reshaped. Specified like \code{ids}.
}
  \item{names}{names of columns to generate, or retrieve column-names from:
        \tabular{lll}{\code{ how }   \tab\tab \emph{ Description }  \cr\cr
                 \code{"longer"}  \tab\tab list of names for the variable and value column in the long format, respectively. If \code{NULL}, the default \code{list("variable", "value")} will be chosen. Alternatively, a named list length 1 or 2 can be provided using "variable"/"value" as keys e.g. \code{list(value = "data_col")}.   \cr\cr
                 \code{ "wider"} \tab\tab column(s) containing names of the new variables. Specified using a vector of column names, indices, a logical vector or selector function e.g. \code{is.character}. Multiple columns will be combined using \code{\link{finteraction}} with \code{"_"} as separator.  \cr\cr
                 \code{ "recast"}   \tab\tab (named) list with the following elements: [[1]]/[["from"]] - column(s) containing names of the new variables, specified as in \code{"wider"}; [[2]]/[["to"]] - name of the variable to generate containing old column names. If \code{NULL}, the default \code{list("variable", "variable")} will be chosen.   \cr
                 }

}
  \item{labels}{
names of columns to generate, or retrieve variable labels from:
        \tabular{lll}{\code{ how }   \tab\tab \emph{ Description }  \cr\cr
                 \code{"longer"}  \tab\tab A string specifying the name of the column to store labels - retrieved from the data using \code{vlabels(values)}. \code{TRUE} will create a column named \code{"label"}. Alternatively, a (named) list with two elements: [[1]]/[["name"]] - the name of the labels column; [[2]]/[["new"]] - a (named) character vector of new labels for the 'variable', 'label' and 'value' columns in the long-format frame. See Examples.   \cr\cr
                 \code{ "wider"} \tab\tab column(s) containing names of the new variables. Specified using a vector of column names, indices, a logical vector or selector function e.g. \code{is.character}. Multiple columns will be combined using \code{\link{finteraction}} with \code{"_"} as separator.  \cr\cr
                 \code{ "recast"}   \tab\tab (named) list with the following elements: [[1]]/[["from"]] - column(s) containing names of the new variables, specified as in \code{"wider"}; [[2]]/[["to"]] - name of the variable to generate containing old column names. If \code{NULL}, the default \code{list("variable", "variable")} will be chosen.   \cr
                 }
}
  \item{how}{
  character. The pivoting method: one of \code{"longer"}, \code{"wider"} or \code{"recast"}. These can be abbreviated by the first letter i.e. \code{"l"/"w"/"r"}.
}
  \item{na.rm}{
logical. \code{TRUE} will remove missing values such that in the reshaped data there is no row missing all data columns - selected through 'values'.
}
  \item{factor}{
character. Whether to generate new 'names' and/or 'labels' columns as factor variables. This is generally recommended as factors are more memory efficient than character vectors and also faster in subsequent filtering and grouping. Internally, this argument is evaluated as \code{factor <- c("names", "labels") \%in\% factor}, so passing anything other than \code{"names"} and/or \code{"labels"} will disable it.
}
  \item{fill}{
%%     ~~Describe \code{fill} here~~
}
  \item{drop}{
%%     ~~Describe \code{drop} here~~
}
  \item{check.dups}{
%%     ~~Describe \code{check.dups} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# -------------------------------- PIVOT LONGER ---------------------------------
# Simple Melting (Reshaping Long)
pivot(mtcars) |> head()
pivot(iris, "Species") |> head()
pivot(wlddev, 1:8, na.rm = TRUE) |> head()
pivot(wlddev, c("iso3c", "year"), c("PCGDP", "LIFEEX"), na.rm = TRUE) |> head()
pivot(GGDC10S, 1:5, na.rm = TRUE) |> head()

# Melting including saving labels
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = TRUE) |> head()
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = "description") |> head()

# Also assigning new labels
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = list("description",
      c("Sector Code", "Sector Description", "Value")) |> namlab()

# Can leave out value column by providing names
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = list("description",
      c(variable = "Sector Code", description = "Sector Description"))) |> namlab()

# This is how I would do it (full code example)
pivot(GGDC10S, 1:5, na.rm = TRUE,
      names = list(variable = "sector_code"),
      labels = list(name = "sector_description",
                    new = c(sector_code = "GGDC10S Sector Code",
                            sector_description = "Long Sector Description",
                            value = "Employment or Value Added"))) |> namlab()

# -------------------------------- PIVOT WIDER ---------------------------------

# -------------------------------- PIVOT RECAST ---------------------------------
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ manip }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
