\name{pivot}
\alias{pivot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast and Easy Data Reshaping
}
\description{
\code{pivot()} is \emph{collapse}'s data reshaping command. It combines longer-, wider-, and recast-pivoting functionality in a single parsimonious API. Notably, it can also accommodate variable labels.
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

      # Only apply if how = "wider" or "recast"
      nthreads = .op[["nthreads"]], # minor gains as grouping remains serial
      fill = NULL,        # value to insert for unbalanced data (default NA/NULL)
      drop = TRUE,        # drop unused levels (=columns) if 'names' is factor
      sort = FALSE,       # "ids": sort 'ids' and/or "names": alphabetic casting

      # Only applies if how = "wider" with multiple long columns ('values')
      transpose = FALSE   # "columns": applies t_list() before flattening, and/or
)                         # "names": sets names nami_colj. default: colj_nami
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{data frame-like object (list of equal-length columns).}
  \item{ids}{identifier columns to keep. Specified using column names, indices, a logical vector or an identifier function e.g. \code{\link{is_categorical}}.}
  \item{values}{columns containing the data to be reshaped. Specified like \code{ids}.
}
  \item{names}{names of columns to generate, or retrieve variable names from:
        \tabular{lll}{\code{ how }   \tab\tab \emph{ Description }  \cr\cr
                 \code{"longer"}  \tab\tab list of names for the variable and value column in the long format, respectively. If \code{NULL}, \code{list("variable", "value")} will be chosen. Alternatively, a named list length 1 or 2 can be provided using "variable"/"value" as keys e.g. \code{list(value = "data_col")}.   \cr\cr
                 \code{ "wider"} \tab\tab column(s) containing names of the new variables. Specified using a vector of column names, indices, a logical vector or selector function e.g. \code{is.character}. Multiple columns will be combined using \code{\link{finteraction}} with \code{"_"} as separator.  \cr\cr
                 \code{ "recast"}   \tab\tab (named) list with the following elements: [[1]]/[["from"]] - column(s) containing names of the new variables, specified as in \code{"wider"}; [[2]]/[["to"]] - name of the variable to generate containing old column names. If \code{NULL}, \code{list("variable", "variable")} will be chosen.   \cr
                 }

}
  \item{labels}{
names of columns to generate, or retrieve variable labels from:
        \tabular{lll}{\code{ how }   \tab\tab \emph{ Description }  \cr\cr
                 \code{"longer"}  \tab\tab A string specifying the name of the column to store labels - retrieved from the data using \code{vlabels(values)}. \code{TRUE} will create a column named \code{"label"}. Alternatively, a (named) list with two elements: [[1]]/[["name"]] - the name of the labels column; [[2]]/[["new"]] - a (named) character vector of new labels for the 'variable', 'label' and 'value' columns in the long-format frame. See Examples.   \cr\cr
                 \code{ "wider"} \tab\tab column(s) containing labels of the new variables. Specified using a vector of column names, indices, a logical vector or selector function e.g. \code{is.character}. Multiple columns will be combined using \code{\link{finteraction}} with \code{" - "} as separator.  \cr\cr
                 \code{ "recast"}   \tab\tab (named) list with the following elements: [[1]]/[["from"]] - column(s) containing labels for the new variables, specified as in \code{"wider"}; [[2]]/[["to"]] - name of the variable to generate containing old labels; [[3]]/[["new"]] - a (named) character vector of new labels for the generated 'variable' and 'label' columns. If [[1]]/[["from"]] is not supplied, this can also include labels for new variables. Omitting one of the elements via a named list or setting it to \code{NULL} in a list of 3 will omit the corresponding operation i.e. either not saving existing labels or not assigning new ones.  \cr
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
  \item{check.dups}{
logical. \code{TRUE} checks for duplicate 'ids'+'names' combinations, and, if 'labels' are specified, also for duplicate 'names'+'labels' combinations. The default \code{FALSE} implies that the algorithm just runs through the data, leading effectively to the last value being chosen i.e. the data is aggregated using the last value. See Details.
}
\item{nthreads}{
  integer. if \code{how = "wider"|"recast"}: number of threads to use with OpenMP (default \code{get_collapse("nthreads")}, initialized to 1). Only the distribution of values to columns with \code{how = "wider"|"recast"} is multithreaded here. Since grouping id columns on a long data frame is expensive and serial, the gains are minor. With \code{how = "long"}, multithreading does not make much sense as the most expensive operation is allocating the long results vectors. The rest is a couple of \code{memset()}'s in C to copy the values.
}
  \item{fill}{if \code{how = "wider"|"recast"}: value to insert for 'ids'-'names' combinations not present in the long format. \code{NULL} uses \code{NA} for atomic vectors and \code{NULL} for lists.
}
  \item{drop}{
logical. if \code{how = "wider"|"recast"} and 'names' is a single factor variable: \code{TRUE} will check for and drop unused levels in that factor, avoiding the generation of empty columns.
}
\item{sort}{
  if \code{how = "wider"|"recast"}: specifying \code{"ids"} applies ordered grouping on the id-columns, returning data sorted by ids. Specifying \code{"names"} sorts the names before casting (unless 'names' is a factor), yielding columns cast in alphabetic order. Both options can be passed as a character vector, or, alternatively, \code{TRUE} can be used to enable both.
}
\item{transpose}{
  if \code{how = "wider"|"recast"} and multiple columns are selected through 'values': specifying \code{"columns"} applies \code{\link{t_list}} to the result before flattening, resulting in a different column order. Specifying \code{"names"} generates names of the form nami_colj, instead of colj_nami. Both options can be passed as a character vector, or, alternatively, \code{TRUE} can be used to enable both.
}

}
\details{
Pivot wider essentially works as follows: compute \code{g_rows = group(ids)} and also \code{g_cols = group(names)} (using \code{\link{group}} if \code{sort = FALSE}). \code{g_rows} gives the row-numbers of the wider data frame and \code{g_cols} the column numbers. Then, a C function generates a wide data frame and runs through each long column ('values'), assigning each value to the corresponding row and column in the wide frame. If there are duplicates, some observations will be assigned to the same slot the new data frame, overwriting the value already set. Since the algorithm runs through the data from first to last row, this amounts to removing duplicates by taking the last value for each set of 'ids' - in first-appearance order. If \code{check.dups = TRUE}, a check of the form  \code{fnunique(list(g_rows, g_cols)) < fnrow(data)} is run, and an informative warning is issued if duplicates are found. Recast pivoting works similarly, and in long pivoting the check simply amounts to \code{fnunique(ids) < fnrow(data)}. Making this check optional ensures greater performance, but also requires the user to exercise discretion i.e. know your data or invoke the check.
}
\value{
A reshaped data frame with the same class and attributes (except for 'names'/'row-names') as the input frame.
}

\note{
Leaving either 'ids' or 'values' empty will assign all other columns (except for \code{"variable"} if \code{how = "wider"|"recast"}) to the non-specified argument. It is also possible to leave both empty, e.g. for complete melting if \code{how = "wider"} or data transposition if \code{how = "recast"} (similar to \code{data.table::transpose} but supporting multiple names columns and variable labels). See Examples.

\code{pivot} currently does not support concurrently melting/pivoting longer to multiple columns. See \code{data.table::melt} or \code{pivot_longer} from \emph{tidyr} or \emph{tidytable} for an efficient alternative with this feature. It is also possible to achieve this with just a little bit of programming. An example is provided below. Currently I don't see a 'complete' (like \code{pivot_longer}) way of including this feature in the API, and in general I don't see this as a very well-defined operation. But I am open to suggestions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{collap}}, \code{\link{vec}}, \code{\link{rowbind}}, \code{\link{unlist2d}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# -------------------------------- PIVOT LONGER ---------------------------------
# Simple Melting (Reshaping Long)
pivot(mtcars) |> head()
pivot(iris, "Species") |> head()
pivot(iris, values = 1:4) |> head() # Same thing

# Using collapse's datasets
head(wlddev)
pivot(wlddev, 1:8, na.rm = TRUE) |> head()
pivot(wlddev, c("iso3c", "year"), c("PCGDP", "LIFEEX"), na.rm = TRUE) |> head()
head(GGDC10S)
pivot(GGDC10S, 1:5, names = list("Sectorcode", "Value"), na.rm = TRUE) |> head()
# Can also set by name: variable and/or value. Note that 'value' here remains lowercase
pivot(GGDC10S, 1:5, names = list(variable = "Sectorcode"), na.rm = TRUE) |> head()

# Melting including saving labels
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = TRUE) |> head()
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = "description") |> head()

# Also assigning new labels
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = list("description",
            c("Sector Code", "Sector Description", "Value"))) |> namlab()

# Can leave out value column by providing named vector of labels
pivot(GGDC10S, 1:5, na.rm = TRUE, labels = list("description",
          c(variable = "Sector Code", description = "Sector Description"))) |> namlab()

# Now here is a nice example that is explicit and respects the dataset naming conventions
pivot(GGDC10S, ids = 1:5, na.rm = TRUE,
      names = list(variable = "Sectorcode",
                   value = "Value"),
      labels = list(name = "Sector",
                    new = c(Sectorcode = "GGDC10S Sector Code",
                            Sector = "Long Sector Description",
                            Value = "Employment or Value Added"))) |>
  namlab(N = TRUE, Nd = TRUE, class = TRUE)

# Note that pivot() currently does not support melting to multiple columns
# But you can tackle the issue with a bit of programming:
wide <- pivot(GGDC10S, c("Country", "Year"), c("AGR", "MAN", "SUM"), "Variable",
              how = "wider", na.rm = TRUE)
head(wide)
library(magrittr)
wide \%>\% {av(pivot(., 1:2, grep("_VA", names(.))), pivot(gvr(., "_EMP")))} \%>\% head()
wide \%>\% {av(av(gv(., 1:2), rm_stub(gvr(., "_VA"), "_VA", pre = FALSE)) \%>\%
                  pivot(1:2, names = list("Sectorcode", "VA"), labels = "Sector"),
             EMP = vec(gvr(., "_EMP")))} \%>\% head()
rm(wide)

# -------------------------------- PIVOT WIDER ---------------------------------
iris_long <- pivot(iris, "Species") # Getting a long frame
head(iris_long)
# If 'names'/'values' not supplied, searches for 'variable' and 'value' columns
pivot(iris_long, how = "wider")
# But here the records are not identified by 'Species': thus aggregation with last value:
pivot(iris_long, how = "wider", check = TRUE) # issues a warning
rm(iris_long)

# This works better, these two are inverse operations
wlddev |> pivot(1:8) |> pivot(how = "w") |> head()
# ...but not perfect, we loose labels
namlab(wlddev)
wlddev |> pivot(1:8) |> pivot(how = "w") |> namlab()
# But pivot() supports labels: these are perfect inverse operations
wlddev |> pivot(1:8, labels = "label") |> print(max = 50) |> # Notice the "label" column
  pivot(how = "w", labels = "label") |> namlab()

# If the data does not have 'variable'/'value' cols: need to specify 'names'/'values'
# Using a single column:
pivot(GGDC10S, c("Country", "Year"), "SUM", "Variable", how = "w") |> head()
SUM_wide <- pivot(GGDC10S, c("Country", "Year"), "SUM", "Variable", how = "w", na.rm = TRUE)
head(SUM_wide) # na.rm = TRUE here removes all new rows completely missing data
tail(SUM_wide) # But there may still be NA's, notice the NA in the final row
# We could use fill to set another value
pivot(GGDC10S, c("Country", "Year"), "SUM", "Variable", how = "w",
      na.rm = TRUE, fill = -9999) |> tail()
# This will keep the label of "SUM", unless we supply a column with new labels
namlab(SUM_wide)
# Such a column is not available here, but we could use "Variable" twice
pivot(GGDC10S, c("Country", "Year"), "SUM", "Variable", "Variable", how = "w",
      na.rm = TRUE) |> namlab()
# Alternatively, can of course relabel ex-post
SUM_wide |> relabel(VA = "Value Added", EMP = "Employment") |> namlab()
rm(SUM_wide)

# Multiple-column pivots
pivot(GGDC10S, c("Country", "Year"), c("AGR", "MAN", "SUM"), "Variable", how = "w",
      na.rm = TRUE) |> head()
# Here we may prefer a transposed column order
pivot(GGDC10S, c("Country", "Year"), c("AGR", "MAN", "SUM"), "Variable", how = "w",
      na.rm = TRUE, transpose = "columns") |> head()
# Can also flip the order of names (independently of columns)
pivot(GGDC10S, c("Country", "Year"), c("AGR", "MAN", "SUM"), "Variable", how = "w",
      na.rm = TRUE, transpose = "names") |> head()
# Can also enable both (complete transposition)
pivot(GGDC10S, c("Country", "Year"), c("AGR", "MAN", "SUM"), "Variable", how = "w",
      na.rm = TRUE, transpose = TRUE) |> head() # or tranpose = c("columns", "names")

# Finally, here is a nice, simple way to reshape the entire dataset.
pivot(GGDC10S, values = 6:16, names = "Variable", na.rm = TRUE, how = "w") |>
  namlab(N = TRUE, Nd = TRUE, class = TRUE)

# -------------------------------- PIVOT RECAST ---------------------------------
# Look at the data again
head(GGDC10S)
# Let's stack the sectors and instead create variable columns
pivot(GGDC10S, .c(Country, Regioncode, Region, Year),
      names = list("Variable", "Sectorcode"), how = "r") |> head()
# Same thing (a bit easier)
pivot(GGDC10S, values = 6:16, names = list("Variable", "Sectorcode"), how = "r") |> head()
# Removing missing values
pivot(GGDC10S, values = 6:16, names = list("Variable", "Sectorcode"), how = "r",
      na.rm = TRUE) |> head()
# Saving Labels
pivot(GGDC10S, values = 6:16, names = list("Variable", "Sectorcode"),
      labels = list(to = "Sector"), how = "r", na.rm = TRUE) |> head()

# Supplying new labels for generated columns: as complete as it gets
pivot(GGDC10S, values = 6:16, names = list("Variable", "Sectorcode"),
      labels = list(to = "Sector",
                    new = c(Sectorcode = "GGDC10S Sector Code",
                            Sector = "Long Sector Description",
                            VA = "Value Added",
                            EMP = "Employment")), how = "r", na.rm = TRUE) |>
  namlab(N = TRUE, Nd = TRUE, class = TRUE)

# Now another (slightly unconventional) use case here is data transposition
# Let's get the data for Botswana
BWA <- GGDC10S |> fsubset(Country == "BWA", Variable, Year, AGR:SUM)
head(BWA)
# By supplying no ids or values, we are simply requesting a transpose operation
pivot(BWA, names = list(from = c("Variable", "Year"), to = "Sectorcode"), how = "r")
# Same with labels
pivot(BWA, names = list(from = c("Variable", "Year"), to = "Sectorcode"),
      labels = list(to = "Sector"), how = "r")
# For simple cases, data.table::transpose() will be more efficient, but with multiple
# columns to generate names and/or variable labels to be saved/assigned, pivot() is handy
rm(BWA)
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
