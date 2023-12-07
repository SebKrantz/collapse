\name{join}
\alias{join}

\title{Fast Table Joins}

\description{
  Join two data frame like objects \code{x} and \code{y} \code{on} columns. Inspired by \emph{polars} and by default uses a vectorized hash join algorithm (workhorse function \code{\link{fmatch}}).
}

\usage{
join(x, y,
     on = NULL,
     how = "left",
     suffix = NULL,
     validate = "m:m",
     multiple = FALSE,
     sort = FALSE,
     keep.col.order = TRUE,
     drop.dup.cols = FALSE,
     verbose = .op[["verbose"]],
     column = NULL,
     attr = NULL,
     \dots
)
}

\arguments{
  \item{x}{a data frame-like object. The result will inherit the attributes of this object. }

  \item{y}{a data frame-like object to join with \code{x}.}

  \item{on}{character. vector of columns to join on. \code{NULL} uses \code{union(names(x), names(y))}. Use a named vector to match columns named differently in \code{x} and \code{y}, e.g. \code{c("x_id" = "y_id")}.}

  \item{how}{character. Join type: \code{"left"}, \code{"right"}, \code{"inner"}, \code{"full"}, \code{"semi"} or \code{"anti"}. The first letter suffices. }

  \item{suffix}{character(1 or 2). Suffix to add to duplicate column names. \code{NULL} renames duplicate \code{y} columns as \code{paste(col, y_name, sep = "_")}, where \code{y_name = as.character(substitute(y))} i.e. the name of the data frame as passed into the function. In general, passing \code{suffix} length 1 will only rename \code{y}, whereas a length 2 suffix will rename both \code{x} and \code{y}, respectively. If \code{verbose > 0} a message will be printed. }

  \item{validate}{character. (Optional) check if join is of specified type. One of \code{"1:1"}, \code{"1:m"}, \code{"m:1"} or \code{"m:m"}. The default \code{"m:m"} does not perform any checks. Checks are done before the actual join step and failure results in an error. \emph{Note} that this argument does not affect the result, it only triggers a check.}

  \item{multiple}{logical. Handling of rows in \code{x} with multiple matches in \code{y}. The default \code{FALSE} takes the first match in \code{y}. \code{TRUE} returns every match in \code{y}, increasing the size of the joined table.}

  \item{sort}{logical. \code{TRUE} implements a sort-merge-join: a completely separate join algorithm that sorts both datasets on the join columns using \code{\link{radixorder}} and then matches the rows without hashing. \emph{Note} that in this case the result will be sorted by the join columns, whereas \code{sort = FALSE} preserves the order of rows in \code{x}.}

  \item{keep.col.order}{logical. Keep order of columns in \code{x}? \code{FALSE} places the \code{on} columns in front.}

  \item{drop.dup.cols}{instead of renaming duplicate columns in \code{x} and \code{y} using \code{suffix}, this option simply drops them: \code{TRUE} or \code{"y"} drops them from \code{y}, \code{"x"} from \code{x}.}

  \item{verbose}{integer. Prints information about the join. One of 0 (off), 1 (default) or 2 (additionally prints the classes of the \code{on} columns). \emph{Note:} \code{verbose > 0} or \code{validate != "m:m"} invoke the \code{count} argument to \code{\link{fmatch}}, so \code{verbose = 0} is slightly more efficient.}

    \item{column}{(optional) name for an extra column to generate in the output indicating which dataset a record came from. \code{TRUE} calls this column \code{".join"} (inspired by STATA's '_merge' column). By default this column is generated as the last column, but, if \code{keep.col.order = FALSE}, it is placed after the 'on' columns. The column is a factor variable with levels corresponding to the dataset names (inferred from the input) or \code{"matched"} for matched records. Alternatively, it is possible to specify a list of 2, where the first element is the column name, and the second a length 3 (!) vector of levels e.g. \code{column = list("joined", c("x", "y", "x_y"))}, where \code{"x_y"} replaces \code{"matched"}. The column has an additional attribute \code{"on.cols"} giving the join columns corresponding to the factor levels. See Examples. }

  \item{attr}{(optional) name for attribute providing information about the join performed (including the output of \code{\link{fmatch}}) to the result. \code{TRUE} calls this attribute \code{"join.match"}.}

    \item{\dots}{further arguments to \code{\link{fmatch}} (if \code{sort = FALSE}). Notably, \code{overid} can bet set to 0 or 2 (default 1) to control the matching process if the join condition more than identifies the records.}
}


\value{
A data frame-like object of the same type and attributes as \code{x}. \code{"row.names"} of \code{x} are only preserved in left-join operations.
}


\examples{
df1 <- data.frame(
  id1 = c(1, 1, 2, 3),
  id2 = c("a", "b", "b", "c"),
  name = c("John", "Jane", "Bob", "Carl"),
  age = c(35, 28, 42, 50)
)
df2 <- data.frame(
  id1 = c(1, 2, 3, 3),
  id2 = c("a", "b", "c", "e"),
  salary = c(60000, 55000, 70000, 80000),
  dept = c("IT", "Marketing", "Sales", "IT")
)

# Different types of joins
for(i in c("l","i","r","f","s","a"))
    join(df1, df2, how = i) |> print()

# Adding join column: useful esp. for full join
join(df1, df2, how = "f", column = TRUE)
# Custom column + rearranging
join(df1, df2, how = "f", column = list("join", c("x", "y", "x_y")), keep = FALSE)

# Attaching match attribute
str(join(df1, df2, attr = TRUE))

}

\seealso{
\code{\link{fmatch}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}

\keyword{manip}
