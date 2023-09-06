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
     keep.col.order = TRUE,
     verbose = 1L,
     overid = 1L,
     column = NULL,
     match.attr = NULL
)
}

\arguments{
  \item{x}{a data frame-like object. The result will inherit the attributes of this object. }

  \item{y}{a data frame-like object to join with \code{x}.}

  \item{on}{character. vector of columns to join on. \code{NULL} uses \code{union(names(x), names(y))}. Use a named vector to match columns named differently in \code{x} and \code{y}, e.g. \code{c("x_id" = "y_id")}.}

  \item{how}{charcater. Join type: \code{"left"}, \code{"right"}, \code{"inner"}, \code{"full"}, \code{"semi"} or \code{"anti"}.}

  \item{suffix}{character(1 or 2). Suffix to add to duplicate column names. \code{NULL} renames duplicate \code{y} columns as \code{paste(col, y_name, sep = "_")}, where \code{y_name = as.character(substitute(y))} i.e. the name of the data frame as passed into the function. In general, passing \code{suffix} length 1 will only rename \code{y}, whereas a length 2 suffix will rename both \code{x} and \code{y}, respectively. If \code{verbose > 0} a message will be printed. }

  \item{validate}{character. (Optional) check if join is of specified type. One of \code{"1:1"}, \code{"1:m"}, \code{"m:1"} or \code{"m:m"}. The default \code{"m:m"} does not perform any checks, first matches in \code{x} and \code{y} are taken. Checks are done before the actual join step and failure results in an error.}

  \item{keep.col.order}{logical. Keep order of columns in \code{x}? \code{FALSE} places the \code{on} columns in front.}

  \item{verbose}{integer. Prints information about the join. One of 0 (off), 1 (default) or 2 (additionally prints the classes of the \code{on} columns). \emph{Note:} \code{verbose > 0} or \code{validate != "m:m"} invoke the \code{count} argument to \code{\link{fmatch}}, so \code{verbose = 0} is slightly more efficient.}

  \item{overid}{integer. One of 0, 1 or 2: passed to \code{\link{fmatch}} and allows checking of (1), and (optional) early termination on (0), overidentified join conditions. Default (1) warns about this but proceeds.}

    \item{column}{(optional) name for an extra column to generate in the output indicating which dataset a record came from. \code{TRUE} calls this column \code{".join"} (inspired by STATA's '_merge' column). By default this column is generated as the last column, but, if \code{keep.col.order = FALSE}, it is placed after the 'on' columns. The column is a factor variable with levels corresponding to the dataset names (inferred from the input) or \code{"matched"} for matched records. The \code{suffix} argument (of length 2) can be used to set specific level-names for \code{x} and \code{y}, respectively.}

  \item{match.attr}{(optional) name for attribute to attach the output of \code{\link{fmatch}} to the result. \code{TRUE} calls this attribute \code{"join.match"}.}
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
}

\seealso{
\code{\link{fmatch}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}

\keyword{manip}
