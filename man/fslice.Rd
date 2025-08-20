\name{fslice}
\alias{fslice}
\alias{fslicev}
\title{
Fast Slicing of Matrix-Like Objects
}
\description{
A fast function to extract rows from a matrix or data frame-like object (by groups).
}
\usage{
fslice(x, ..., n = 1, how = "first", order.by = NULL,
       na.rm = .op[["na.rm"]], sort = FALSE, with.ties = FALSE)

fslicev(x, cols = NULL, n = 1, how = "first", order.by = NULL,
        na.rm = .op[["na.rm"]], sort = FALSE, with.ties = FALSE, ...)
}
\arguments{
  \item{x}{a matrix, data frame or list-like object, including 'grouped_df'.}
  \item{\dots}{for \code{fslice}: names or sequences of columns to group by - passed to \code{\link{fselect}}. If \code{x} is a matrix: atomic vectors to group \code{x}. Can be empty to operate on (un)grouped data. For \code{fslicev}: further arguments passed to \code{\link{GRP}} (such as \code{decreasing}, \code{na.last}, \code{method}). }
  \item{cols}{select columns to group by, using column names, indices, a logical vector or a selector function (e.g. \code{is_categorical}). It can also be a list of vectors, or, if \code{x} is a matrix, a single vector.}
  \item{n}{integer or proportion (if < 1). Number of rows to select from each group. If a proportion is provided, it is converted to the equivalent number of rows using \code{max(1, round(n * nrow(x)))} or \code{max(1, round(n * nrow(x) / N.groups))} for grouped data.}
  \item{how}{character. Method to select rows. One of:
    \itemize{
      \item \code{"first"}: select first \code{n} rows
      \item \code{"last"}: select last \code{n} rows
      \item \code{"min"}: select \code{n} rows with minimum values of \code{order.by}
      \item \code{"max"}: select \code{n} rows with maximum values of \code{order.by}
    }
  }
  \item{order.by}{vector or column name to order by when \code{how} is \code{"min"} or \code{"max"}. Must be same length as rows in \code{x}. In \code{fslice} it must not be quoted.}
  \item{na.rm}{logical. If \code{TRUE}, missing values in \code{order.by} are removed before selecting rows.}
  \item{sort}{logical. If \code{TRUE}, sort selected rows on the grouping columns. \code{FALSE} uses first-appearance order (including grouping columns if \code{how} is \code{"first"} or \code{"last"}) - fastest.}
  \item{with.ties}{logical. If \code{TRUE} and \code{how} is \code{"min"} or \code{"max"}, returns all rows with the extreme value. Currently only supported for \code{n = 1} and \code{sort = FALSE}.}
}
\value{
A subset of \code{x} containing the selected rows.
}
\seealso{
\code{\link{fsubset}}, \code{\link{fcount}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# Basic usage
fslice(mtcars, n = 3)                    # First 3 rows
fslice(mtcars, n = 3, how = "last")      # Last 3 rows
fslice(mtcars, n = 0.1)                  # First 10\% of rows

# Using order.by
fslice(mtcars, n = 3, how = "min", order.by = mpg)  # 3 cars with lowest mpg
fslice(mtcars, n = 3, how = "max", order.by = mpg)  # 3 cars with highest mpg

# With grouping
mtcars |> fslice(cyl, n = 2)                        # First 2 cars per cylinder
mtcars |> fslice(cyl, n = 2, sort = TRUE)           # with sorting (slightly less efficient)
mtcars |> fslice(cyl, n = 2, how = "min", order.by = mpg)  # 2 lowest mpg cars per cylinder

# Using with.ties
mtcars |> fslice(cyl, n = 1, how = "min", order.by = mpg, with.ties = TRUE)

# With grouped data
mtcars |>
  fgroup_by(cyl) |>
  fslice(n = 2, how = "max", order.by = mpg)        # 2 highest mpg cars per cylinder
}
\keyword{manip}
