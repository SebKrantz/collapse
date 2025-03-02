\name{fslice}
\alias{fslice}
\title{
Fast Slicing of Matrix-Like Objects
}
\description{
A fast function to extract rows from a matrix or data frame, with support for grouped operations.
}
\usage{
fslice(x, ..., n = 1, how = "first", order.by = NULL,
       na.rm = FALSE, with.ties = FALSE)
}
\arguments{
  \item{x}{a matrix, data frame or list-like object, including 'grouped_df'.}
  \item{\dots}{names or sequences of columns to group by - passed to \code{\link{fselect}}. If \code{x} is a matrix: atomic vectors to group \code{x}. Can be empty to operate on (un)grouped data.}
  \item{n}{integer or proportion (if < 1). Number of rows to select from each group. If a proportion is provided, it is converted to the equivalent number of rows.}
  \item{how}{character. Method to select rows. One of:
    \itemize{
      \item "first": select first n rows
      \item "last": select last n rows
      \item "min": select n rows with minimum values of order.by
      \item "max": select n rows with maximum values of order.by
    }
  }
  \item{order.by}{vector or column name to order by when \code{how} is "min" or "max". Must be same length as rows in x.}
  \item{na.rm}{logical. If \code{TRUE}, missing values in \code{order.by} are removed before selecting rows.}
  \item{sort}{logical. If \code{TRUE}, sort selected rows on the grouping columns.}
  \item{with.ties}{logical. If \code{TRUE} and \code{how} is "min" or "max", returns all rows that tie with the extreme value. Only applies when n = 1.}
}
\value{
A subset of the input matrix/data frame containing the selected rows. The output preserves the input's attributes and class.
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
