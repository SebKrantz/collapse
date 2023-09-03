\name{fmatch}
\alias{fmatch}

\title{Fast Matching}

\description{
  Fast matching of elements/rows in \code{x} to elements/rows in \code{table}.

  This is a much faster replacement for \code{\link[base]{match}} that works
  with atomic vectors and data frames / lists of equal-length vectors. It is the workhorse function of \code{\link{join}}.
}

\usage{
fmatch(x, table, nomatch = NA_integer_,
       count = FALSE, overid = 1L)

# Infix operators based on fmatch():
\%!in\%   # Opposite of \%in\%
\%iin\%   # = which(x \%in\% table), but more efficient
\%!iin\%  # = which(x \%!in\% table), but more efficient
# Use set_collapse(mask = "\%in\%") to replace \%in\% with
# a much faster version based on fmatch()
}

\arguments{
  \item{x}{a vector, list or data frame whose elements are matched against \code{table}. If a list/data frame, matches are found by comparing rows, unlike \code{\link{match}} which compares columns. }
  \item{table}{a vector, list or data frame to match against.}
  \item{nomatch}{integer. Value to be returned in the case when no match is found. Default is \code{NA_integer_}.}
  \item{count}{logical. Counts number of (unique) matches and attaches 3 attributes:
   \itemize{
   \item \code{"n_nomatch"}: The number of elements in \code{x} not matched \code{= sum(result == nomatch)}.
   \item \code{"table_size"}: The size of the table \code{ = NROW(table)}.
   \item \code{"n_distinct"}: The number of unique matches \code{ = fndistinct(result[result != nomatch])}.
   }
   \emph{Note} that computing these attributes requires an extra pass through the matching vector. Also note that these attributes contain no direct information about whether either \code{x} or \code{table} are unique. This can be checked with \code{\link{any_duplicated}}.
  }
  \item{overid}{integer. If \code{x/table} are lists/data frames, \code{fmatch} compares the rows incrementally, starting with the first two columns, and matching further columns as necessary (see Details). Overidentification corresponds to the case when a subset of the columns uniquely identify the data. In this case this argument controls the behavior:
   \itemize{
   \item \code{0}: Early termination: stop matching additional columns. Most efficient.
   \item \code{1}: Continue matching columns and issue a warning that the data is overidentified.
   \item \code{2}: Continue matching columns without warning.
   }

  }
}

\value{
  Integer vector containing the positions of first matches of \code{x} in \code{table}. \code{nomatch} is returned for elements of \code{x} that have no match in \code{table}.

  If \code{count = TRUE}, the result has additional attributes \code{"n_nomatch"} counting the number of nomatches, \code{"table_size"} containing the number of elements in \code{table}, and \code{"n_distinct"} counting the number of distinct matches.
}

\details{
  With data frames / lists, \code{fmatch} compares the rows but moves through the data on a column-by-column basis (like a vectorized hash join algorithm). With two or more columns, the first two colums are hased simultaneously for speed. Further columns can be added to this match. It is likely that the first 2, 3, 4 etc. columns of a data frame fully identify the data. After each column \code{fmatch()} internally checks whether the \code{table} rows that are still eligible for matching (eliminating \code{nomatch} rows from earlier columns) are unique. If this is the case and \code{overid = 0}, \code{fmatch()} terminates early without considering further columns. This is efficient but may give undesirable/wrong results if considering further columns would turn some additional elements of the result vector into \code{nomatch} values.
}

\seealso{
\code{\link{ckmatch}}, \code{\link{funique}}, \code{\link{group}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}

\examples{

x <- c("b", "c", "a")
table <- c("a", "b", "c")

fmatch(x, table)
#> [1] 3 1 2

fmatch(x, table, nomatch = 0)
#> [1] 2 1 3

fmatch(x, table, count = TRUE)
#> [1] 3 1 2
#> attr(,"n_nomatch")
#> [1] 0
#> attr(,"table_size")
#> [1] 3
#> attr(,"n_distinct")
#> [1] 3

}



\keyword{manip}
