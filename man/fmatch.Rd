\name{fmatch}
\alias{fmatch}
\alias{ckmatch}
\alias{\%!in\%}
\alias{\%!iin\%}
\alias{\%iin\%}

\title{Fast Matching}

\description{
  Fast matching of elements/rows in \code{x} to elements/rows in \code{table}.

  This is a much faster replacement for \code{\link[base]{match}} that works
  with atomic vectors and data frames / lists of equal-length vectors. It is the workhorse function of \code{\link{join}}.
}

\usage{
fmatch(x, table, nomatch = NA_integer_,
       count = FALSE, overid = 1L)

# Check match: throws an informative error for non-matched elements
# Default message reflects frequent internal use to check data frame columns
ckmatch(x, table, e = "Unknown columns:", \dots)

# Infix operators based on fmatch():
x \%!in\% table   # Opposite of \%in\%
x \%iin\% table   # = which(x \%in\% table), but more efficient
x \%!iin\% table  # = which(x \%!in\% table), but more efficient
# Use set_collapse(mask = "\%in\%") to replace \%in\% with
# a much faster version based on fmatch()
}

\arguments{
  \item{x}{a vector, list or data frame whose elements are matched against \code{table}. If a list/data frame, matches are found by comparing rows, unlike \code{\link{match}} which compares columns. }
  \item{table}{a vector, list or data frame to match against.}
  \item{nomatch}{integer. Value to be returned in the case when no match is found. Default is \code{NA_integer_}.}
  \item{count}{logical. Counts number of (unique) matches and attaches 4 attributes:
   \itemize{
   \item \code{"N.nomatch"}: The number of elements in \code{x} not matched \code{= sum(result == nomatch)}.
   \item \code{"N.groups"}: The size of the table \code{ = NROW(table)}.
   \item \code{"N.distinct"}: The number of unique matches \code{ = fndistinct(result[result != nomatch])}.
   \item \code{"class"}: The \code{\link[=qG]{"qG"}} class: needed for optimized computations on the results object (e.g. \code{funique(result)}, which is needed for a full join).
   }
   \emph{Note} that computing these attributes requires an extra pass through the matching vector. Also note that these attributes contain no general information about whether either \code{x} or \code{table} are unique, except for two special cases when N.groups = N.distinct (table is unique) or length(result) = N.distinct (x is unique). Otherwise use \code{\link{any_duplicated}} to check x/table.
  }
  \item{overid}{integer. If \code{x/table} are lists/data frames, \code{fmatch} compares the rows incrementally, starting with the first two columns, and matching further columns as necessary (see Details). Overidentification corresponds to the case when a subset of the columns uniquely identify the data. In this case this argument controls the behavior:
   \itemize{
   \item \code{0}: Early termination: stop matching additional columns. Most efficient.
   \item \code{1}: Continue matching columns and issue a warning that the data is overidentified.
   \item \code{2}: Continue matching columns without warning.
   }
  }
\item{e}{the error message thrown by \code{ckmatch} for non-matched elements. The message is followed by the comma-separated non-matched elements.}
\item{\dots}{further arguments to \code{fmatch}.}
}

\value{
  Integer vector containing the positions of first matches of \code{x} in \code{table}. \code{nomatch} is returned for elements of \code{x} that have no match in \code{table}. If \code{count = TRUE}, the result has additional attributes and a class \code{\link[=qG]{"qG"}}.
}

\details{
  With data frames / lists, \code{fmatch} compares the rows but moves through the data on a column-by-column basis (like a vectorized hash join algorithm). With two or more columns, the first two columns are hashed simultaneously for speed. Further columns can be added to this match. It is likely that the first 2, 3, 4 etc. columns of a data frame fully identify the data. After each column \code{fmatch()} internally checks whether the \code{table} rows that are still eligible for matching (eliminating \code{nomatch} rows from earlier columns) are unique. If this is the case and \code{overid = 0}, \code{fmatch()} terminates early without considering further columns. This is efficient but may give undesirable/wrong results if considering further columns would turn some additional elements of the result vector into \code{nomatch} values.
}

\seealso{
\code{\link{join}}, \code{\link{funique}}, \code{\link{group}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}

\examples{
x <- c("b", "c", "a", "e", "f", "ff")
fmatch(x, letters)
fmatch(x, letters, nomatch = 0)
fmatch(x, letters, count = TRUE)

# Table 1
df1 <- data.frame(
  id1 = c(1, 1, 2, 3),
  id2 = c("a", "b", "b", "c"),
  name = c("John", "Bob", "Jane", "Carl")
)
head(df1)
# Table 2
df2 <- data.frame(
  id1 = c(1, 2, 3, 3),
  id2 = c("a", "b", "c", "e"),
  name = c("John", "Janne", "Carl", "Lynne")
)
head(df2)

# This gives an overidentification warning: columns 1:2 identify the data
if(FALSE) fmatch(df1, df2)
# This just runs through without warning
fmatch(df1, df2, overid = 2)
# This terminates computation after first 2 columns
fmatch(df1, df2, overid = 0)
fmatch(df1[1:2], df2[1:2])  # Same thing!
# -> note that here we get an additional match based on the unique ids,
# which we didn't get before because "Jane" != "Janne"
}


\keyword{manip}
