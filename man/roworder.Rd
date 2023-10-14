\name{roworder}
\alias{roworder}
\alias{roworderv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Reordering of Data Frame Rows
}
\description{
A fast substitute for \code{dplyr::arrange}. It returns a sorted copy of the data frame, unless the data is already sorted in which case no copy is made. In addition, rows can be manually re-ordered. Use \code{data.table::setorder} to sort a data frame without creating a copy. %\code{roworder} also does not support grouped tibbles or pdata.frame's, i.e. every data frame is treated the same.
}
\usage{
roworder(X, \dots, na.last = TRUE, verbose = .op[["verbose"]])

roworderv(X, cols = NULL, neworder = NULL, decreasing = FALSE,
          na.last = TRUE, pos = "front", verbose = .op[["verbose"]])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data frame or list of equal-length columns. }
  \item{\dots}{comma-separated columns of \code{X} to sort by e.g. \code{var1, var2}. Negatives i.e. \code{-var1, var2} can be used to sort in decreasing order of \code{var1}. Internally all expressions are turned into strings and \code{startsWith(expr, "-")} is used to detect this, thus it does not negate the actual values (which may as well be strings), and you cannot apply any other functions to columns inside \code{roworder()} to induce different sorting behavior.}
  \item{cols}{select columns to sort by using a function, column names, indices or a logical vector. The default \code{NULL} sorts by all columns in order of occurrence (from left to right). }
  \item{na.last}{logical. If \code{TRUE}, missing values in the sorting columns are placed last; if \code{FALSE}, they are placed first; if \code{NA} they are removed (argument passed to \code{\link{radixorderv}}).}
  \item{decreasing}{logical. Should the sort order be increasing or decreasing? Can also be a vector of length equal to the number of arguments in \code{cols} (argument passed to \code{\link{radixorderv}}).}
  \item{neworder}{an ordering vector, can be \code{< nrow(X)}. if \code{pos = "front"} or \code{pos = "end"}, a logical vector can also be supplied. This argument overwrites \code{cols}.}
    \item{pos}{integer or character. Different arrangement options if \code{!is.null(neworder) && length(neworder) < nrow(X)}.
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "front"   \tab\tab move rows in \code{neworder} to the front (top) of \code{X} (the default). \cr
                 2 \tab\tab "end" \tab\tab move rows in \code{neworder} to the end (bottom) of \code{X}. \cr
                 3 \tab\tab "exchange"   \tab\tab just exchange the order of rows in \code{neworder}, other rows remain in the same position. \cr
                 4 \tab\tab "after"  \tab\tab place all further selected rows behind the first selected row. \cr

  }
  }
  \item{verbose}{logical. \code{TRUE} (default) prints a message when ordering a grouped or indexed frame, indicating that this is not efficient and encouraging reordering the data prior to the grouping/indexing step.}

}

\value{
A copy of \code{X} with rows reordered. If \code{X} is already sorted, \code{X} is simply returned.
}
\note{
If you don't require a copy of the data, use \code{data.table::setorder} (you can also use it in a piped call as it invisibly returns the data).

\code{roworder(v)} has internal facilities to deal with \link[=indexing]{indexed data}.
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{colorder}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
head(roworder(airquality, Month, -Ozone))
head(roworder(airquality, Month, -Ozone, na.last = NA))  # Removes the missing values in Ozone

## Same in standard evaluation
head(roworderv(airquality, c("Month", "Ozone"), decreasing = c(FALSE, TRUE)))
head(roworderv(airquality, c("Month", "Ozone"), decreasing = c(FALSE, TRUE), na.last = NA))

## Custom reordering
head(roworderv(mtcars, neworder = 3:4))               # Bring rows 3 and 4 to the front
head(roworderv(mtcars, neworder = 3:4, pos = "end"))  # Bring them to the end
head(roworderv(mtcars, neworder = mtcars$vs == 1))    # Bring rows with vs == 1 to the top
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
