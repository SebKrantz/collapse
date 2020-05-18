\name{ftransform}
\alias{ftransform}
\alias{tfm}
\alias{settransform}
\alias{settfm}
\alias{fcompute}
\title{Fast Transform and Compute Columns on a Data Frame}
\description{
  \code{ftransform} is a much faster update of \code{base::transform} for data frames. It returns the data frame with new columns computed and/or existing columns modified or deleted. \code{settransform} does all of that by reference i.e. it modifies the data frame in the global environment. \code{fcompute} can be used to compute new columns from the columns in a data frame and returns only the computed columns.
}

\usage{
# Modify and return 'data.frame'
ftransform(X, ...)
tfm(X, ...)            # Shortcut for ftransform

# Modify 'data.frame' by reference
settransform(X, ...)
settfm(X, ...)         # Shortcut for settransform

# Compute and return new 'data.frame' from existing one
fcompute(X, ...)
}
\arguments{
    \item{X}{a data.frame.}
    \item{\dots}{further arguments of the form \code{column = value}. The \code{value} can be a combination of other columns, a scalar value, or \code{NULL}, which deletes \code{column}. }
}
\details{
  The \code{\dots} arguments to \code{ftransform} are tagged
  vector expressions, which are evaluated in the data frame
  \code{X}.  The tags are matched against \code{names(X)}, and for
  those that match, the value replace the corresponding variable in
  \code{X}, and the others are appended to \code{X}. It is also possible to delete columns by assigning \code{NULL} to them, i.e. \code{ftransform(data, column = NULL)} removes \code{column} from the data.

  The function \code{settransform} does all of that by reference, but uses base-R's copy-on modify semantics, which is equivalent to replacing the data with \code{<-} (thus it is still memory efficient but the data will have a different memory address after each call of \code{settransform}).

  Finally, the function \code{fcompute} functions just like \code{ftransform}, but returns only the changed / computed columns without modifying or appending the data in \code{X}.
}
\value{
  The modified data.frame \code{X}, or, for \code{fcompute}, a new data.frame with the columns computed on \code{X}. All attributes of \code{X} are preserved.
}
\seealso{\code{\link{with}}, \code{\link{within}}, \link[=A3-data-frame-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## ftransform modifies and returns a data.frame
ftransform(airquality, Ozone = -Ozone)
ftransform(airquality, new = -Ozone, Temp = (Temp-32)/1.8)
ftransform(airquality, new = -Ozone, new2 = 1, Temp = NULL)  # Deleting Temp
ftransform(airquality, Ozone = NULL, Temp = NULL)            # Deleting columns

## settransform modifies a data.frame in the global environment
airquality_c <- airquality
settransform(airquality_c, Ratio = Ozone / Temp, Ozone = NULL, Temp = NULL)
head(airquality_c)
rm(airquality_c)

## fcompute only returns the modified / computed data
fcompute(airquality, Ozone = -Ozone)
fcompute(airquality, new = -Ozone, Temp = (Temp-32)/1.8)
fcompute(airquality, new = -Ozone, new2 = 1)

}
\keyword{manip}
