\name{ftransform}
\alias{ftransform}
\alias{tfm}
\alias{ftransformv}
\alias{tfmv}
\alias{ftransform<-}
\alias{tfm<-}
\alias{settransform}
\alias{settfm}
\alias{settransformv}
\alias{settfmv}
\alias{fcompute}
\title{Fast Transform and Compute Columns on a Data Frame}
\description{
  \code{ftransform} is a much faster version of \code{base::transform} and \code{dplyr::mutate} for data frames. It returns the data frame with new columns computed and/or existing columns modified or deleted. \code{settransform} does all of that by reference i.e. it modifies the data frame in the global environment. \code{fcompute} can be used to compute new columns from the columns in a data frame and returns only the computed columns.
}

\usage{
# Modify and return data frame
ftransform(X, \dots)
ftransformv(X, vars, FUN, \dots, apply = TRUE)
tfm(X, \dots)               # Shortcut for ftransform
tfmv(X, vars, FUN, \dots, apply = TRUE)

# Modify data frame by reference
settransform(X, \dots)
settransformv(X, vars, FUN, \dots, apply = TRUE)
settfm(X, \dots)            # Shortcut for settransform
settfmv(X, vars, FUN, \dots, apply = TRUE)

# Replace modified columns in a data frame
ftransform(X) <- value
tfm(X) <- value           # Shortcut for ftransform<-

# Compute and return new data frame from existing one
fcompute(X, \dots)
}
\arguments{
    \item{X}{a data frame or named list of columns.}
    \item{\dots}{further arguments of the form \code{column = value}. The \code{value} can be a combination of other columns, a scalar value, or \code{NULL}, which deletes \code{column}. Alternatively it is also possible to place a single list here, which will be treated the same as a list of \code{column = value} arguments. For \code{ftransformv}, \code{\dots} can be used to pass further arguments to \code{FUN} (arguments are evaluated within the data frame environment). }
    \item{vars}{variables to be transformed by applying \code{FUN} to them: select using names, indices, a logical vector or a selector function (i.e. \code{is.numeric}).}
    \item{apply}{logical. \code{TRUE} (default) will apply \code{FUN} to each column selected in \code{vars}; \code{FALSE} will apply \code{FUN} to the subsetted data frame i.e. \code{FUN(get_vars(X, vars), ...)}. The latter is useful for \emph{collapse} functions with data frame or grouped tibble / panel-data frame methods, yielding performance gains and enabling grouped mutations.}
    \item{value}{a named list of replacements, it will be treated like an evaluated list of \code{column = value} arguments.}
}
\details{
  The \code{\dots} arguments to \code{ftransform} are tagged
  vector expressions, which are evaluated in the data frame
  \code{X}.  The tags are matched against \code{names(X)}, and for
  those that match, the values replace the corresponding variable in
  \code{X}, whereas the others are appended to \code{X}. It is also possible to delete columns by assigning \code{NULL} to them, i.e. \code{ftransform(data, colk = NULL)} removes \code{colk} from the data.

  Since \emph{collapse} v1.3.0, is is also possible to pass a single named list to \code{\dots}, i.e. \code{ftransform(data, newdata)} or \code{ftransform(data, fmean(list(col1mean = col1, col2mean = col2), drop = FALSE))} etc.. This list will be treated the same as a list of tagged vector expressions.

  Note that such arguments should not be named i.e. \code{ftransform(data, newcol = as.list(col1))} creates a list column (which also works fine), \code{ftransform(data, as.list(col1))} gives an error because an unnamed list is passed, but \code{ftransform(data, setNames(as.list(col1), col1))} will work and add the values of \code{col1} as separate columns. See Examples.

 The function \code{ftransformv} added in v1.3.2 provides a fast replacement for the functions \code{dplyr::mutate_at} and \code{dplyr::mutate_if} facilitating mutations of groups of columns (\code{dplyr::mutate_all} is accounted for by \code{\link{dapply}}). For example \code{ftransformv(data, 1:3, log)} is the same as \code{ftransform(data, lapply(get_vars(data, 1:3), log))}, and \code{ftransformv(data, 1:3, log, apply = FALSE)} is the same as \code{ftransform(data, log(get_vars(data, 1:3)))}. See Examples.

  The function \code{settransform} does all of that by reference, but uses base-R's copy-on modify semantics, which is equivalent to replacing the data with \code{<-} (thus it is still memory efficient but the data will have a different memory address after each call of \code{settransform}).

  Finally, the function \code{fcompute} functions just like \code{ftransform}, but returns only the changed / computed columns without modifying or appending the data in \code{X}.
}
\value{
  The modified data frame \code{X}, or, for \code{fcompute}, a new data frame with the columns computed on \code{X}. All attributes of \code{X} are preserved.
}
\seealso{\code{\link{with}}, \code{\link{within}}, \link[=A3-fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## ftransform modifies and returns a data.frame
head(ftransform(airquality, Ozone = -Ozone))
head(ftransform(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
head(ftransform(airquality, new = -Ozone, new2 = 1, Temp = NULL))  # Deleting Temp
head(ftransform(airquality, Ozone = NULL, Temp = NULL))            # Deleting columns

# This computes the median and standard-deviation of Ozone in each month
head(ftransform(airquality,
                   Ozone_Month_median = fmedian(Ozone, Month, TRA = "replace_fill"),
                   Ozone_Month_sd = fsd(Ozone, Month, TRA = "replace_fill")))

# Grouping by month and above/below average temperature in each month
head(ftransform(airquality, Ozone_Month_high_median =
     fmedian(Ozone, list(Month, Temp > fbetween(Temp, Month)), TRA = "replace_fill")))

## Since v1.3.0 one can pass a list of columns, and there is a replacement method
head(ftransform(airquality, STD(airquality, cols = 1:3)))  # Could use magrittr::`\%<>\%`
ftransform(airquality) <- fscale(get_vars(airquality, 1:3))
rm(airquality)

# This feature also allows to flexibly do grouped operations creating multiple new columns
head(ftransform(airquality,
   fmedian(list(Wind_Month_median = Wind,
                Ozone_Month_median = Ozone), Month, TRA = "replace_fill")))

# This performs 2 different multi-column grouped operations (need c() to make it one list)
head(ftransform(airquality, c(fmedian(list(Wind_Day_median = Wind,
                                      Ozone_Day_median = Ozone), Day, TRA = "replace_fill"),
                         fsd(list(Wind_Month_sd = Wind,
                                  Ozone_Month_sd = Ozone), Month, TRA = "replace_fill"))))

## settransform works like ftransform but modifies a data frame in the global environment..
airquality_c <- airquality
settransform(airquality_c, Ratio = Ozone / Temp, Ozone = NULL, Temp = NULL)
settransform(airquality_c, STD(airquality_c, cols = 1:2))    # This also works..
head(airquality_c)
rm(airquality_c)

## fcompute only returns the modified / computed columns
head(fcompute(airquality, Ozone = -Ozone))
head(fcompute(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
head(fcompute(airquality, new = -Ozone, new2 = 1))

}
\keyword{manip}
