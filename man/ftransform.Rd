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
\alias{fcomputev}
\title{Fast Transform and Compute Columns on a Data Frame}
\description{
  \code{ftransform} is a much faster version of \code{\link{transform}} and \code{dplyr::mutate} for data frames. It returns the data frame with new columns computed and/or existing columns modified or deleted. \code{settransform} does all of that by reference i.e. it modifies the data frame in the global environment. \code{fcompute} can be used to compute new columns from the columns in a data frame and returns only the computed columns.
}

\usage{
# Modify and return data frame
ftransform(.data, \dots)
ftransformv(.data, vars, FUN, \dots, apply = TRUE)
tfm(.data, \dots)               # Shortcut for ftransform
tfmv(.data, vars, FUN, \dots, apply = TRUE)

# Modify data frame by reference
settransform(.data, \dots)
settransformv(.data, vars, FUN, \dots, apply = TRUE)
settfm(.data, \dots)            # Shortcut for settransform
settfmv(.data, vars, FUN, \dots, apply = TRUE)

# Replace/add modified columns in/to a data frame
ftransform(.data) <- value
tfm(.data) <- value           # Shortcut for ftransform<-

# Compute columns, returned as a new data frame
fcompute(.data, \dots, keep = NULL)
fcomputev(.data, vars, FUN, \dots, apply = TRUE, keep = NULL)
}
\arguments{
    \item{.data}{a data frame or named list of columns.}
    \item{\dots}{further arguments of the form \code{column = value}. The \code{value} can be a combination of other columns, a scalar value, or \code{NULL}, which deletes \code{column}. Alternatively it is also possible to place a single list here, which will be treated like a list of \code{column = value} arguments. For \code{ftransformv} and \code{fcomputev}, \code{\dots} can be used to pass further arguments to \code{FUN}. \emph{Note:} The ellipsis (\code{\dots}) is always evaluated within the data frame (\code{.data}) environment. See Examples.}
    \item{vars}{variables to be transformed by applying \code{FUN} to them: select using names, indices, a logical vector or a selector function (e.g. \code{is.numeric}).}
    \item{FUN}{a single function yielding a result of length \code{NROW(.data)} or 1. See also \code{apply}.}
    \item{apply}{logical. \code{TRUE} (default) will apply \code{FUN} to each column selected in \code{vars}; \code{FALSE} will apply \code{FUN} to the subsetted data frame i.e. \code{FUN(get_vars(.data, vars), ...)}. The latter is useful for \emph{collapse} functions with data frame or grouped / panel data frame methods, yielding performance gains and enabling grouped transformations. See Examples.}
    \item{value}{a named list of replacements, it will be treated like an evaluated list of \code{column = value} arguments.}
    \item{keep}{select columns to preserve using column names, indices or a function (e.g. \code{is.numeric}). By default computed columns are added after the preserved ones, unless they are assigned the same name in which case the preserved columns will be replaced in order. See Examples.}
}
\details{
  The \code{\dots} arguments to \code{ftransform} are tagged
  vector expressions, which are evaluated in the data frame
  \code{.data}.  The tags are matched against \code{names(.data)}, and for
  those that match, the values replace the corresponding variable in
  \code{.data}, whereas the others are appended to \code{.data}. It is also possible to delete columns by assigning \code{NULL} to them, i.e. \code{ftransform(data, colk = NULL)} removes \code{colk} from the data. \emph{Note} that \code{names(.data)} and the names of the \code{...} arguments are checked for uniqueness beforehand, yielding an error if this is not the case.

  Since \emph{collapse} v1.3.0, is is also possible to pass a single named list to \code{\dots}, i.e. \code{ftransform(data, newdata)}. This list will be treated like a list of tagged vector expressions. \emph{Note} the different behavior: \code{ftransform(data, list(newcol = col1))} is the same as \code{ftransform(data, newcol = col1)}, whereas \code{ftransform(data, newcol = as.list(col1))} creates a list column. Something like \code{ftransform(data, as.list(col1))} gives an error because the list is not named. See Examples.   % and \code{ftransform(data, as.list(col1))} gives an error because an unnamed list is passed. % , but \code{ftransform(data, setNames(as.list(col1), col1))} will work and add the values of \code{col1} as separate columns.

  % \code{ftransform(data, fmean(list(col1mean = col1, col2mean = col2), drop = FALSE))} etc.

  % For example \code{ftransformv(data, 1:3, log)} is the same as \code{ftransform(data, lapply(get_vars(data, 1:3), log))}, and \code{ftransformv(data, 1:3, log, apply = FALSE)} is the same as \code{ftransform(data, log(get_vars(data, 1:3)))}.

 The function \code{ftransformv} added in v1.3.2 provides a fast replacement for the functions \code{dplyr::mutate_at} and \code{dplyr::mutate_if} facilitating mutations of groups of columns (\code{dplyr::mutate_all} is already accounted for by \code{\link{dapply}}). See Examples.

  The function \code{settransform} does all of that by reference, but uses base-R's copy-on modify semantics, which is equivalent to replacing the data with \code{<-} (thus it is still memory efficient but the data will have a different memory address afterwards).

The function \code{fcompute(v)} works just like \code{ftransform(v)}, but returns only the changed / computed columns without modifying or appending the data in \code{.data}. See Examples.
}
\note{
Please note that unlike \code{dplyr::mutate}, \code{ftransform} does not do anything per se with a grouped data frame. This is not a limitation as it affords greater flexibility and performance in programming with the \code{.FAST_FUN}, through which \emph{collapse} supports \emph{various kinds} of fully vectorized grouped transformations (see \code{\link{TRA}} for a list of available transformations). See Examples or the answer provided \href{https://stackoverflow.com/questions/67349744/using-ftransform-along-with-fgroup-by-from-collapse-r-package}{here}.
}
\value{
  The modified data frame \code{.data}, or, for \code{fcompute}, a new data frame with the columns computed on \code{.data}. All attributes of \code{.data} are preserved.
}
\seealso{\code{\link{with}}, \code{\link{within}}, \code{\link{add_vars}}, \link[=A3-fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}


\examples{
## ftransform modifies and returns a data.frame
head(ftransform(airquality, Ozone = -Ozone))
head(ftransform(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
head(ftransform(airquality, new = -Ozone, new2 = 1, Temp = NULL))  # Deleting Temp
head(ftransform(airquality, Ozone = NULL, Temp = NULL))            # Deleting columns

# With collapse's grouped and weighted functions, complex operations are done on the fly
head(ftransform(airquality, # Grouped operations by month:
                Ozone_Month_median = fmedian(Ozone, Month, TRA = "replace_fill"),
                Ozone_Month_sd = fsd(Ozone, Month, TRA = "replace"),
                Ozone_Month_centered = fwithin(Ozone, Month)))

# Grouping by month and above/below average temperature in each month
head(ftransform(airquality, Ozone_Month_high_median =
                  fmedian(Ozone, list(Month, Temp > fbetween(Temp, Month)), TRA = "replace_fill")))

## ftransformv can be used to modify multiple columns using a function
head(ftransformv(airquality, 1:3, log))
head(`[<-`(airquality, 1:3, value = lapply(airquality[1:3], log))) # Same thing in base R

head(ftransformv(airquality, 1:3, log, apply = FALSE))
head(`[<-`(airquality, 1:3, value = log(airquality[1:3])))         # Same thing in base R

# Using apply = FALSE yields meaningful performance gains with collapse functions
# This calls fwithin.default, and repeates the grouping by month 3 times:
head(ftransformv(airquality, 1:3, fwithin, Month))

# This calls fwithin.data.frame, and only groups one time -> 5x faster!
head(ftransformv(airquality, 1:3, fwithin, Month, apply = FALSE))

library(magrittr) # Pipe operators
# This also works for grouped and panel data frames (calling fwithin.grouped_df)
airquality \%>\% fgroup_by(Month) \%>\%
  ftransformv(1:3, fwithin, apply = FALSE) \%>\% head

# But this gives the WRONG result (calling fwithin.default). Need option apply = FALSE!!
airquality \%>\% fgroup_by(Month) \%>\%
  ftransformv(1:3, fwithin) \%>\% head

# For grouped modification of single columns in a grouped dataset, we can use GRP():
airquality \%>\% fgroup_by(Month) \%>\%
  ftransform(W_Ozone = fwithin(Ozone, GRP(.)),                 # Grouped centering
             sd_Ozone_m = fsd(Ozone, GRP(.), TRA = "replace"), # In-Month standard deviation
             sd_Ozone = fsd(Ozone, TRA = "replace"),           # Overall standard deviation
             sd_Ozone2 = fsd(Ozone, TRA = "replace_fill"),     # Same, overwriting NA's
             sd_Ozone3 = fsd(Ozone)) \%>\% head                  # Same thing (calling alloc())

rm(airquality)

## For more complex mutations we can use ftransform with compound pipes
airquality \%>\% fgroup_by(Month) \%>\%
  ftransform(get_vars(., 1:3) \%>\% fwithin \%>\% flag(0:2)) \%>\% head

airquality \%>\% ftransform(STD(., cols = 1:3) \%>\% replace_NA(0)) \%>\% head

# The list argument feature also allows flexible operations creating multiple new columns
airquality \%>\% # The variance of Wind and Ozone, by month, weighted by temperature:
  ftransform(fvar(list(Wind_var = Wind, Ozone_var = Ozone), Month, Temp, "replace")) \%>\% head

# Same as above using a grouped data frame (a bit more complex)
airquality \%>\% fgroup_by(Month) \%>\%
  ftransform(fselect(., Wind, Ozone) \%>\% fvar(Temp, "replace") \%>\% add_stub("_var", FALSE)) \%>\%
  fungroup \%>\% head

# This performs 2 different multi-column grouped operations (need c() to make it one list)
ftransform(airquality, c(fmedian(list(Wind_Day_median = Wind,
                                      Ozone_Day_median = Ozone), Day, TRA = "replace"),
                         fsd(list(Wind_Month_sd = Wind,
                                  Ozone_Month_sd = Ozone), Month, TRA = "replace"))) \%>\% head

## settransform(v) works like ftransform(v) but modifies a data frame in the global environment..
settransform(airquality, Ratio = Ozone / Temp, Ozone = NULL, Temp = NULL)
head(airquality)
rm(airquality)

# Grouped and weighted centering
settransformv(airquality, 1:3, fwithin, Month, Temp, apply = FALSE)
head(airquality)
rm(airquality)

# Suitably lagged first-differences
settransform(airquality, get_vars(airquality, 1:3) \%>\% fdiff \%>\% flag(0:2))
head(airquality)
rm(airquality)

# Same as above using magrittr::`\%<>\%`
airquality \%<>\% ftransform(get_vars(., 1:3) \%>\% fdiff \%>\% flag(0:2))
head(airquality)
rm(airquality)

# It is also possible to achieve the same thing via a replacement method (if needed)
ftransform(airquality) <- get_vars(airquality, 1:3) \%>\% fdiff \%>\% flag(0:2)
head(airquality)
rm(airquality)

## fcompute only returns the modified / computed columns
head(fcompute(airquality, Ozone = -Ozone))
head(fcompute(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
head(fcompute(airquality, new = -Ozone, new2 = 1))

# Can preserve existing columns, computed ones are added to the right if names are different
head(fcompute(airquality, new = -Ozone, new2 = 1, keep = 1:3))

# If given same name as preserved columns, preserved columns are replaced in order...
head(fcompute(airquality, Ozone = -Ozone, new = 1, keep = 1:3))

# Same holds for fcomputev
head(fcomputev(iris, is.numeric, log)) # Same as:
iris \%>\% get_vars(is.numeric) \%>\% dapply(log) \%>\% head()

head(fcomputev(iris, is.numeric, log, keep = "Species"))   # Adds in front
head(fcomputev(iris, is.numeric, log, keep = names(iris))) # Preserve order

# Keep a subset of the data, add standardized columns
head(fcomputev(iris, 3:4, STD, apply = FALSE, keep = names(iris)[3:5]))
}
\keyword{manip}
