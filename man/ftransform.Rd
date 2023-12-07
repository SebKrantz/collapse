\name{ftransform}
\alias{fmutate}
\alias{mtt}
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
  \code{ftransform} is a much faster version of \code{\link{transform}} for data frames. It returns the data frame with new columns computed and/or existing columns modified or deleted. \code{settransform} does all of that by reference. \code{fcompute} computes and returns new columns. These functions evaluate all arguments simultaneously, allow list-input (nested pipelines) and disregard grouped data.

Catering to the \emph{tidyverse} user, v1.7.0 introduced \code{fmutate}, providing familiar functionality i.e. arguments are evaluated sequentially, computation on grouped data is done by groups, and functions can be applied to multiple columns using \code{\link{across}}. See also the Details.
}

\usage{
# dplyr-style mutate (sequential evaluation + across() feature)
fmutate(.data, ..., .keep = "all", .cols = NULL)
mtt(.data, ..., .keep = "all", .cols = NULL) # Shorthand for fmutate

# Modify and return data frame
ftransform(.data, \dots)
ftransformv(.data, vars, FUN, \dots, apply = TRUE)
tfm(.data, \dots)               # Shorthand for ftransform
tfmv(.data, vars, FUN, \dots, apply = TRUE)

# Modify data frame by reference
settransform(.data, \dots)
settransformv(.data, \dots)     # Same arguments as ftransformv
settfm(.data, \dots)            # Shorthand for settransform
settfmv(.data, \dots)

# Replace/add modified columns in/to a data frame
ftransform(.data) <- value
tfm(.data) <- value           # Shorthand for ftransform<-

# Compute columns, returned as a new data frame
fcompute(.data, \dots, keep = NULL)
fcomputev(.data, vars, FUN, \dots, apply = TRUE, keep = NULL)
}
\arguments{
    \item{.data}{a data frame or named list of columns.}
    \item{\dots}{further arguments of the form \code{column = value}. The \code{value} can be a combination of other columns, a scalar value, or \code{NULL}, which deletes \code{column}. Alternatively it is also possible to place a single list here, which will be treated like a list of \code{column = value} arguments. For \code{ftransformv} and \code{fcomputev}, \code{\dots} can be used to pass further arguments to \code{FUN}. The ellipsis (\code{\dots}) is always evaluated within the data frame (\code{.data}) environment. See Examples.
\code{fmutate} additionally supports \code{\link{across}} statements, and evaluates tagged vector expressions sequentially. With grouped execution, \code{dots} can also contain arbitrary expressions that result in a list of data-length columns. See Examples.}
    \item{vars}{variables to be transformed by applying \code{FUN} to them: select using names, indices, a logical vector or a selector function (e.g. \code{is.numeric}). Since v1.7 \code{vars} is evaluated within the \code{.data} environment, permitting expressions on columns e.g. \code{c(col1, col3:coln)}.}
    \item{FUN}{a single function yielding a result of length \code{NROW(.data)} or 1. See also \code{apply}.}
    \item{apply}{logical. \code{TRUE} (default) will apply \code{FUN} to each column selected in \code{vars}; \code{FALSE} will apply \code{FUN} to the subsetted data frame i.e. \code{FUN(get_vars(.data, vars), ...)}. The latter is useful for \emph{collapse} functions with data frame or grouped / panel data frame methods, yielding performance gains and enabling grouped transformations. See Examples.}
    \item{value}{a named list of replacements, it will be treated like an evaluated list of \code{column = value} arguments.}
    \item{keep}{select columns to preserve using column names, indices or a function (e.g. \code{is.numeric}). By default computed columns are added after the preserved ones, unless they are assigned the same name in which case the preserved columns will be replaced in order.}
   \item{.keep}{either one of \code{"all", "used", "unused"} or \code{"none"} (see \code{\link[dplyr]{mutate}}), or columns names/indices/function as \code{keep}. \emph{Note} that this does not work well with \code{across()} or other expressions supported since v1.9.0. The only sensible option you have there is to supply a character vector of all columns in the final dataset that you want to keep.   }

   \item{.cols}{for expressions involving \code{.data}, \code{.cols} can be used to subset columns, e.g. \code{mtcars |> gby(cyl) |> mtt(broom::augment(lm(mpg ~., .data)), .cols = 1:7)}. Can pass column names, indices, a logical vector or a selector function (e.g. \code{is.numericr}).}
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

 The function \code{ftransformv} added in v1.3.2 provides a fast replacement for the functions \code{dplyr::mutate_at} and \code{dplyr::mutate_if} (without the grouping feature) facilitating mutations of groups of columns (\code{dplyr::mutate_all} is already accounted for by \code{\link{dapply}}). See Examples.

  The function \code{settransform} does all of that by reference, but uses base-R's copy-on modify semantics, which is equivalent to replacing the data with \code{<-} (thus it is still memory efficient but the data will have a different memory address afterwards).

The function \code{fcompute(v)} works just like \code{ftransform(v)}, but returns only the changed / computed columns without modifying or appending the data in \code{.data}. See Examples.

 The function \code{fmutate} added in v1.7.0, provides functionality familiar from \emph{dplyr} 1.0.0 and higher. It evaluates tagged vector expressions sequentially and does operations by groups on a grouped frame (thus it is slower than \code{ftransform} if you have many tagged expressions or a grouped data frame). Note however that \emph{collapse} does not depend on \emph{rlang}, so things like lambda expressions are not available. \emph{Note also} that \code{fmutate} operates differently on grouped data whether you use \code{.FAST_FUN} or base R functions / functions from other packages. With \code{.FAST_FUN} (including \code{.OPERATOR_FUN}, excluding \code{fhdbetween} / \code{fhdwithin} / \code{HDW} / \code{HDB}), \code{fmutate} performs an efficient vectorized execution, i.e. the grouping object from the grouped data frame is passed to the \code{g} argument of these functions, and for \code{.FAST_STAT_FUN} also \code{TRA = "replace_fill"} is set (if not overwritten by the user), yielding internal grouped computation by these functions without the need for splitting the data by groups. For base R and other functions, \code{fmutate} performs classical split-apply combine computing i.e. the relevant columns of the data are selected and split into groups, the expression is evaluated for each group, and the result is recombined and suitably expanded to match the original data frame. \bold{Note} that it is not possible to mix vectorized and standard execution in the same expression!! Vectorized execution is performed if \bold{any} \code{.FAST_FUN} or \code{.OPERATOR_FUN} is part of the expression, thus a code like \code{mtcars |> gby(cyl) |> fmutate(new = fmin(mpg) / min(mpg))} will be expanded to something like \code{mtcars \%>\% gby(cyl) \%>\% ftransform(new = fmin(mpg, g = GRP(.), TRA = "replace_fill") / min(mpg))} and then executed, i.e. \code{fmin(mpg)} will be executed in a vectorized way, and \code{min(mpg)} will not be executed by groups at all.
}
\note{
\code{ftransform} ignores grouped data. This is on purpose as it allows non-grouped transformation inside a pipeline on grouped data, and affords greater flexibility and performance in programming with the \code{.FAST_FUN}. In particular, you can run a nested pipeline inside \code{ftransform}, and decide which expressions should be grouped, and you can use the ad-hoc grouping functionality of the \code{.FAST_FUN}, allowing operations where different groupings are applied simultaneously in an expression. See Examples or the answer provided \href{https://stackoverflow.com/questions/67349744/using-ftransform-along-with-fgroup-by-from-collapse-r-package}{here}.

\code{fmutate} on the other hand supports grouped operations just like \code{dplyr::mutate}, but works in two different ways depending on whether you use \code{.FAST_FUN} in an expression or other functions. See the Examples.
}
\value{
  The modified data frame \code{.data}, or, for \code{fcompute}, a new data frame with the columns computed on \code{.data}. All attributes of \code{.data} are preserved.
}

\seealso{
\code{\link{across}}, \code{\link{fsummarise}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}

\examples{
library(magrittr)

## fmutate() examples ---------------------------------------------------------------

# Please note that expressions are vectorized whenever they contain 'ANY' fast function
mtcars \%>\%
  fgroup_by(cyl, vs, am) \%>\%
  fmutate(mean_mpg = fmean(mpg),                     # Vectorized
          mean_mpg_base = mean(mpg),                 # Non-vectorized
          mpg_cumpr = fcumsum(mpg) / fsum(mpg),      # Vectorized
          mpg_cumpr_base = cumsum(mpg) / sum(mpg),   # Non-vectorized
          mpg_cumpr_mixed = fcumsum(mpg) / sum(mpg)) # Vectorized: division by overall sum

# Using across: here fmean() gets vectorized across both groups and columns (requiring a single
# call to fmean.data.frame which goes to C), whereas weighted.mean needs to be called many times.
mtcars \%>\% fgroup_by(cyl, vs, am) \%>\%
  fmutate(across(disp:qsec, list(mu = fmean, mu2 = weighted.mean), w = wt, .names = "flip"))

# Can do more complex things...
mtcars \%>\% fgroup_by(cyl) \%>\%
  fmutate(res = resid(lm(mpg ~ carb + hp, weights = wt)))

# Since v1.9.0: supports arbitrary expressions returning suitable lists
\dontrun{
mtcars \%>\% fgroup_by(cyl) \%>\%
  fmutate(broom::augment(lm(mpg ~ carb + hp, weights = wt)))

# Same thing using across() (supported before 1.9.0)
modelfun <- function(data) broom::augment(lm(mpg ~ carb + hp, data, weights = wt))
mtcars \%>\% fgroup_by(cyl) \%>\%
  fmutate(across(c(mpg, carb, hp, wt), modelfun, .apply = FALSE))
}


## ftransform() / fcompute() examples: ----------------------------------------------

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
\donttest{ % No code relying on suggested package
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

airquality \%>\% ftransform(STD(., cols = 1:3) \%>\% replace_na(0)) \%>\% head

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
}
## fcompute only returns the modified / computed columns
head(fcompute(airquality, Ozone = -Ozone))
head(fcompute(airquality, new = -Ozone, Temp = (Temp-32)/1.8))
head(fcompute(airquality, new = -Ozone, new2 = 1))

# Can preserve existing columns, computed ones are added to the right if names are different
head(fcompute(airquality, new = -Ozone, new2 = 1, keep = 1:3))

# If given same name as preserved columns, preserved columns are replaced in order...
head(fcompute(airquality, Ozone = -Ozone, new = 1, keep = 1:3))
\donttest{ % No code relying on suggested package
# Same holds for fcomputev
head(fcomputev(iris, is.numeric, log)) # Same as:
iris \%>\% get_vars(is.numeric) \%>\% dapply(log) \%>\% head()
}
head(fcomputev(iris, is.numeric, log, keep = "Species"))   # Adds in front
head(fcomputev(iris, is.numeric, log, keep = names(iris))) # Preserve order

# Keep a subset of the data, add standardized columns
head(fcomputev(iris, 3:4, STD, apply = FALSE, keep = names(iris)[3:5]))
}
\keyword{manip}
