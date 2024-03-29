\name{fsummarise}
\alias{fsummarise}
\alias{fsummarize}
\alias{smr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Summarise
}
\description{
\code{fsummarise} is a much faster version of \code{dplyr::summarise}, when used together with the \link[=fast-statistical-functions]{Fast Statistical Functions}.

\code{fsummarize} and \code{fsummarise} are synonyms.
}
\usage{
fsummarise(.data, ..., keep.group_vars = TRUE, .cols = NULL)
fsummarize(.data, ..., keep.group_vars = TRUE, .cols = NULL)
smr(.data, ..., keep.group_vars = TRUE, .cols = NULL)        # Shorthand
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{
 a (grouped) data frame or named list of columns. Grouped data can be created with \code{\link{fgroup_by}} or \code{dplyr::group_by}.
}
  \item{\dots}{
  name-value pairs of summary functions, \code{\link{across}} statements, or arbitrary expressions resulting in a list. See Examples. For fast performance use the \link[=fast-statistical-functions]{Fast Statistical Functions}. % The name will be the name of the variable in the result. Functions when applied to a vector need to return a scalar.
}
  \item{keep.group_vars}{
logical. \code{FALSE} removes grouping variables after computation.
}
  \item{.cols}{
for expressions involving \code{.data}, \code{.cols} can be used to subset columns, e.g. \code{mtcars |> gby(cyl) |> smr(mctl(cor(.data), TRUE), .cols = 5:7)}. Can pass column names, indices, a logical vector or a selector function (e.g. \code{is.numericr}).
}
}
\value{
If \code{.data} is grouped by \code{\link{fgroup_by}} or \code{dplyr::group_by}, the result is a data frame of the same class and attributes with rows reduced to the number of groups. If \code{.data} is not grouped, the result is a data frame of the same class and attributes with 1 row.
}
\note{
Since v1.7, \code{fsummarise} is fully featured, allowing expressions using functions and columns of the data as well as external scalar values (just like \code{dplyr::summarise}). \bold{NOTE} however that once a \link[=fast-statistical-functions]{Fast Statistical Function} is used, the execution will be vectorized instead of split-apply-combine computing over groups. Please see the first Example.
}
\seealso{
\code{\link{across}}, \code{\link{collap}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}

}
\examples{
## Since v1.7, fsummarise supports arbitrary expressions, and expressions
## containing fast statistical functions receive vectorized execution:

# (a) This is an expression using base R functions which is executed by groups
mtcars |> fgroup_by(cyl) |> fsummarise(res = mean(mpg) + min(qsec))

# (b) Here, the use of fmean causes the whole expression to be executed
# in a vectorized way i.e. the expression is translated to something like
# fmean(mpg, g = cyl) + min(mpg) and executed, thus the result is different
# from (a), because the minimum is calculated over the entire sample
mtcars |> fgroup_by(cyl) |> fsummarise(mpg = fmean(mpg) + min(qsec))

# (c) For fully vectorized execution, use fmin. This yields the same as (a)
mtcars |> fgroup_by(cyl) |> fsummarise(mpg = fmean(mpg) + fmin(qsec))

# More advanced use: vectorized grouped regression slopes: mpg ~ carb
mtcars |>
  fgroup_by(cyl) |>
  fmutate(dm_carb = fwithin(carb)) |>
  fsummarise(beta = fsum(mpg, dm_carb) \%/=\% fsum(dm_carb^2))


# In across() statements it is fine to mix different functions, each will
# be executed on its own terms (i.e. vectorized for fmean and standard for sum)
mtcars |> fgroup_by(cyl) |> fsummarise(across(mpg:hp, list(fmean, sum)))

# Note that this still detects fmean as a fast function, the names of the list
# are irrelevant, but the function name must be typed or passed as a character vector,
# Otherwise functions will be executed by groups e.g. function(x) fmean(x) won't vectorize
mtcars |> fgroup_by(cyl) |> fsummarise(across(mpg:hp, list(mu = fmean, sum = sum)))

# We can force none-vectorized execution by setting .apply = TRUE
mtcars |> fgroup_by(cyl) |> fsummarise(across(mpg:hp, list(mu = fmean, sum = sum), .apply = TRUE))

# Another argument of across(): Order the result first by function, then by column
mtcars |> fgroup_by(cyl) |>
     fsummarise(across(mpg:hp, list(mu = fmean, sum = sum), .transpose = FALSE))


# Since v1.9.0, can also evaluate arbitrary expressions
mtcars |> fgroup_by(cyl, vs, am) |>
   fsummarise(mctl(cor(cbind(mpg, wt, carb)), names = TRUE))

# This can also be achieved using across():
corfun <- function(x) mctl(cor(x), names = TRUE)
mtcars |> fgroup_by(cyl, vs, am) |>
   fsummarise(across(c(mpg, wt, carb), corfun, .apply = FALSE))

#----------------------------------------------------------------------------
# Examples that also work for pre 1.7 versions

# Simple use
fsummarise(mtcars, mean_mpg = fmean(mpg),
                   sd_mpg = fsd(mpg))

# Using base functions (not a big difference without groups)
fsummarise(mtcars, mean_mpg = mean(mpg),
                   sd_mpg = sd(mpg))

# Grouped use
mtcars |> fgroup_by(cyl) |>
  fsummarise(mean_mpg = fmean(mpg),
             sd_mpg = fsd(mpg))

# This is still efficient but quite a bit slower on large data (many groups)
mtcars |> fgroup_by(cyl) |>
  fsummarise(mean_mpg = mean(mpg),
             sd_mpg = sd(mpg))

# Weighted aggregation
mtcars |> fgroup_by(cyl) |>
  fsummarise(w_mean_mpg = fmean(mpg, wt),
             w_sd_mpg = fsd(mpg, wt))

\donttest{ % The tidyverse regularly causes havoc to CRAN tests in other packages, therefore this is not tested
## Can also group with dplyr::group_by, but at a conversion cost, see ?GRP
library(dplyr)
mtcars |> group_by(cyl) |>
  fsummarise(mean_mpg = fmean(mpg),
             sd_mpg = fsd(mpg))

# Again less efficient...
mtcars |> group_by(cyl) |>
  fsummarise(mean_mpg = mean(mpg),
             sd_mpg = sd(mpg))

}
}
\keyword{manip}
