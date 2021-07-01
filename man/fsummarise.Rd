\name{fsummarise}
\alias{fsummarise}
\alias{smr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Summarise
}
\description{
\code{fsummarize} is a much faster version of \code{dplyr::summarise}, when used together with the \link[=A1-fast-statistical-functions]{Fast Statistical Functions}.
}
\usage{
fsummarise(.data, ..., keep.group_vars = TRUE)
smr(.data, ..., keep.group_vars = TRUE)        # Shortcut
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.data}{
 a (grouped) data frame or named list of columns. Grouped data can be created with \code{\link{fgroup_by}} or \code{dplyr::group_by}.
}
  \item{\dots}{
  name-value pairs of summary functions. The name will be the name of the variable in the result. Functions when applied to a vector need to return a scalar. For fast performance and weighted aggregation use the \link[=A1-fast-statistical-functions]{Fast Statistical Functions}.
}
  \item{keep.group_vars}{
logical. \code{FALSE} removes grouping variables after computation.
}
}
\value{
If \code{.data} is grouped by \code{\link{fgroup_by}} or \code{dplyr::group_by}, the result is a data frame of the same class and attributes with rows reduced to the number of groups. If \code{.data} is not grouped, the result is a data frame of the same class and attributes with 1 row.
}
\note{
Functions other than the \link[=A1-fast-statistical-functions]{Fast Statistical Functions} cannot be passed weight vectors as \code{fsummarize} only splits the \code{.data}, not any arguments passed to functions. If you want to summarize data with a weighted function you created yourself, use \emph{data.table}.
}
\seealso{
\code{\link{collap}}, \link[=A3-fast-data-manipulation]{Data Frame Manipulation}, \link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}

}
\examples{
# Simple use
fsummarise(mtcars, mean_mpg = fmean(mpg),
                   sd_mpg = fsd(mpg))

# Using base functions (not a big difference without groups)
fsummarise(mtcars, mean_mpg = mean(mpg),
                   sd_mpg = sd(mpg))
\donttest{ % No code relying on suggested package or base Pipe
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
