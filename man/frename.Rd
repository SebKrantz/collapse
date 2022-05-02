\name{frename}
\alias{rnm}
\alias{frename}
\alias{setrename}
\alias{relabel}
\alias{setrelabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Renaming and Relabelling Objects
}
\description{
A fast substitute for \code{dplyr::rename}. \code{setrename} renames objects by reference. These functions also work with objects other than data frames that have a 'names' attribute. \code{relabel} and \code{setrelabel} do that same for labels attached to list elements / data frame columns.
}
\usage{
frename(.x, \dots, cols = NULL)
rnm(.x, \dots, cols = NULL)     # Shorthand for frename()

setrename(.x, \dots, cols = NULL)

relabel(.x, \dots, cols = NULL, attrn = "label")

setrelabel(.x, \dots, cols = NULL, attrn = "label")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.x}{for \code{(f/set)rename}: an R object with a \code{"names"} attribute. For \code{(set)relabel}: a named list.
}
  \item{\dots}{either tagged vector expressions of the form \code{name = newname} / \code{name = newlabel}, or a single function (+ optional arguments to the function) applied to all names/labels (of columns/elements selected in \code{cols}).
}
  \item{cols}{If \code{\dots} is a function, select a subset of columns/elements to rename/relabel using names, indices, a logical vector or a function applied to the columns if \code{.x} is a list (e.g. \code{is.numeric}).}

  \item{attrn}{character. Name of attribute to store labels or retrieve labels from.}
}
\value{
\code{.x} renamed / relabelled. \code{setrename} and \code{setrelabel} return \code{.x} invisibly.
}
\note{
Note that both \code{relabel} and \code{setrelabel} modify \code{.x} by reference. This is because labels are attached to columns themselves, making it impossible to avoid permanent modification by taking a shallow copy of the encompassing list / data.frame. On the other hand \code{frename} makes a shallow copy whereas \code{setrename} also modifies by reference.
}
\seealso{
\link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
## Using tagged expressions
head(frename(iris, Sepal.Length = SL, Sepal.Width = SW,
                   Petal.Length = PL, Petal.Width = PW))
head(frename(iris, Sepal.Length = "S L", Sepal.Width = "S W",
                   Petal.Length = "P L", Petal.Width = "P W"))

## Using a function
head(frename(iris, tolower))
head(frename(iris, tolower, cols = 1:2))
head(frename(iris, tolower, cols = is.numeric))
head(frename(iris, paste, "new", sep = "_", cols = 1:2))

## Renaming by reference
# setrename(iris, tolower)
# head(iris)
# rm(iris)

## Relabelling (by reference)
# namlab(relabel(wlddev, PCGDP = "GDP per Capita", LIFEEX = "Life Expectancy"))
# namlab(relabel(wlddev, toupper))


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
