\name{rsplit}
\alias{rsplit}
\alias{rsplit.default}
\alias{rsplit.data.frame}
\title{
Recursive Splitting
}
\description{
\code{rsplit} (recursively) splits a vector or data frame into subsets according to combinations of (multiple) vectors / factors and returns a (nested) list. If \code{flatten = TRUE}, the list is flattened yielding the same result as \code{\link{split}}. \code{rsplit} is implemented as a wrapper around \code{\link{gsplit}}, and significantly faster than \code{\link{split}}.
}
\usage{
rsplit(x, \dots)

\method{rsplit}{default}(x, fl, drop = TRUE, flatten = FALSE, use.names = TRUE, \dots)

\method{rsplit}{data.frame}(x, by, drop = TRUE, flatten = FALSE, cols = NULL,
       keep.by = FALSE, simplify = TRUE, use.names = TRUE, \dots)
}
\arguments{
  \item{x}{a vector, data.frame or list.}
  \item{fl}{a \code{\link{GRP}} object, or a (list of) vector(s) / factor(s) (internally converted to a \code{\link{GRP}} object(s)) used to split \code{x}.}
  \item{by}{\emph{data.frame method}: Same as \code{fl}, but also allows one- or two-sided formulas i.e. \code{~ group1} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{drop}{logical. \code{TRUE} removes unused levels or combinations of levels from factors before splitting; \code{FALSE} retains those combinations yielding empty list elements in the output.}
  \item{flatten}{logical. If \code{fl} is a list of vectors / factors, \code{TRUE} calls \code{\link{GRP}} on the list, creating a single grouping used for splitting; \code{FALSE} yields recursive splitting.}
  \item{use.names}{logical. \code{TRUE} returns a named list (like \code{\link{split}}); \code{FALSE} returns a plain list.}
  \item{cols}{\emph{data.frame method}: Select columns to split using a function, column names, indices or a logical vector. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
  \item{keep.by}{logical. If a formula is passed to \code{by}, then \code{TRUE} preserves the splitting (right-hand-side) variables in the data frame.}
  \item{simplify}{\emph{data.frame method}: Logical. \code{TRUE} calls \code{rsplit.default} if a single column is split e.g. \code{rsplit(data, col1 ~ group1)} becomes the same as \code{rsplit(data$col1, data$group1)}.}
\item{\dots}{further arguments passed to \code{\link{GRP}}. Sensible choices would be \code{sort = FALSE}, \code{decreasing = TRUE} or \code{na.last = FALSE}. Note that these options only apply if \code{fl} is not already a (list of) factor(s).}
}
\value{
a (nested) list containing the subsets of \code{x}.
}
\seealso{
\code{\link{gsplit}}, \code{\link{rapply2d}}, \code{\link{unlist2d}}, \link[=list-processing]{List Processing}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
rsplit(mtcars$mpg, mtcars$cyl)
rsplit(mtcars, mtcars$cyl)

rsplit(mtcars, mtcars[.c(cyl, vs, am)])
rsplit(mtcars, ~ cyl + vs + am, keep.by = TRUE)  # Same thing
rsplit(mtcars, ~ cyl + vs + am)

rsplit(mtcars, ~ cyl + vs + am, flatten = TRUE)

rsplit(mtcars, mpg ~ cyl)
rsplit(mtcars, mpg ~ cyl, simplify = FALSE)
rsplit(mtcars, mpg + hp ~ cyl + vs + am)
rsplit(mtcars, mpg + hp ~ cyl + vs + am, keep.by = TRUE)

# Split this sectoral data, first by Variable (Emloyment and Value Added), then by Country
GGDCspl <- rsplit(GGDC10S, ~ Variable + Country, cols = 6:16)
str(GGDCspl)

# The nested list can be reassembled using unlist2d()
head(unlist2d(GGDCspl, idcols = .c(Variable, Country)))
rm(GGDCspl)

# Another example with mtcars (not as clean because of row.names)
nl <- rsplit(mtcars, mpg + hp ~ cyl + vs + am)
str(nl)
unlist2d(nl, idcols = .c(cyl, vs, am), row.names = "car")
rm(nl)
}
\keyword{manip}
