\name{qF, qG, finteraction}
\alias{qF}
\alias{qG}
\alias{is_qG}
\alias{as_factor_qG}
\alias{finteraction}
\title{
Fast Factor Generation, Interactions and Vector Grouping
}
\description{
\code{qF}, shorthand for 'quick-factor' implements very fast (ordered) factor generation from atomic vectors using either radix ordering or index hashing.

\code{qG}, shorthand for 'quick-group', generates a kind of factor-light without the levels attribute but instead an attribute providing the number of levels. Optionally the levels / groups can be attached, but without converting them to character. Objects have a class 'qG'.

\code{finteraction} generates a factor by interacting multiple vectors or factors. In that process missing values are always replaced with a level and unused levels are always dropped.
}
\usage{
qF(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE, drop = FALSE,
   keep.attr = TRUE, method = c("auto", "radix", "hash"))

qG(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE,
   return.groups = FALSE, method = c("auto", "radix", "hash"))

is_qG(x)

as_factor_qG(x, ordered = FALSE, na.exclude = TRUE)

finteraction(\dots, ordered = FALSE, sort = TRUE)
}
\arguments{
  \item{x}{a atomic vector, factor or quick-group.}
  \item{ordered}{logical. Adds a class 'ordered'.}
  \item{na.exclude}{logical. \code{TRUE} preserves missing values (i.e. no level is generated for \code{NA}).}
  \item{sort}{logical. \code{TRUE} sorts the levels in ascending order (like \code{\link{factor}}); \code{FALSE} provides the levels in an arbitrary order (determined by the hash map if \code{method = "hash"}), which can be slightly faster.}
  \item{drop}{logical. If \code{x} is a factor, \code{TRUE} efficiently drops unused factor levels beforehand using \code{\link{fdroplevels}}.}
  \item{keep.attr}{logical. If \code{TRUE} and \code{x} has additional attributes apart from 'levels' and 'class', these are preserved in the conversion to factor.}
  \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto"   \tab\tab automatic selection: hash for character, logical or if \code{length(x) < 500}, else radix. \cr
                 2 \tab\tab "radix" \tab\tab use radix ordering to generate factors. See Details. \cr
                 3 \tab\tab "hash"   \tab\tab use index hashing to generate factors. See Details. \cr
  }
}
\item{return.groups}{logical. \code{TRUE} returns the unique elements / groups / levels of \code{x} in an attribute called 'groups'. Unlike \code{qF}, they are not converted to character.}
\item{\dots}{multiple atomic vectors or factors, or a single list of equal-length vectors or factors. See Details. }
}
\details{
These functions are quite important. Whenever a vector is passed to a \emph{collapse} function such as \code{fmean(mtcars, mtcars$cyl)}, is is grouped using \code{qF} or \code{qG}.

\code{qF} is a combination of \code{as.factor} and \code{factor}. Applying it to a vector i.e. \code{qF(x)} gives the same result as \code{as.factor(x)}. \code{qF(x, ordered = TRUE)} generates and ordered factor (same as \code{factor(x, ordered = TRUE)}), and \code{qF(x, na.exclude = FALSE)} generates a level for missing values (same as \code{factor(x, exclude = NULL)}). An important addition is that \code{qF(x, na.exclude = FALSE)} also adds a class 'na.included'. This prevents \emph{collapse} functions from checking missing values in the factor, and is thus computationally more efficient. Therefore factors used in grouped operations should preferably be generated using \code{qF(x, na.exclude = FALSE)}. Setting \code{sort = FALSE} gathers the levels in a random order (unless \code{method = "radix"} and \code{x} is numeric, in which case the levels are always sorted). This can provide a speed improvement, particularly for character data. % for non-numeric \code{x}.

There are two methods of computation: radix ordering and index hashing. Radix ordering is done through combining the functions \code{\link{radixorder}} and \code{\link{groupid}}. It is generally faster than index hashing for large numeric data (although there are exceptions). Index hashing is done using \code{Rcpp::sugar::sort_unique} and \code{Rcpp::sugar::match}. It is generally faster for character data. For logical data, a super fast one-pass method was written which is subsumed in the hash method. Regarding speed: In general \code{qF} is around 5x faster than \code{as.factor} on character data and about 30x faster on numeric data. Automatic method dispatch typically does a good job delivering optimal performance.

\code{qG} is in the first place a programmers function. It generates a factor-'light' consisting of only an integer grouping vector and an attribute providing the number of groups. It is faster and more memory efficient than \code{\link{GRP}} for grouping atomic vectors, which is the main reason it exists. The fact that it (optionally) returns the unique groups / levels without converting them to character is an added bonus (this also provides a small performance gain compared to \code{qF}).

\code{finteraction} is simply a wrapper around \code{as_factor_GRP(GRP.default(X, sort = TRUE))}, where X is replaced by the arguments in '\dots' combined in a list. See \code{\link{GRP}} for computational details. In general: All vectors, factors, or lists of vectors / factors passed can be interacted. Interactions always create a level for missing values and always drop any unused levels.
}
\value{
\code{qF} returns an (ordered) factor. \code{qG} returns an object of class 'qG': an integer grouping vector with an attribute 'N.groups' indicating the number of groups, and, if \code{return.groups = TRUE}, an attribute 'groups' containing the vector of unique groups / elements in \code{x} corresponding to the integer-id.
}
\note{
Neither \code{qF} nor \code{qG} can reorder groups / factor levels. These objects can however be converted into one another using \code{qF/qG} or the direct method \code{as_factor_qG}, and it is also possible to add a class 'ordered' (\code{ordered = TRUE}) and to create am extra level / integer for missing values (\code{na.exclude = FALSE}). % Apart from that \code{qF} and \code{qG} don't do much to each others objects.
}
\seealso{
\code{\link{groupid}}, \code{\link{GRP}}, \link[=A2-fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
cylF <- qF(mtcars$cyl)     # Factor from atomic vector
cylG <- qG(mtcars$cyl)     # Quick-group from atomic vector
cylG                       # See the simple structure of this object

cf  <- qF(wlddev$country)  # Bigger data
cf2 <- qF(wlddev$country, na.exclude = FALSE)  # With na.included class
dat <- num_vars(wlddev)
\donttest{ % No code relying on suggested package
# cf2 is faster in grouped operations because no missing value check is performed
library(microbenchmark)
microbenchmark(fmax(dat, cf), fmax(dat, cf2))
}
finteraction(mtcars$cyl, mtcars$vs)  # Interacting two variables (can be factors)
head(finteraction(mtcars))           # A more crude example..

}
\keyword{manip}
