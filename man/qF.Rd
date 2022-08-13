\name{qF-qG-finteraction}
\alias{qF}
\alias{qG}
\alias{is_qG}
\alias{as_factor_qG}
\alias{finteraction}
\alias{itn}
\title{
Fast Factor Generation, Interactions and Vector Grouping
}
\description{
\code{qF}, shorthand for 'quick-factor' implements very fast factor generation from atomic vectors using either radix ordering or index hashing followed by sorting.

\code{qG}, shorthand for 'quick-group', generates a kind of factor-light without the levels attribute but instead an attribute providing the number of levels. Optionally the levels / groups can be attached, but without converting them to character (which can have large performance implications). Objects have a class 'qG'.

\code{finteraction} generates a factor or 'qG' object by interacting multiple vectors or factors. In that process missing values are always replaced with a level and unused levels/combinations are always dropped.

\emph{collapse} internally makes optimal use of factors and 'qG' objects when passed as grouping vectors to statistical functions (\code{g/by}, or \code{t} arguments) i.e. typically no further grouping or ordering is performed and objects are used directly by statistcial C/C++ code.
}
\usage{
qF(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE, drop = FALSE,
   keep.attr = TRUE, method = "auto")

qG(x, ordered = FALSE, na.exclude = TRUE, sort = TRUE,
   return.groups = FALSE, method = "auto")

is_qG(x)

as_factor_qG(x, ordered = FALSE, na.exclude = TRUE)

finteraction(\dots, factor = TRUE, ordered = FALSE, sort = factor,
             method = "auto")
itn(\dots) # Shorthand for finteraction
}
\arguments{
  \item{x}{a atomic vector, factor or quick-group.}
  \item{ordered}{logical. Adds a class 'ordered'.}
  \item{na.exclude}{logical. \code{TRUE} preserves missing values (i.e. no level is generated for \code{NA}). \code{FALSE} attaches an additional class \code{"na.included"} which is used to skip missing value checks performed before sending objects to C/C++. See Details. }
  \item{sort}{logical. \code{TRUE} sorts the levels in ascending order (like \code{\link{factor}}); \code{FALSE} provides the levels in order of first appearance, which can be significantly faster. Note that if a factor is passed as input, only \code{sort = FALSE} takes effect and unused levels will be dropped (as factors usually have sorted levels and checking sortedness can be expensive).}
  \item{drop}{logical. If \code{x} is a factor, \code{TRUE} efficiently drops unused factor levels beforehand using \code{\link{fdroplevels}}.}
  \item{keep.attr}{logical. If \code{TRUE} and \code{x} has additional attributes apart from 'levels' and 'class', these are preserved in the conversion to factor.}
  \item{method}{an integer or character string specifying the method of computation:
        \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 1 \tab\tab "auto"   \tab\tab automatic selection: \code{if(is.double(x) && sort) "radix" else if(sort && length(x) < 1e5) "rcpp_hash" else "hash"}. \cr
                 2 \tab\tab "radix" \tab\tab use radix ordering to generate factors. Supports \code{sort = FALSE} only for character vectors. See Details. \cr
                 3 \tab\tab "hash"   \tab\tab use hashing to generate factors. Since v1.8.3 this is a fast hybrid implementation using \code{\link{group}} and radix ordering applied to the unique elements. See Details. \cr
                 4 \tab\tab "rcpp_hash" \tab\tab the previos "hash" algorithm prior to v1.8.3: uses \code{Rcpp::sugar::sort_unique} and \code{Rcpp::sugar::match}. Only supports \code{sort = TRUE}. \cr
  }
  Note that for \code{finteraction}, \code{method = "hash"} is always unsorted and \code{method = "rcpp_hash"} is not available.
}
\item{return.groups}{logical. \code{TRUE} returns the unique elements / groups / levels of \code{x} in an attribute called \code{"groups"}. Unlike \code{qF}, they are not converted to character.}
\item{factor}{logical. \code{TRUE} returns an factor, \code{FALSE} returns a 'qG' object. }
\item{\dots}{multiple atomic vectors or factors, or a single list of equal-length vectors or factors. See Details. }
}
\details{
Whenever a vector is passed to a \link[=fast-statistical-functions]{Fast Statistical Function} such as \code{fmean(mtcars, mtcars$cyl)}, is is grouped using \code{qF}, or \code{qG} if \code{use.g.names = FALSE}.

\code{qF} is a combination of \code{as.factor} and \code{factor}. Applying it to a vector i.e. \code{qF(x)} gives the same result as \code{as.factor(x)}. \code{qF(x, ordered = TRUE)} generates an ordered factor (same as \code{factor(x, ordered = TRUE)}), and \code{qF(x, na.exclude = FALSE)} generates a level for missing values (same as \code{factor(x, exclude = NULL)}). An important addition is that \code{qF(x, na.exclude = FALSE)} also adds a class 'na.included'. This prevents \emph{collapse} functions from checking missing values in the factor, and is thus computationally more efficient. Therefore factors used in grouped operations should preferably be generated using \code{qF(x, na.exclude = FALSE)}. Setting \code{sort = FALSE} gathers the levels in first-appearance order (unless \code{method = "radix"} and \code{x} is numeric, in which case the levels are always sorted). This often gives a noticeable speed improvement. % for non-numeric \code{x}.

There are 3 internal methods of computation: radix ordering, hashing, and Rcpp sugar hashing. Radix ordering is done by combining the functions \code{\link{radixorder}} and \code{\link{groupid}}. It is generally faster than hashing for large numeric data and pre-sorted data (although there are exceptions). Hashing uses \code{\link{group}}, followed by \code{\link{radixorder}} on the unique elements if \code{sort = TRUE}. It is generally fastest for character data. Rcpp hashing uses \code{Rcpp::sugar::sort_unique} and \code{Rcpp::sugar::match}. This is often less efficient than the former on large data, but the sorting properties (relying on \code{std::sort}) may be superior in borderline cases where \code{\link{radixorder}} fails to deliver exact lexicographic ordering of factor levels.  % If \code{sort = FALSE}, \code{\link{group}} is used which is generally very fast. % The hashing methods have very fast For logical data, a super fast one-pass method was written which is subsumed in the hash method.

Regarding speed: In general \code{qF} is around 5x faster than \code{as.factor} on character data and about 30x faster on numeric data. Automatic method dispatch typically does a good job delivering optimal performance.

\code{qG} is in the first place a programmers function. It generates a factor-'light' class 'qG' consisting of only an integer grouping vector and an attribute providing the number of groups. It is slightly faster and more memory efficient than \code{\link{GRP}} for grouping atomic vectors, and also convenient as it can be stored in a data frame column, which are the main reasons for its existance. %The fact that it (optionally) returns the unique groups / levels without converting them to character is an added bonus (this also provides a small performance gain compared to \code{qF}). Since v1.7, you can also call a C-level function \code{\link{group}} directly, which works for multivariate data as well, but does not sort the data and does not preserve missing values.

\code{finteraction} is simply a wrapper around \code{as_factor_GRP(GRP.default(X))}, where X is replaced by the arguments in '\dots' combined in a list (so its not really an interaction function but just a multivariate grouping converted to factor, see \code{\link{GRP}} for computational details). In general: All vectors, factors, or lists of vectors / factors passed can be interacted. Interactions always create a level for missing values and always drop unused levels.
}
\value{
\code{qF} returns an (ordered) factor. \code{qG} returns an object of class 'qG': an integer grouping vector with an attribute \code{"N.groups"} indicating the number of groups, and, if \code{return.groups = TRUE}, an attribute \code{"groups"} containing the vector of unique groups / elements in \code{x} corresponding to the integer-id. \code{finteraction} can return either.
}
\note{
An efficient alternative for character vectors with multithreading support is provided by \code{kit::charToFact}.

\code{qG(x, sort = FALSE, na.exclude = FALSE, method = "hash")} internally calls \code{\link[=group]{group(x)}} which can also be used directly and also supports multivariate groupings where \code{x} can be a list of vectors.

Neither \code{qF} nor \code{qG} reorder groups / factor levels. An exception was added in v1.7, when calling \code{qF(f, sort = FALSE)} on a factor \code{f}, the levels are recast in first appearance order. These objects can however be converted into one another using \code{qF/qG} or the direct method \code{as_factor_qG} (called inside \code{qF}). It is also possible to add a class 'ordered' (\code{ordered = TRUE}) and to create am extra level / integer for missing values (\code{na.exclude = FALSE}) if factors or 'qG' objects are passed to \code{qF} or \code{qG}. % Apart from that \code{qF} and \code{qG} don't do much to each others objects.
}
\seealso{
\code{\link{group}}, \code{\link{groupid}}, \code{\link{GRP}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}
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

finteraction(mtcars$cyl, mtcars$vs, factor = FALSE) # Returns 'qG', by default unsorted
group(mtcars[c("cyl", "vs")]) # Same thing. Use whatever syntax is more convenient

}
\keyword{manip}
