\name{TRA}
\alias{TRA}
\alias{TRA.default}
\alias{TRA.matrix}
\alias{TRA.data.frame}
\alias{TRA.grouped_df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Data by (Groupwise) Replacing or Sweeping out Statistics
}
\description{
\code{TRA} is an S3 generic that efficiently transforms data by either (column-wise) replacing data values with supplied statistics or sweeping the statistics out of the data. \code{TRA} supports grouped operations and data frame's, and is thus a generalization of \code{base::sweep}.
}
\usage{
TRA(x, STATS, FUN = "-", ...)

\method{TRA}{default}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{matrix}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{data.frame}(x, STATS, FUN = "-", g = NULL, ...)

\method{TRA}{grouped_df}(x, STATS, FUN = "-", keep.group_vars = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a atomic vector, matrix, data frame or grouped tibble (\code{dplyr::grouped_df}).}
  \item{STATS}{a matching set of summary statistics computed on \code{x}. If \code{g = NULL} (no groups), all methods support an atomic vector of statistics of length \code{NCOL(x)}. The matrix and data.frame methods also support a 1-row matrix or 1-row data.frame/list, repsectively. If groups are supplied to \code{g}, \code{STATS} needs to be of the same type as \code{x} and of appropriate dimensions (such that \code{NCOL(x) == NCOL(STATS)} and \code{NROW(STATS)} matches the number of groups supplied to \code{g} i.e. the number of levels if \code{g} is a factor, with the first row of \code{STATS} corresponding to the first level of \code{g} etc...)}
  \item{FUN}{an integer or character string indicating the operation to perform. There are 8 supported operations:
    \tabular{lll}{\emph{ Int. } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ String } \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
                 1 \tab "replace_fill" \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab replace and overwrite missing values \cr
                 2 \tab "replace" \tab replace but preserve missing values \cr
                 3 \tab "-"   \tab subtract (center) \cr
                 4 \tab "-+"  \tab subtract group-statistics but add group-frequency weighted average of group statistics (center on global average statistic) \cr
                 5 \tab "/"   \tab divide (scale) \cr
                 6 \tab "\%"     \tab compute percentages (divide and multiply by 100) \cr
                 7 \tab "+" \tab add \cr
                 8 \tab "*" \tab multiply
  }
  }
\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to ordered factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}. Number of groups must match rows of \code{STATS}. See \code{STATS} and Details.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. Remove grouping variables after computation. In contrat to the other methods, \code{TRA.grouped_df} matches column names exactly, thus \code{STATS} can be any subset of aggregated columns in \code{x} in any order, with or without grouping columns. \code{TRA.grouped_df} will transform the columns in \code{x} with their aggregated versions matched from \code{STATS} (ignoring grouping columns found in \code{x} or \code{STATS} and columns in \code{x} not found in \code{STATS}), and return \code{x} again. If \code{keep.group_vars = FALSE}, \code{x} is retured again without grouping columns. See Details and Examples.}

\item{...}{arguments to be passed to or from other methods.}
}

\details{
Without groups (\code{g = NULL}), \code{TRA} is nothing more than a column based version of \code{base::sweep}, albeit 4-times more efficient on matrices and many times more efficient on data frames. \code{TRA} always preserves all attributes of \code{x}.

With groups passed to \code{g}, \code{TRA} expects (and checks for) a set of statistics such that \code{NROW(STATS)} equals the number of groups. If this condition is satisfied, \code{TRA} will (without further checks) assume that the first row of \code{STATS} is the set of statistics computed on the first group of \code{g}, the second row on the second group etc. and do groupwise replacing or sweeping out accordingly.

For example Let \code{x = c(1.2, 4.6, 2.5, 9.1, 8.7, 3.3)}, g is an integer vector in 3 groups \code{g = c(1,3,3,2,1,2)} and \code{STATS = fmean(x,g) = c(4.95, 6.20, 3.55)}. Then \code{out = TRA(x,fmean(x,g),"-",g) = c(-3.75,  1.05, -1.05,  2.90,  3.75, -2.90)} (same as \code{fmean(x, g, TRA = "-")}) does the equavalent to the following for-loop: \code{for(i in 1:6) out[i] = x[i] - fmean(x,g)[g[i]]}.

Correct computation requires that \code{g} as used in \code{fmean} and \code{g} passed to \code{TRA} are exactly the same vector. Using \code{g = c(1,3,3,2,1,2)} for \code{fmean} and \code{g = c(3,1,1,2,3,2)} for \code{TRA} will not give the right result. The safest way of programming with \code{TRA} is thus to repeatedly employ the same factor or \code{\link{GRP}} object for all grouped computations. Atomic vectors passed to \code{g} will be converted to ordered factors (see \code{\link{qF}}) and lists will be converted to ordered \code{\link{GRP}} objects. This is also done by all \link[=fast-statistical-functions]{Fast Statistical Functions} and by default by \code{\link{BY}}, thus together with these functions, \code{TRA} can also safely be used with atomic- or list-groups. Problems may arise if other functions internally convert atomic vectors or lists to groups in a non-sorted way. \emph{Note}: \code{as.factor} conversions are ok as this also involves sorting.

If \code{x} is a grouped tibble (\code{grouped_df}), \code{TRA} matches the columns of \code{x} and \code{STATS} and also checks for grouping columns \code{names(attr(x, "groups"))} in \code{x} and \code{STATS}. \code{TRA.grouped_df} will then only transform those columns in \code{x} for which matching counterparts were found in \code{STATS}, excempting grouping columns, and returns \code{x} again (with columns in the name order). If \code{keep.group_vars = FALSE}, the grouping columns are dropped after computation, however the "groups" attribute is not dropped (it can be removed using \code{dplyr::ungroup()}).
}
\value{
\code{x} with columns replaced or swept out using \code{STATS}, grouped by \code{g}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
\note{
I have tried to make \code{TRA} as redundant as possible by adding a TRA-argument to all \link[=fast-statistical-functions]{Fast Statistical Functions} (ensuring that the exact same grouping vector is used for aggregation and transformation), and by creating the \code{\link{B}} (between-transform) and \code{\link{W}} (within-transform) operators as well as \code{\link[=fscale]{fscale/STD}} as safer and faster solutions for frequent scaling, centering and averaging tasks. \code{TRA} is really a programmers function for cases when both aggregate statistics and transformed data need to be retained, or to work with more complex statistics (i.e. together with \code{\link{dapply}} or \code{\link{BY}}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dapply}}, \code{\link{BY}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}
}
\examples{
v <- iris$Sepal.Length          # A numeric vector
f <- iris$Species               # A factor
dat <- num_vars(iris)           # Numeric columns
m <- qM(dat)                    # Matrix of numeric data

head(TRA(v, fmean(v)))                # Simple (column) centering     [same as fmean(v, TRA = "-") or W(v)]
head(TRA(m, fmean(m)))                # same as sweep(m, 2, fmean(m)) [same as fmean(m, TRA = "-") or W(m)]
head(TRA(dat, fmean(dat)))            # [same as fmean(dat, TRA = "-") or W(dat)]
head(TRA(v, fmean(v), "replace"))     # Simple (column) replacing     [same as fmean(v, TRA = "replace") or B(v)]
head(TRA(m, fmean(m), "replace"))     # same as sweep(m, 2, fmean(m)) [same as fmean(m, TRA = "replace") or B(m)]
head(TRA(dat, fmean(dat), "replace")) # [same as fmean(dat, TRA = "replace") or B(dat)]
head(TRA(m, fsd(m), "/"))             # Simple (column) scaling... [same as fsd(m, TRA = "/")]...

# Note: All grouped examples also apply for v and dat...
head(TRA(m, fmean(m, f), "-", f))       # Grouped (column) centering... [same as fmean(m, f, TRA = "-") or W(m, f)]...
head(TRA(m, fmean(m, f), "replace", f)) # Grouped (column) replacing... [same as fmean(m, f, TRA = "replace") or B(m, f)]...
head(TRA(m, fsd(m, f), "/", f))         # Grouped (column) scaling...   [same as fsd(m, f, TRA = "/")]...

head(TRA(m, fmean(m, f), "-+", f))      # Grouped (column) centering on the overall mean... [same as fmean(m, f, TRA = "-+") or W(m, f, add.global.mean = TRUE)]
head(TRA(TRA(m, fmean(m, f), "-", f),   # Also the same thing done manually !!
     fmean(m), "+"))

# grouped tibble method
library(dplyr)
iris \%>\%  group_by(Species) \%>\% TRA(.,fmean(.))
iris \%>\%  group_by(Species) \%>\% fmean(TRA = "-")        # Same thing
iris \%>\%  group_by(Species) \%>\% TRA(.,fmean(.)[c(2,4)]) # Only transforming 2 columns
iris \%>\%  group_by(Species) \%>\% TRA(.,fmean(.)[c(2,4)], # Dropping species column
                            keep.group_vars = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
