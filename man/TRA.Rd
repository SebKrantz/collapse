\name{TRA}
\alias{TRA}
\alias{setTRA}
\alias{TRA.default}
\alias{TRA.matrix}
\alias{TRA.data.frame}
\alias{TRA.grouped_df}
\title{
Transform Data by (Grouped) Replacing or Sweeping out Statistics
}
\description{
\code{TRA} is an S3 generic that efficiently transforms data by either (column-wise) replacing data values with supplied statistics or sweeping the statistics out of the data. \code{TRA} supports grouped operations and data transformation by reference, and is thus a generalization of \code{\link{sweep}}.
}
\usage{
   TRA(x, STATS, FUN = "-", ...)
setTRA(x, STATS, FUN = "-", ...) # Shorthand for invisible(TRA(..., set = TRUE))

\method{TRA}{default}(x, STATS, FUN = "-", g = NULL, set = FALSE, ...)

\method{TRA}{matrix}(x, STATS, FUN = "-", g = NULL, set = FALSE, ...)

\method{TRA}{data.frame}(x, STATS, FUN = "-", g = NULL, set = FALSE, ...)

\method{TRA}{grouped_df}(x, STATS, FUN = "-", keep.group_vars = TRUE, set = FALSE, ...)
}
\arguments{
  \item{x}{a atomic vector, matrix, data frame or grouped data frame (class 'grouped_df').}
  \item{STATS}{a matching set of summary statistics. See Details and Examples.}
  \item{FUN}{an integer or character string indicating the operation to perform. There are 11 supported operations:
    \tabular{lllll}{\emph{ Int. }   \tab\tab \emph{ String }   \tab\tab \emph{ Description }  \cr
                 0 \tab\tab "na" or "replace_na"   \tab\tab replace missing values in \code{x} \cr
                 1 \tab\tab "fill" or "replace_fill"   \tab\tab replace data and missing values in \code{x} \cr
                 2 \tab\tab "replace" \tab\tab replace data but preserve missing values in \code{x} \cr
                 3 \tab\tab "-"   \tab\tab subtract (i.e. center) \cr
                 4 \tab\tab "-+"  \tab\tab subtract group-statistics but add group-frequency weighted average of group statistics (i.e. center on overall average statistic) \cr
                 5 \tab\tab "/"   \tab\tab divide (i.e. scale. For mean-preserving scaling see also \code{\link{fscale}}) \cr
                 6 \tab\tab "\%"     \tab\tab compute percentages (i.e. divide and multiply by 100) \cr
                 7 \tab\tab "+" \tab\tab add \cr
                 8 \tab\tab "*" \tab\tab multiply \cr
                 9 \tab\tab "\%\%" \tab\tab modulus (i.e. remainder from division by  \code{STATS}) \cr
                 10 \tab\tab "-\%\%" \tab\tab subtract modulus (i.e. make data divisible by \code{STATS})
  }
  }
\item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}. Number of groups must match rows of \code{STATS}. See Details.}

\item{set}{logical. \code{TRUE} transforms data by reference i.e. performs in-place modification of the data without creating a copy.}

\item{keep.group_vars}{\emph{grouped_df method:} Logical. \code{FALSE} removes grouping variables after computation. See Details and Examples.}

\item{...}{arguments to be passed to or from other methods.}
}

\details{
Without groups (\code{g = NULL}), \code{TRA} is little more than a column based version of \code{\link{sweep}}, albeit many times more efficient. In this case all methods support an atomic vector of statistics of length \code{NCOL(x)} passed to \code{STATS}. The matrix and data frame methods also support a 1-row matrix or 1-row data frame / list, respectively. \code{TRA} always preserves all attributes of \code{x}.

With groups passed to \code{g}, \code{STATS} needs to be of the same type as \code{x} and of appropriate dimensions [such that \code{NCOL(x) == NCOL(STATS)} and \code{NROW(STATS)} equals the number of groups (i.e. the number of levels if \code{g} is a factor)]. If this condition is satisfied, \code{TRA} will assume that the first row of \code{STATS} is the set of statistics computed on the first group/level of \code{g}, the second row on the second group/level etc. and do groupwise replacing or sweeping out accordingly.

For example Let \code{x = c(1.2, 4.6, 2.5, 9.1, 8.7, 3.3)}, g is an integer vector in 3 groups \code{g = c(1,3,3,2,1,2)} and \code{STATS = fmean(x,g) = c(4.95, 6.20, 3.55)}. Then \code{out = TRA(x,STATS,"-",g) = c(-3.75,  1.05, -1.05,  2.90,  3.75, -2.90)} [same as \code{fmean(x, g, TRA = "-")}] does the equivalent of the following for-loop: \code{for(i in 1:6) out[i] = x[i] - STATS[g[i]]}.

Correct computation requires that \code{g} as used in \code{fmean} and \code{g} passed to \code{TRA} are exactly the same vector. Using \code{g = c(1,3,3,2,1,2)} for \code{fmean} and \code{g = c(3,1,1,2,3,2)} for \code{TRA} will not give the right result. The safest way of programming with \code{TRA} is thus to repeatedly employ the same factor or \code{\link{GRP}} object for all grouped computations. Atomic vectors passed to \code{g} will be converted to factors (see \code{\link{qF}}) and lists will be converted to \code{\link{GRP}} objects. This is also done by all \link[=fast-statistical-functions]{Fast Statistical Functions} and \code{\link{BY}}, thus together with these functions, \code{TRA} can also safely be used with atomic- or list-groups (as long as all functions apply sorted grouping, which is the default in \emph{collapse}). %Problems may arise if functions from other packages internally group atomic vectors or lists in a non-sorted way. [\emph{Note}: \code{as.factor} conversions are ok as this also involves sorting.]

%In contrast to the other methods, \code{TRA.grouped_df} matches column names exactly, thus \code{STATS} can be any subset of aggregated columns in \code{x} in any order, with or without grouping columns. \code{TRA.grouped_df} will transform the columns in \code{x} with their aggregated versions matched from \code{STATS} (ignoring grouping columns found in \code{x} or \code{STATS} and columns in \code{x} not found in \code{STATS}), and return \code{x} again.

If \code{x} is a grouped data frame ('grouped_df'), \code{TRA} matches the columns of \code{x} and \code{STATS} and also checks for grouping columns in \code{x} and \code{STATS}. \code{TRA.grouped_df} will then only transform those columns in \code{x} for which matching counterparts were found in \code{STATS} (exempting grouping columns) and return \code{x} again (with columns in the same order). If \code{keep.group_vars = FALSE}, the grouping columns are dropped after computation, however the "groups" attribute is not dropped (it can be removed using \code{\link[=fungroup]{fungroup()}} or \code{dplyr::ungroup()}).
}
\value{
\code{x} with columns replaced or swept out using \code{STATS}, (optionally) grouped by \code{g}.
}
\note{
In most cases there is no need to call the \code{TRA()} function, because of the TRA-argument to all \link[=fast-statistical-functions]{Fast Statistical Functions} (ensuring that the exact same grouping vector is used for computing statistics and subsequent transformation). In addition the functions \code{\link[=fbetween]{fbetween/B}} and \code{\link[=fwithin]{fwithin/W}} and \code{\link[=fscale]{fscale/STD}} provide optimized solutions for frequent scaling, centering and averaging tasks. %\code{TRA} is really a programmers function for cases when both aggregate statistics and transformed data need to be retained, or to work with more complex statistics (i.e. together with \code{\link{dapply}} or \code{\link{BY}}).
}

\seealso{
\code{\link{sweep}}, \link[=fast-statistical-functions]{Fast Statistical Functions}, \link[=data-transformations]{Data Transformations}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
v <- iris$Sepal.Length          # A numeric vector
f <- iris$Species               # A factor
dat <- num_vars(iris)           # Numeric columns
m <- qM(dat)                    # Matrix of numeric data

head(TRA(v, fmean(v)))                # Simple centering [same as fmean(v, TRA = "-") or W(v)]
head(TRA(m, fmean(m)))                # [same as sweep(m, 2, fmean(m)), fmean(m, TRA = "-") or W(m)]
head(TRA(dat, fmean(dat)))            # [same as fmean(dat, TRA = "-") or W(dat)]
head(TRA(v, fmean(v), "replace"))     # Simple replacing [same as fmean(v, TRA = "replace") or B(v)]
head(TRA(m, fmean(m), "replace"))     # [same as sweep(m, 2, fmean(m)), fmean(m, TRA = 1L) or B(m)]
head(TRA(dat, fmean(dat), "replace")) # [same as fmean(dat, TRA = "replace") or B(dat)]
head(TRA(m, fsd(m), "/"))             # Simple scaling... [same as fsd(m, TRA = "/")]...

# Note: All grouped examples also apply for v and dat...
head(TRA(m, fmean(m, f), "-", f))       # Centering [same as fmean(m, f, TRA = "-") or W(m, f)]
head(TRA(m, fmean(m, f), "replace", f)) # Replacing [same fmean(m, f, TRA = "replace") or B(m, f)]
head(TRA(m, fsd(m, f), "/", f))         # Scaling [same as fsd(m, f, TRA = "/")]

head(TRA(m, fmean(m, f), "-+", f))      # Centering on the overall mean ...
                                        # [same as fmean(m, f, TRA = "-+") or
                                        #           W(m, f, mean = "overall.mean")]
head(TRA(TRA(m, fmean(m, f), "-", f),   # Also the same thing done manually !!
     fmean(m), "+"))

# Grouped data method
library(magrittr)
iris \%>\% fgroup_by(Species) \%>\% TRA(fmean(.))
iris \%>\% fgroup_by(Species) \%>\% fmean(TRA = "-")        # Same thing
iris \%>\% fgroup_by(Species) \%>\% TRA(fmean(.)[c(2,4)])   # Only transforming 2 columns
iris \%>\% fgroup_by(Species) \%>\% TRA(fmean(.)[c(2,4)],   # Dropping species column
                                        keep.group_vars = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line
