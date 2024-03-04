\name{seqid}
\alias{seqid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Group-Id from Integer Sequences
}
\description{
\code{seqid} can be used to group sequences of integers in a vector, e.g. \code{seqid(c(1:3, 5:7))} becomes \code{c(rep(1,3), rep(2,3))}. It also supports increments \code{> 1}, unordered sequences, and missing values in the sequence.

Some applications are to facilitate identification of, and grouped operations on, (irregular) time series and panels.
}
\usage{
seqid(x, o = NULL, del = 1L, start = 1L, na.skip = FALSE,
      skip.seq = FALSE, check.o = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a factor or integer vector. Numeric vectors will be converted to integer i.e. rounded downwards.}

  \item{o}{an (optional) integer ordering vector specifying the order by which to pass through \code{x}.}

  \item{del}{integer. The integer deliminating two consecutive points in a sequence. \code{del = 1} lets \code{seqid} track sequences of the form \code{c(1,2,3,..)}, \code{del = 2} tracks sequences \code{c(1,3,5,..)} etc.}

  \item{start}{integer. The starting value of the resulting sequence id. Default is starting from 1. }

  \item{na.skip}{logical. \code{TRUE} skips missing values in the sequence. The default behavior is skipping such that \code{seqid(c(1, NA, 2))} is regarded as one sequence and coded as \code{c(1, NA, 1)}.}

 \item{skip.seq}{logical. If \code{na.skip = TRUE}, this changes the behavior such that missing values are viewed as part of the sequence, i.e. \code{seqid(c(1, NA, 3))} is regarded as one sequence and coded as \code{c(1, NA, 1)}.}

 \item{check.o}{logical. Programmers option: \code{FALSE} prevents checking that each element of \code{o} is in the range \code{[1, length(x)]}, it only checks the length of \code{o}. This gives some extra speed, but will terminate R if any element of \code{o} is too large or too small. }

}
\details{
\code{seqid} was created primarily as a workaround to deal with problems of computing lagged values, differences and growth rates on irregularly spaced time series and panels before \emph{collapse} version 1.5.0 (\href{https://github.com/SebKrantz/collapse/issues/26}{#26}). Now \code{flag}, \code{fdiff} and \code{fgrowth} natively support irregular data so this workaround is superfluous, except for iterated differencing which is not yet supported with irregular data.

% panels because they do not pre-compute an ordering of the data but directly compute the ordering from the supplied id and time variables while providing errors for gaps and repeated time values. see \code{\link{flag}} for computational details.

The theory of the workaround was to express an irregular time series or panel series as a regular panel series with a group-id created such that the time-periods within each group are consecutive. \code{seqid} makes this very easy: For an irregular panel with some gaps or repeated values in the time variable, an appropriate id variable can be generated using \code{settransform(data, newid = seqid(time, radixorder(id, time)))}. Lags can then be computed using \code{L(data, 1, ~newid, ~time)} etc.

%A simple solution to applying existing functionality (\code{flag}, \code{fdiff} and \code{fgrowth}) to irregular time series and panels is thus to create a group-id that fully identifies the data together with the time variable.
% This way \emph{collapse} maintains a balance between offering very fast computations on regular time series and panels (which may be unbalanced but where observations for each entity are consecutive in time), and flexibility of application.

In general, for any regularly spaced panel the identity given by \code{identical(groupid(id, order(id, time)), seqid(time, order(id, time)))} should hold.

% Regularly spaced panels with gaps in time (such as a panel-survey with measurements every 2 years) can be handled either by \code{seqid(\dots, del = gap)} or, in most cases, simply by converting the time variable to factor using \code{\link{qF}}, which will make observations consecutive.

% \enumerate{
% \item Sort the data in ascending order (e.g. using \code{data.table::setorder(data, time)} for time series and \code{data.table::setorder(data, id, time)} for panels)
% \item Generate a new id variable using \code{seqid} (e.g. \code{settransform(data, newid = seqid(time))})
% \item Use the new id to identify the data together with the time variable (e.g. compute a panel-lag using \code{L(data, 1, ~newid, ~time)} or create a panel data frame: \code{pdata <- plm::pdata.frame(data, index = c("newid", "time")); L(pdata)})
% }
%There are potentially other more analytical applications for \code{seqid}\dots

For the opposite operation of creating a new time-variable that is consecutive in each group, see \code{data.table::rowid}.


}
\value{
An integer vector of class 'qG'. See \code{\link{qG}}.
}

\seealso{
\code{\link{timeid}}, \code{\link{groupid}}, \code{\link{qG}}, \link[=fast-grouping-ordering]{Fast Grouping and Ordering}, \link[=collapse-documentation]{Collapse Overview}

}
\examples{
## This creates an irregularly spaced panel, with a gap in time for id = 2
data <- data.frame(id = rep(1:3, each = 4),
                   time = c(1:4, 1:2, 4:5, 1:4),
                   value = rnorm(12))
data

## This gave a gaps in time error previous to collapse 1.5.0
L(data, 1, value ~ id, ~time)

## Generating new id variable (here seqid(time) would suffice as data is sorted)
settransform(data, newid = seqid(time, order(id, time)))
data

## Lag the panel this way
L(data, 1, value ~ newid, ~time)

## A different possibility: Creating a consecutive time variable
settransform(data, newtime = data.table::rowid(id))
data
L(data, 1, value ~ id, ~newtime)

## With sorted data, the time variable can also just be omitted..
L(data, 1, value ~ id)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
\keyword{manip}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
