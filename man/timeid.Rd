\name{timeid}
\alias{timeid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Integer-Id From Time/Date Sequences
}
\description{
\code{timeid} groups time vectors in a way that preserves the temporal structure. It generate an integer id where unit steps represent the greatest common divisor in the original sequence e.g \code{c(4, 6, 10) -> c(1, 2, 4)} or \code{c(0.25, 0.75, 1) -> c(1, 3, 4)}.
}
\usage{
timeid(x, factor = FALSE, ordered = factor, extra = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric time object such as a \code{Date}, \code{POSIXct} or other integer or double vector representing time.}
  \item{factor}{logical. \code{TRUE} returns an (ordered) factor with levels corresponding to the full sequence (without irregular gaps) of time. This is useful for inclusion in the \link[=findex]{index} but might be computationally expensive for long sequences, see Details. \code{FALSE} returns a simpler object of class '\code{\link{qG}}'. }
  \item{ordered}{logical. \code{TRUE} adds a class 'ordered'. }
  \item{extra}{logical. \code{TRUE} attaches a set of 4 diagnostic items as attributes to the result:
   \itemize{
   \item \code{"unique_ints"}: \code{unique(unattrib(timeid(x)))} - the unique integer time steps in first-appearance order. This can be useful to check the size of gaps in the sequence. %The \code{\link{seqid}} function can help in the exploration of this attribute, e.g. \code{seqid(attr(timeid(x, extra = TRUE), "unique"))} shows the number and position of the dicontinuities.
   \item \code{"sort_unique_x"}: \code{sort(unique(x))}.
   \item \code{"range_x"}: \code{range(x)}.
   \item \code{"step_x"}: \code{vgcd(sort(unique(diff(sort(unique(x))))))} - the greatest common divisor.
   }
   \emph{Note} that returning these attributes does not incur additional computations.
  }
}
\details{
Let \code{range_x} and \code{step_x} be the like-named attributes returned when \code{extra = TRUE}, then, if \code{factor = TRUE}, a complete sequence of levels is generated as \code{seq(range_x[1], range_x[2], by = step_x) |> copyMostAttrib(x) |> as.character()}. If \code{factor = FALSE}, the number of timesteps recorded in the \code{"N.groups"} attribute is computed as \code{(range_x[2]-range_x[1])/step_x + 1}, which is equal to the number of factor levels. In both cases the underlying integer id is the same and preserves gaps in time. Large gaps (strong irregularity) can result in many unused factor levels, the generation of which can become expensive. Using \code{factor = FALSE} (the default) is thus more efficient.
}
\value{
A factor or '\code{\link{qG}}' object, optionally with additional attributes attached.
}

\seealso{
\code{\link{seqid}}, \link[=indexing]{Indexing}, \link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# A normal use case
timeid(wlddev$decade)
timeid(wlddev$decade, factor = TRUE)
timeid(wlddev$decade, extra = TRUE)

# Here a large number of levels is generated, which is expensive
timeid(wlddev$date, factor = TRUE)
tid <- timeid(wlddev$date, extra = TRUE) # Much faster
str(tid)

# The reason for step = 1 are leap years with 366 days every 4 years
diff(attr(tid, "unique"))

# So in this case simple factor generation gives a better result
qF(wlddev$date, ordered = TRUE, na.exclude = FALSE)

# The best way to deal with this data would be to convert it
# to zoo::yearmon and then use timeid:
timeid(zoo::as.yearmon(wlddev$date), factor = TRUE, extra = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
