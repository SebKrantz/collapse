\name{rowbind}
\alias{rowbind}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Row-Bind Lists / Data Frame-Like Objects
}
\description{
\emph{collapse}'s version of \code{data.table::rbindlist} and \code{rbind.data.frame}. The core code is copied from \emph{data.table}, which deserves all credit for the implementation. \code{rowbind} only binds lists/data.frame's. For a more flexible recursive version see \code{\link{unlist2d}}. To combine lists column-wise see \code{\link{add_vars}} or \code{\link{ftransform}} (with replacement).
}
\usage{
rowbind(\dots, idcol = NULL, row.names = FALSE,
        use.names = TRUE, fill = FALSE, id.factor = "auto",
        return = c("as.first", "data.frame", "data.table", "tibble", "list"))
}
\arguments{
  \item{\dots}{a single list of list-like objects (data.frames) or comma separated objects (internally assembled using \code{list(\dots)}). Names can be supplied if \code{!is.null(idcol)}.}

  \item{idcol}{character. The name of an id-column to be generated identifying the source of rows in the final object. Using \code{idcol = TRUE} will set the name to \code{".id"}. If the input list has names, these will form the content of the id column, otherwise integers are used. To save memory, it is advised to keep \code{id.factor = TRUE}.}

    \item{row.names}{\code{TRUE} extracts row names from all the objects in \code{l} and adds them to the output in a column named \code{"row.names"}. Alternatively, a column name i.e. \code{row.names = "variable"} can be supplied. }

  \item{use.names}{logical. \code{TRUE} binds by matching column name, \code{FALSE} by position. }

  \item{fill}{logical. \code{TRUE} fills missing columns with NAs. When \code{TRUE}, \code{use.names} is set to \code{TRUE}.}

  \item{id.factor}{if \code{TRUE} and \code{!isFALSE(idcols)}, create id column as factor instead of character or integer vector. It is also possible to specify \code{"ordered"} to generate an ordered factor id. \code{"auto"} uses \code{TRUE} if \code{!is.null(names(l))} where \code{l} is the input list (because factors are much more memory efficient than character vectors). }

\item{return}{an integer or string specifying what to return. \code{1 - "as.first"} preserves the attributes of the first element of the list, \code{2/3/4 - "data.frame"/"data.table"/"tibble"} coerces to specific objects, and \code{5 - "list"} returns a (named) list. }

}

\value{
a long list or data frame-like object formed by combining the rows / elements of the input objects. The \code{return} argument controls the exact format of the output.
}


\seealso{
\code{\link{unlist2d}}, \code{\link{add_vars}}, \code{\link{ftransform}}, \link[=fast-data-manipulation]{Data Frame Manipulation}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
# These are the same
rowbind(mtcars, mtcars)
rowbind(list(mtcars, mtcars))

# With id column
rowbind(mtcars, mtcars, idcol = "id")
rowbind(a = mtcars, b = mtcars, idcol = "id")

# With saving row-names
rowbind(mtcars, mtcars, row.names = "cars")
rowbind(a = mtcars, b = mtcars, idcol = "id", row.names = "cars")

# Filling up columns
rowbind(mtcars, mtcars[2:8], fill = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
