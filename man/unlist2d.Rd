\name{unlist2d}
\alias{unlist2d}
\title{
Recursive Row-Binding / Unlisting in 2-Dimensions - to Data Frame
}
\description{
\code{unlist2d} efficiently unlists lists of regular R objects (objects built up from atomic elements) and creates a data.frame representation of the list. It is a faithful 2-dimensional generalization of \code{base::unlist}, and can also be understood as a recursive generalization of \code{do.call(rbind, l)}, for lists of vectors, data.frames, arrays or heterogenous objects.
}
\usage{
unlist2d(l, idcols = ".id", row.names = FALSE, recursive = TRUE,
         id.factor = FALSE, DT = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{a unlistable list, see \code{\link{is.unlistable}}.}
  \item{idcols}{a character stub or a vector of names for id-columns automatically added - one for each level of nesting in \code{l}. By default the stub is '.id', so columns will be of the form '.id.1', '.id.2', etc... . if \code{idcols = TRUE}, the stub is also set to '.id'. If \code{idcols = FALSE}, id-columns are omitted. The content of the id columns are the list names, or (if missing) integers for the list elements. Missing elements in asymmetric nested structures are filled up with \code{NA}. See examples. }
  \item{row.names}{\code{TRUE} extracts row names from all the objects in \code{l} (where available) and adds them to the output in a column named 'row.names'. Alternatively, a column name i.e. row.names = "file" can be supplied.}
  \item{recursive}{if \code{FALSE}, only process the lowest (deepest) level of \code{l}.}
  \item{id.factor}{if \code{TRUE} and \code{idcols != FALSE}, create id columns as ordered factors instead of character or integer vectors. This is useful if id's are used for further analysis e.g. as inputs to \code{ggplot2}. }
  \item{DT}{if \code{TRUE}, return a \code{data.table}, not a \code{data.frame}.}
}
\details{
The data.frame representation created by \code{unlist2d} is built as follows:
\itemize{
  \item Recurse down to the lowest level of the list-tree, data.frames are exempted and treated as a final elements.
  \item Check out the objects, if they are vectors, matrices or arrays convert them to data.frame (in the case of atomic vectors each element becomes a column).
  \item Row-bind these data.frame's using data.table's \code{Rbindlist} function. Columns are matched by name. If the number of columns differ, fill empty spaces with NA's. Create an id-column on the left, filled with the object names or indices (if unnamed). If row.names = TRUE, store row.names of the objects (if available) in a separate column.
  \item Move up to the next higher level of the list-tree and repeat: Convert atomic objects to data.frame and row-bind while matching all columns and filling unmatched ones with NA's. Create another id-column for each level of nesting passed through. If the list-tree is asymmetric, fill empty spaces in lower-level id columns with NA's.
}
The result of this iterative procedure is a single data.frame containing on the left side id-columns for each level of nesting (from higher to lower level), followed by a column containing all the row.names of the objects if row.names = TRUE, followed by the object columns, matched at each level of recursion. Optimal results are of course obtained with symmetric lists of arrays, matrices or data.frames, which \code{unlist2d} nicely converts to a beautiful data.frame ready for plotting or further analysis.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
