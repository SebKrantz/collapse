\name{unlist2d}
\alias{unlist2d}
\title{
Recursive Row-Binding / Unlisting in 2D - to Data Frame
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
  \item{idcols}{a character stub or a vector of names for id-columns automatically added - one for each level of nesting in \code{l}. By default the stub is \code{".id"}, so columns will be of the form \code{".id.1", ".id.2",} etc... . if \code{idcols = TRUE}, the stub is also set to \code{".id"}. If \code{idcols = FALSE}, id-columns are omitted. The content of the id columns are the list names, or (if missing) integers for the list elements. Missing elements in asymmetric nested structures are filled up with \code{NA}. See examples. }
  \item{row.names}{\code{TRUE} extracts row names from all the objects in \code{l} (where available) and adds them to the output in a column named \code{"row.names"}. Alternatively, a column name i.e. \code{row.names = "file"} can be supplied.}
  \item{recursive}{if \code{FALSE}, only process the lowest (deepest) level of \code{l}.}
  \item{id.factor}{if \code{TRUE} and \code{idcols != FALSE}, create id columns as ordered factors instead of character or integer vectors. This is useful if id's are used for further analysis e.g. as inputs to \code{ggplot2}. }
  \item{DT}{if \code{TRUE}, return a \code{data.table}, not a \code{data.frame}.}
}
\details{
The data.frame representation created by \code{unlist2d} is built as follows:
\itemize{
  \item Recurse down to the lowest level of the list-tree, data.frames are exempted and treated as a final elements.
  \item Check out the objects, if they are vectors, matrices or arrays convert them to data.frame (in the case of atomic vectors each element becomes a column).
  \item Row-bind these data.frame's using data.table's \code{Rbindlist} function. Columns are matched by name. If the number of columns differ, fill empty spaces with NA's. Create an id-column on the left, filled with the object names or indices (if unnamed). If \code{row.names = TRUE}, store row.names of the objects (if available) in a separate column.
  \item Move up to the next higher level of the list-tree and repeat: Convert atomic objects to data.frame and row-bind while matching all columns and filling unmatched ones with NA's. Create another id-column for each level of nesting passed through. If the list-tree is asymmetric, fill empty spaces in lower-level id columns with NA's.
}
The result of this iterative procedure is a single \code{data.frame} containing on the left side id-columns for each level of nesting (from higher to lower level), followed by a column containing all the row.names of the objects if \code{row.names = TRUE}, followed by the object columns, matched at each level of recursion. Optimal results are of course obtained with symmetric lists of arrays, matrices or data.frames, which \code{unlist2d} nicely converts to a beautiful data.frame ready for plotting or further analysis.
}
\value{
A \code{data.frame} or (if \code{DT = TRUE}) \code{data.table}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
% \note{
%%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rapply2d}}, \link[=list-processing]{List Processing}
}
\examples{
## basic examples:
l <- list(mtcars, list(mtcars, mtcars))
unlist2d(l)
unlist2d(rapply2d(l, fmean))
l = list(a = qM(mtcars[1:8]),
         b = list(c = mtcars[4:11], d = list(e = mtcars[2:10], f = mtcars)))
unlist2d(l, row.names = TRUE)
unlist2d(rapply2d(l, fmean))
unlist2d(rapply2d(l, fmean), recursive = FALSE)

## Groningen Growth and Development Center 10-Sector Database
head(GGDC10S) # See ?GGDC10S
namlab(GGDC10S, class = TRUE)

# Panel-Summarize this data by Variable (Emloyment and Value Added), output as list (instead of 4D array)
l <- qsu(GGDC10S, by = ~Variable,
         xt = ~ Variable + Country,
         cols = 6:16, array = FALSE)
str(l)                                         # A list of 2-levels containing matrices of statistics
head(unlist2d(l))                              # Default output, but we are missing the variables (row-names of the matrices)
head(unlist2d(l, row.names = TRUE))            # Here we go, but this is still not very nice
head(unlist2d(l, idcols = c("Sector","Trans"), # Now this is looking pretty good
              row.names = "Variable"))

dat <- unlist2d(l, c("Sector","Trans"),        # Id-columns can also be generated as ordered factors
                "Variable", id.factor = TRUE)
str(dat)

# Split this sectoral data, first by Variable (Emloyment and Value Added), then by Country
sdat <- rapply2d(split(GGDC10S[c(1,6:16)], GGDC10S$Variable), function(x) split(x[-1],x[[1]]))

# Compute pairwise correlations between sectors and recombine:
dat <- unlist2d(rapply2d(sdat, pwcor),
                idcols = c("Variable","Country"),
                row.names = "Sector")
head(dat)
plot(hclust(as.dist(1-pwcor(dat[-(1:3)]))))    # Using these correlations as a distance metric to cluster sectors

# Together with other functions like psmat, unlist2d can also effectively help reshape data:
head(unlist2d(psmat(subset(GGDC10S, Variable == "VA"), ~Country, ~Year, cols = 6:16, array = FALSE),
              idcols = "Sector", row.names = "Country"))


}

% # We can also examine the correlations of Growth rates of VA in each sector across countries
% dat <- G(subset(GGDC10S, Variable == "VA"),1,1, ~ Country, ~Year, cols = 6:16)
% dat <- psmat(dat, ~ Country, ~Year)
% plot(dat, legend = TRUE)
% dat[dat > 100] = NA       # remove outliers
% plot(dat, legend = TRUE)
% sort(apply(dat, 3, function(x) fmean.default(pwcor(x))))

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
