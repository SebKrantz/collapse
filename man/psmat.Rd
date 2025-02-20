\name{psmat}
\alias{psmat}
\alias{psmat.default}
\alias{psmat.pseries}
\alias{psmat.data.frame}
\alias{psmat.pdata.frame}
\alias{plot.psmat}
\alias{aperm.psmat}
\alias{[.psmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix / Array from Panel Series
}
\description{
\code{psmat} efficiently expands a panel-vector or 'indexed_series' ('pseries') into a matrix. If a data frame or 'indexed_frame' ('pdata.frame') is passed, \code{psmat} returns a 3D array or a list of matrices.  % By default the matrix is created such that group-identifiers constitute the rows and time the columns.
}
\usage{
psmat(x, \dots)

\method{psmat}{default}(x, g, t = NULL, transpose = FALSE, fill = NULL, \dots)

\method{psmat}{data.frame}(x, by, t = NULL, cols = NULL, transpose = FALSE, fill = NULL, array = TRUE, \dots)

# Methods for indexed data / compatibility with plm:

\method{psmat}{pseries}(x, transpose = FALSE, fill = NULL, drop.index.levels = "none", \dots)

\method{psmat}{pdata.frame}(x, cols = NULL, transpose = FALSE, fill = NULL, array = TRUE,
      drop.index.levels = "none", \dots)


\method{plot}{psmat}(x, legend = FALSE, colours = legend, labs = NULL, grid = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, indexed series 'indexed_series' ('pseries'), data frame or 'indexed_frame' ('pdata.frame').}
  \item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}. If the panel is balanced an integer indicating the number of groups can also be supplied. See Examples.}
\item{by}{\emph{data.frame method}: Same input as \code{g}, but also allows one- or two-sided formulas using the variables in \code{x}, i.e. \code{~ idvar} or \code{var1 + var2 ~ idvar1 + idvar2}.}
\item{t}{same inputs as \code{g/by}, to indicate the time-variable(s) or second identifier(s). \code{g} and \code{t} together should fully identify the panel. If \code{t = NULL}, the data is assumed sorted and \code{seq_col} is used to generate rownames for the output matrix.}
\item{cols}{\emph{data.frame method}: Select columns using a function, column names, indices or a logical vector. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
\item{transpose}{logical. \code{TRUE} generates the matrix such that \code{g/by -> columns, t -> rows}. Default is \code{g/by -> rows, t -> columns}.}
\item{fill}{element to fill empty slots of matrix / array if panel is unbalanced. \code{NULL} will generate a \code{NA} of the right type.}
\item{array}{\emph{data.frame / pdata.frame methods}: logical. \code{TRUE} returns a 3D array (if just one column is selected a matrix is returned). \code{FALSE} returns a list of matrices.}
 \item{drop.index.levels}{character. Either \code{"id"}, \code{"time"}, \code{"all"} or \code{"none"}. See \link{indexing}.}
  \item{\dots}{arguments to be passed to or from other methods, or for the plot method additional arguments passed to \code{\link{ts.plot}}.}

\item{legend}{logical. Automatically create a legend of panel-groups.}
\item{colours}{either \code{TRUE} to automatically colour by panel-groups using \code{\link{rainbow}} or a character vector of colours matching the number of panel-groups (series).}
\item{labs}{character. Provide a character-vector of variable labels / series titles when plotting an array.}
\item{grid}{logical. Calls \code{\link{grid}} to draw gridlines on the plot.}
}
\details{
If n > 2 index variables are attached to an indexed series or frame, the first n-1 variables in the index are interacted.
}
\value{
A matrix or 3D array containing the data in \code{x}, where by default the rows constitute the groups-ids (\code{g/by}) and the columns the time variable or individual ids (\code{t}). 3D arrays contain the variables in the 3rd dimension. The objects have a class 'psmat', and also a 'transpose' attribute indicating whether \code{transpose = TRUE}.
}
% \references{
%% ~put references to the literature/web site here ~
% }
% \author{
%%  ~~who you are~~
% }
\note{
The \code{pdata.frame} method only works for properly subsetted objects of class 'pdata.frame'. A list of 'pseries' won't work. There also exist simple \code{aperm} and \code{[} (subset) methods for 'psmat' objects. These differ from the default methods only by keeping the class and the 'transpose' attribute.
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\link[=time-series-panel-series]{Time Series and Panel Series}, \link[=collapse-documentation]{Collapse Overview} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## World Development Panel Data
head(wlddev)                                                    # View data
qsu(wlddev, pid = ~ iso3c, cols = 9:12, vlabels = TRUE)         # Sumarizing data
str(psmat(wlddev$PCGDP, wlddev$iso3c, wlddev$year))             # Generating matrix of GDP
r <- psmat(wlddev, PCGDP ~ iso3c, ~ year)                       # Same thing using data.frame method
plot(r, main = vlabels(wlddev)[9], xlab = "Year")               # Plot the matrix
str(r)                                                          # See srructure
str(psmat(wlddev$PCGDP, wlddev$iso3c))                          # The Data is sorted, could omit t
str(psmat(wlddev$PCGDP, 216))                                   # This panel is also balanced, so
# ..indicating the number of groups would be sufficient to obtain a matrix

ar <- psmat(wlddev, ~ iso3c, ~ year, 9:12)                      # Get array of transposed matrices
str(ar)
plot(ar)
plot(ar, legend = TRUE)
plot(psmat(collap(wlddev, ~region+year, cols = 9:12),           # More legible and fancy plot
           ~region, ~year), legend = TRUE,
     labs = vlabels(wlddev)[9:12])

psml <- psmat(wlddev, ~ iso3c, ~ year, 9:12, array = FALSE)     # This gives list of ps-matrices
head(unlist2d(psml, "Variable", "Country", id.factor = TRUE),2) # Using unlist2d, can generate DF

## Indexing simplifies things
wldi <- findex_by(wlddev, iso3c, year)  # Creating an indexed frame
PCGDP <- wldi$PCGDP                     # An indexed_series of GDP per Capita
head(psmat(PCGDP), 2)                   # Same as above, more parsimonious
plot(psmat(PCGDP))
plot(psmat(wldi[9:12]))
plot(psmat(G(wldi[9:12])))              # Here plotting panel-growth rates

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{array}
\keyword{ts}
