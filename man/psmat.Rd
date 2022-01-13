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
\code{psmat} efficiently expands a panel-vector or \code{plm::pseries} into a matrix. If a data frame or \code{plm::pdata.frame} is passed, \code{psmat} returns (default) a 3D array or a list of matrices.  % By default the matrix is created such that group-identifiers constitute the rows and time the columns.
}
\usage{
psmat(x, \dots)

\method{psmat}{default}(x, g, t = NULL, transpose = FALSE, \dots)

\method{psmat}{pseries}(x, transpose = FALSE, \dots)

\method{psmat}{data.frame}(x, by, t = NULL, cols = NULL, transpose = FALSE, array = TRUE, \dots)

\method{psmat}{pdata.frame}(x, cols = NULL, transpose = FALSE, array = TRUE, \dots)


\method{plot}{psmat}(x, legend = FALSE, colours = legend, labs = NULL, grid = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, panel series (\code{plm::pseries}), data frame or panel data frame (\code{plm::pdata.frame}).}
  \item{g}{a factor, \code{GRP} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP} object) used to group \code{x}. If the panel is balanced an integer indicating the number of groups can also be supplied. See Examples.}
\item{by}{\emph{data.frame method}: Same input as \code{g}, but also allows one- or two-sided formulas using the variables in \code{x}, i.e. \code{~ idvar} or \code{var1 + var2 ~ idvar1 + idvar2}.}
\item{t}{same inputs as \code{g}, to indicate the time-variable(s) or second identifier(s). \code{g} and \code{t} together should fully identify the panel. If \code{t = NULL}, the data is assumed sorted and \code{seq_col} is used to generate rownames for the output matrix.}
\item{cols}{\emph{data.frame method}: Select columns using a function, column names, indices or a logical vector. \emph{Note}: \code{cols} is ignored if a two-sided formula is passed to \code{by}.}
\item{transpose}{logical. \code{TRUE} generates the matrix such that \code{g/by -> columns, t -> rows}. Default is \code{g/by -> rows, t -> columns}.}
\item{array}{\emph{data.frame / pdata.frame methods}: logical. \code{TRUE} returns a 3D array (if just one column is selected a matrix is returned). \code{FALSE} returns a list of matrices.}
  \item{\dots}{arguments to be passed to or from other methods, or for the plot method additional arguments passed to \code{ts.plot}.}

\item{legend}{logical. Automatically create a legend of panel-groups.}
\item{colours}{either \code{TRUE} to automatically colour by panel-groups using \code{\link{rainbow}} or a character vector of colours matching the number of panel-groups (series).}
\item{labs}{character. Provide a character-vector of variable labels / series titles when plotting an array.}
\item{grid}{logical. Calls \code{\link{grid}} to draw gridlines on the plot.}
}
\details{
For \code{plm::pseries}, the first \code{index} variable is taken to be the group-id and the second the time variable. If more than 2 index variables are attached to \code{plm::pseries}, the last one is taken as the time variable and the others are taken as group-id's and interacted.
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
\donttest{ % No code relying on suggested package
## Using plm simplifies things
pwlddev <- plm::pdata.frame(wlddev, index = c("iso3c","year"))  # Creating a Panel Data Frame
PCGDP <- pwlddev$PCGDP                                          # A panel-Series of GDP per Capita
head(psmat(PCGDP), 2)                                           # Same as above, more parsimonious
plot(psmat(PCGDP))
plot(psmat(pwlddev[9:12]))
plot(psmat(G(pwlddev[9:12])))                                   # Here plotting panel- growth rates
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
\keyword{array}
\keyword{ts}
