\name{psmat}
\alias{psmat}
\alias{psmat.default}
\alias{psmat.pseries}
\alias{psmat.data.frame}
\alias{psmat.pdata.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix / Array from Panel-Series
}
\description{
\code{psmat} efficiently expands a panel-vector or \code{plm::pseries} into a matrix. If a data frame or \code{plm::pdata.frame} is passed, \code{psmat} returns a list of such matrices or a 3D array.  % By default the matrix is created such that group-identifiers constitute the rows and time the columns.
}
\usage{
\method{psmat}{default}(x, g, t = NULL, transpose = FALSE)

\method{psmat}{pseries}(x, transpose = FALSE)

\method{psmat}{data.frame}(x, by, t = NULL, cols = NULL, transpose = FALSE, simplify = FALSE)

\method{psmat}{pdata.frame}(x, transpose = FALSE, simplify = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, panel-series (\code{plm::pseries}), data.frame or panel-data.frame (\code{plm::pdata.frame}).}
  \item{g}{a factor, \code{GRP()} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{GRP()} object and then interacted) used to group \code{x}. If the panel is balanced an integer indicating the number of groups can also be supplied. See Examples.}
\item{by}{\emph{data frame method}: Same input as \code{g}, but also allows one- or two-sided formulas using the variables in \code{x}, i.e. \code{~ idvar} or \code{var1 + var2 ~ idvar1 + idvar2}.}
  \item{t}{same inputs as \code{g}, to indicate the time-variable or second identifier(s). \code{g} and \code{t} together should fully identify the panel. If \code{t = NULL}, the data is assumed sorted and \code{seq_col} is used to generate rownames.}
    \item{cols}{\emph{data.frame method}: Select columns using a function, column names or indices. \emph{Note}: \code{cols} is ignored if a two-sided formula is passsed to \code{by}.}
  \item{transpose}{logiccal. \code{TRUE} generates the matrix such that \code{g -> columns, t -> rows}. Default is \code{g -> rows, t -> columns}.}
  \item{simplify}{\emph{data.frame / pdata.frame methods}: logical. \code{TRUE} returns a 3D array (if just one column is selected a matrix is returned). Default is to always return a list of matrices.}
}
\details{
For \code{plm::pseries}, the first \code{index} variable is taken to be the group-id and the second the time variable. If more than 2 index variables are attached to \code{plm::pseries}, the last one is taken as the time variable and the others are taken as group-id's and interacted.
}
\value{
a matrix containing the data in \code{x}, where by default the rows constitute the groups-ids (\code{g}) and the columns the time variable or individual ids (\code{t}).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
The \code{pdata.frame} method only works for properly subsetted objects of class 'pdata.frame'. A list of 'pseries' won't work.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## World Development Panel Data
View(wlddev)                                                    # see also help(wlddev)
View(psmat(wlddev$PCGDP, wlddev$country, wlddev$year))          # Generating a matrix of GDP per Capita
psmat(wlddev, PCGDP ~ country, ~ year, simplify = TRUE)         # Same thing using data.frame method
View(psmat(wlddev$PCGDP, wlddev$country))                       # The Data is sorted, so we could leave out t
View(psmat(wlddev$PCGDP, 216))                                  # This panel is also balanced, so indicating the number of groups would be sufficient to obtain a matrix
ts.plot(psmat(wlddev$PCGDP, 216, transpose = TRUE))             # To plot this data, we need the transpose

tsml <- psmat(wlddev, ~ country, ~ year, 9:12, TRUE)            # Get list of transposed matrices from all 4 series in this dataset
par(mfrow = c(2,2))                                             # Plot the 4 matrices
for(n in names(tsml)) ts.plot(ts(tsml[[n]], 1960, 2018), main = n)
par(mfrow = c(1,1))
str(simplify2array(tsml))                                       # This generates 3D array of the data, same as setting option simplify = TRUE
View(unlist2d(tsml, "Variable", "Year", id.factor = TRUE))      # Using unlist2d, we can also generate a nice data.frame from this data

## Using plm:
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))# Creating a Panel-Data Frame of this data
PCGDP <- pwlddev$PCGDP                                          # A panel-Series of GDP per Capita
psmat(PCGDP)                                                    # Same as above, but more parsimonious and much faster than plm::as.matrix.pseries
ts.plot(psmat(PCGDP, TRUE))
View(unlist2d(psmat(pwlddev[9:12]), "Variable", "Country"))     # Again creating a data.frame, this time from untransposed matrices


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
