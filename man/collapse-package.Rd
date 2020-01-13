\name{collapse-package}
\alias{collapse-package}
\alias{collapse}
\docType{package}
\title{ %   \code{collapse:}
Advanced and Fast Data Transformation % and Exploration % in R % and Fast  % Advanced Data Aggregation, Transformation and Exploration in R % Fast Advanced Data Transformation and Exploration % Advanced and Fast % in R % Advanced Data Aggregation, Transformation and Exploration in R % or: Advanced and Fast Data Manipulation and Exploration. Or: Advanced and Fast Data Transformation and Exploration.. Or: Advanced Data Manipulation and Exploration. Or: Advanced and Fast Data Aggregation, Transformation and Exploration. Or: Fast Advanced Data Manipulation and Exploration. (Nah, 2 adjectives is bad, and you have to earn the title Fast !!) % Advanced and Fast Data Manipulation and Exploration (Manipulation sounds so meticulous, but your package moves in broad strokes)
  % Fazit: The title is Great !!
}
\description{
\code{collapse} is a C/C++ based package for data manipulation in R. It's aims are
\enumerate{
\item to facilitate complex data transformation and exploration tasks in R
\item to help make R code fast, flexible, parsimonious and programmer friendly. % \code{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, it's not a misc package
}
}
\details{
% Put this below bullet points ?? % \code{collapse} provides a carefully conceived % \code{collapse} provides a compact set of functions % organized roughly into several topics
\code{collapse} provides an integrated set of functions organized into several topics (see \link[=collapse-documentation]{Collapse Overview}). Many are S3 generic with methods for vectors, matrices and data.frames. Inputs are quickly passed to compiled C/C++ code, enabling flexible and parsimonious coding at extreme speeds.

The package avoids non-standard evaluation and exports core methods for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data.frames, ...) % high computation  %(aggregation and transformations ~10x \code{data.table} on data <1 Mio obs.).
Smart attribute handling and additional methods ensure compatibility and support for \code{dplyr}, \code{data.table} and the \code{plm} approach to panel-data. \code{collapse} comes with a built-in hierarchical \link[=collapse-documentation]{documentation} facilitating the use of the package. % and \code{base/stats}. % extra methods warrant, provide ?

%Beyond speed, flexibility and parsimony in coding, a central objective of \code{collapse} is to facilitate advanced / complex operations on data.
\bold{Key Features:}     (More details in \link[=collapse-documentation]{Collapse Overview}) % Key functionality: % Key areas Key topics addressed by \code{collapse} are: % where \code{collapse} offers innovative solutions are:

% \tabular{lll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
\enumerate{  % (1) \tab\tab ... \cr
                \item  \emph{Advanced data programming}: A full set of fast statistical functions supporting grouped and weighted computations on vectors, matrices and data.frames. Quick (ordered) and reusable grouping, quick data-conversions, and quick select and replace variables.
                 \item \emph{Advanced aggregation}: Fast and easy multi-data-type, multi-function, weighted, parallelized and fully customized data aggregation.
                 \item \emph{Advanced transformations}: Efficient (grouped, weighted) sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out.
                 \item \emph{Advanced time-computations}: Efficient (sequences of) lags/leads, and (iterated) differences and growth rates on (unordered) time-series and panel data. Multivariate auto, partial and cross-correlation functions for panel data. Panel data to (ts-)array conversions.
                 \item \emph{List Processing}: (Recursive) list-identification, extraction / subsetting, apply, and row-binding / unlisting in 2D.
                 \item \emph{Advanced data exploration}: Fast (grouped, weighted, panel-decomposed) summary statistics for cross-sectional and complex multilevel/panel data.  % Elegant (multivariate, grouped, weighted and panel-decomposed) histograms, density, scatter and correlation plots.  \cr
% }
}
\code{collapse} is built using \code{Rcpp} and imports \code{C} functions from \code{data.table}, \code{lfe} and \code{stats}. %, and uses \code{ggplot2} visualizations.

}
\section{Getting Started}{
Please see \link[=collapse-documentation]{Collapse Documentation & Overview}, or one of the Vignettes. Some examples are also provided below.
}
\section{Author(s)}{ % \author{
\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}

Other contributors:
\itemize{
\item Matt Dowle, Arun Srinivasan and contributors worldwide (data.table)
\item Simen Gaure (lfe)
\item Dirk Eddelbuettel and contributors worldwide (Rcpp)
\item R Core Team and contributors worldwide (stats)
}
I thank Ralf Stubner, Joseph Wood and Dirk Eddelbuettel for helpful answers on Stackoverflow,
and Joris Meys on R-Devel for encouraging me and helping to set up the github page for \code{collapse}.
}
\section{Developing / Feature Requests / Bug Reporting}{
\itemize{
\item If you are interested in extending or optimizing this package, see the source code at \url{https://github.com/SebKrantz/collapse/tree/master}, fork and send pull-requests, or e-mail me.
\item Please send feature requests via e-mail.
\item Please report issues at \url{https://github.com/SebKrantz/collapse/issues} or e-mail me.
}
}

%\references{
%  This optional section can contain literature or other references for
%  background information.
%}
\keyword{ package }
%\seealso{
%  Optional links to other man pages
%}
\examples{
# World Bank World Development Data: 216 countries, 59 years, 4 series (columns 9-12)
head(wlddev)

# See variable names, labels and classes
namlab(wlddev, class = TRUE)

# Panel-summarize columns 9 though 12 of this data (within and between countries)
qsu(wlddev, pid = ~ country, cols = 9:12, vlabels = TRUE)

# Do all of that by region and also compute higher moments -> returns a 4D array
qsu(wlddev, ~ region, ~ country, cols = 9:12, higher = TRUE)

# Return as nested list of statistics-matrices instead
suml <- qsu(wlddev, ~ region, ~ country,
            cols = 9:12, higher = TRUE, array = FALSE)
str(suml)

# Create data.frame from this list with 3 identifier columns
unlist2d(suml, idcols = c("Variable","Trans"), row.names = "Region")

# Compute the means of all the regions and create a simpler data.frame instead
unlist2d(rapply2d(suml, fmean), idcols = c("Variable","Trans"))

# Select columns from wlddev: same as wlddev[9:12] but 2x faster and works with data.tables etc.
series <- get_vars(wlddev, 9:12)

# Replace columns, 3x faster than wlddev[9:12] <- series and also replaces names
get_vars(wlddev, 9:12) <- series

# Fast subset by fast conversion to data.table with qDT
\dontrun{
library(data.table)
qDT(wlddev)[country == "Ireland"]
}

# Calculating fast colum-wise statistics
fNobs(series)                   # Number of non-missing values
fmean(series)                   # means of series
fmedian(series)                 # medians of series
fmin(series)                    # mins of series

# Fast grouped statistics
fNobs(series, wlddev$region)    # regional number of obs
fmean(series, wlddev$region)    # regional means
fmedian(series, wlddev$region)  # regional medians
fsd(series, wlddev$region)      # regional standard-deviations

# Means by region and income
fmean(series, get_vars(wlddev, c("region","income")))

# Same using GRP objects:
g <- GRP(wlddev, ~ region + income)
print(g)
plot(g)

# GRP objects are extremely efficient inputs to fast functions
fmean(series, g)
fmedian(series, g)
fsd(series, g)

# Fast aggregations using dplyr:
\dontrun{
library(dplyr)
wlddev \%>\% group_by(region,income) \%>\% select(PCGDP,LIFEEX) \%>\% fmean
}

# qM -> quickly convert data to matrix, qDF/qDT do the reverse
fmean(rowSums(qM(series), na.rm = TRUE))

# Data-Apply to columns
head(dapply(series, log))
dapply(series, quantile, na.rm = TRUE)

# Data-Apply to rows
head(dapply(series, sum, MARGIN = 1, na.rm = TRUE))

# Split-apply combine computing on columns
BY(series, wlddev$region, sum, na.rm = TRUE)  # Please use: fsum(series, wlddev$region) -> faster
BY(series, wlddev$region, quantile, na.rm = TRUE)
BY(series, wlddev$region, quantile, na.rm = TRUE, expand.wide = TRUE)

# Convert panel-data to array
psar <- psmat(wlddev, ~country, ~year, cols = 9:12)
str(psar)
psar["Ireland",,]                              # Fast data access
psar["Ireland",,"PCGDP"]
psar[,"2016",]
qDF(psar[,"2016",], row.names.col = "Country") # Convert to data.frame
plot(psar)                                     # Visualize
plot(psar, colour = TRUE, labs = vlabels(wlddev)[9:12])
plot(psar[c("Brazil","India","South Africa","Russian Federation","China"),,
          c("PCGDP","LIFEEX","ODA")], legend = TRUE, labs = vlabels(wlddev)[c(9:10,12)])
plot(ts(psar["Brazil",,], 1960, 2018), main = "Brazil, 1960-2018")

# Aggregate this data by country and decade: Numeric columns with mean, categorical with mode
head(collap(wlddev, ~ country + decade, fmean, fmode))

# Multi-function aggregation of certain columns
head(collap(wlddev, ~ country + decade,
       list(fmean, fmedian, fsd),
       list(ffirst, flast), cols = c(3,9:12)))

# Customized Aggregation: Assign columns to functions
head(collap(wlddev, ~ country + decade,
       custom = list(fmean = 9:10, fsd = 9:12, flast = 3, ffirst = 6:8)))

# Fast functions can also do grouped transformations:
head(fsd(series, g, TRA = "/"))     # Scale series by region and income
head(fsum(series, g, TRA = "\%"))    # Percentages by region and income
head(fmean(series, g, TRA = "-"))   # Demean / center by region and income
head(fmedian(series, g, TRA = "-")) # De-median by region and income
gmeds <- fmedian(series, g)         # Same thing in 2 steps
head(TRA(series, gmeds, "-", g))

# Fast transformations using dplyr:
\dontrun{
library(dplyr) # (here we are demeaning PCGDP and LIFEEX using weighted means, weighted by ODA)
wlddev \%>\% group_by(region,income) \%>\% select(PCGDP,LIFEEX,ODA) \%>\% fmean(ODA, "-")
}

## But there are also specialized transformation operators for common jobs:
# Centering (Within-transforming) the 4 series by country
head(W(wlddev, ~ country, cols = 9:12))

# Same but adding overall mean back after subtracting out group means
head(W(wlddev, ~ country, cols = 9:12, add.global.mean = TRUE))

# Partialling out country and year fixed effects from 2 series (qF = quick-factor)
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year)))

# Same, adding ODA as continuous regressor
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year) + ODA))

# Standardizing (scaling and centering) by country
head(STD(wlddev, ~ country, cols = 9:12))

# Computing 1 lead and 3 lags of the 4 series: Panel-computations efficient and exactly identified
head(L(wlddev, -1:3, ~ country, ~year, cols = 9:12))

# Computing the 1- and 10-year first differences of the 4 series
head(D(wlddev, c(1,10), 1, ~ country, ~year, cols = 9:12))
head(D(wlddev, c(1,10), 1:2, ~ country, ~year, cols = 9:12))  # first and second differences
head(D(wlddev, -1:1, 1, ~ country, ~year, cols = 9:12))       # 1-year lagged and leaded FD

# Computing the 1- and 10-year growth rates of the 4 series (also keeping the level series)
head(G(wlddev, c(0,1,10), 1, ~ country, ~year, cols = 9:12))

# Adding exactly identified growth rates using data.table
\dontrun{
library(data.table)
setDT(wlddev)[, paste0("G.", names(wlddev)[9:12]) := fgrowth(.SD,1,1,iso3c,year), .SDcols = 9:12]
}

# Computing the 1- and 10-year log-differences of GDP per capita and Life-Expectancy
head(G(wlddev, c(0,1,10), 1, PCGDP + LIFEEX ~ country, ~year, logdiff = TRUE))

# Same transformations using plm package:
\dontrun{
library(plm)
pwlddev <- pdata.frame(wlddev, index = c("country","year"))
W(pwlddev$PCGDP)                      # Country-demeaning
W(pwlddev, cols = 9:12)
W(pwlddev$PCGDP, effect = 2)          # Time-demeaning
W(pwlddev, effect = 2, cols = 9:12)
HDW(pwlddev$PCGDP)                    # Country- and time-demeaning
HDW(pwlddev, cols = 9:12)
STD(pwlddev$PCGDP)                    # Standardizing by country
STD(pwlddev, cols = 9:12)
L(pwlddev$PCGDP, -1:3)                # Panel-lags
L(pwlddev, -1:3, 9:12)
G(pwlddev$PCGDP)                      # Panel-Growth rates
G(pwlddev, 1, 1, 9:12)
}
}
\keyword{package}
\keyword{manip}
