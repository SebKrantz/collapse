\name{collapse-package}
\alias{collapse-package}
\alias{collapse}
\docType{package}
\title{ %   \code{collapse:}
Advanced and Fast Data Transformation
}
\description{
\emph{collapse} is a C/C++ based package for data manipulation in R. It's aims are
\itemize{
\item to facilitate complex data transformation and exploration tasks in R
\item to help make R code fast, flexible, parsimonious and programmer friendly. % \emph{collapse} is a fast  %to facilitate (advanced) data manipulation in R   % To achieve the latter,
% collapse provides a broad set.. -> Nah, it's not a misc package
}
It is compatible with \emph{dplyr}, \emph{data.table} and the \emph{plm} approach to panel-data.

% \cr
\bold{Key Features:}    % (In more detail in \link[=collapse-documentation]{Collapse Overview}) % Key functionality: % Key areas Key topics addressed by \emph{collapse} are: % where \emph{collapse} offers innovative solutions are:

% \tabular{lll}{ % \Sexpr{"\u200B"} \Sexpr{"\u200B"}  \tab \emph{ Description }  \cr
\enumerate{  % (1) \tab\tab ... \cr

                \item  \emph{Advanced data programming}: A full set of fast statistical functions
        supporting grouped and weighted computations on vectors, matrices and
        data frames. Fast (ordered) and programmable grouping, and fast factor
        generation, manipulation of data frames, and data object conversions. %  select, subset, transform, replace, add and delete data.frame columns.
                 \item \emph{Advanced aggregation}: Fast and easy multi-data-type, multi-function,
        weighted, parallelized and fully customized data aggregation.
                 \item \emph{Advanced transformations}: Fast (grouped, weighted) replacing and
        sweeping out of statistics, scaling / standardizing, centering (i.e.
        between and within transformations), higher-dimensional centering
        (i.e. multiple fixed effects transformations), linear
        prediction and partialling-out.
                 \item \emph{Advanced time-computations}: Fast (sequences of) lags / leads, and
        (lagged / leaded, iterated, quasi-, log-) differences and growth
        rates on (unordered) time-series and panel data. Multivariate auto,
        partial and cross-correlation functions for panel data.
        Panel data to (ts-)array conversions.
                 \item \emph{List processing}: (Recursive) list search / identification, extraction /
        subsetting, data-apply, and generalized row-binding / unlisting in 2D.
                 \item \emph{Advanced data exploration}: Fast (grouped, weighted, panel-decomposed)
        summary statistics for complex multilevel / panel data.
% }
}
}
\section{Getting Started}{
Please see \link[=collapse-documentation]{Collapse Documentation & Overview}, or the introductory vignette. A compact but non-exhaustive set of examples is also provided below.
}
\section{Details}{
% Put this below bullet points ?? % \emph{collapse} provides a carefully conceived % \emph{collapse} provides a compact set of functions % organized roughly into several topics
\emph{collapse} provides an integrated set of functions organized into several topics (see \link[=collapse-documentation]{Collapse Overview}). Many functions are S3 generic with core methods for vectors, matrices and data.frames. Inputs are quickly passed to
compiled C/C++ code, enabling flexible and parsimonious programming in R at extreme speeds.

Broad areas of use are fast grouped programming and data manipulation to implement complex statistical techniques, and fast data transformation and exploration code (i.e. for \code{shiny} apps). Applications include fast panel data estimators and techniques, fast weighted programming (i.e. for survey techniques), fast programming with and aggregation of categorical data, fast programming with time-series and panel-series data, and programming with lists of data objects.

The package largely avoids non-standard evaluation and exports core methods for maximum programmability.  % Most are S3 generic with methods for common \code{R} objects (vectors, matrices, data.frames, ...) % high computation  %(aggregation and transformations ~10x \emph{data.table} on data <1 Mio obs.).
Smart attribute handling and additional (not-exported) methods ensure compatibility and support for \emph{dplyr}, \emph{data.table} and the \emph{plm} approach to panel-data. \emph{collapse} comes with a built-in hierarchical \link[=collapse-documentation]{documentation} facilitating the use of the package. % and \code{base/stats}. % extra methods warrant, provide ?

% Beyond speed, flexibility and parsimony in coding, a central objective of \emph{collapse} is to facilitate advanced / complex operations on data.
\emph{collapse} is mainly coded in C++ and built with \code{Rcpp}, but also uses C functions from \emph{data.table} (grouping, subsetting, row-binding), \emph{lfe} (centering on multiple factors) and \emph{stats} (ACF and PACF).
% \emph{collapse} is built with \code{Rcpp} and imports \code{C} functions from \emph{data.table}, \emph{lfe} and \emph{stats}. %, and uses \code{ggplot2} visualizations.

}
\section{Author(s)}{ % \author{
\bold{Maintainer}: Sebastian Krantz \email{sebastian.krantz@graduateinstitute.ch}

Other contributors from packages \emph{collapse} utilizes:
\itemize{
\item Matt Dowle, Arun Srinivasan and contributors worldwide (data.table)
\item Simen Gaure (lfe)
\item Dirk Eddelbuettel and contributors worldwide (Rcpp)
\item R Core Team and contributors worldwide (stats)
}
I also thank Ralf Stubner, Joseph Wood and Dirk Eddelbuettel for helpful answers on Stackoverflow,
and Joris Meys on R-Devel for encouraging me and helping to set up the \href{https://github.com/SebKrantz/collapse}{github repository} for \emph{collapse}.
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
%\seealso{
%  Optional links to other man pages
%}
\examples{
# World Bank World Development Data: 216 countries, 59 years, 4 series (columns 9-12)
head(wlddev)

# Describe data
descr(wlddev)

# Pairwise correlations with p-value
pwcor(num_vars(wlddev), P = TRUE)

# Panel-summarize columns 9 though 12 of this data (within and between countries)
qsu(wlddev, pid = ~ country, cols = 9:12, vlabels = TRUE)

# Do all of that by region and also compute higher moments -> returns a 4D array
qsu(wlddev, ~ region, ~ country, cols = 9:12, higher = TRUE)

# Return as nested list of statistics-matrices instead
suml <- qsu(wlddev, ~ region, ~ country,
            cols = 9:12, higher = TRUE, array = FALSE)
str(suml)

# Create data.frame from this list with 3 identifier columns
head(unlist2d(suml, idcols = c("Variable","Trans"), row.names = "Region"))
% # Compute the means of all the regions and create a simpler data.frame instead % unlist2d(rapply2d(suml, fmean), idcols = c("Variable","Trans"))
# Select columns from wlddev
series <- get_vars(wlddev, 9:12)     # same as wlddev[9:12] but 2x faster and works with data.tables
series <- fselect(wlddev, PCGDP:ODA) # Same thing: > 100x faster t. dplyr::select(wlddev, PCGDP:ODA)

# Replace columns, 8x faster than wlddev[9:12] <- series and also replaces names
get_vars(wlddev, 9:12) <- series

# Fast subsetting
head(fsubset(wlddev, country == "Ireland", -country, -iso3c))
head(fsubset(wlddev, country == "Ireland" & year > 1990, year, PCGDP:ODA))
ss(wlddev, 1:10, 1:10) # This is an order of magnitude faster than wlddev[1:10, 1:10]

# Fast transforming
head(ftransform(wlddev, ODA_GDP = ODA / PCGDP, ODA_LIFEEX = sqrt(ODA) / LIFEEX))
head(ftransform(wlddev, ODA_GDP = ODA / PCGDP, PCGDP = NULL, ODA = NULL, GINI_sum = fsum(GINI)))

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
fmean(series, fselect(wlddev, region, income))

# Same using GRP objects:
g <- GRP(wlddev, ~ region + income)
print(g)
plot(g)

# GRP objects are extremely efficient inputs to fast functions
fmean(series, g)
fmedian(series, g)
fsd(series, g)

# Another option is creating a grouped_df, using dplyr::group_by or the faster fgroup_by
gseries <- fgroup_by(fselect(wlddev, region, income, PCGDP:ODA), region, income)
str(gseries)
fmean(gseries)          # grouped mean
fmean(gseries, w = ODA) # weighted grouped mean, weighted by ODA
fsd(gseries, w = ODA)   # Weighted group standard deviation
\donttest{
# Faster aggregations with dplyr:
library(dplyr) # This is already a lot faster than summarize_all(mean)
wlddev \%>\% group_by(region,income) \%>\% select(PCGDP,LIFEEX) \%>\% fmean
# Now this is getting fast, apart from the pipe which still slows things down...
wlddev \%>\% fgroup_by(region,income) \%>\% fselect(PCGDP,LIFEEX) \%>\% fmean
}
# Data-Apply to columns
head(dapply(series, log))
dapply(series, quantile, na.rm = TRUE)

# Data-Apply to rows (for sum use rowSums(qM(series), na.rm = TRUE), same for rowMeans ...)
head(dapply(mtcars, max, MARGIN = 1, na.rm = TRUE))
head(dapply(mtcars, quantile, MARGIN = 1))

# qM -> quickly convert data to matrix, qDF/qDT do the reverse
fmean(rowSums(qM(series), na.rm = TRUE))

# Split-apply combine computing on columns
BY(series, wlddev$region, sum, na.rm = TRUE)  # Please use: fsum(series, wlddev$region) -> faster
BY(series, wlddev$region, quantile, na.rm = TRUE)
BY(series, wlddev$region, quantile, na.rm = TRUE, expand.wide = TRUE)

# Convert panel-data to array
psar <- psmat(wlddev, ~country, ~year, cols = 9:12)
str(psar)
psar["Ireland",,]                                       # Fast data access
psar["Ireland",,"PCGDP"]
psar[,"2016",]
qDF(psar[,"2016",], row.names.col = "Country")          # Convert to data.frame
plot(psar, colour = TRUE, labs = vlabels(wlddev)[9:12]) # Visualize
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
\donttest{
# Faster transformations with dplyr:
wlddev \%>\% fgroup_by(region,income) \%>\% fselect(PCGDP,LIFEEX,ODA) \%>\%
fwithin(ODA) # Centering using weighted means, weighted by ODA
}
## But there are also tidy transformation operators for common jobs:
# Centering (within-transforming) the 4 series by country
head(W(wlddev, ~ country, cols = 9:12))

# Same but adding overall mean back after subtracting out group means
head(W(wlddev, ~ country, cols = 9:12, mean = "overall.mean"))

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
\donttest{
# Adding exactly identified growth rates using data.table
library(data.table)
setDT(wlddev)[, paste0("G.", names(wlddev)[9:12]) := fgrowth(.SD,1,1,iso3c,year), .SDcols = 9:12]
}
# Deleting again and doing the dame thing with add_vars
get_vars(wlddev, "G1.", regex = TRUE) <- NULL
add_vars(wlddev) <- fgrowth(gv(wlddev, 9:12), 1, 1, wlddev$iso3c, wlddev$year)
get_vars(wlddev, "G1.", regex = TRUE) <- NULL

# Computing the 1- and 10-year log-differences of GDP per capita and Life-Expectancy
head(G(wlddev, c(0,1,10), 1, PCGDP + LIFEEX ~ country, ~year, logdiff = TRUE))

# Same transformations using plm package:
\donttest{
library(plm)
pwlddev <- pdata.frame(wlddev, index = c("country","year"))
head(W(pwlddev$PCGDP))                      # Country-demeaning
head(W(pwlddev, cols = 9:12))
head(W(pwlddev$PCGDP, effect = 2))          # Time-demeaning
head(W(pwlddev, effect = 2, cols = 9:12))
head(HDW(pwlddev$PCGDP))                    # Country- and time-demeaning
head(HDW(pwlddev, cols = 9:12))
head(STD(pwlddev$PCGDP))                    # Standardizing by country
head(STD(pwlddev, cols = 9:12))
head(L(pwlddev$PCGDP, -1:3))                # Panel-lags
head(L(pwlddev, -1:3, 9:12))
head(G(pwlddev$PCGDP))                      # Panel-Growth rates
head(G(pwlddev, 1, 1, 9:12))
}
}
\keyword{package}
\keyword{manip}
