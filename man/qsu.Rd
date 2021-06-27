\name{qsu}
\alias{qsu}
\alias{qsu.default}
\alias{qsu.matrix}
\alias{qsu.data.frame}
\alias{qsu.pseries}
\alias{qsu.pdata.frame}
\alias{qsu.sf}
\alias{print.qsu}
% - Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast (Grouped, Weighted) Summary Statistics for Cross-Sectional and Panel Data
}
\description{
\code{qsu}, shorthand for quick-summary, is an extremely fast summary command inspired by the (xt)summarize command in the STATA statistical software.

It computes a set of 7 statistics (nobs, mean, sd, min, max, skewness and kurtosis) using a numerically stable one-pass method generalized from Welford's Algorithm. Statistics can be computed weighted, by groups, and also within-and between entities (for panel data, see Details).
}
\usage{
qsu(x, \dots)

\method{qsu}{default}(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, \dots)

\method{qsu}{matrix}(x, g = NULL, pid = NULL, w = NULL, higher = FALSE, array = TRUE, \dots)

\method{qsu}{data.frame}(x, by = NULL, pid = NULL, w = NULL, cols = NULL,
    higher = FALSE, array = TRUE, vlabels = FALSE, \dots)

# Methods for compatibility with plm:

\method{qsu}{pseries}(x, g = NULL, w = NULL, effect = 1L, higher = FALSE, array = TRUE, \dots)

\method{qsu}{pdata.frame}(x, by = NULL, w = NULL, cols = NULL, effect = 1L,
    higher = FALSE, array = TRUE, vlabels = FALSE, \dots)

# Methods for compatibility with sf:

\method{qsu}{sf}(x, by = NULL, pid = NULL, w = NULL, cols = NULL,
    higher = FALSE, array = TRUE, vlabels = FALSE, \dots)


\method{print}{qsu}(x, digits = 4, nonsci.digits = 9, na.print = "-",
      return = FALSE, print.gap = 2, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector, matrix, data frame, panel series (\code{plm::pseries}) or panel data frame (\code{plm::pdata.frame}).}
  \item{g}{a factor, \code{\link{GRP}} object, atomic vector (internally converted to factor) or a list of vectors / factors (internally converted to a \code{\link{GRP}} object) used to group \code{x}.}
  \item{by}{\emph{(p)data.frame method}: Same as \code{g}, but also allows one- or two-sided formulas i.e. \code{~ group1 + group2} or \code{var1 + var2 ~ group1 + group2}. See Examples.}
  \item{pid}{same input as \code{g/by}: Specify a panel-identifier to also compute statistics on between- and within- transformed data. Data frame method also supports one- or two-sided formulas. Transformations are taken independently from grouping with \code{g/by} (grouped statistics are computed on the transformed data if \code{g/by} is also used). However, passing any LHS variables to \code{pid} will overwrite any \code{LHS} variables passed to \code{by}.}
  \item{w}{a vector of (non-negative) weights. Adding weights will compute the weighted mean, sd, skewness and kurtosis, and transform the data using weighted individual means if \code{pid} is used.}
  \item{cols}{select columns to summarize using column names, indices, a logical vector or a function (e.g. \code{is.numeric}). Two-sided formulas passed to \code{by} or \code{pid} overwrite \code{cols}.}

  \item{higher}{logical. Add higher moments (skewness and kurtosis).}
  \item{array}{logical. If computations have more than 2 dimensions (up to a maximum of 4D: variables, statistics, groups and panel-decomposition) output to array, else output (nested) list of matrices.}
  \item{vlabels}{logical. Use variable labels in the summary. See \code{\link{vlabels}}.}
  \item{effect}{\emph{plm} methods: Select which panel identifier should be used for between and within transformations of the data. 1L takes the first variable in the \code{plm::index}, 2L the second etc.. Index variables can also be called by name using a character string. More than one variable can be supplied. }
  \item{\dots}{arguments to be passed to or from other methods.}
  \item{digits}{the number of digits to print after the comma/dot.}
  \item{nonsci.digits}{the number of digits to print before resorting to scientific notation (default is to print out numbers with up to 9 digits and print larger numbers scientifically).}
  \item{na.print}{character string to substitute for missing values.}
  \item{return}{logical. Don't print but instead return the formatted object.}
  \item{print.gap}{integer. Spacing between printed columns. Passed to \code{print.default}.}
}
\details{
The algorithm used to compute statistics is well described \href{https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance}{here} [see sections \emph{Welford's online algorithm}, \emph{Weighted incremental algorithm} and \emph{Higher-order statistics}. Skewness and kurtosis are calculated as described in \emph{Higher-order statistics} and are mathematically identical to those implemented in the \emph{moments} package. Just note that \code{qsu} computes the kurtosis (like \code{momens::kurtosis}), not the excess-kurtosis (=  kurtosis - 3) defined in \emph{Higher-order statistics}. The \emph{Weighted incremental algorithm} described can easily be generalized to higher-order statistics].

Grouped computations specified with \code{g/by} are carried out extremely efficiently as in \code{fsum} (in a single pass, without splitting the data).

If \code{pid} is used, \code{qsu} performs a panel-decomposition of each variable and computes 3 sets of statistics: Statistics computed on the 'Overall' (raw) data, statistics computed on the 'Between' - transformed (pid - averaged) data, and statistics computed on the 'Within' - transformed (pid - demeaned) data.

More formally, let \bold{\code{x}} (bold) be a panel vector of data for \code{N} individuals indexed by \code{i}, recorded for \code{T} periods, indexed by \code{t}. \code{xit} then denotes a single data-point belonging to individual \code{i} in time-period \code{t} (\code{t/T} must not represent time). Then \code{xi.} denotes the average of all values for individual \code{i} (averaged over \code{t}), and by extension \bold{\code{xN.}} is the vector (length \code{N}) of such averages for all individuals. If no groups are supplied to \code{g/by}, the 'Between' statistics are computed on \bold{\code{xN.}}, the vector of individual averages. (This means that for a non-balanced panel or in the presence of missing values, the 'Overall' mean computed on \bold{\code{x}} can be slightly different than the 'Between' mean computed on \bold{\code{xN.}}). If groups are supplied to \code{g/by}, \bold{\code{xN.}} is expanded to the vector \bold{\code{xi.}} (length \code{N x T}) by replacing each value \code{xit} in \bold{\code{x}} with \code{xi.}, while preserving missing values in \bold{\code{x}}. Grouped Between-statistics are then computed on \bold{\code{xi.}}, with the only difference that the number of observations ('Between-N') reported for each group is the number of distinct non-missing values of \bold{\code{xi.}} in each group (not the total number of non-missing values of \bold{\code{xi.}} in each group, which is already reported in 'Overall-N').

'Within' statistics are always computed on the vector \bold{\code{x - xi. + x..}}, where \bold{\code{x..}} is simply the 'Overall' mean computed from \bold{\code{x}}, which is added back to preserve the level of the data. The 'Within' mean computed on this data will always be identical to the 'Overall' mean. In the summary output, \code{qsu} reports not 'N', which would be identical to the 'Overall-N', but 'T', the average number of time-periods of data available for each individual obtained as 'T' = 'Overall-N / 'Between-N'. See Examples.

Apart from 'N/T' and the extrema, the standard-deviations ('SD') computed on between- and within- transformed data are extremely valuable because they indicate how much of the variation in a panel-variable is between-individuals and how much of the variation is within-individuals (over time). At the extremes, variables that have common values across individuals (such as the time-variable(s) 't' in a balanced panel), can readily be identified as individual-invariant because the 'Between-SD' on this variable is 0 and the 'Within-SD' is equal to the 'Overall-SD'. Analogous, time-invariant individual characteristics (such as the individual-id 'i') have a 0 'Within-SD' and a 'Between-SD' equal to the 'Overall-SD'.

\code{qsu} comes with it's own print method which by default writes out up to 9 digits at 4 decimal places. Larger numbers are printed in scientific format. for numbers between 7 and 9 digits, a comma ',' is placed after the 6th digit to designate the millions. Missing values are printed using '-'.

The \emph{sf} method simply ignores the geometry column.
}
\value{
A vector, matrix, array or list of matrices of summary statistics. All matrices and arrays have a class 'qsu' and a class 'table' attached.
}
\references{
Welford, B. P. (1962). Note on a method for calculating corrected sums of squares and products. \emph{Technometrics}. 4 (3): 419-420. doi:10.2307/1266577.
}
% \author{
%%  ~~who you are~~
% }
\note{
If weights \code{w} are used together with \code{pid}, transformed data is computed using weighted individual means i.e. weighted \bold{\code{xi.}} and weighted \bold{\code{x..}}. Weighted statistics are subsequently computed on this weighted-transformed data.
}

%% ~Make other sections like Warning with \section{Warning }{\dots.} ~

\seealso{
\code{\link{descr}}, \code{\link{pwcor}}, \link[=A1-fast-statistical-functions]{Fast Statistical Functions}, \link[=collapse-documentation]{Collapse Overview}  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## World Development Panel Data
# Simple Summaries -------------------------
qsu(wlddev)                                 # Simple summary
qsu(wlddev, vlabels = TRUE)                 # Display variable labels
qsu(wlddev, higher = TRUE)                  # Add skewness and kurtosis

# Grouped Summaries ------------------------
qsu(wlddev, ~ region, vlabels = TRUE)       # Statistics by World Bank Region
qsu(wlddev, PCGDP + LIFEEX ~ income)        # Summarize GDP per Capita and Life Expectancy by
stats <- qsu(wlddev, ~ region + income,     # World Bank Income Level
             cols = 9:10, higher = TRUE)    # Same variables, by both region and income
aperm(stats)                                # A different perspective on the same stats

# Panel Data Summaries ---------------------
qsu(wlddev, pid = ~ iso3c, vlabels = TRUE)   # Adding between and within countries statistics
# -> They show amongst other things that year and decade are individual-invariant,
# that we have GINI-data on only 161 countries, with only 8.42 observations per country on average,
# and that GDP, LIFEEX and GINI vary more between-countries, but ODA received varies more within
# countries over time.
\donttest{ % No code relying on suggested package
# Using plm:
pwlddev <- plm::pdata.frame(wlddev,         # Creating a Panel Data Frame frame from this data
                  index = c("iso3c","year"))
qsu(pwlddev)                                # Summary for pdata.frame -> qsu(wlddev, pid = ~ iso3c)
qsu(pwlddev$PCGDP)                          # Default summary for Panel Series  (class pseries)
qsu(G(pwlddev$PCGDP))                       # Summarizing GDP growth, see also ?G

# Grouped Panel Data Summaries -------------
qsu(wlddev, ~ region, ~ iso3c, cols = 9:12) # Panel-Statistics by region
psr <- qsu(pwlddev, ~ region, cols = 9:12)  # Same on plm pdata.frame
psr                                         # -> Gives a 4D array
print.qsu(psr[,"N/T",,])                    # Checking out the number of observations:
# In North america we only have 3 countries, for the GINI we only have 3.91 observations on average
# for 45 Sub-Saharan-African countries, etc..
print.qsu(psr[,"SD",,])                     # Considering only standard deviations
# -> In all regions variations in inequality (GINI) between countries are greater than variations
# in inequality within countries. The opposite is true for Life-Expectancy in all regions apart
# from Europe, etc..
}
psrl <- qsu(wlddev, ~ region, ~ iso3c,      # Same, but output as nested list
            cols = 9:12, array = FALSE)
psrl                                        # We can use unlist2d to create a tidy data.frame
head(unlist2d(psrl, c("Variable","Trans"),
                    row.names = "Region"))

# Weighted Summaries -----------------------
n <- nrow(wlddev)
weights <- abs(rnorm(n))                    # Generate random weights
qsu(wlddev, w = weights, higher = TRUE)     # Computed weighted mean, SD, skewness and kurtosis
weightsNA <- weights                        # Weights may contain missing values.. inserting 1000
weightsNA[sample.int(n, 1000)] <- NA
qsu(wlddev, w = weightsNA, higher = TRUE)   # But now these values are removed from all variables

# Grouped and panel-summaries can also be weighted in the same manor
}
% View(psrdat)
%                                            # We've gotten this far, let's give it a ggplot2 finish:
% psrdat <- reshape2::melt(psrdat, 1:3,
%           variable.name = "Statistic")      # Looks freakin rediculous, but still a nice demonstation
% library(ggplot2)
% ggplot(psrdat, aes(x = Trans, y = value, fill = Region)) +
%        geom_bar(stat = "identity", position = position_dodge()) +
%        facet_wrap(Statistic ~ Variable, scales = "free", ncol = 4)


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar} % __ONLY ONE__ keyword per line % use one of  RShowDoc("KEYWORDS")
