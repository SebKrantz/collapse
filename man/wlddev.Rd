\name{wlddev}
\alias{wlddev}
\docType{data}
\title{
World Development Dataset
}
\description{
This dataset contains 5 indicators from the World Bank's World Development Indicators (WDI) database: (1) GDP per capita, (2) Life expectancy at birth, (3) GINI index, (4) Net ODA and official aid received and (5) Population. The panel data is balanced and covers 216 present and historic countries from 1960-2020 (World Bank aggregates and regional entities are excluded).

Apart from the indicators the data contains a number of identifiers (character country name, factor ISO3 country code, World Bank region and income level, numeric year and decade) and 2 generated variables: A logical variable indicating whether the country is an OECD member, and a fictitious variable stating the date the data was recorded. These variables were added so that all common data-types are represented in this dataset, making it an ideal test-dataset for certain \emph{collapse} functions.
}
\usage{data("wlddev")}
\format{
  A data frame with 13176 observations on the following 13 variables. All variables are labeled e.g. have a 'label' attribute.
  \describe{
    \item{\code{country}}{\emph{chr} Country Name}
    \item{\code{iso3c}}{\emph{fct} Country Code}
    \item{\code{date}}{\emph{date} Date Recorded (Fictitious)}
    \item{\code{year}}{\emph{int} Year}
    \item{\code{decade}}{\emph{int} Decade}
    \item{\code{region}}{\emph{fct} World Bank Region}
    \item{\code{income}}{\emph{fct} World Bank Income Level}
    \item{\code{OECD}}{\emph{log} Is OECD Member Country?}
    \item{\code{PCGDP}}{\emph{num} GDP per capita (constant 2010 US$)}
    \item{\code{LIFEEX}}{\emph{num} Life expectancy at birth, total (years)}
    \item{\code{GINI}}{\emph{num} GINI index (World Bank estimate)}
    \item{\code{ODA}}{\emph{num} Net official development assistance and official aid received (constant 2018 US$)}
    \item{\code{POP}}{\emph{num} Population, total}
  }
}
% \details{
%%  ~~ If necessary, more details than the __description__ above ~~
% }
\source{
\url{https://data.worldbank.org/}, accessed via the \code{WDI} package. The codes for the series are \code{c("NY.GDP.PCAP.KD", "SP.DYN.LE00.IN", "SI.POV.GINI", "DT.ODA.ALLD.KD", "SP.POP.TOTL")}.
}
% \references{
%%  ~~ possibly secondary sources and usages ~~
% }
\seealso{
\code{\link{GGDC10S}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
data(wlddev)

# Panel-summarizing the 5 series
qsu(wlddev, pid = ~iso3c, cols = 9:13, vlabels = TRUE)

# By Region
qsu(wlddev, by = ~region, cols = 9:13, vlabels = TRUE)

# Panel-summary by region
qsu(wlddev, by = ~region, pid = ~iso3c, cols = 9:13, vlabels = TRUE)

# Pairwise correlations: Ovarall
print(pwcor(get_vars(wlddev, 9:13), N = TRUE, P = TRUE), show = "lower.tri")

# Pairwise correlations: Between Countries
print(pwcor(fmean(get_vars(wlddev, 9:13), wlddev$iso3c), N = TRUE, P = TRUE), show = "lower.tri")

# Pairwise correlations: Within Countries
print(pwcor(fwithin(get_vars(wlddev, 9:13), wlddev$iso3c), N = TRUE, P = TRUE), show = "lower.tri")

}


\keyword{datasets}
