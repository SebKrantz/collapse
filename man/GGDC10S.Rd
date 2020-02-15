\name{GGDC10S}
\alias{GGDC10S}
\docType{data}
\title{
Groningen Growth and Development Centre 10-Sector Database
}
\description{
The GGDC 10-Sector Database provides a long-run internationally comparable dataset on sectoral productivity performance in Africa, Asia, and Latin America. Variables covered in the data set are annual series of value added (in local currency), and persons employed for 10 broad sectors.
}
\usage{data("GGDC10S")}
\format{
  A data frame with 5027 observations on the following 16 variables.
  \describe{
    \item{\code{Country}}{\emph{char}: Country (43 countries)}
    \item{\code{Regioncode}}{\emph{char}: ISO3 Region code}
    \item{\code{Region}}{\emph{char}: Region (6 World Regions)}
    \item{\code{Variable}}{\emph{char}: Variable (Value Added or Employment)}
    \item{\code{Year}}{\emph{num}: Year (67 Years, 1947-2013)}
    \item{\code{AGR}}{\emph{num}: Agriculture}
    \item{\code{MIN}}{\emph{num}: Mining}
    \item{\code{MAN}}{\emph{num}: Manufacturing}
    \item{\code{PU}}{\emph{num}: Utilities}
    \item{\code{CON}}{\emph{num}: Construction}
    \item{\code{WRT}}{\emph{num}: Trade, restaurants and hotels}
    \item{\code{TRA}}{\emph{num}: Transport, storage and communication}
    \item{\code{FIRE}}{\emph{num}: Finance, insurance, real estate and business services}
    \item{\code{GOV}}{\emph{num}: Government services}
    \item{\code{OTH}}{\emph{num}: Community, social and personal services}
    \item{\code{SUM}}{\emph{num}: Summation of sector GDP}
  }
}
% \details{
%%  ~~ If necessary, more details than the __description__ above ~~
% }
\source{
\url{https://www.rug.nl/ggdc/productivity/10-sector/}
}
\references{
Timmer, M. P., de Vries, G. J., & de Vries, K. (2015). "Patterns of Structural Change in Developing Countries." . In J. Weiss, & M. Tribe (Eds.), \emph{Routledge Handbook of Industry and Development.} (pp. 65-83). Routledge.
}
\seealso{
\code{\link{wlddev}}, \link[=collapse-documentation]{Collapse Overview}
}
\examples{
namlab(GGDC10S, class = TRUE)
qsu(GGDC10S, ~ Variable, ~ Variable + Country, vlabels = TRUE)

library(data.table)
library(ggplot2)

\dontrun{
## World Regions Structural Change Plot

dat <- GGDC10S
get_vars(dat, 6:15) <- replace_outliers(sweep(get_vars(dat, 6:15), 1, dat$SUM, "/"), 0,
                                        single.limit = "min")
names(dat) <- stringr::str_wrap(vlabels(dat), width = 20)
dat$Variable <- Recode(dat$Variable,"VA"="Value Added Share","EMP"="Employment Share")
dat <- collap(dat, ~ Variable + Region + Year, cols = 6:15)
dat <- melt(qDT(dat), 1:3, variable.name = "Sector")

ggplot(aes(x = Year, y = value, fill = Sector), data = dat) +
  geom_area(position = "fill", alpha = 0.9) + labs(x = NULL, y = NULL) +
  theme_linedraw(base_size = 14) + facet_grid(Variable ~ Region, scales = "free_x") +
  scale_fill_manual(values = sub("#00FF66FF", "#00CC66", rainbow(10))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), expand = c(0, 0))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0, 0),
                     labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 315, hjust = 0, margin = ggplot2::margin(t = 0)),
        strip.background = element_rect(colour = "grey30", fill = "grey30"))
}


# A function to plot the structural change of an arbitrary country

plotGGDC <- function(ctry) {
dat <- qDT(GGDC10S)[Country == ctry]
dat <- cbind(get_vars(dat, c("Variable","Year")),
             replace_outliers(sweep(get_vars(dat, 6:15), 1, dat$SUM, "/"), 0, NA, "min"))
dat$Variable <- Recode(dat$Variable,"VA"="Value Added Share","EMP"="Employment Share")
dat <- melt(dat, 1:2, variable.name = "Sector")

ggplot(aes(x = Year, y = value, fill = Sector), data = dat) +
  geom_area(position = "fill", alpha = 0.9) + labs(x = NULL, y = NULL) +
  theme_linedraw(base_size = 14) + facet_wrap( ~ Variable) +
  scale_fill_manual(values = sub("#00FF66FF", "#00CC66", rainbow(10))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0, 0),
                     labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 315, hjust = 0, margin = ggplot2::margin(t = 0)),
        strip.background = element_rect(colour = "grey20", fill = "grey20"),
        strip.text = element_text(face = "bold"))
}
plotGGDC("TZA")

}
\keyword{datasets}
