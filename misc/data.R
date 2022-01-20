PRIO <- haven::read_dta("C:/Users/Sebastian Krantz/Documents/Geneva Graduate Institute/S411006 - Environment and Development Economics/My Course 2017/Paper Material/Data and Analysis/Prio GRID/v2.0/Priogrid2.0.dta")
WDI <- haven::read_dta("C:/Users/Sebastian Krantz/Documents/University College Roosevelt/SSC 400 - Honors Thesis/Data/World BANK/WDI/The World Bank Development Indicators (No duplicate series).dta")
GGDC <- haven::read_dta("C:/Users/Sebastian Krantz/Documents/Geneva Graduate Institute/Dissertation/Data/GGDC_10sd_jan15_2014/10SD_jan15.dta")
source("R/collapse R/GRP.R")

gid <- as.integer(PRIO$gid)[order(rnorm(nrow(PRIO)))]
gwno <- as.integer(PRIO$gwno)[order(rnorm(nrow(PRIO)))]
rm(PRIO)
gc()

library(haven)
library(magrittr)
DHSBR <- read_dta("C:/Users/Sebastian Krantz/Google Drive (keine Synchronisierung)/6. Data and Data Systems/DHS2016/Data/UGBR7BDT - Births Recode/UGBR7BFL.DTA") %>%
  get_vars(fNobs(.) > 0) %>% as_factor

library(sf)
# Level 4: (Census) ---------------------------------------------------------------------------------
UGA_sf <- st_read("C:/Users/Sebastian Krantz/Documents/ODI Fellowship/MoFPED/Census Mapping App/shapefiles/parish other/UG_Parish_Level_Census_Data_GIS_Layer.gpkg") # st_read("C:/Users/Sebastian Krantz/Downloads/gadm36_UGA.gpkg")


NWDI = WDI[sapply(WDI,is.numeric)]
WDIM = as.matrix(NWDI)
g = GRP(WDI$countrycode)
t = GRP(WDI$year)

NGGDC = GGDC[sapply(GGDC,is.numeric)]
GGDCM = as.matrix(NGGDC)
g = GRP(GGDC, ~ Variable + Country)
t = GRP(GGDC$Year)

# Update wlddev:
library(collapse)
codes <- c("NY.GDP.PCAP.KD", "SP.DYN.LE00.IN", "SI.POV.GINI", "DT.ODA.ALLD.KD", "SP.POP.TOTL")
data <- WDI::WDI(country = as.character(unique(wlddev$iso3c)), # API for the world bank development indicators:, see ?WDI::WDI. For more indicators see als install.packages("wbstats")
                 indicator = codes, extra = TRUE)
OECD <- c("AUS", "AUT", "BEL", "CAN", "CHL", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "ISL", "IRL", "ISR", "ITA", "JPN", "KOR", "LVA", "LTU", "LUX", "MEX", "NLD", "NZL", "NOR", "POL", "PRT", "SVK", "SVN", "ESP", "SWE", "CHE", "TUR", "GBR", "USA")
data <- data |>
  ftransform(date = as.Date(paste0(year + 1L, "-01-01")),
             decade = as.integer(floor(year / 10L) * 10L),
             OECD = iso3c %in% OECD,
             iso3c = qF(iso3c),
             region = qF(region),
             income = qF(income)) |>
        get_vars(c("country", "iso3c", "date", "year", "decade", "region", "income", "OECD", codes)) |>
        roworder(country, year)
vlabels(data)[1:8] <- vlabels(wlddev)[1:8]
names(data)[9:13] <- .c(PCGDP, LIFEEX, GINI, ODA, POP)
namlab(data)
str(data)
vclasses(data)[1:12] == vclasses(wlddev)
rm(OECD)

wlddev <- data
usethis::use_data(wlddev, overwrite = TRUE)

# Something more random: regresing agricultural VA on mining and manufacturing VA with country-standardized data
mylm <- lm(STD.AGR ~ STD.MIN + STD.MAN,
           data = STD(subset(GGDC10S, Variable == "VA"), ~Country))
summary(mylm)


# panel-lag test-data:
load("~/R/WDItest.RData")
wlddev <- data
rm(data)
pwlddev <- plm::pdata.frame(wlddev, index = c("country","year"))

dt = subset(data, year >= 1990 & year <= 2000 & iso3c %in% c("FRA","DEU","USA"))[c(2,4,9:10)]
dt = dt[do.call(order,dt[1:2]),]
dtuo = dt[order(rnorm(33)),]
o = do.call(order,dtuo[1:2])

# panel test data:
y = 1:10
yd = list(y,y)
ng = 3
g = c(1,1,1,2,2,2,3,3,3,3)
gs = c(3,3,4)
t = c(1:3,1:3,1:4)
t2 = c(1,2,2,1:3,1:4)
t3 = c(1:3,c(3,4,6),1:4)

# longer:
y = 1:30
ng = 3
g = rep(1:3, each = 10)
gs = rep(10,3)
t = c(1:10,1:10,1:10)
ord = order(rnorm(30))
yuo = y[ord]
guo = g[ord]
tuo = t[ord]
o = order(ord)

# large data !!
y = 1:1e6
ng = 1e5
g = rep(1:ng, each = 10)
gs = rep(10, ng)
t = data.table::rowidv(g)
yuo = y[order(rnorm(1e6))]
guo = g[yuo]
tuo = t[yuo]
o = order(yuo)

# very large data !!
y = 1:1e7
ng = 1e6
g = rep(1:ng, each = 10)
gs = rep(10, ng)
t = data.table::rowidv(g)
yuo = y[order(rnorm(1e7))]
guo = g[yuo]
tuo = t[yuo]
o = order(yuo)


sourceCpp("R/C++/fdiff.cpp", rebuild = TRUE)
r1 = fdiffCpp(y,-3:3,1:3,NA,ng,g,gs)
r2 = fdiffCpp(y,-3:3,1:3,NA,ng,g)
r3 = fdiffCpp(yuo,-3:3,1:3,NA,ng,guo,gs,tuo)[o,]
all.identical(r1,r2,r3)

identical(fdiffCpp(yuo,-3:3,1:3,NA,ng,guo,gs,tuo)[o,],fdiffCpp(y,-3:3,1:3,NA,ng,g))

fdiffCpp(yuo,1:5,1,NA,ng,guo,gs,tuo)[o,] # positive lags works !!
fdiffCpp(yuo,1:5,1:2,NA,ng,guo,gs,tuo)[o,] # positive lags and differences also works !!
fdiffCpp(yuo,-3:3,1,NA,ng,guo,gs,tuo)[o,] # negative lags also works !!
fdiffCpp(yuo,-3:3,1:3,NA,ng,guo,gs,tuo)[o,] # negative lags also works !!

fdiffCpp(y,1,1,NA,ng,g,gs,t)
fdiffCpp(yuo,1,1,NA,ng,guo,gs,tuo)


identical(fdiffCpp(yuo,-3:3,1:3, t = yuo)[o,],fdiffCpp(y,-3:3,1:3, t = y))



# still check looping through groups option.. -> Nah, slower !!

# Manually trying omap !!
ng = 3
g = rep(1:3, each = 10)
gs = rep(10,3)
t = c(1:10,1:10,1:10)
y = t
o = order(rnorm(30))
y = y[o]
g = g[o]
t = t[o]
o = order(o)
omap = list(numeric(10),numeric(10),numeric(10))
for(i in 1:30) omap[[g[i]]][t[i]] = i
outp = numeric(30)
for(i in 1:30) if(t[i]>1) outp[i] = y[i] - y[omap[[g[i]]][t[i]-1]] else outp[i] = NA
outp
outp[o]
for(i in 1:30) if(t[i]>2) outp[i] = outp[i] - outp[omap[[g[i]]][t[i]-1]] else outp[i] = NA
outp2 = outp
for(i in 1:30) if(t[i]>2) outp[i] = outp2[i] - outp2[omap[[g[i]]][t[i]-1]] else outp[i] = NA







# testing lag function:
sourceCpp("R/C++/Experiment/test3.cpp", rebuild = TRUE)
sourceCpp("R/C++/Experiment/test5.cpp", rebuild = TRUE)
sourceCpp("R/C++/Experiment/test6.cpp", rebuild = TRUE)
sourceCpp("R/C++/flag.cpp", rebuild = TRUE)


flagleadCpp2(yuo, t = yuo)[o]

flagleadCpp2(yuo,1,NA,ng,guo,gs,tuo)[o]
flagleadCpp1(yuo,1,NA,ng,guo,gs,tuo)[o]
flagleadCpp(yuo,1,NA,ng,guo,gs,tuo)[o]

# In general: Subtracting 1 (i.e. from g or ord) does not cost much time !!!!

# The 2D map is faster for multiple lags, but the 1D map generates slightly faster !!
microbenchmark(flagleadCpp2(yuo,1:3,NA,ng,guo,gs,tuo),flagleadCpp(yuo,1:3,NA,ng,guo,gs,tuo))

# howver if you access 1D map in the classical way, they are the same speed !!






m = as.matrix(mtcars)
mtcNA = missForest::prodNA(mtcars)
mNA = as.matrix(mtcNA)

x = rnorm(10000000)
x = 100000000 + abs(rnorm(10000000, sd = 0.0001))
xNA = x
xNA[sample.int(10000000, 1000000)] = NA
xNA[1] = NA
xNA[10000000] = NA
gx = GRP(sample.int(1e6,1e7,TRUE))
w = rep(1,10000000)

v = iris$Sepal.Length
vNA = v
vNA[sample.int(150,30)] = NA
vNA[1] = NA
vNA[150] = NA
gv = GRP(iris$Species)
wv = rep(1,150)

datNA = as.data.frame(list(xNA,xNA,xNA,xNA,xNA))
