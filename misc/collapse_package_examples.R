
## Let's start with some statistical programming ---------------------------------------
v <- iris$Sepal.Length
d <- num_vars(iris)    # Saving numeric variables
f <- iris$Species      # Factor

# Simple statistics
fmean(v)               # vector
fmean(qM(d))           # matrix (qM is a faster as.matrix)
fmean(d)               # data.frame

# Preserving data structure
fmean(qM(d), drop = FALSE)     # Still a matrix
fmean(d, drop = FALSE)         # Still a data.frame

# Weighted statistics, supported by most functions...
w <- abs(rnorm(fnrow(iris)))
fmean(d, w = w)

# Grouped statistics...
fmean(d, f)

# Groupwise-weighted statistics...
fmean(d, f, w)

# Simple Transformations...
head(fmode(d, TRA = "replace"))    # Replacing values with the mode
head(fmedian(d, TRA = "-"))        # Subtracting the median
head(fsum(d, TRA = "%"))           # Computing percentages
head(fsd(d, TRA = "/"))            # Dividing by the standard-deviation (scaling), etc...

# Weighted Transformations...
head(fnth(d, 0.75, w = w, TRA = "replace"))  # Replacing by the weighted 3rd quartile

# Grouped Transformations...
head(fvar(d, f, TRA = "replace"))  # Replacing values with the group variance
head(fsd(d, f, TRA = "/"))         # Grouped scaling
head(fmin(d, f, TRA = "-"))        # Setting the minimum value in each species to 0
head(fsum(d, f, TRA = "/"))        # Dividing by the sum (proportions)
head(fmedian(d, f, TRA = "-"))     # Groupwise de-median
head(ffirst(d, f, TRA = "%%"))     # Taking modulus of first group-value, etc. ...

# Grouped and weighted transformations...
head(fsd(d, f, w, "/"), 3)         # weighted scaling
head(fmedian(d, f, w, "-"), 3)     # subtracting the weighted group-median
head(fmode(d, f, w, "replace"), 3) # replace with weighted statistical mode

## Some more advanced transformations...
fbetween(d)                             # Averaging (same as fmean(d, TRA = "replace"), but a bit faster..)
fwithin(d)                              # Centering (same as fmean(d, TRA = "-"))
fwithin(d, f, w)                        # Grouped and weighted Centering (same as fmean(d, f, w, "-"))
fwithin(d, f, w, mean = 5)              # Setting a custom mean
fwithin(d, f, w, theta = 0.76)          # Quasi-centering i.e. d - theta*fbetween(d, f, w)
fwithin(d, f, w, mean = "overall.mean") # Preserving the overall mean of the data after grouped centering
fscale(d)                               # Scaling and centering
fscale(d, mean = 5, sd = 3)             # Custom scaling and centering
fscale(d, mean = FALSE, sd = 3)         # Mean preserving scaling
fscale(d, f, w)                         # Grouped and weighted scaling and centering
fscale(d, f, w, mean = 5, sd = 3)       # Custom grouped and weighted scaling and centering
fscale(d, f, w, mean = FALSE,           # Preserving group means
       sd = "within.sd")                # and setting the group-sd to fsd(fwithin(d, f, w), w = w)
fscale(d, f, w, mean = "overall.mean",  # Full harmonization of group means and variances,
       sd = "within.sd")                # while preserving the level and scale of the data.

get_vars(iris, 1:2)                      # Use get_vars for fast selecting data.frame columns, gv is a shortcut
fhdbetween(gv(iris, 1:2), gv(iris, 3:5)) # Linear prediction with factors and continuous covariates
fhdwithin(gv(iris, 1:2), gv(iris, 3:5))  # Linear partialling out factors and continuous covariates
ss(iris, 1:10, 1:2)                      # Similarly fsubset/ss is super fast for subsetting rows of data.frame.

# Simple Time-Computations..
flag(AirPassengers, -1:3)                # One lead and three lags
fdiff(EuStockMarkets,                    # Suitably lagged first and second differences
      c(1, frequency(EuStockMarkets)), diff = 1:2)
fdiff(EuStockMarkets, rho = 0.87)        # Quasi-differences (x_t - rho*x_t-1)
fdiff(EuStockMarkets, log = TRUE)        # Log-differences
fgrowth(EuStockMarkets)                  # Exact growth rates (percentage change)
fgrowth(EuStockMarkets, logdiff = TRUE)  # Log-difference growth rates (percentage change)

# Note that it is not necessary to use factors for grouping.
fmean(gv(mtcars, -c(2,8:9)), mtcars$cyl)            # Can also use a vector (internally converted to factor using qF())
fmean(gv(mtcars, -c(2,8:9)), gv(mtcars, c(2,8:9)))  # or a list of vector (internally grouped using GRP())
g <- GRP(mtcars, ~ cyl + vs + am)                   # It is also possible to create grouping objects [or GRP(mtcars, c(2,8:9))]
print(g)                                            # These are instructive to learn about the grouping,
plot(g)                                             # and are directly handed down to C++ code by collapse functions.
fmean(gv(mtcars, -c(2,8:9)), g)                     # This can speed up multiple computations over the same groups..
fsd(gv(mtcars, -c(2,8:9)), g)

# Factors can efficiently be created using qF
f1 <- qF(mtcars$cyl)                            # Unlike GRP objects, factors have to be checked for missing values
f2 <- qF(mtcars$cyl, na.exclude = FALSE)        # This can however be avoided through this option
class(f2)                                       # Note the added class
library(microbenchmark)
microbenchmark(fmean(mtcars, f1), fmean(mtcars, f2)) # A minor difference, but worthwhile on larger data

with(mtcars, finteraction(cyl, vs, am))  # This creates efficient interactions of vectors and/or factors
finteraction(gv(mtcars, c(2,8:9)))       # .. or lists of vectors/factors

# For simple row-or columnwise computations on matrices or data.frames with custom functions, use dapply
dapply(mtcars, quantile)             # column quantiles
dapply(mtcars, quantile, MARGIN = 1) # row-quantiles
# dapply preserves the data structure of any matrices / data.frames passed. Some fast matrix row/column functions provided by matrixStats.
# Similarly, BY performs grouped comptations
BY(mtcars, f2, quantile)
BY(mtcars, f2, quantile, expand.wide = TRUE)
# For efficient (grouped) replacing and sweeping out computed statistics, use TRA()
sds <- fsd(mtcars)
TRA(mtcars, sds, "/")     # Simple scaling (if sd's are not needed, use fsd(mtcars, TRA = "/"))
microbenchmark(TRA(mtcars, sds, "/"), sweep(mtcars, 2, sds, "/")) # A remarkable performance gain..
sds <- fsd(mtcars, f2)
TRA(mtcars, sds, "/", f2) # Groupd scaling (if sd's are not needed, use fsd(mtcars, f2, TRA = "/"))

# All functions above perserve the structure of matrices/data.frame's. If conversions are required, use these efficient functions:
mtcarsM <- qM(mtcars)                # Matrix from data.frame
qDF(mtcarsM)                         # data.frame from matrix columns
mrtl(mtcarsM, TRUE, "data.frame")    # data.frame from matrix rows, etc..
qDT(mtcarsM, "cars")                 # Saving row.names when converting matrix to data.table
qDT(mtcars, "cars")                  # Same use a data.frame


## Now let's get some real data and see how we can use this power for data manipulation ------------------
library(magrittr)
head(wlddev) # World Bank World Development Data: 216 countries, 59 years, 4 series (columns 9-12)

# Starting with some discriptive tools...
namlab(wlddev, class = TRUE)           # Show variable names, labels and classes
fnobs(wlddev)                          # Observation count
pwnobs(wlddev)                         # Pairwise observation count
fnobs(wlddev, wlddev$country)          # Grouped observation count
fndistinct(wlddev)                     # Distinct values
descr(wlddev)                          # Describe data
varying(wlddev, ~ country)             # Show which variables vary within countries
qsu(wlddev, pid = ~ country,           # Panel-summarize columns 9 though 12 of this data
    cols = 9:12, vlabels = TRUE)       # (between and within countries)
qsu(wlddev, ~ region, ~ country,       # Do all of that by region and also compute higher moments
    cols = 9:12, higher = TRUE)        # -> returns a 4D array
qsu(wlddev, ~ region, ~ country, cols = 9:12,
    higher = TRUE, array = FALSE) %>%                       # Return a list of statistics matrices..
  unlist2d(c("Variable","Trans"), row.names = "Region")     # and turn this into a tidy data.frame
pwcor(num_vars(wlddev), P = TRUE)                           # Pairwise correlations with p-value
pwcor(fmean(num_vars(wlddev), wlddev$country), P = TRUE)    # Correlating country means
pwcor(fwithin(num_vars(wlddev), wlddev$country), P = TRUE)  # Within-country correlations
psacf(wlddev, ~country, ~year, cols = 9:12)                 # Panel-data Autocorrelation function
pspacf(wlddev, ~country, ~year, cols = 9:12)                # Partial panel-autocorrelations
psmat(wlddev, ~iso3c, ~year, cols = 9:12) %>% plot          # Convert panel data to 3D array and plot

## collapse offers a few very efficent functions for data manipulation:
# Fast selecting and replacing columns
series <- get_vars(wlddev, 9:12)     # Same as wlddev[9:12] but 2x faster and works with data.table's
series <- fselect(wlddev, PCGDP:ODA) # Same thing: > 100x faster than dplyr::select(wlddev, PCGDP:ODA)
get_vars(wlddev, 9:12) <- series     # Replace columns, 8x faster than wlddev[9:12] <- series and also replaces names
fselect(wlddev, PCGDP:ODA) <- series # Same thing
# Fast subsetting
head(fsubset(wlddev, country == "Ireland", -country, -iso3c))
head(fsubset(wlddev, country == "Ireland" & year > 1990, year, PCGDP:ODA))
ss(wlddev, 1:10, 1:10) # This is an order of magnitude faster than wlddev[1:10, 1:10]
# Fast transforming
head(ftransform(wlddev, ODA_GDP = ODA / PCGDP, ODA_LIFEEX = sqrt(ODA) / LIFEEX))
settransform(wlddev, ODA_GDP = ODA / PCGDP, ODA_LIFEEX = sqrt(ODA) / LIFEEX) # transfor a data.frame by reference
head(ftransform(wlddev, PCGDP = NULL, ODA = NULL, GINI_sum = fsum(GINI)))
ftransform(wlddev, fscale(gv(wlddev, 9:12)))   # Can also transform with lists of columns
settransform(wlddev, fscale(gv(wlddev, 9:12))) # Changing the data by reference
ftransform(wlddev) <- fscale(gv(wlddev, 9:12)) # Same thing
wlddev %<>% ftransform(fscale(gv(., 9:12)))    # Same thing, using magrittr
# Fast reordering
roworder(wlddev, -country, year)
colorder(wlddev, country, year)
# Fast renaming
frename(wlddev, country = Ctry, year = Yr)
setrename(wlddev, country = Ctry, year = Yr)  # By reference
frename(wlddev, tolower, cols = 9:12)
# Fast grouping
fgroup_by(wlddev, Ctry, decade) %>% fgroup_vars  # fgroup_by is a lot faster than dplyr::group_by, but only works with collapse functions
rm(wlddev)

## Now lets start putting things together
wlddev[["weights"]] <- abs(rnorm(fnrow(wlddev))) # Adding some weights

wlddev %>% fsubset(year > 1990, region, income, PCGDP:ODA) %>%
  fgroup_by(region, income) %>% fmean            # Fast aggregation using the mean

# Same thing using dplyr manipulation verbs
library(dplyr)
wlddev %>% filter(year > 1990) %>% select(region, income, PCGDP:ODA) %>%
  group_by(region,income) %>% fmean              # This is already a lot faster than summarize_all(mean)


wlddev %>% fsubset(year > 1990, region, income, PCGDP:weights) %>%
  fgroup_by(region, income) %>% fmean(weights)   # Weighted group means

wlddev %>% fsubset(year > 1990, region, income, PCGDP:weights) %>%
  fgroup_by(region, income) %>% fsd(weights)     # Weighted group standard deviations

wlddev %>% fgroup_by(region, income) %>% fselect(PCGDP:weights) %>%
   fnth(0.75, weights)                           # Weighted group third quartile

wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA) %>% fwithin              # Within transformation
wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA) %>% fmedian(TRA = "-")   # Grouped centering using the median
wlddev %>% fgroup_by(country) %>% fselect(country, year, PCGDP:weights) %>%   # Replacing data points by the weighted first quartile
  ftransform(., fselect(., -country, -year) %>% fnth(0.25, weights, "replace_fill"))
wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA) %>% fscale               # Standardizing
wlddev %>% fgroup_by(country) %>% fselect(PCGDP:weights) %>% fscale(weights)  # Weigted standardizing

wlddev %>% fselect(country, year, PCGDP:ODA) %>%  # Adding 1 lead and 2 lags of each variable
  fgroup_by(country) %>% flag(-1:2, year)
wlddev %>% fselect(country, year, PCGDP:ODA) %>%  # Adding 1 lead and 10-year growth rates
  fgroup_by(country) %>% fgrowth(c(0:1,10), 1, year)

# etc...

# Aggregation with multiple functions
wlddev %>% fsubset(year > 1990, region, income, PCGDP:ODA) %>%
  fgroup_by(region, income) %>% {
    add_vars(fgroup_vars(., "unique"),
             fmedian(., keep.group_vars = FALSE) %>% add_stub("median_"),
             fmean(., keep.group_vars = FALSE) %>% add_stub("mean_"),
             fsd(., keep.group_vars = FALSE) %>% add_stub("sd_"))
  }

# Transformation with multiple functions
wlddev %>% fselect(country, year, PCGDP:ODA) %>%
  fgroup_by(country) %>% {
    add_vars(fdiff(., c(1,10), 1, year) %>% flag(0:2, year),  # Sequence of lagged 1- and 10-year differences
             ftransform(., fselect(., PCGDP:ODA) %>% fwithin %>% add_stub("W.")) %>%
               flag(0:2, year, keep.ids = FALSE))             # Sequence of lagged demeaned variables
  }

# With ftransform, can also easily do one or more grouped mutations on the fly..
settransform(wlddev, median_ODA = fmedian(ODA, list(region, income), TRA = "replace_fill"))

settransform(wlddev, sd_ODA = fsd(ODA, list(region, income), TRA = "replace_fill"),
                     mean_GDP = fmean(PCGDP, country, TRA = "replace_fill"))

wlddev %<>% ftransform(fmedian(list(median_ODA = ODA, median_GDP = PCGDP),
                               list(region, income), TRA = "replace_fill"))
rm(wlddev)

## For multi-type data aggregation, the function collap offers ease and flexibility
# Aggregate this data by country and decade: Numeric columns with mean, categorical with mode
head(collap(wlddev, ~ country + decade, fmean, fmode))

wlddev[["weights"]] <- abs(rnorm(fnrow(wlddev))) # Adding some weights: taking weighted mean and weighted mode
head(collap(wlddev, ~ country + decade, fmean, fmode, w = ~ weights, wFUN = fsum))

# Multi-function aggregation of certain columns
head(collap(wlddev, ~ country + decade,
            list(fmean, fmedian, fsd),
            list(ffirst, flast), cols = c(3,9:12)))

# Customized Aggregation: Assign columns to functions
head(collap(wlddev, ~ country + decade,
            custom = list(fmean = 9:10, fsd = 9:12, flast = 3, ffirst = 6:8)))

# For grouped tibbles use collapg
wlddev %>% fsubset(year > 1990, country, region, income, PCGDP:ODA) %>%
  fgroup_by(country) %>% collapg(fmean, ffirst) %>%
  ftransform(AMGDP = PCGDP > fmedian(PCGDP, list(region, income), TRA = "replace_fill"),
             AMODA = ODA > fmedian(ODA, income, TRA = "replace_fill"))

## Additional flexibility for data transformation tasks is offerend by tidy transformation operators:
head(W(wlddev, ~ country, cols = 9:12, mean = "overall.mean"))   # Within-transformation (centering on overall mean)
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year)))       # Partialling out country and year fixed effects
head(HDW(wlddev, PCGDP + LIFEEX ~ qF(country) + qF(year) + ODA)) # Same, adding ODA as continuous regressor
head(STD(wlddev, ~ country, cols = 9:12))                        # Standardizing (scaling and centering) by country
head(L(wlddev, -1:3, ~ country, ~year, cols = 9:12))             # Computing 1 lead and 3 lags of the 4 series
head(D(wlddev, c(1,10), 1, ~ country, ~year, cols = 9:12))       # Computing the 1- and 10-year first differences
head(D(wlddev, c(1,10), 1:2, ~ country, ~year, cols = 9:12))     # ..first and second differences
head(G(wlddev, c(1,10), 1, ~ country, ~year, cols = 9:12))       # Computing the 1- and 10-year growth rates
# Adding growth rate variables to dataset
add_vars(wlddev) <- G(wlddev, c(1, 10), 1, ~ country, ~year, cols = 9:12, keep.ids = FALSE)
get_vars(wlddev, "G1.", regex = TRUE) <- NULL # Deleting again

# These operators can conveniently be used in regression formulas:
lm(LIFEEX ~ log(PCGDP) + OECD + B(log(PCGDP), country),          # Using a Mundlak (1978) procedure to estimate the effect
   wlddev %>% fselect(country, OECD, PCGDP, LIFEEX) %>% na_omit) # of OECD membership on LIFEEX, controlling for PCGDP

# Adding 10-year lagged life-expectancy to allow for some convergence effects (dynamic panel model)
lm(LIFEEX ~ L(LIFEEX, 10, country) + log(PCGDP) + OECD + B(log(PCGDP), country),
   wlddev %>% fselect(country, OECD, PCGDP, LIFEEX) %>% na_omit)

# Tranformation functions and operators also support plm panel data classes:
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

# Remove all objects used in this example section
rm(v, d, w, f, f1, f2, g, mtcarsM, pwlddev, sds, series, wlddev)
# wlddev %>% fselect(country, year, PCGDP:ODA) %>%  # Transformation with multiple functions
#   fgroup_by(country) %>% {
#     add_vars(W(.) %>% flag(0:2),
#              B(.) %>% flag(0:2, keep.ids = FALSE))
#   }
