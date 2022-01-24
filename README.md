# collapse <img src='misc/figures/collapse_logo_small.png' width="150px" align="right" />


<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/collapse)](https://cran.r-project.org/package=collapse) 
[![cran checks](https://cranchecks.info/badges/worst/collapse)](https://cran.r-project.org/web/checks/check_results_collapse.html)
[![R build status](https://github.com/SebKrantz/collapse/workflows/R-CMD-check/badge.svg)](https://github.com/SebKrantz/collapse/actions)
[![Codecov test coverage](https://codecov.io/gh/SebKrantz/collapse/branch/master/graph/badge.svg)](https://codecov.io/gh/SebKrantz/collapse?branch=master)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-2.10-6666ff.svg)](https://cran.r-project.org/)
[![status](https://tinyverse.netlify.com/badge/collapse)](https://CRAN.R-project.org/package=collapse)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
![downloads per month](http://cranlogs.r-pkg.org/badges/collapse?color=blue)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/collapse?color=blue)
<!-- [![Travis build status](https://travis-ci.com/SebKrantz/collapse.svg?branch=master)](https://travis-ci.com/SebKrantz/collapse) -->
<!-- badges: end -->

*collapse* is a C/C++ based package for data transformation and statistical computing in R. It's aims are:

* To facilitate complex data transformation, exploration and computing tasks in R.
* To help make R code fast, flexible, parsimonious and programmer friendly. 

It further implements a class-agnostic approach to data manipulation in R, supporting base R, *dplyr* (*tibble*), *data.table*, *sf*, *plm* classes for panel data ('pseries' and 'pdata.frame'), and non-destructively handling other matrix or data frame based classes (including most time series classes such as 'ts', 'xts' / 'zoo', 'timeSeries', 'tsibble', 'tibbletime', etc.). 

<!-- *collapse* thus provides a robust, flexible, class-agnostic and computationally advanced toolkit for data manipulation in R. -->

<!-- Core functions are implicit-type generic and attribute preserving, supporting other matrix or data frame based classes e.g. time series (*ts*, *xts* / *zoo*, *timeSeries* etc.), *sf* data frames etc. -->


**Key Features:**

*  **Advanced statistical programming**: A full set of fast statistical functions 
        supporting grouped and weighted computations on vectors, matrices and 
        data frames. Fast and programmable grouping, ordering, unique values / rows, 
        factor generation and interactions. Fast and flexible functions for data 
        manipulation, data object conversions, and memory efficient R programming.

*  **Advanced aggregation**: Fast and easy multi-data-type, multi-function, 
        weighted, parallelized and fully custom data aggregation.

*  **Advanced transformations**: Fast row / column arithmetic, (grouped) replacing 
        and sweeping out of statistics, (grouped, weighted) scaling / standardizing, 
        between (averaging) and (quasi-)within (demeaning) transformations, 
        higher-dimensional centering (i.e. multiple fixed effects or polynomials), 
        linear prediction, model fitting and testing exclusion restrictions.

*  **Advanced time-computations**: Fast (sequences of) lags / leads, and 
        (lagged / leaded, iterated, quasi-, log-) differences and (compounded) 
        growth rates on (irregular) time series and panel data. 
        Multivariate auto-, partial- and cross-correlation functions for panel data. 
        Panel data to (ts-)array conversions.

*  **List processing**: (Recursive) list search, splitting, 
        extraction / subsetting, data-apply, and generalized recursive 
        row-binding / unlisting in 2D.

* **Advanced data exploration**: Fast (grouped, weighted, panel-decomposed) 
        summary statistics for complex multilevel / panel data.

*collapse* utilizes both C and C++ via *Rcpp*, and also uses C/C++ functions from *data.table*, *kit*, *fixest*, *weights*, *RcppArmadillo*, *RcppEigen* and *stats*. Currently no low-level parallelism is implemented. Effort has been expended to minimize the execution speed of R code employed. 

## Installation

``` r
# Install the current version on CRAN
install.packages("collapse")

# Install previous versions from the CRAN Archive (Requires Rtools)
install.packages("https://cran.r-project.org/src/contrib/Archive/collapse/collapse_1.6.5.tar.gz", repos = NULL, type = "source")
install.packages("https://cran.r-project.org/src/contrib/Archive/collapse/collapse_1.5.3.tar.gz", repos = NULL, type = "source")

# Install a stable development version from Github (Requires Rtools)
remotes::install_github("SebKrantz/collapse")
```
<!--
# install the development version
devtools::install_github("SebKrantz/collapse")
-->

## Documentation
*collapse* installs with a built-in structured [documentation](<https://sebkrantz.github.io/collapse/reference/index.html>), implemented via a set of separate help pages. Calling `help('collapse-documentation')` from the R console brings up the the top-level documentation page, which provides an overview of the entire functionality of the package and links to all other documentation pages. 

In addition, *collapse* provides 5 vignettes (available online):

* [**Introduction to *collapse***](<https://sebkrantz.github.io/collapse/articles/collapse_intro.html>): Introduces all main features of the package in a structured way.

* [***collapse* and *dplyr***](<https://sebkrantz.github.io/collapse/articles/collapse_and_dplyr.html>): Demonstrates the integration of *collapse* with *dplyr* / *tidyverse* workflows and associated performance improvements.

* [***collapse* and *plm***](<https://sebkrantz.github.io/collapse/articles/collapse_and_plm.html>): Demonstrates the integration of *collapse* with the *plm* package and provides examples of fast and easy programming with panel data. 

* [***collapse* and *data.table***](<https://sebkrantz.github.io/collapse/articles/collapse_and_data.table.html>): Shows how *collapse* and *data.table* may be used together in a harmonious way. 

* [***collapse* and *sf***](<https://sebkrantz.github.io/collapse/articles/collapse_and_sf.html>): Shows how collapse can be used to efficiently manipulate *sf* data frames.

### Cheatsheet

<a href="https://raw.githubusercontent.com/SebKrantz/cheatsheets/master/collapse.pdf"><img src="https://raw.githubusercontent.com/SebKrantz/cheatsheets/master/pngs/collapse.png" width="330" height="227"/></a> <!-- 294 -->

## Example Usage
This provides a simple set of examples introducing some important features of *collapse*. It should be easy to follow for readers familiar with R. 


<details>
  <summary><b><a style="cursor: pointer;">Click here to expand </a></b> </summary>
  
``` r
library(collapse)
data("iris")            # iris dataset in base R
v <- iris$Sepal.Length  # Vector
d <- num_vars(iris)     # Saving numeric variables (could also be a matrix, statistical functions are S3 generic)
g <- iris$Species       # Grouping variable (could also be a list of variables)

## Advanced Statistical Programming -----------------------------------------------------------------------------

# Simple (column-wise) statistics...
fmedian(v)                       # Vector
fsd(qM(d))                       # Matrix (qM is a faster as.matrix)
fmode(d)                         # data.frame
fmean(qM(d), drop = FALSE)       # Still a matrix
fmax(d, drop = FALSE)            # Still a data.frame

# Fast grouped and/or weighted statistics
w <- abs(rnorm(fnrow(iris)))
fmedian(d, w = w)                 # Simple weighted statistics
fnth(d, 0.75, g)                  # Grouped statistics (grouped third quartile)
fmedian(d, g, w)                  # Groupwise-weighted statistics
fsd(v, g, w)                      # Similarly for vectors
fmode(qM(d), g, w, ties = "max")  # Or matrices (grouped and weighted maximum mode) ...

# A fast set of data manipulation functions allows complex piped programming at high speeds
library(magrittr)                            # Pipe operators
iris %>% fgroup_by(Species) %>% fndistinct   # Grouped distinct value counts
iris %>% fgroup_by(Species) %>% fmedian(w)   # Weighted group medians 
iris %>% add_vars(w) %>%                     # Adding weight vector to dataset
  fsubset(Sepal.Length < fmean(Sepal.Length), Species, Sepal.Width:w) %>% # Fast selecting and subsetting
  fgroup_by(Species) %>%                     # Grouping (efficiently creates a grouped tibble)
  fvar(w) %>%                                # Frequency-weighted group-variance, default (keep.w = TRUE)  
  roworder(sum.w)                            # also saves group weights in a column called 'sum.w'

# Can also use dplyr (but dplyr manipulation verbs are a lot slower)
library(dplyr)
iris %>% add_vars(w) %>% 
  filter(Sepal.Length < fmean(Sepal.Length)) %>% 
  select(Species, Sepal.Width:w) %>% 
  group_by(Species) %>% 
  fvar(w) %>% arrange(sum.w)

## Advanced Aggregation -----------------------------------------------------------------------------------------

collap(iris, Sepal.Length + Sepal.Width ~ Species, fmean)  # Simple aggregation using the mean..
collap(iris, ~ Species, list(fmean, fmedian, fmode))       # Multiple functions applied to each column
add_vars(iris) <- w                                        # Adding weights, return in long format..
collap(iris, ~ Species, list(fmean, fmedian, fmode), w = ~ w, return = "long")

# Generate some additional logical data
settransform(iris, AWMSL = Sepal.Length > fmedian(Sepal.Length, w = w), 
                   AWMSW = Sepal.Width > fmedian(Sepal.Width, w = w))

# Multi-type data aggregation: catFUN applies to all categorical columns (here AMWSW)
collap(iris, ~ Species + AWMSL, list(fmean, fmedian, fmode), 
       catFUN = fmode, w = ~ w, return = "long")

# Custom aggregation gives the greatest possible flexibility: directly mapping functions to columns
collap(iris, ~ Species + AWMSL, 
       custom = list(fmean = 2:3, fsd = 3:4, fmode = "AWMSL"), w = ~ w, 
       wFUN = list(fsum, fmin, fmax), # Here also aggregating the weight vector with 3 different functions
       keep.col.order = FALSE)        # Column order not maintained -> grouping and weight variables first

# Can also use grouped tibble: weighted median for numeric, weighted mode for categorical columns
iris %>% fgroup_by(Species, AWMSL) %>% collapg(fmedian, fmode, w = w)

## Advanced Transformations -------------------------------------------------------------------------------------

# All Fast Statistical Functions have a TRA argument, supporting 10 different replacing and sweeping operations
fmode(d, TRA = "replace")     # Replacing values with the mode
fsd(v, TRA = "/")             # dividing by the overall standard deviation (scaling)
fsum(d, TRA = "%")            # Computing percentages
fsd(d, g, TRA = "/")          # Grouped scaling
fmin(d, g, TRA = "-")         # Setting the minimum value in each species to 0
ffirst(d, g, TRA = "%%")      # Taking modulus of first value in each species
fmedian(d, g, w, "-")         # Groupwise centering by the weighted median
fnth(d, 0.95, g, w, "%")      # Expressing data in percentages of the weighted species-wise 95th percentile
fmode(d, g, w, "replace",     # Replacing data by the species-wise weighted minimum-mode
      ties = "min")

# TRA() can also be called directly to replace or sweep with a matching set of computed statistics
TRA(v, sd(v), "/")                       # Same as fsd(v, TRA = "/")
TRA(d, fmedian(d, g, w), "-", g)         # Same as fmedian(d, g, w, "-")
TRA(d, BY(d, g, quantile, 0.95), "%", g) # Same as fnth(d, 0.95, g, TRA = "%") (apart from quantile algorithm)

# For common uses, there are some faster and more advanced functions
fbetween(d, g)                           # Grouped averaging [same as fmean(d, g, TRA = "replace") but faster]
fwithin(d, g)                            # Grouped centering [same as fmean(d, g, TRA = "-") but faster]
fwithin(d, g, w)                         # Grouped and weighted centering [same as fmean(d, g, w, "-")]
fwithin(d, g, w, theta = 0.76)           # Quasi-centering i.e. d - theta*fbetween(d, g, w)
fwithin(d, g, w, mean = "overall.mean")  # Preserving the overall weighted mean of the data

fscale(d)                                # Scaling and centering (default mean = 0, sd = 1)
fscale(d, mean = 5, sd = 3)              # Custom scaling and centering
fscale(d, mean = FALSE, sd = 3)          # Mean preserving scaling
fscale(d, g, w)                          # Grouped and weighted scaling and centering
fscale(d, g, w, mean = "overall.mean",   # Setting group means to overall weighted mean,
       sd = "within.sd")                 # and group sd's to fsd(fwithin(d, g, w), w = w)

get_vars(iris, 1:2)                      # Use get_vars for fast selecting data.frame columns, gv is shortcut
fhdbetween(gv(iris, 1:2), gv(iris, 3:5)) # Linear prediction with factors and continuous covariates
fhdwithin(gv(iris, 1:2), gv(iris, 3:5))  # Linear partialling out factors and continuous covariates

# This again opens up new possibilities for data manipulation...
iris %>%  
  ftransform(ASWMSL = Sepal.Length > fmedian(Sepal.Length, Species, w, "replace")) %>%
  fgroup_by(ASWMSL) %>% collapg(w = w, keep.col.order = FALSE)

iris %>% fgroup_by(Species) %>% num_vars %>% fwithin(w)  # Weighted demeaning


## Time Series and Panel Series ---------------------------------------------------------------------------------

flag(AirPassengers, -1:3)                      # A sequence of lags and leads
EuStockMarkets %>%                             # A sequence of first and second seasonal differences
  fdiff(0:1 * frequency(.), 1:2)  
fdiff(EuStockMarkets, rho = 0.95)              # Quasi-difference [x - rho*flag(x)]
fdiff(EuStockMarkets, log = TRUE)              # Log-difference [log(x/flag(x))]
EuStockMarkets %>% fgrowth(c(1, frequency(.))) # Ordinary and seasonal growth rate
EuStockMarkets %>% fgrowth(logdiff = TRUE)     # Log-difference growth rate [log(x/flag(x))*100]

# Creating panel data
pdata <- EuStockMarkets %>% list(`A` = ., `B` = .) %>% 
         unlist2d(idcols = "Id", row.names = "Time")  

L(pdata, -1:3, ~Id, ~Time)                   # Sequence of fully identified panel-lags (L is operator for flag) 
pdata %>% fgroup_by(Id) %>% flag(-1:3, Time) # Same thing..

# collapse supports pseries and pdata.frame's, provided by the plm package
pdata <- plm::pdata.frame(pdata, index = c("Id", "Time"))         
L(pdata, -1:3)          # Same as above, ...
psacf(pdata)            # Multivariate panel-ACF
psmat(pdata) %>% plot   # 3D-array of time series from panel data + plotting

HDW(pdata)              # This projects out id and time fixed effects.. (HDW is operator for fhdwithin)
W(pdata, effect = "Id") # Only Id effects.. (W is operator for fwithin)

## List Processing ----------------------------------------------------------------------------------------------

# Some nested list of heterogenous data objects..
l <- list(a = qM(mtcars[1:8]),                                   # Matrix
          b = list(c = mtcars[4:11],                             # data.frame
                   d = list(e = mtcars[2:10], 
                            f = fsd(mtcars))))                   # Vector

ldepth(l)                       # List has 4 levels of nesting (considering that mtcars is a data.frame)
is_unlistable(l)                # Can be unlisted
has_elem(l, "f")                # Contains an element by the name of "f"
has_elem(l, is.matrix)          # Contains a matrix

get_elem(l, "f")                # Recursive extraction of elements..
get_elem(l, c("c","f"))         
get_elem(l, c("c","f"), keep.tree = TRUE)
unlist2d(l, row.names = TRUE)   # Intelligent recursive row-binding to data.frame   
rapply2d(l, fmean) %>% unlist2d # Taking the mean of all elements and repeating

# Application: extracting and tidying results from (potentially nested) lists of model objects
list(mod1 = lm(mpg ~ carb, mtcars), 
     mod2 = lm(mpg ~ carb + hp, mtcars)) %>%
  lapply(summary) %>% 
  get_elem("coef", regex = TRUE) %>%   # Regular expression search and extraction
  unlist2d(idcols = "Model", row.names = "Predictor")

## Summary Statistics -------------------------------------------------------------------------------------------

irisNA <- na_insert(iris, prop = 0.15)  # Randmonly set 15% missing
fnobs(irisNA)                           # Observation count
pwnobs(irisNA)                          # Pairwise observation count
fnobs(irisNA, g)                        # Grouped observation count
fndistinct(irisNA)                      # Same with distinct values... (default na.rm = TRUE skips NA's)
fndistinct(irisNA, g)  

descr(iris)                                   # Detailed statistical description of data

varying(iris, ~ Species)                      # Show which variables vary within Species
varying(pdata)                                # Which are time-varying ? 
qsu(iris, w = ~ w)                            # Fast (one-pass) summary (with weights)
qsu(iris, ~ Species, w = ~ w, higher = TRUE)  # Grouped summary + higher moments
qsu(pdata, higher = TRUE)                     # Panel-data summary (between and within entities)
pwcor(num_vars(irisNA), N = TRUE, P = TRUE)   # Pairwise correlations with p-value and observations
pwcor(W(pdata, keep.ids = FALSE), P = TRUE)   # Within-correlations

```

</details>
<p> </p>

Evaluated and more extensive sets of examples are provided on the [package page](<https://sebkrantz.github.io/collapse/reference/collapse-package.html>) (also accessible from R by calling `example('collapse-package')`), and further in the [vignettes](<https://sebkrantz.github.io/collapse/articles/index.html>) and  [documentation](<https://sebkrantz.github.io/collapse/reference/index.html>).

## Additional Notes
### Regarding Performance 
Some simple benchmarks against *dplyr*, *data.table* and *plm* are provided in [this](<https://sebkrantz.github.io/Rblog/2020/08/31/welcome-to-collapse/>) blog post and in the [vignettes](<https://sebkrantz.github.io/collapse/articles/index.html>). In general:

<!-- using functions that *data.table* also GeForce optimizes, -->

* For simple aggregations of large data (~ 10 mio. obs) the performance is comparable to *data.table* (e.g. see [here](<https://sebkrantz.github.io/collapse/reference/fast-statistical-functions.html#benchmark>) and [here](<https://sebkrantz.github.io/Rblog/2020/08/31/welcome-to-collapse/>))^[Collapse has quite efficient algorithms but no low-level parallelism. Thus huge aggregations with simple functions like `mean` or `sum` and meaningful parallel processing power are faster on *data.table*, whereas *collapse* can still be faster on 2-core machines.].

* For more complex categorical or weighed aggregations and for data transformations, *collapse* can be ~10x faster than *data.table*. Notable are very fast algorithms for (grouped) statistical mode and distinct value counts, variance, various weighted statistics, scaling, centering, panel-lags, differences and growth rates.

* Due to its highly optimized R code, *collapse* is very efficient for programming. On smaller data a *collapse* implementation will execute within microseconds, whereas packages like *dplyr* or *data.table* will typically evaluate in the millisecond domain (up to ~100x slower).

* This performance extends to grouped and weighted computations on vectors and matrices (*collapse* provides separate vector, matrix and data.frame methods written in C++, the performance in matrix computations is comparable to *Rfast* and *matrixStats*).

### Regarding the Integration with *dplyr*, *plm*, *data.table*, *sf* and Other Classes

* ***collapse*** **and** ***dplyr***: The [Fast Statistical Functions](<https://sebkrantz.github.io/collapse/reference/fast-statistical-functions.html>) and [transformation functions and operators](<https://sebkrantz.github.io/collapse/reference/data-transformations.html>) provided by *collapse* have a *grouped_df* method, allowing them to be seamlessly integrated into *dplyr* / *tidyverse* workflows. Doing so facilitates advanced operations in *dplyr* and provides remarkable performance improvements. In addition, *collapse* provides some faster replacements for common base R / *dplyr* verbs (`fselect`/`get_vars`, `fgroup_by`, `fsubset`, `fmutate`, `fsummarise`, `across`, `roworder`, `colorder`, `frename`, `frelabel`, `funique`, `na_omit`, etc.). `options(collapse_mask = "manip")` can be used to export copies of these functions named `select`, `group_by`, `summarise`, `mutate`, `rename` etc. so that *dplyr* codes can be translated and optimized without much change of syntax. See [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_dplyr.html>) for further examples and benchmarks. 

<!-- 
, providing further performance improvements for programming with piped expressions and non-standard evaluation
(bringing *dplyr* close to *data.table* on large data aggregations, and making it faster than *data.table* for advanced transformations) -->

* ***collapse*** **and** ***plm***: The fast [transformation functions and operators](<https://sebkrantz.github.io/collapse/reference/data-transformations.html>) provided by *collapse* also have *pseries* (panel-series) and *pdata.frame* (panel-data.frame) methods. This integrates them seamlessly into *plm* workflows and facilitates the manipulation of panel data. For typical panel data operations like between- and within-transformations or panel lags / leads / differences, *collapse* functions are 20-100x faster than *plm* equivalents, and provide greater versatility. See also [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_plm.html>).

<!-- (e.g. for applying transformations to multiple variables in a *pdata.frame*) -->

* ***collapse*** **and** ***data.table***: All collapse functions can be applied to *data.table*'s and they will also return a *data.table* again. The C/C++ programming of *collapse* was inspired by *data.table* and directly relies on some *data.table* C source code (e.g. for grouping and row-binding). The function `qDT` efficiently converts various R objects to *data.table*, and several functions (`mrtl`, `mctl`, `unlist2d`, ...) have an option to return a *data.table*. See also [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_data.table.html>).

* ***collapse*** **and** ***sf***: *collapse* now directly supports *sf* data frames through functions like `fselect`, `fsubset`, `num_vars`, `qsu`, `descr`, `varying`, `funique`, `roworder`, `rsplit`, `fcompute` etc., which will take along the geometry column even if it is not explicitly selected (mirroring *dplyr* methods for *sf* data frames). See also [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_sf.html>).

* **Time series and other classes**: Besides explicit support for *dplyr* / *tibble*, *data.table*, *sf* and *plm* panel data classes, *collapse*'s statistical and transformation functions are S3 generic, with 'default', 'matrix' and 'data.frame' methods which dispatch on the implicit data type. Furthermore, these methods intelligently preserve the attributes of the objects passed. Therefore *collapse* can handle many other matrix or data frame based classes, including *ts*, *xts* / *zoo*, *timeSeries*, *tsibble* and *tibbletime*. Compatibility is of course limited if manipulating an object requires further actions besides preservation of the attributes and suitable modification of 'names', 'dim', 'dimnames' and 'row.names'. 

<!--

fndistinct(wlddev)
fndistinct(wlddev, wlddev$iso3c)

wlddev %>% fgroup_by(iso3c) %>% fndistinct

collap(wlddev, ~ country + decade, fmean, fmode)

fscale(num_vars(wlddev), wlddec$iso3c)
fwithin(num_vars(wlddev), wlddec$iso3c)

wlddev %>% fgroup_by(country, decade) %>% fselect(PCGDP:ODA) %>% fwithin(ODA)

L(wlddev, -1:1, ~iso3c, ~year, cols = 9:12)

wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA, year) %>% flag(-1:1, year)

wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA, year) %>% fdiff(-1:1, 1:2, year)
wlddev %>% fgroup_by(country) %>% fselect(PCGDP:ODA, year) %>% fgrowth(1, 1, year)
-->
<!--
## Contributing 
If you want to contribute, please fork and create a pull request for merging with the **development** branch. Presently I am particularly interested in fast algorithms to compute weighted medians and (weighted) quantiles. 


settransform(Species.AWMSL = finteraction(Species, AWMSL))



<!-- *collapse* is not limited to programming with data.frames and it is class-secure and attribute-preserving (thus it can be applied to data.table's, tibbles, grouped tibbles etc. and also to special atomic objects like time-series and time-series matrices etc.). -->



