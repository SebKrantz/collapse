# collapse <img src='misc/figures/collapse_logo_small.png' width="150px" align="right" />


<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/collapse)](https://cran.r-project.org/package=collapse)
[![cran checks](https://cranchecks.info/badges/worst/collapse)](https://cran.r-project.org/web/checks/check_results_collapse.html)
[![Travis build status](https://travis-ci.com/SebKrantz/collapse.svg?branch=master)](https://travis-ci.com/SebKrantz/collapse)
[![Codecov test coverage](https://codecov.io/gh/SebKrantz/collapse/branch/master/graph/badge.svg)](https://codecov.io/gh/SebKrantz/collapse?branch=master)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![status](https://tinyverse.netlify.com/badge/collapse)](https://CRAN.R-project.org/package=collapse)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
![downloads per month](http://cranlogs.r-pkg.org/badges/collapse?color=blue)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/collapse?color=blue)
<!-- badges: end -->

*collapse* is a C/C++ based package for data transformation and statistical computing in R. It's aims are:

* To facilitate complex data transformation, exploration and computing tasks in R.
* To help make R code fast, flexible, parsimonious and programmer friendly.

It is made compatible with *dplyr*, *data.table* and the *plm* approach to panel data.

**Key Features:**

*  **Advanced statistical programming**: A full set of fast statistical functions 
        supporting grouped and weighted computations on vectors, matrices and 
        data frames. Fast and programmable grouping, ordering, unique values / rows, 
        factor generation and interactions. Fast and flexible functions for data 
        manipulation and data object conversions.

*  **Advanced aggregation**: Fast and easy multi-data-type, multi-function, 
        weighted, parallelized and fully customized data aggregation.

*  **Advanced transformations**: Fast (grouped) replacing and sweeping out of 
        statistics, and (grouped, weighted) scaling / standardizing, between 
        (averaging) and (quasi-)within (centering / demeaning) transformations, 
        higher-dimensional centering (i.e. multiple fixed effects transformations), 
        linear prediction and partialling-out. 

*  **Advanced time-computations**: Fast (sequences of) lags / leads, and 
        (lagged / leaded, iterated, quasi-, log-) differences and growth 
        rates on (unordered) time series and panel data. Multivariate auto-, 
        partial- and cross-correlation functions for panel data. 
        Panel data to (ts-)array conversions. 

*  **List processing**: (Recursive) list search / identification, extraction / 
        subsetting, data-apply, and generalized row-binding / unlisting in 2D.

* **Advanced data exploration**: Fast (grouped, weighted, panel-decomposed) 
        summary statistics for complex multilevel / panel data.

*collapse* is mainly coded in C++ and built with *Rcpp*, but also uses C functions from *data.table*, *lfe* and *stats*. Effort has been expended to minimize the execution speed of R code employed. 

## Installation

``` r
# From the R console call
install.packages("collapse")
```
<!--
# install the development version
devtools::install_github("SebKrantz/collapse")
-->

## Documentation
*collapse* installs with a built-in structured [documentation](<https://sebkrantz.github.io/collapse/reference/index.html>), implemented via a set of separate help pages. Calling `help('collapse-documentation')` from the R console brings up the the top-level documentation page, which provides an overview of the entire functionality of the package and links to all other documentation pages. 

In addition, *collapse* provides 3 vignettes:

* [Introduction to *collapse*](<https://sebkrantz.github.io/collapse/articles/collapse_intro.html>): Introduces all main features of the package in a structured way.

* [*collapse* and *dplyr*](<https://sebkrantz.github.io/collapse/articles/collapse_and_dplyr.html>): Demonstrates the integration of *collapse* with *dplyr* / *tidyverse* workflows and associated performance improvements.

* [*collapse* and *plm*](<https://sebkrantz.github.io/collapse/articles/collapse_and_plm.html>): Demonstrates the integration of *collapse* with the *plm* package and provides examples of fast and easy programming with panel data. 

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
wt <- abs(rnorm(fnrow(iris)))
fmedian(d, w = wt)                # Simple weighted statistics
fnth(d, 0.75, g)                  # Grouped statistics (grouped third quartile)
fmedian(d, g, wt)                 # Groupwise-weighted statistics
fsd(v, g, wt)                     # Similarly for vectors
fmode(qM(d), g, wt, ties = "max") # Or matrices (grouped and weighted maximum mode) ...

# A fast set of data manipulation functions allows complex piped programming at high speeds
library(magrittr)                            # Pipe operators
iris %>% fgroup_by(Species) %>% fNdistinct   # Grouped distinct value counts
iris %>% fgroup_by(Species) %>% fmedian(wt)  # Weighted group medians 
iris %>% add_vars(wt) %>%                    # Adding weight vector to dataset
  fsubset(Sepal.Length < fmean(Sepal.Length), Species, Sepal.Width:wt) %>% # Fast selecting and subsetting
  fgroup_by(Species) %>%                     # Grouping (efficiently creates a grouped tibble)
  fvar(wt) %>%                               # Frequency-weighted group-variance, default (keep.w = TRUE)  
  roworder(sum.wt)                           # also saves group weights in a column called 'sum.wt'

# Can also use dplyr (but dplyr manipulation verbs are a lot slower)
library(dplyr)
iris %>% add_vars(wt) %>% 
  filter(Sepal.Length < fmean(Sepal.Length)) %>% 
  select(Species, Sepal.Width:wt) %>% 
  group_by(Species) %>% 
  fvar(wt) %>% arrange(sum.wt)

## Advanced Aggregation -----------------------------------------------------------------------------------------

collap(iris, Sepal.Length + Sepal.Width ~ Species, fmean)  # Simple aggregation using the mean..
collap(iris, ~ Species, list(fmean, fmedian, fmode))       # Multiple functions applied to each column
add_vars(iris) <- wt                                       # Adding weights, return in long format..
collap(iris, ~ Species, list(fmean, fmedian, fmode), w = ~ wt, return = "long")

# Generate some additional logical data
settransform(iris, AWMSL = Sepal.Length > fmedian(Sepal.Length, w = wt), 
                   AWMSW = Sepal.Width > fmedian(Sepal.Width, w = wt))

# Multi-type data aggregation: catFUN applies to all categorical columns (here AMWSW)
collap(iris, ~ Species + AWMSL, list(fmean, fmedian, fmode), 
       catFUN = fmode, w = ~ wt, return = "long")

# Custom aggregation gives the greatest possible flexibility: directly mapping functions to columns
collap(iris, ~ Species + AWMSL, 
       custom = list(fmean = 2:3, fsd = 3:4, fmode = "AWMSL"), w = ~ wt, 
       wFUN = list(fsum, fmin, fmax), # Here also aggregating the weight vector with 3 different functions
       keep.col.order = FALSE)        # Column order not maintained -> grouping and weight variables first

# Can also use grouped tibble: weighted median for numeric, weighted mode for categorical columns
iris %>% fgroup_by(Species, AWMSL) %>% collapg(fmedian, fmode, w = wt)

## Advanced Transformations -------------------------------------------------------------------------------------

# All Fast Statistical Functions have a TRA argument, supporting 10 different replacing and sweeping operations
fmode(d, TRA = "replace")     # Replacing values with the mode
fsd(v, TRA = "/")             # dividing by the overall standard deviation (scaling)
fsum(d, TRA = "%")            # Computing percentages
fsd(d, g, TRA = "/")          # Grouped scaling
fmin(d, g, TRA = "-")         # Setting the minimum value in each species to 0
ffirst(d, g, TRA = "%%")      # Taking modulus of first value in each species
fmedian(d, g, wt, "-")        # Groupwise centering by the weighted median
fnth(d, 0.95, g, wt, "%")     # Expressing data in percentages of the weighted species-wise 95th percentile
fmode(d, g, wt, "replace",    # Replacing data by the species-wise weighted minimum-mode
      ties = "min")

# TRA() can also be called directly to replace or sweep with a matching set of computed statistics
TRA(v, sd(v), "/")                       # Same as fsd(v, TRA = "/")
TRA(d, fmedian(d, g, wt), "-", g)        # Same as fmedian(d, g, wt, "-")
TRA(d, BY(d, g, quantile, 0.95), "%", g) # Same as fnth(d, 0.95, g, TRA = "%") (apart from quantile algorithm)

# For common uses, there are some faster and more advanced functions
fbetween(d, g)                           # Grouped averaging [same as fmean(d, g, TRA = "replace") but faster]
fwithin(d, g)                            # Grouped centering [same as fmean(d, g, TRA = "-") but faster]
fwithin(d, g, wt)                        # Grouped and weighted centering [same as fmean(d, g, wt, "-")]
fwithin(d, g, wt, theta = 0.76)          # Quasi-centering i.e. d - theta*fbetween(d, g, wt)
fwithin(d, g, wt, mean = "overall.mean") # Preserving the overall weighted mean of the data

fscale(d)                                # Scaling and centering (default mean = 0, sd = 1)
fscale(d, mean = 5, sd = 3)              # Custom scaling and centering
fscale(d, mean = FALSE, sd = 3)          # Mean preserving scaling
fscale(d, g, wt)                         # Grouped and weighted scaling and centering
fscale(d, g, wt, mean = "overall.mean",  # Setting group means to overall weighted mean,
       sd = "within.sd")                 # and group sd's to fsd(fwithin(d, g, wt), w = wt)

get_vars(iris, 1:2)                      # Use get_vars for fast selecting data.frame columns, gv is shortcut
fHDbetween(gv(iris, 1:2), gv(iris, 3:5)) # Linear prediction with factors and continuous covariates
fHDwithin(gv(iris, 1:2), gv(iris, 3:5))  # Linear partialling out factors and continuous covariates

# This again opens up new possibilities for data manipulation...
iris %>%  
  ftransform(ASWMSL = Sepal.Length > fmedian(Sepal.Length, Species, wt, "replace")) %>%
  fgroup_by(ASWMSL) %>% collapg(w = wt, keep.col.order = FALSE)

iris %>% fgroup_by(Species) %>% num_vars %>% fwithin(wt)  # Weighted demeaning


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

HDW(pdata)              # This projects out id and time fixed effects.. (HDW is operator for fHDwithin)
W(pdata, effect = "Id") # Only Id effects.. (W is operator for fwithin)

## List Processing ----------------------------------------------------------------------------------------------

# Some nested list of heterogenous data objects..
l <- list(a = qM(mtcars[1:8]),                                   # Matrix
          b = list(c = mtcars[4:11],                             # data.frame
                   d = list(e = mtcars[2:10], 
                            f = fsd(mtcars))))                   # Vector

ldepth(l)                       # List has 4 levels of nesting (considering that mtcars is a data.frame)
is.unlistable(l)                # Can be unlisted
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
fNobs(irisNA)                           # Observation count
pwNobs(irisNA)                          # Pairwise observation count
fNobs(irisNA, g)                        # Grouped observation count
fNdistinct(irisNA)                      # Same with distinct values... (default na.rm = TRUE skips NA's)
fNdistinct(irisNA, g)  

descr(iris)                                   # Detailed statistical description of data

varying(iris, ~ Species)                      # Show which variables vary within Species
varying(pdata)                                # Which are time-varying ? 
qsu(iris, w = ~ wt)                           # Fast (one-pass) summary (with weights)
qsu(iris, ~ Species, w = ~ wt, higher = TRUE) # Grouped summary + higher moments
qsu(pdata, higher = TRUE)                     # Panel-data summary (between and within entities)
pwcor(num_vars(irisNA), N = TRUE, P = TRUE)   # Pairwise correlations with p-value and observations
pwcor(W(pdata, keep.ids = FALSE), P = TRUE)   # Within-correlations

```

</details>
<p> </p>

Evaluated and more extensive sets of examples are provided on the [package page](<https://sebkrantz.github.io/collapse/reference/collapse-package.html>) (also accessible from R by calling `example('collapse-package')`), and further in the [vignettes](<https://sebkrantz.github.io/collapse/articles/index.html>) and  [documentation](<https://sebkrantz.github.io/collapse/reference/index.html>).

## Additional Notes
### Regarding Performance 
Some simple benchmarks are provided in the [vignettes](<https://sebkrantz.github.io/collapse/articles/index.html>). In general:

* For simple aggregations of large data (~ 10 mio. obs) the performance is comparable to *data.table* (using functions that *data.table* also GeForce optimizes), e.g. see [here](<https://sebkrantz.github.io/collapse/reference/fast-statistical-functions.html#benchmark>).

* For more complex categorical or weighed aggregations, and for nearly all transformations like grouped replacing and sweeping out statistics, scaling, centering, panel-lags or differences etc. *collapse* is ~10x faster than *data.table*. 

* Due to its highly optimized R code, *collapse* is very efficient for programming. On smaller data a *collapse* implementation will execute within microseconds, whereas packages like *dplyr* or *data.table* will typically evaluate in the millisecond domain.

* This performance extends to grouped and weighted computations on vectors and matrices (no internal conversions, vector and matrix methods are also written in C++). With matrices *collapse* performs similar to fast packages like *Rfast* or *matrixStats*.

### Regarding the Integration with *dplyr*, *plm* and *data.table* 

* ***collapse*** **and** ***dplyr***: The [Fast Statistical Functions](<https://sebkrantz.github.io/collapse/reference/fast-statistical-functions.html>) and [transformation functions and operators](<https://sebkrantz.github.io/collapse/reference/data-transformations.html>) provided by *collapse* have a *grouped_df* method, allowing them to be seamlessly integrated into *dplyr* / *tidyverse* workflows. Doing so facilitates advanced operations in *dplyr* and provides remarkable performance improvements. In addition, *collapse* provides some faster replacements for common base R / *dplyr* verbs (`fselect`/`get_vars`, `fgroup_by`, `fsubset`, `ftransform`/`TRA`, `roworder`, `colorder`, `frename`, `funique`, `na_omit`, etc.). See also [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_dplyr.html>). 

<!-- 
, providing further performance improvements for programming with piped expressions and non-standard evaluation
(bringing *dplyr* close to *data.table* on large data aggregations, and making it faster than *data.table* for advanced transformations) -->

* ***collapse*** **and** ***plm***: The fast [transformation functions and operators](<https://sebkrantz.github.io/collapse/reference/data-transformations.html>) provided by *collapse* also have *pseries* (panel-series) and *pdata.frame* (panel-data.frame) methods. This integrates them seamlessly into *plm* workflows and facilitates the manipulation of panel data. For typical panel data operations like between- and within-transformations or panel lags / leads / differences, *collapse* functions are 20-100x faster than *plm* equivalents, and provide greater versatility. See also [this vignette](<https://sebkrantz.github.io/collapse/articles/collapse_and_plm.html>).

<!-- (e.g. for applying transformations to multiple variables in a *pdata.frame*) -->

* ***collapse*** **and** ***data.table***: All collapse functions can be applied to *data.table*'s and they will also return a *data.table* again. The C/C++ programming of *collapse* was inspired by *data.table* and directly relies on some *data.table* C source code (e.g. for grouping and row-binding). The function `qDT` efficiently converts various R objects to *data.table*, and several functions (`mrtl`, `mctl`, `unlist2d`, ...) have an option to return a *data.table*. 

<!--

fNdistinct(wlddev)
fNdistinct(wlddev, wlddev$iso3c)

wlddev %>% fgroup_by(iso3c) %>% fNdistinct

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



