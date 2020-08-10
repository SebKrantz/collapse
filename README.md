# collapse <img src='misc/figures/collapse_logo_small.png' width="150px" align="right" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/collapse)](https://cran.r-project.org/package=collapse)
[![Travis build status](https://travis-ci.com/SebKrantz/collapse.svg?branch=master)](https://travis-ci.com/SebKrantz/collapse)
![](http://cranlogs.r-pkg.org/badges/collapse?color=blue)
![](http://cranlogs.r-pkg.org/badges/grand-total/collapse?color=blue)
<!-- badges: end -->

*collapse* is a C/C++ based package for data manipulation in R. It's aims are

* to facilitate complex data transformation and exploration tasks in R
* to help make R code fast, flexible, parsimonious and programmer friendly 

It is compatible with *dplyr*, *data.table* and the *plm* approach to panel-data.

**Key Features:**

*  **Advanced data programming**: A full set of fast statistical functions 
        supporting grouped and weighted computations on vectors, matrices and 
        data frames. Fast (ordered) and programmable grouping, factor 
        generation, manipulation of data frames and data object conversions.

*  **Advanced aggregation**: Fast and easy multi-data-type, multi-function, 
        weighted, parallelized and fully customized data aggregation.

*  **Advanced transformations**: Fast (grouped, weighted) replacing and 
        sweeping out of statistics, scaling / standardizing, centering (i.e. 
        between and within transformations), higher-dimensional centering 
        (i.e. multiple fixed effects transformations), linear 
        prediction and partialling-out. 

*  **Advanced time-computations**: Fast (sequences of) lags / leads, and 
        (lagged / leaded, iterated, quasi-, log-) differences and growth 
        rates on (unordered) time-series and panel data. Multivariate auto, 
        partial and cross-correlation functions for panel data. 
        Panel data to (ts-)array conversions. 

*  **List processing**: (Recursive) list search / identification, extraction / 
        subsetting, data-apply, and generalized row-binding / unlisting in 2D.

* **Advanced data exploration**: Fast (grouped, weighted, panel-decomposed) 
        summary statistics for complex multilevel / panel data. 

*collapse* is mainly coded in C++ and built with *Rcpp*, but also uses C functions from *data.table*, *lfe* and *stats*. Effort has been expended to minimize the 
execution speed of R code employed. 

## Installation

```{r}
install.packages("collapse")

# install the development version
devtools::install_github("SebKrantz/collapse")
```

## Package Documentation
*collapse* installs with a built-in hierarchically structured documentation, implemented via a set of separate help pages. The top-level documentation page provides a quick overview of the entire functionality of the package and links to all other documentation pages. It can be accessed from the R console by calling `help('collapse-documentation')`. 

In addition, *collapse* provides 3 vignettes:

* [Introduction to *collapse*](<https://cran.r-project.org/web/packages/collapse/vignettes/collapse_intro.html>): Introduces all main features of the package in a structured way.

* [*collapse* and *dplyr*](<https://cran.r-project.org/web/packages/collapse/vignettes/collapse_and_dplyr.html>): Demonstrates the integration of *collapse* with *dplyr* and the *tidyverse*.

* [*collapse* and *plm*](<https://cran.r-project.org/web/packages/collapse/vignettes/collapse_and_plm.html>): Demonstrates the integration of *collapse* with the *plm* package and provides examples of fast and easy programming with panel data. 

## Demonstration

```{r}
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

```

<!--
## Contributing 
If you want to contribute, please fork and create a pull request for merging with the **development** branch. Presently I am particularly interested in fast algorithms to compute weighted medians and (weighted) quantiles. -->


***

### Notes on Performance 
Some simple benchmarks are provided in the vignettes. In general:

* On simple aggregations of large data (~ 10 mio. obs) the performance is comparable to *data.table* (when using functions that *data.table* GeForce optimizes. The C/C++ programming principles applied and the grouping mechanisms of *collapse* are similar to *data.table*). On very large data (100 mio. obs +), *data.table*'s thread parallelization will let it run faster on a multicore machine. 

* For more complex categorical or weighed aggregations, and for nearly all transformations like grouped scaling, centering, panel-differences etc. *collapse* is ~10x faster than *data.table*. 

* Due to its minimized R overhead and avoidance of non-standard evaluation, *collapse* is very efficient for programming. On smaller data a *collapse* implementation will execute within the microsecond domain, whereas packages like *dplyr* or *data.table* will typically evaluate in the millisecond domain (~10x slower).

* This performance extends to grouped and weighted computations on vectors and matrices (no internal conversions, vector and matrix methods are also written in C++). 

<!-- *collapse* is not limited to programming with data.frames and it is class-secure and attribute-preserving (thus it can be applied to data.table's, tibbles, grouped tibbles etc. and also to special atomic objects like time-series and time-series matrices etc.). -->

### Notes on the Integration with *dplyr*, *plm* and *data.table* 

* ***collapse*** **and** ***dplyr***: The *Fast Statistical Functions* and transformation functions and operators provided by *collapse* all have a *grouped_df* method, allowing them to be seamlessly integrated into *dplyr* / *tidyverse* workflows. Doing so facilitates advanced operations in *dplyr* and provides remarkable performance improvements (bringing *dplyr* close to *data.table* on large data aggregations, and making it faster than *data.table* for advanced transformations). In addition, *collapse* provides some simpler and faster replacements for common *dplyr* verbs (`fselect`, `fgroup_by`, `fsubset` (faster `dplyr::filter`), `ftransform` (faster `dplyr::mutate`) and `TRA` (faster `dplyr::mutate` for grouped replacing and sweeping out statistics)), providing further performance improvements for programming with piped expressions and non-standard evaluation. See also [this vignette](<https://cran.r-project.org/web/packages/collapse/vignettes/collapse_and_dplyr.html>). 

* ***collapse*** **and** ***plm***: Fast transformation functions and transformation operators provided by *collapse* also have *pseries* (panel-series) and *pdata.frame* (panel-data.frame) methods. This integrates them seamlessly into *plm* workflows and facilitates the manipulation of panel data. For typical panel-data operations like between- and within-transformations or panel lags / leads / differences, *collapse* functions are 20-100x faster than *plm* equivalents, and provide greater versatility (i.e. for applying transformations to multiple variables in a *pdata.frame*). See also [this vignette](<https://cran.r-project.org/web/packages/collapse/vignettes/collapse_and_plm.html>).

* ***collapse*** **and** ***data.table***: All collapse functions can be applied to *data.table*'s and they will also return a *data.table* again. The C/C++ programming of *collapse* was inspired by *data.table* and directly relies on some *data.table* source code (i.e. for grouping and row-binding). The function `qDT` also exists to efficiently convert various objects to *data.table*, and various functions (`mrtl`, `mctl`, `unlist2d`, ...) have an option to output a *data.table*. 


