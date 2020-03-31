# *collapse*: Advanced and Fast Data Transformation in R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/collapse)](https://cran.r-project.org/package=collapse)
![](http://cranlogs.r-pkg.org/badges/collapse?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/collapse?color=brightgreen)
<!-- badges: end -->

*collapse* is a C/C++ based package for data manipulation in R. It's aims are

* to facilitate complex data transformation and exploration tasks in R
* to help make R code fast, flexible, parsimonious and programmer friendly 

It is compatible with *dplyr*, *data.table* and the *plm* approach to panel-data.

**Key Features:**

*  *Advanced data programming*: A full set of fast statistical functions supporting grouped and/or weighted computations on vectors, matrices and data.frames. Fast (ordered) and reusable grouping, quick data conversions, and quick select, replace or add data.frame columns. 

*  *Advanced aggregation*: Fast and easy multi-data-type, multi-function, weighted, parallelized and fully customized data aggregation. 

*  *Advanced transformations*: Fast (grouped, weighted) replacing and sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out. 

*  *Advanced time-computations*: Fast (sequences of) lags / leads, and (lagged / leaded, iterated) differences and growth rates on (unordered) time-series and panel data. Multivariate auto, partial and cross-correlation functions for panel data. Panel data to (ts-)array conversions. 

*  *List Processing*: (Recursive) list search / identification, extraction / subsetting, apply, and row-binding / unlisting in 2D. 

* *Advanced data exploration*: Fast (grouped, weighted, panel-decomposed) summary statistics for cross-sectional and complex multilevel / panel data. 

*collapse* is mainly coded in C++ and built with *Rcpp*, but also uses C functions from *data.table*, *lfe* and *stats*. Effort has been expended to minimize the 
execution speed of R code employed. 

## Installation

The package can be installed in R using the following code:

remotes::install_github("SebKrantz/collapse")

It is also available on CRAN. 

## Contributing

If you want to contribute, please fork and create a pull request for merging with the `development` branch.

## Package Documentation
*collapse* installs with a built-in hierarchically structured documentation, implemented via a set of separate help pages. The top-level documentation page provides a quick overview of the entire functionality of the package and links to all other documentation pages. It can be accessed from the R console by calling `help('collapse-documentation')`. 

In addition, *collapse* provides 3 vignettes:

* 'Introduction to *collapse*': Introduces all main features of the package in a structured way.

* '*collapse* and *dplyr*': Demonstrates the integration of *collapse* with *dplyr* and the *tidyverse*.

* '*collapse* and *plm*': Demonstrates the integration of *collapse* with the *plm* package and provides examples of fast and easy programming with panel data. 

### Notes on Performance 
Simple benchmarks are provided in the vignettes. In general:

* For simple aggregations of large data (<= 10 mio. obs) the performance is identical to *data.table* (when using functions that *data.table* internally optimizes. The C/C++ programming principles applied and the grouping mechanism of *collapse* is the same as *data.table*). On very large data (100 mio. obs +), *data.table*'s thread parallelization will let it run faster on a multicore machine. 

* For more complex categorical or weighed aggregations, or for data transformations like grouped scaling, centering or panel-differences, *collapse* is ~10x faster than *data.table* in nearly all applications. 

* Due to its minimized R overhead and a complete avoidance of non-standard evaluation, *collapse* is very efficient and easy to use for advanced programming purposes. On smaller data a *collapse* implementation will execute within the microsecond domain, whereas packages like *dplyr* or *data.table* will typically evaluate in the millisecond domain (~10x slower).

* This performance extends to grouped and weighted computations on vectors and matrices (no internal conversions, vector and matrix methods are also written in C++). *collapse* is not limited to programming with data.frames and it is class-secure and attribute-preserving (thus it can be applied to data.table's, tibbles, grouped tibbles etc. and also to special atomic objects like time-series and time-series matrices etc.).

### Notes on the Integration with *dplyr*, *plm* and *data.table* 

* *collapse* and *dplyr*: The *Fast Statistical Functions* and transformation functions and operators provided by *collapse* all have a *grouped_df* method, allowing them to be seamlessly integrated into *dplyr* / *tidyverse* workflows. Doing so facilitates advanced operations in *dplyr* and provides stunning performance improvements (bringing *dplyr* close to *data.table* on large data aggregations, and making it faster than *data.table* for advanced transformations). This integration is discussed and demonstrated in a separate vignette. 

* *collapse* and *plm*: Fast transformation functions and transformation operators provided by collapse also have *pseries* (panel-series) and *pdata.frame* (panel-data.frame) methods. This integrates them seamlessly into *plm* workflows and facilitates the manipulation of panel data. For typical panel-data operations like between- and within-transformations or panel lags / leads / differences, *collapse* functions are 20-100x faster than *plm* equivalents, and provide greater versatility (i.e. for applying transformations to multiple variables in a *pdata.frame*). This integration is also discussed and demonstrated in a separate vignette. 

* *collapse* and *data.table*: All collapse functions can be applied to *data.table*'s and they will also return a *data.table* again. The C/C++ programming of *collapse* was inspired by *data.table* and directly relies on some *data.table* source code (i.e. for grouping and row-binding). 


