# collapse R Package: Advanced and Fast Data Transformation

*collapse* is a C/C++ based package for data manipulation in R. It's aims are

* to facilitate complex data transformation and exploration tasks in R
* to help make R code fast, flexible, parsimonious and programmer friendly 

It is compatible with *dplyr*, *data.table* and the *plm* approach to panel-data.

**Key Features:**

*  *Advanced data programming*: A full set of fast statistical functions supporting grouped and weighted computations on vectors, matrices and data.frames. Quick (ordered) and reusable grouping, quick data-converions, and quick select and replace variables. 
*  *Advanced aggregation*: Fast and easy multi-data-type, multi-function, weighted, parallelized and fully customized data aggregation. 
*  *Advanced transformations*: Efficient (grouped, weighted) sweeping out of statistics, scaling, centering, higher-dimensional centering, complex linear prediction and partialling-out. 
*  *Advanced time-computations*: Efficient (sequences of) lags/leads, and (iterated) differences and growth rates on (unordered) time-series and panel data. Multivariate auto, partial and cross-correlation functions for panel data. Panel data to (ts-)array conversions. 
*  *List Processing*: (Recursive) list-identification, extraction / subsetting, apply, and row-binding / unlisting in 2D. 
* *Advanced data exploration*: Fast (grouped, weighted, panel-decomposed) summary statistics for cross-sectional and complex multilevel/panel data. 

*collapse* is built using *Rcpp*, uses C functions from *data.table* and *stats*, and imports from *lfe*.

## Installation

The package can be installed in R using the following code:

remotes::install_github("SebKrantz/collapse")

## Contributing

If you want to contribute, please fork and create a pull request for merging with the `development` branch.

