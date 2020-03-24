# collapse 1.1.0
collapse 1.1.0 released end of March 2020: Some small fixes, improvements, and addressing test failures on some operating systems:

* Fixed gcc10 compiler issues, and added some more tests.

* Fixed the issue that supplying an unnamed list to `GRP()`, i.e. `GRP(list(v1, v2))` would give an error. Unnamed lists are now automatically named 'Group.1', 'Group.2', etc...

* Fixed an issue where aggregating by a single id in `collap()` (i.e. `collap(data, ~ id1)`), the id would be coded as factor in the aggregated data.frame. All variables including id's now retain their class and attributes in the aggregated data.

* Added an argument `mean = 0` to `fwithin / W`. This allows simple and grouped centering on an arbitrary mean, `0` being the default. For grouped centering `mean = "overall.mean"` can be specified, which will center data on the overall mean of the data. The logical argument `add.global.mean = TRUE` used to toggle this in *collapse* 1.0.0 is therefore depreciated. 

* Added arguments `mean = 0` (the default) and `sd = 1` (the default) to `fscale / STD`. These arguments now allow to (group) scale and center data to an arbitrary mean and standard deviation. Setting `mean = FALSE` will just scale data while preserving the mean(s). Special options for grouped scaling are `mean = "overall.mean"` (same as `fwithin / W`), and `sd = "within.sd"`, which will scale the data such that the standard deviation of each group is equal to the within- standard deviation (= the standard deviation computed on the group-centered data). Thus group scaling a panel-dataset with `mean = "overall.mean"` and `sd = "within.sd"` harmonizes the data across all groups in terms of both mean and variance. The fast algorithm for variance calculation toggled with `stable.algo = FALSE` was removed from `fscale`. Welford's numerically stable algorithm used by default is fast enough for all practical purposes. The fast algorithm is still available for `fvar` and `fsd`. 

* Added the modulus (`%%`) and subtract modulus (`-%%`) operations to `TRA()`. 

* Added the function `finteraction`, for fast interactions, and `as.character_factor` to coerce a factor, or all factors in a list, to character (analogous to `as.numeric_factor`). Also exported the function `ckmatch`, for matching with error message showing non-matched elements.


# collapse 1.0.0 and earlier

* First version of the package featuring only the functions `collap` and `qsu` based on code shared by Sebastian Martin Krantz on R-devel, February 2019.

* Major rework of the package using Rcpp and data.table internals, introduction of fast statistical functions and operators and expansion of the scope of the package to a broad set of data transformation and exploration tasks. Several iterations of enhancing speed of R code used. Seamless integration of *collapse* with *dplyr*, *plm* and *data.table*. CRAN release of *collapse* 1.0.0 on 19th March 2020. 

