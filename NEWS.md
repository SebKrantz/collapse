# collapse 1.3.2
collapse 1.3.2, released mid September 2020, is a minor update:

* Fixed a small bug in `fNdistinct` for grouped distinct value counts on logical vectors. 

* Additional security for `ftransform`, which now efficiently checks the names of the data and replacement arguments for uniqueness, and also allows computing and transforming list-columns.  

* Added function `ftransformv` to facilitate transforming selected columns with function - a very efficient replacement for `dplyr::mutate_if` and `dplyr::mutate_at`. 

* `frename` now allows additional arguments to be passed to a renaming function.  

# collapse 1.3.1
collapse 1.3.1, released end of August 2020, is a patch for v1.3.0 that takes care of some unit test failures on certain operating systems (mostly because of numeric precision issues). It provides no changes to the code or functionality.

# collapse 1.3.0
collapse 1.3.0, released mid August 2020, is another major update:

### Changes to Functionality

* `dapply` and `BY` now drop all unnecessary attributes if `return = "matrix"` or `return = "data.frame"` are explicitly requested (the default `return = "same"` still seeks to preserve the input data structure).

* `unlist2d` now saves integer rownames if `row.names = TRUE` and a list of matrices without rownames is passed, and `id.factor = TRUE` generates a normal factor not an ordered factor. It is however possible to write `id.factor = "ordered"` to get an ordered factor id.  

* `fdiff` argument `logdiff` renamed to `log`, and taking logs is now done in R (reduces size of C++ code and does not generate as many NaN's). `logdiff` may still be used, but it may be deactivated in the future. Also in the matrix and data.frame methods for `flag`, `fdiff` and `fgrowth`, columns are only stub-renamed if more than one lag/difference/growth rate is computed. 

### Additions

* Added `fnth` for fast (grouped, weighted) n'th element/quantile computations.

* Added `roworder(v)` and `colorder(v)` for fast row and column reordering.  

* Added `frename` and `setrename` for fast and flexible renaming (by reference).  

* Added function `fungroup`, as replacement for `dplyr::ungroup`, intended for use with `fgroup_by`. 

* The shortcut `gvr` was created for `get_vars(..., regex = TRUE)`. Also a helper `.c` was introduced for non-standard concatenation (i.e. `.c(a, b) == c("a", "b")`). 


### Improvements

* `fmedian` now supports weights, computing a decently fast (grouped) weighted median based on radix ordering. 

* `fmode` now has the option to compute min and max mode, the default is still simply the first mode. 

* `fwithin` now supports quasi-demeaning (added argument `theta`) and can thus be used to manually estimate random-effects models. 

* `fmode` and `fNdistinct` have become a bit faster.

* `fgroup_by` now preserves *data.table*'s.

* `funique` is now generic with a default vector and data.frame method, providing fast unique values and rows of data. The default was changed to `sort = FALSE`.   

* `ftransform` now also supports a data.frame as replacement argument, which automatically replaces matching columns and adds unmatched ones. Also `ftransform<-` was created as a more formal replacement method for this feature.

* `collap` columns selected through `cols` argument are returned in the order selected if `keep.col.order = FALSE`. Argument `sort.row` is depreciated, and replace by argument `sort`. In addition the `decreasing` and `na.last` arguments were added and handed down to `GRP.default`. 

* `radixorder` 'sorted' attribute is now always attached.

* `stats::D` which is masked when collapse is attached, is now preserved through methods `D.expression` and `D.call`. 

* `GRP` option `call = FALSE` to omit a call to `match.call` -> minor performance improvement.

* Several small performance improvements through rewriting some internal helper functions in C and reworking some R code. 

* Performance improvements for some helper functions, `setRownames` / `setColnames`, `na_insert` etc.  

* Increased scope of testing statistical functions. The functionality of the package is now secured by 7700 unit tests covering all central bits and pieces. 


# collapse 1.2.1
collapse 1.2.1, released end of May 2020, is a patch for v1.2.0:

* Minor fixes for 1.2.0 issues that prevented correct installation on Mac OS X and a vignette rebuilding error on solaris.

* *fmode.grouped_df* with groups and weights now saves the sum of the weights instead of the max (this makes more sense as the max only applies if all elements are unique). 

# collapse 1.2.0
collapse 1.2.0, released mid May 2020, is a major update of the package - changes and additions:

### Changes to Functionality
* *grouped_df* methods for fast statistical functions now always attach the grouping variables to the output in aggregations, unless argument `keep.group_vars = FALSE`. (formerly grouping variables were only attached if also present in the data. Code hinged on this feature should be adjusted)

* `qF` `ordered` argument default was changed to `ordered = FALSE`, and the `NA` level is only added if `na.exclude = FALSE`. Thus `qF` now behaves exactly like `as.factor`. 

* `Recode` is depreciated in favor of `recode_num` and `recode_char`, it will be removed soon. Similarly `replace_non_finite` was renamed to `replace_Inf`. 

* In `mrtl` and `mctl` the argument `ret` was renamed `return` and now takes descriptive character arguments (the previous version was a direct C++ export and unsafe, code written with these functions should be adjusted). 

* `GRP` argument `order` is depreciated in favor of argument `decreasing`. `order` can still be used but will be removed at some point. 

### Bug Fixes
* Fixed a bug in `flag` where unused factor levels caused a group size error. 

<!-- It is still recommended to remove unused factor levels when programming with factors, some functions check for them, others not. For example `fmean(data, f)` will simply generate a missing row for each unused factor level. If in doubt, use safer `GRP` objects for grouped programming. A general level check for all functions will not be implemented as this requires an additional pass in some cases. -->

### Additions

* Added a suite of functions for fast data manipulation: 
  + `fselect` selects variables from a data frame and is equivalent but much faster than `dplyr::select`.
  + `fsubset` is a much faster version of `base::subset` to subset vectors, matrices and data.frames. The function `ss` was also added as a faster alternative to `[.data.frame`. 
  + `ftransform` is a much faster update of `base::transform`, to transform data frames by adding, modifying or deleting columns. The function `settransform` does all of that by reference.
  + `fcompute` is equivalent to `ftransform` but returns a new data frame containing only the columns computed from an existing one. 
  + `na_omit` is a much faster and enhanced version of `base::na.omit`. 
  + `replace_NA` efficiently replaces missing values in multi-type data. 
  
  
* Added function `fgroup_by` as a much faster version of `dplyr::group_by` based on *collapse* grouping. It attaches a 'GRP' object to a data frame, but only works with *collapse*'s fast functions. This allows *dplyr* like manipulations that are fully *collapse* based and thus significantly faster, i.e. `data %>% fgroup_by(g1,g2) %>% fselect(cola,colb) %>% fmean`. Note that `data %>% dplyr::group_by(g1,g2) %>% dplyr::select(cola,colb) %>% fmean` still works, in which case the *dplyr* 'group' object is converted to 'GRP' as before. However `data %>% fgroup_by(g1,g2) %>% dplyr::summarize(...)` does not work.

* Added function `varying` to efficiently check the variation of multi-type data over a dimension or within groups.

* Added function `radixorder`, same as `base::order(..., method = "radix")` but more accessible and with built-in grouping features. 

* Added functions `seqid` and `groupid` for generalized run-length type id variable generation from grouping and time variables. `seqid` in particular strongly facilitates lagging / differencing irregularly spaced panels using `flag`, `fdiff` etc. 


### Improvements
* Faster grouping with `GRP` and faster factor generation with added radix method + automatic dispatch between hash and radix method. `qF` is now ~ 5x faster than `as.factor` on character and around 30x faster on numeric data. Also `qG` was enhanced. 

* Further slight speed tweaks here and there. 

* `collap` now provides more control for weighted aggregations with additional arguments `w`, `keep.w` and `wFUN` to aggregate the weights as well. The defaults are `keep.w = TRUE` and `wFUN = fsum`. A specialty of `collap` remains that `keep.by` and `keep.w` also work for external objects passed, so code of the form `collap(data, by, FUN, catFUN, w = data$weights)` will now have an aggregated `weights` vector in the first column. 

<!-- In such cases use `keep.w = FALSE` to omit the weights or `collap(data, by, FUN, catFUN, w = ~ weights)` to keep the column order. -->

* `qsu` now also allows weights to be passed in formula i.e. `qsu(data, by = ~ group, pid = ~ panelid, w = ~ weights)`. 

* `fdiff` now supports quasi-differences i.e. $x_t - \rho x_{t-1}$ and quasi-log differences i.e. $log(x_t) - \rho log(x_{t-1})$. an arbitrary $\rho$ can be supplied.

* Added a `Dlog` operator for faster access to log-differences. 

* `fgrowth` has a `scale` argument, the default is `scale = 100` which provides growth rates in percentage terms (as before), but this may now be changed. 

* All statistical and transformation functions now have a hidden list method, so they can be applied to unclassed list-objects as well. An error is however provided in grouped operations with unequal-length columns. 



# collapse 1.1.0
collapse 1.1.0 released early April 2020 - some small fixes and additions:

* Fixed remaining gcc10, LTO and valgrind issues in C/C++ code, and added some more tests (there are now ~ 5300 tests ensuring that *collapse* statistical functions perform as expected).

* Fixed the issue that supplying an unnamed list to `GRP()`, i.e. `GRP(list(v1, v2))` would give an error. Unnamed lists are now automatically named 'Group.1', 'Group.2', etc...

* Fixed an issue where aggregating by a single id in `collap()` (i.e. `collap(data, ~ id1)`), the id would be coded as factor in the aggregated data.frame. All variables including id's now retain their class and attributes in the aggregated data.

* Added weights (`w`) argument to `fsum` and `fprod`. *Note*: `fmedian` will also support weights as soon as I am able to implement a sufficiently fast (i.e. linear time) algorithm. I also hope to introduce (weighted) quantiles. I am happy for any help with these features. 

* Added an argument `mean = 0` to `fwithin / W`. This allows simple and grouped centering on an arbitrary mean, `0` being the default. For grouped centering `mean = "overall.mean"` can be specified, which will center data on the overall mean of the data. The logical argument `add.global.mean = TRUE` used to toggle this in *collapse* 1.0.0 is therefore depreciated. 

* Added arguments `mean = 0` (the default) and `sd = 1` (the default) to `fscale / STD`. These arguments now allow to (group) scale and center data to an arbitrary mean and standard deviation. Setting `mean = FALSE` will just scale data while preserving the mean(s). Special options for grouped scaling are `mean = "overall.mean"` (same as `fwithin / W`), and `sd = "within.sd"`, which will scale the data such that the standard deviation of each group is equal to the within- standard deviation (= the standard deviation computed on the group-centered data). Thus group scaling a panel-dataset with `mean = "overall.mean"` and `sd = "within.sd"` harmonizes the data across all groups in terms of both mean and variance. The fast algorithm for variance calculation toggled with `stable.algo = FALSE` was removed from `fscale`. Welford's numerically stable algorithm used by default is fast enough for all practical purposes. The fast algorithm is still available for `fvar` and `fsd`. 

* Added the modulus (`%%`) and subtract modulus (`-%%`) operations to `TRA()`. 

* Added the function `finteraction`, for fast interactions, and `as.character_factor` to coerce a factor, or all factors in a list, to character (analogous to `as.numeric_factor`). Also exported the function `ckmatch`, for matching with error message showing non-matched elements.


# collapse 1.0.0 and earlier

* First version of the package featuring only the functions `collap` and `qsu` based on code shared by Sebastian Martin Krantz on R-devel, February 2019.

* Major rework of the package using Rcpp and data.table internals, introduction of fast statistical functions and operators and expansion of the scope of the package to a broad set of data transformation and exploration tasks. Several iterations of enhancing speed of R code used. Seamless integration of *collapse* with *dplyr*, *plm* and *data.table*. CRAN release of *collapse* 1.0.0 on 19th March 2020. 

