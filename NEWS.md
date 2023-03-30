# collapse 1.9.4

* Improvements in `get_elem()/has_elem()`: Option `invert = TRUE` is implemented more robustly, and a function passed to `get_elem()/has_elem()` is now applied to all elements in the list, including elements that are themselves list-like. This enables the use of `inherits` to find list-like objects inside a broader list structure e.g. `get_elem(l, inherits, what = "lm")` fetches all linear model objects inside `l`. 

* Fixed a small bug in `descr()` introduced in v1.9.0, producing an error if a data frame contained no numeric columns - because an internal function was not defined in that case. Also, POSIXct columns are handled better in print - preserving the time zone (thanks @cdignam-chwy #392).

* `fmean()` and `fsum()` with `g = NULL`, as well as `TRA()`, `setop()`, and related operators `%r+%`, `%+=%` etc., `setv()` and `fdist()` now utilize Single Instruction Multiple Data (SIMD) vectorization by default (if OpenMP is enabled), enabling potentially very fast computing speeds. Whether these instructions are utilized during compilation depends on your system. In general, if you want to max out *collapse* on your system, consider compiling from source with `CFLAGS += -O3 -march=native -fopenmp` and `CXXFLAGS += -O3 -march=native` in your [`.R/Makevars`](https://cran.r-project.org/doc/manuals/r-devel/R-admin.html#Customizing-package-compilation).

# collapse 1.9.3

* Added functions `fduplicated()` and `any_duplicated()`, for vectors and lists / data frames. Thanks @NicChr (#373)

* `sort` option added to `set_collapse()` to be able to set unordered grouping as a default. E.g. setting `set_collapse(sort = FALSE)` will affect `collap()`, `BY()`, `GRP()`, `fgroup_by()`, `qF()`, `qG()`, `finteraction()`, `qtab()` and internal use of these functions for ad-hoc grouping in fast statistical functions. Other uses of `sort`, for example in `funique()` where the default is `sort = FALSE`, are not affected by the global default setting. 

* Fixed a small bug in `group()` / `funique()` resulting in an unnecessary memory allocation error in rare cases. Thanks @NicChr (#381).

# collapse 1.9.2

* Further fix to an Address Sanitizer issue as required by CRAN (eliminating an unused out of bounds access at the end of a loop). 

* `qsu()` finally has a grouped_df method. 

* Added options `option("collapse_nthreads")` and `option("collapse_na.rm")`, which allow you to load *collapse* with different defaults e.g. through an `.Rprofile` or `.fastverse` configuration file. Once *collapse* is loaded, these options take no effect, and users need to use `set_collapse()` to change `.op[["nthreads"]]` and `.op[["na.rm"]]` interactively. 

* Exported method `plot.psmat()` (can be useful to plot time series matrices).


# collapse 1.9.1

* Fixed minor C/C++ issues flagged by CRAN's detailed checks.

* Added functions `set_collapse()` and `get_collapse()`, allowing you to globally set defaults for the `nthreads` and `na.rm` arguments to all functions in the package. E.g. `set_collapse(nthreads = 4, na.rm = FALSE)` could be a suitable setting for larger data without missing values. This is implemented using an internal environment by the name of `.op`, such that these defaults are received using e.g. `.op[["nthreads"]]`, at the computational cost of a few nanoseconds (8-10x faster than `getOption("nthreads")` which would take about 1 microsecond). `.op` is not accessible by the user, so function `get_collapse()` can be used to retrieve settings. Exempt from this are functions `.quantile`, and a new function `.range` (alias of `frange`), which go directly to C for maximum performance in repeated executions, and are not affected by these global settings. Function `descr()`, which internally calls a bunch of statistical functions, is also not affected by these settings. 

* Further improvements in thread safety for `fsum()` and `fmean()` in grouped computations across data frame columns. All OpenMP enabled functions in *collapse* can now be considered thread safe i.e. they pass the full battery of tests in multithreaded mode. 

# collapse 1.9.0

*collapse* 1.9.0 released mid of January 2023, provides improvements in performance and versatility in many areas, as well as greater statistical capabilities, most notably efficient (grouped, weighted) estimation of sample quantiles. 

### Changes to functionality

* All functions renamed in *collapse* 1.6.0 are now depreciated, to be removed end of 2023. These functions had already been giving messages since v1.6.0. See `help("collapse-renamed")`.

* The lead operator `F()` is not exported anymore from the package namespace, to avoid clashes with `base::F` flagged by multiple people. The operator is still part of the package and can be accessed using `collapse:::F`. I have also added an option `"collapse_export_F"`, such that setting `options(collapse_export_F = TRUE)` before loading the package exports the operator as before. Thanks @matthewross07 (#100), @edrubin (#194), and @arthurgailes (#347). 

* Function `fnth()` has a new default `ties = "q7"`, which gives the same result as `quantile(..., type = 7)` (R's default). More details below. 

### Bug Fixes

* `fmode()` gave wrong results for singleton groups (groups of size 1) on *unsorted* data. I had optimized `fmode()` for singleton groups to directly return the corresponding element, but it did not access the element through the (internal) ordering vector, so the first element/row of the entire vector/data was taken. The same mistake occurred for `fndistinct` if singleton groups were `NA`, which were counted as `1` instead of `0` under the `na.rm = TRUE` default (provided the first element of the vector/data was not `NA`). The mistake did not occur with data sorted by the groups, because here the data pointer already pointed to the first element of the group. (My apologies for this bug, it took me more than half a year to discover it, using *collapse* on a daily basis, and it escaped 700 unit tests as well).

* Function `groupid(x, na.skip = TRUE)` returned uninitialized first elements if the first values in `x` where `NA`. Thanks for reporting @Henrik-P (#335). 

* Fixed a bug in the `.names` argument to `across()`. Passing a naming function such as `.names = function(c, f) paste0(c, "-", f)` now works as intended i.e. the function is applied to all combinations of columns (c) and functions (f) using `outer()`. Previously this was just internally evaluated as `.names(cols, funs)`, which did not work if there were multiple cols and multiple funs. There is also now a possibility to set `.names = "flip"`, which names columns `f_c` instead of `c_f`. 

* `fnrow()` was rewritten in C and also supports data frames with 0 columns. Similarly for `seq_row()`. Thanks @NicChr (#344). 

### Additions

* Added functions `fcount()` and `fcountv()`: a versatile and blazing fast alternative to `dplyr::count`. It also works with vectors, matrices, as well as grouped and indexed data. 

* Added function `fquantile()`: Fast (weighted) continuous quantile estimation (methods 5-9 following Hyndman and Fan (1996)), implemented fully in C based on quickselect and radixsort algorithms, and also supports an ordering vector as optional input to speed up the process. It is up to 2x faster than `stats::quantile` on larger vectors, but also especially fast on smaller data, where the R overhead of `stats::quantile` becomes burdensome. For maximum performance during repeated executions, a programmers version `.quantile()` with different defaults is also provided. 

* Added function `fdist()`: A fast and versatile replacement for `stats::dist`. It computes a full euclidian distance matrix around 4x faster than `stats::dist` in serial mode, with additional gains possible through multithreading along the distance matrix columns (decreasing thread loads as the matrix is lower triangular). It also supports computing the distance of a matrix with a single row-vector, or simply between two vectors. E.g. `fdist(mat, mat[1, ])` is the same as `sqrt(colSums((t(mat) - mat[1, ])^2)))`, but about 20x faster in serial mode, and `fdist(x, y)` is the same as `sqrt(sum((x-y)^2))`, about 3x faster in serial mode. In both cases (sub-column level) multithreading is available. *Note* that `fdist` does not skip missing values i.e. `NA`'s will result in `NA` distances. There is also no internal implementation for integers or data frames. Such inputs will be coerced to numeric matrices. 

* Added function `GRPid()` to easily fetch the group id from a grouping object, especially inside grouped `fmutate()` calls. This addition was warranted especially by the new improved `fnth.default()` method which allows orderings to be supplied for performance improvements. See commends on `fnth()` and the example provided below. 

* `fsummarize()` was added as a synonym to `fsummarise`. Thanks @arthurgailes for the PR. 

* **C API**: *collapse* exports around 40 C functions that provide functionality that is either convenient or rather complicated to implement from scratch. The exported functions can be found at the bottom of `src/ExportSymbols.c`. The API does not include the *Fast Statistical Functions*, which I thought are too closely related to how *collapse* works internally to be of much use to a C programmer (e.g. they expect grouping objects or certain kinds of integer vectors). But you are free to request the export of additional functions, including C++ functions. 

### Improvements

* `fnth()` and `fmedian()` were rewritten in C, with significant gains in performance and versatility. Notably, `fnth()` now supports (grouped, weighted) continuous quantile estimation like `fquantile()` (`fmedian()`, which is a wrapper around `fnth()`, can also estimate various quantile based weighted medians). The new default for `fnth()` is `ties = "q7"`, which gives the same result as `(f)quantile(..., type = 7)` (R's default). OpenMP multithreading across groups is also much more effective in both the weighted and unweighted case. Finally, `fnth.default` gained an additional argument `o` to pass an ordering vector, which can dramatically speed up repeated invocations of the function on the dame data: 

  ```r
  # Estimating multiple weighted-grouped quantiles on mpg: pre-computing an ordering provides extra speed. 
  mtcars %>% fgroup_by(cyl, vs, am) %>%
      fmutate(o = radixorder(GRPid(), mpg)) %>% # On grouped data, need to account for GRPid()
      fsummarise(mpg_Q1 = fnth(mpg, 0.25, o = o, w = wt),
                 mpg_median = fmedian(mpg, o = o, w = wt),
                 mpg_Q3 = fnth(mpg, 0.75, o = o, w = wt))
  # Note that without weights this is not always faster. Quickselect can be very efficient, so it depends 
  # on the data, the number of groups, whether they are sorted (which speeds up radixorder), etc...
  ```

* `BY` now supports data-length arguments to be passed e.g. `BY(mtcars, mtcars$cyl, fquantile, w = mtcars$wt)`, making it effectively a generic grouped `mapply` function as well. Furthermore, the grouped_df method now also expands grouping columns for output length > 1. 

* `collap()`, which internally uses `BY` with non-*Fast Statistical Functions*, now also supports arbitrary further arguments passed down to functions to be split by groups. Thus users can also apply custom weighted functions with `collap()`. Furthermore, the parsing of the `FUN`, `catFUN` and `wFUN` arguments was improved and brought in-line with the parsing of `.fns` in `across()`. The main benefit of this is that *Fast Statistical Functions* are now also detected and optimizations carried out when passed in a list providing a new name e.g. `collap(data, ~ id, list(mean = fmean))` is now optimized! Thanks @ttrodrigz (#358) for requesting this.    

* `descr()`, by virtue of `fquantile` and the improvements to `BY`, supports full-blown grouped and weighted descriptions of data. This is implemented through additional `by` and `w` arguments. The function has also been turned into an S3 generic, with a default and a 'grouped_df' method. The 'descr' methods `as.data.frame` and `print` also feature various improvements, and a new `compact` argument to `print.descr`, allowing a more compact printout. Users will also notice improved performance, mainly due to `fquantile`: on the M1 `descr(wlddev)` is now 2x faster than `summary(wlddev)`, and 41x faster than `Hmisc::describe(wlddev)`. Thanks @statzhero for the request (#355).


* `radixorder` is about 25% faster on characters and doubles. This also benefits grouping performance. Note that `group()` may still be substantially faster on unsorted data, so if performance is critical try the `sort = FALSE` argument to functions like `fgroup_by` and compare. 

* Most list processing functions are noticeably faster, as checking the data types of elements in a list is now also done in C, and I have made some improvements to *collapse*'s version of `rbindlist()` (used in `unlist2d()`, and various other places). 

* `fsummarise` and `fmutate` gained an ability to evaluate arbitrary expressions that result in lists / data frames without the need to use `across()`. For example: `mtcars |> fgroup_by(cyl, vs, am) |> fsummarise(mctl(cor(cbind(mpg, wt, carb)), names = TRUE))` or `mtcars |> fgroup_by(cyl) |> fsummarise(mctl(lmtest::coeftest(lm(mpg ~ wt + carb)), names = TRUE))`. There is also the possibility to compute expressions using `.data` e.g. `mtcars |> fgroup_by(cyl) |> fsummarise(mctl(lmtest::coeftest(lm(mpg ~ wt + carb, .data)), names = TRUE))` yields the same thing, but is less efficient because the whole dataset (including 'cyl') is split by groups. For greater efficiency and convenience, you can pre-select columns using a global `.cols` argument, e.g. `mtcars |> fgroup_by(cyl, vs, am) |> fsummarise(mctl(cor(.data), names = TRUE), .cols = .c(mpg, wt, carb))` gives the same as above. Three *Notes* about this:

  + No grouped vectorizations for fast statistical functions i.e. the entire expression is evaluated for each group. (Let me know if there are applications where vectorization would be possible and beneficial. I can't think of any.)  
  + All elements in the result list need to have the same length, or, for `fmutate`, have the same length as the data (in each group). 
  + If `.data` is used, the entire expression (`expr`) will be turned into a function of `.data` (`function(.data) expr`), which means columns are only available when accessed through `.data` e.g. `.data$col1`. 
  
* `fsummarise` supports computations with mixed result lengths e.g. `mtcars |> fgroup_by(cyl) |> fsummarise(N = GRPN(), mean_mpg = fmean(mpg), quantile_mpg = fquantile(mpg))`, as long as all computations result in either length 1 or length k vectors, where k is the maximum result length (e.g. for `fquantile` with default settings k = 5).   

* List extraction function `get_elem()` now has an option `invert = TRUE` (default `FALSE`) to remove matching elements from a (nested) list. Also the functionality of argument `keep.class = TRUE` is implemented in a better way, such that the default `keep.class = FALSE` toggles classes from (non-matched) list-like objects inside the list to be removed. 

* `num_vars()` has become a bit smarter: columns of class 'ts' and 'units' are now also recognized as numeric. In general, users should be aware that `num_vars()` does not regard any R methods defined for `is.numeric()`, it is implemented in C and simply checks whether objects are of type integer or double, and do not have a class. The addition of these two exceptions now guards against two common cases where `num_vars()` may give undesirable outcomes. Note that `num_vars()`  is also called in `collap()` to distinguish between numeric (`FUN`) and non-numeric (`catFUN`) columns. 

* Improvements to `setv()` and `copyv()`, making them more robust to borderline cases: `integer(0)` passed to `v` does nothing (instead of error), and it is also possible to pass a single real index if `vind1 = TRUE` i.e. passing `1` instead of `1L` does not produce an error. 

* `alloc()` now works with all types of objects i.e. it can replicate any object. If the input is non-atomic, atomic with length > 1 or `NULL`, the output is a list of these objects, e.g. `alloc(NULL, 10)` gives a length 10 list of `NULL` objects, or `alloc(mtcars, 10)` gives a list of `mtcars` datasets. Note that in the latter case the datasets are not deep-copied, so no additional memory is consumed. 

* `missing_cases()` and `na_omit()` have gained an argument `prop = 0`, indicating the proportion of values missing for the case to be considered missing/to be omitted. The default value of `0` indicates that at least 1 value must be missing. Of course setting `prop = 1` indicates that all values must be missing. For data frames/lists the checking is done efficiently in C. For matrices this is currently still implemented using `rowSums(is.na(X)) >= max(as.integer(prop * ncol(X)), 1L)`, so the performance is less than optimal. 

* `missing_cases()` has an extra argument `count = FALSE`. Setting `count = TRUE` returns the case-wise missing value count (by `cols`).  

* Functions `frename()` and `setrename()` have an additional argument `.nse = TRUE`, conforming to the default non-standard evaluation of tagged vector expressions e.g. `frename(mtcars, mpg = newname)` is the same as `frename(mtcars, mpg = "newname")`. Setting `.nse = FALSE` allows `newname` to be a variable holding a name e.g. `newname = "othername"; frename(mtcars, mpg = newname, .nse = FALSE)`. Another use of the argument is that a (named) character vector can now be passed to the function to rename a (subset of) columns e.g. `cvec = letters[1:3]; frename(mtcars, cvec, cols = 4:6, .nse = FALSE)` (this works even with `.nse = TRUE`), and `names(cvec) = c("cyl", "vs", "am"); frename(mtcars, cvec, .nse = FALSE)`. Furthermore, `setrename()` now also returns the renamed data invisibly, and `relabel()` and `setrelabel()` have also gained similar flexibility to allow (named) lists or vectors of variable labels to be passed. *Note* that these function have no NSE capabilities, so they work essentially like `frename(..., .nse = FALSE)`.

* Function `add_vars()` became a bit more flexible and also allows single vectors to be added with tags e.g. `add_vars(mtcars, log_mpg = log(mtcars$mpg), STD(mtcars))`, similar to `cbind`. However `add_vars()` continues to not replicate length 1 inputs. 
* Safer multithreading: OpenMP multithreading over parts of the R API is minimized, reducing errors that occurred especially when multithreading across data frame columns. Also the number of threads supplied by the user to all OpenMP enabled functions is ensured to not exceed either of `omp_get_num_procs()`, `omp_get_thread_limit()`, and `omp_get_max_threads()`. 


# collapse 1.8.9

* Fixed some warnings on rchk and newer C compilers (LLVM clang 10+). 

* `.pseries` / `.indexed_series` methods also change the implicit class of the vector (attached after `"pseries"`), if the data type changed. e.g. calling a function like `fgrowth` on an integer pseries changed the data type to double, but the "integer" class was still attached after "pseries".

* Fixed bad testing for SE inputs in `fgroup_by()` and `findex_by()`. See #320. 

* Added `rsplit.matrix` method. 

* `descr()` now by default also reports 10% and 90% quantiles for numeric variables (in line with STATA's detailed summary statistics), and can also be applied to 'pseries' / 'indexed_series'. Furthermore, `descr()` itself now has an argument `stepwise` such that `descr(big_data, stepwise = TRUE)` yields computation of summary statistics on a variable-by-variable basis (and the finished 'descr' object is returned invisibly). The printed result is thus identical to `print(descr(big_data), stepwise = TRUE)`, with the difference that the latter first does the entire computation whereas the former computes statistics on demand.   

<!-- * Added method `funique.grouped_df`. ???? -->

* Function `ss()` has a new argument `check = TRUE`. Setting `check = FALSE` allows subsetting data frames / lists with positive integers without checking whether integers are positive or in-range. For programmers. 

* Function `get_vars()` has a new argument `rename` allowing select-renaming of columns in standard evaluation programming, e.g. `get_vars(mtcars, c(newname = "cyl", "vs", "am"), rename = TRUE)`. The default is `rename = FALSE`, to warrant full backwards compatibility. See #327. 

* Added helper function `setattrib()`, to set a new attribute list for an object by reference + invisible return. This is different from the existing function `setAttrib()` (note the capital A), which takes a shallow copy of list-like objects and returns the result.   


# collapse 1.8.8

* `flm` and `fFtest` are now internal generic with an added formula method e.g. `flm(mpg ~ hp + carb, mtcars, weights = wt)` or `fFtest(mpg ~ hp + carb | vs + am, mtcars, weights = wt)` in addition to the programming interface. Thanks to Grant McDermott for suggesting. 

* Added method `as.data.frame.qsu`, to efficiently turn the default array outputs from `qsu()` into tidy data frames. 

* Major improvements to `setv` and `copyv`, generalizing the scope of operations that can be performed to all common cases. This means that even simple base R operations such as `X[v] <- R` can now be done significantly faster using `setv(X, v, R)`. 

* `n` and `qtab` can now be added to `options("collapse_mask")` e.g. `options(collapse_mask = c("manip", "helper", "n", "qtab"))`. This will export a function `n()` to get the (group) count in `fsummarise` and `fmutate` (which can also always be done using `GRPN()` but `n()` is more familiar to *dplyr* users), and will mask `table()` with `qtab()`, which is principally a fast drop-in replacement, but with some different further arguments. 

* Added C-level helper function `all_funs`, which fetches all the functions called in an expression, similar to `setdiff(all.names(x), all.vars(x))` but better because it takes account of the syntax. For example let `x = quote(sum(sum))` i.e. we are summing a column named `sum`. Then `all.names(x) = c("sum", "sum")` and `all.vars(x) = "sum"` so that the difference is `character(0)`, whereas `all_funs(x)` returns `"sum"`. This function makes *collapse* smarter when parsing expressions in `fsummarise` and `fmutate` and deciding which ones to vectorize. 

<!--
For example, take the `.OPERATOR_FUN` such as `W()` for within-transforming/centering data. These receive vectorized execution in `fmutate` e.g. `mtcars |> gby(cyl) |> fmutate(mpg = W(mpg))` executes `W(mpg, g = GRP(.))`. Now suppose you have a column called `W` 
-->

# collapse 1.8.7

* Fixed a bug in `fscale.pdata.frame` where the default C++ method was being called instead of the list method (i.e. the method didn't work at all).

* Fixed 2 minor rchk issues (the remaining ones are spurious).

* `fsum` has an additional argument `fill = TRUE` (default `FALSE`) that initializes the result vector with `0` instead of `NA` when `na.rm = TRUE`, so that `fsum(NA, fill = TRUE)` gives `0` like `base::sum(NA, na.rm = TRUE)`. 

* Slight performance increase in `fmean` with groups if `na.rm = TRUE` (the default). 

* Significant performance improvement when using base R expressions involving multiple functions and one column e.g. `mid_col = (min(col) + max(col)) / 2` or `lorentz_col = cumsum(sort(col)) / sum(col)` etc. inside `fsummarise` and `fmutate`. Instead of evaluating such expressions on a data subset of one column for each group, they are now turned into a function e.g. `function(x) cumsum(sort(x)) / sum(x)` which is applied to a single vector split by groups. 
* `fsummarise` now also adds groupings to transformation functions and operators, which allows full vectorization of more complex tasks involving transformations which are subsequently aggregated. A prime example is grouped bivariate linear model fitting, which can now be done using `mtcars |> fgroup_by(cyl) |> fsummarise(slope = fsum(W(mpg), hp) / fsum(W(mpg)^2))`. Before 1.8.7 it was necessary to do a mutate step first e.g. `mtcars |> fgroup_by(cyl) |> fmutate(dm_mpg = W(mpg)) |> fsummarise(slope = fsum(dm_mpg, hp) / fsum(dm_mpg^2))`, because `fsummarise` did not add groupings to transformation functions like `fwithin/W`. Thanks to Brodie Gaslam for making me aware of this. 

* Argument `return.groups` from `GRP.default` is now also available in `fgroup_by`, allowing grouped data frames without materializing the unique grouping columns. This allows more efficient mutate-only operations e.g. `mtcars |> fgroup_by(cyl, return.groups = FALSE) |> fmutate(across(hp:carb, fscale))`. Similarly for aggregation with dropping of grouping columns `mtcars |> fgroup_by(cyl, return.groups = FALSE) |> fmean()` is equivalent and faster than `mtcars |> fgroup_by(cyl) |> fmean(keep.group_vars = FALSE)`.

# collapse 1.8.6

* Fixed further minor issues: 
  - some inline functions in TRA.c needed to be declared 'static' to be local in scope (#275)
  - timeid.Rd now uses *zoo* package conditionally and limits size of printout

# collapse 1.8.5

* Fixed some issues flagged by CRAN:
  - Installation on some linux distributions failed because omp.h was included after Rinternals.h
  - Some signed integer overflows while running tests caused UBSAN warnings. (This happened inside a hash function where overflows are not a problem. I changed to unsigned int to avoid the UBSAN warning.)
  - Ensured that package passes R CMD Check without suggested packages
  
# collapse 1.8.4

* Makevars text substitution hack to have CRAN accept a package that combines C, C++ and OpenMP. Thanks also to @MichaelChirico for pointing me in the right direction. 

# collapse 1.8.3

* Significant speed improvement in `qF/qG` (factor-generation) for character vectors with more than 100,000 obs and many levels if `sort = TRUE` (the default). For details see the `method` argument of `?qF`. 

* Optimizations in `fmode` and `fndistinct` for singleton groups. 

# collapse 1.8.2

* Fixed some rchk issues found by Thomas Kalibera from CRAN. 

* faster `funique.default` method. 

* `group` now also internally optimizes on 'qG' objects. 

# collapse 1.8.1

* Added function `fnunique` (yet another alternative to `data.table::uniqueN`, `kit::uniqLen` or `dplyr::n_distinct`, and principally a simple wrapper for `attr(group(x), "N.groups")`). At present `fnunique` generally outperforms the others on data frames.

* `finteraction` has an additional argument `factor = TRUE`. Setting `factor = FALSE` returns a 'qG' object, which is more efficient if just an integer id but no factor object itself is required.  

* Operators (see `.OPERATOR_FUN`) have been improved a bit such that id-variables selected in the `.data.frame` (`by`, `w` or `t` arguments) or `.pdata.frame` methods (variables in the index) are not computed upon even if they are numeric (since the default is `cols = is.numeric`). In general, if `cols` is a function used to select columns of a certain data type, id variables are excluded from computation even if they are of that data type. It is still possible to compute on id variables by explicitly selecting them using names or indices passed to `cols`, or including them in the lhs of a formula passed to `by`. 

* Further efforts to facilitate adding the group-count in `fsummarise` and `fmutate`: 
  - if `options(collapse_mask = "all")` before loading the package, an additional function `n()` is exported that works just like `dplyr:::n()`. 
  - otherwise the same can now always be done using `GRPN()`. The previous uses of `GRPN` are unaltered i.e. `GRPN` can also:
    + fetch group sizes directly grouping object or grouped data frame i.e. `data |> gby(id) |> GRPN()` or `data %>% gby(id) %>% ftransform(N = GRPN(.))` (note the dot). 
    + compute group sizes on the fly, for example `fsubset(data, GRPN(id) > 10L)` or `fsubset(data, GRPN(list(id1, id2)) > 10L)` or `GRPN(data, by = ~ id1 + id2)`.

# collapse 1.8.0

*collapse* 1.8.0, released mid of May 2022, brings enhanced support for indexed computations on time series and panel data by introducing flexible 'indexed_frame' and 'indexed_series' classes and surrounding infrastructure, sets a modest start to OpenMP multithreading as well as data transformation by reference in statistical functions, and enhances the packages descriptive statistics toolset. 

### Changes to functionality

* Functions `Recode`, `replace_non_finite`, depreciated since *collapse* v1.1.0 and `is.regular`, depreciated since *collapse* v1.5.1 and clashing with a more important function in the *zoo* package, are now removed.

* *Fast Statistical Functions* operating on numeric data (such as `fmean`, `fmedian`, `fsum`, `fmin`, `fmax`, ...) now preserve attributes in more cases. Previously these functions did not preserve attributes for simple computations using the default method, and only preserved attributes in grouped computations if `!is.object(x)` (see NEWS section for collapse 1.4.0). This meant that `fmin` and `fmax` did not preserve the attributes of Date or POSIXct objects, and none of these functions preserved 'units' objects (used a lot by the *sf* package). Now, attributes are preserved if `!inherits(x, "ts")`, that is the new default of these functions is to generally keep attributes, except for 'ts' objects where doing so obviously causes an unwanted error (note that 'xts' and others are handled by the matrix or data.frame method where other principles apply, see NEWS for 1.4.0). An exception are the functions `fnobs` and `fndistinct` where the previous default is kept. 

* *Time Series Functions* `flag`, `fdiff`, `fgrowth` and `psacf/pspacf/psccf` (and the operators `L/F/D/Dlog/G`) now internally process time objects passed to the `t` argument (where `is.object(t) && is.numeric(unclass(t))`) via a new function called `timeid` which turns them into integer vectors based on the greatest common divisor (GCD) (see below). Previously such objects were converted to factor. This can change behavior of code e.g. a 'Date' variable representing monthly data may be regular when converted to factor, but is now irregular and regarded as daily data (with a GCD of 1) because of the different day counts of the months. Users should fix such code by either by calling `qG` on the time variable (for grouping / factor-conversion) or using appropriate classes e.g. `zoo::yearmon`. Note that plain numeric vectors where `!is.object(t)` are still used directly for indexation without passing them through `timeid` (which can still be applied manually if desired).

* `BY` now has an argument `reorder = TRUE`, which casts elements in the original order if `NROW(result) == NROW(x)` (like `fmutate`). Previously the result was just in order of the groups, regardless of the length of the output. To obtain the former outcome users need to set `reorder = FALSE`. 

* `options("collapse_DT_alloccol")` was removed, the default is now fixed at 100. The reason is that *data.table* automatically expands the range of overallocated columns if required (so the option is not really necessary), and calling R options from C slows down C code and can cause problems in parallel code.

### Bug Fixes

* Fixed a bug in `fcumsum` that caused a segfault during grouped operations on larger data, due to flawed internal memory allocation. Thanks @Gulde91 for reporting #237. 

* Fixed a bug in `across` caused by two `function(x)` statements being passed in a list e.g. `mtcars |> fsummarise(acr(mpg, list(ssdd = function(x) sd(x), mu = function(x) mean(x))))`. Thanks @trang1618 for reporting #233. 

* Fixed an issue in `across()` when logical vectors were used to select column on grouped data e.g. `mtcars %>% gby(vs, am) %>% smr(acr(startsWith(names(.), "c"), fmean))` now works without error. 

* `qsu` gives proper output for length 1 vectors e.g. `qsu(1)`. 

* *collapse* depends on R > 3.3.0, due to the use of newer C-level macros introduced then. The earlier indication of R > 2.1.0 was only based on R-level code and misleading. Thanks @ben-schwen for reporting #236. I will try to maintain this dependency for as long as possible, without being too restrained by development in R's C API and the ALTREP system in particular, which *collapse* might utilize in the future.

### Additions

* Introduction of 'indexed_frame','indexed_series' and 'index_df' classes: fast and flexible indexed time series and panel data classes that inherit from *plm*'s 'pdata.frame', 'pseries' and 'pindex' classes. These classes take full advantage of *collapse*'s computational infrastructure, are class-agnostic i.e. they can be superimposed upon any data frame or vector/matrix like object while maintaining most of the functionality of that object, support both time series and panel data, natively handle irregularity, and supports ad-hoc computations inside arbitrary data masking functions and model formulas. This infrastructure comprises of additional functions and methods, and modification of some existing functions and 'pdata.frame' / 'pseries' methods. 

  - New functions: `findex_by/iby`, `findex/ix`, `unindex`, `reindex`, `is_irregular`, `to_plm`. 
  
  - New methods: `[.indexed_series`, `[.indexed_frame`, `[<-.indexed_frame`, `$.indexed_frame`, 
  `$<-.indexed_frame`, `[[.indexed_frame`, `[[<-.indexed_frame`, `[.index_df`, `fsubset.pseries`, `fsubset.pdata.frame`, `funique.pseries`, `funique.pdata.frame`, `roworder(v)` (internal) `na_omit` (internal), `print.indexed_series`, `print.indexed_frame`, `print.index_df`, `Math.indexed_series`, `Ops.indexed_series`.  
  
  - Modification of 'pseries' and 'pdata.frame' methods for functions `flag/L/F`, `fdiff/D/Dlog`, `fgrowth/G`, `fcumsum`, `psmat`, `psacf/pspacf/psccf`, `fscale/STD`, `fbetween/B`, `fwithin/W`, `fhdbetween/HDB`, `fhdwithin/HDW`, `qsu` and `varying` to take advantage of 'indexed_frame' and 'indexed_series' while continuing to work as before with 'pdata.frame' and 'pseries'.
  
  For more information and details see `help("indexing")`.
  
* Added function `timeid`: Generation of an integer-id/time-factor from time or date sequences represented by integer of double vectors (such as 'Date', 'POSIXct', 'ts', 'yearmon', 'yearquarter' or plain integers / doubles) by a numerically quite robust greatest common divisor method (see below). This function is used internally in `findex_by`, `reindex` and also in evaluation of the `t` argument to functions like `flag`/`fdiff`/`fgrowth` whenever `is.object(t) && is.numeric(unclass(t))` (see also note above). 

* Programming helper function `vgcd` to efficiently compute the greatest common divisor from a vector or positive integer or double values (which should ideally be unique and sorted as well, `timeid` uses `vgcd(sort(unique(diff(sort(unique(na_rm(x)))))))`). Precision for doubles is up to 6 digits. 

* Programming helper function `frange`: A significantly faster alternative to `base::range`, which calls both `min` and `max`. Note that `frange` inherits *collapse*'s global `na.rm = TRUE` default. 

* Added function `qtab/qtable`: A versatile and computationally more efficient alternative to `base::table`. Notably, it also supports tabulations with frequency weights, and computation of a statistic over combinations of variables. Objects are of class 'qtab' that inherits from 'table'. Thus all 'table' methods apply to it.

* `TRA` was rewritten in C, and now has an additional argument `set = TRUE` which toggles data transformation by reference. The function `setTRA` was added as a shortcut which additionally returns the result invisibly. Since `TRA` is usually accessed internally through the like-named argument to *Fast Statistical Functions*, passing `set = TRUE` to those functions yields an internal call to `setTRA`. For example `fmedian(num_vars(iris), g = iris$Species, TRA = "-", set = TRUE)` subtracts the species-wise median from the numeric variables in the iris dataset, modifying the data in place and returning the result invisibly. Similarly the argument can be added in other workflows such as `iris |> fgroup_by(Species) |> fmutate(across(1:2, fmedian, set = TRUE))` or  `mtcars |> ftransform(mpg = mpg %+=% hp, wt = fsd(wt, cyl, TRA = "replace_fill", set = TRUE))`. Note that such chains must be ended by `invisible()` if no printout is wanted. 

* Exported helper function `greorder`, the companion to `gsplit` to reorder output in `fmutate` (and now also in `BY`): let `g` be a 'GRP' object (or something coercible such as a vector) and `x` a vector, then `greorder` orders data in `y = unlist(gsplit(x, g))` such that `identical(greorder(y, g), x)`. 

### Improvements

* `fmean`, `fprod`, `fmode` and `fndistinct` were rewritten in C, providing performance improvements, particularly in `fmode` and `fndistinct`, and improvements for integers in `fmean` and `fprod`. 

* OpenMP multithreading in `fsum`, `fmean`, `fmedian`, `fnth`, `fmode` and `fndistinct`, implemented via an additional `nthreads` argument. The default is to use 1 thread, which internally calls a serial version of the code in `fsum` and `fmean` (thus no change in the default behavior). The plan is to slowly roll this out over all statistical functions and then introduce options to set alternative global defaults. Multi-threading internally works different for different functions, see the `nthreads` argument documentation of a particular function. Unfortunately I currently cannot guarantee thread safety, as parallelization of complex loops entails some tricky bugs and I have limited time to sort these out. So please report bugs, and if you happen to have experience with OpenMP please consider examining the code and making some suggestions. 

* `TRA` has an additional option `"replace_NA"`, e.g. `wlddev |> fgroup_by(iso3c) |> fmutate(across(PCGDP:POP, fmedian, TRA = "replace_NA"))` performs median value imputation of missing values. Similarly for a matrix `X <- matrix(na_insert(rnorm(1e7)), ncol = 100)`, `fmedian(X, TRA = "replace_NA", set = TRUE)` (column-wise median imputation by reference).

* All *Fast Statistical Functions* support zero group sizes (e.g. grouping with a factor that has unused levels will always produce an output of length `nlevels(x)` with `0` or `NA` elements for the unused levels). Previously this produced an error message with counting/ordinal functions `fmode`, `fndistinct`, `fnth` and `fmedian`. 

* 'GRP' objects now also contain a 'group.starts' item in the 8'th slot that gives the first positions of the unique groups, and is returned alongside the groups whenever `return.groups = TRUE`. This now benefits `ffirst` when invoked with `na.rm = FALSE`, e.g. `wlddev %>% fgroup_by(country) %>% ffirst(na.rm = FALSE)` is now just as efficient as `funique(wlddev, cols = "country")`. Note that no additional computing cost is incurred by preserving the 'group.starts' information. 

* Conversion methods `GRP.factor`, `GRP.qG`, `GRP.pseries`, `GRP.pdata.frame` and `GRP.grouped_df` now also efficiently check if grouping vectors are sorted (the information is stored in the "ordered" element of 'GRP' objects). This leads to performance improvements in `gsplit` / `greorder` and dependent functions such as `BY` and `rsplit` if factors are sorted. 

* `descr()` received some performance improvements (up to 2x for categorical data), and has an additional argument `sort.table`, allowing frequency tables for categorical variables to be sorted by frequency (`"freq"`) or by table values (`"value"`). The new default is (`"freq"`), which presents tables in decreasing order of frequency. A method `[.descr` was added allowing 'descr' objects to be subset like a list. The print method was also enhanced, and by default now prints 14 values with the highest frequency and groups the remaining values into a single `... %s Others` category. Furthermore, if there are any missing values in the column, the percentage of values missing is now printed behind `Statistics `. Additional arguments `reverse` and `stepwise` allow printing in reverse order and/or one variable at a time. 

* `whichv` (and operators `%==%`, `%!=%`) now also support comparisons of equal-length arguments e.g. `1:3 %==% 1:3`. Note that this should not be used to compare 2 factors. 

* Added some code to the `.onLoad` function that checks for the existence of a `.fastverse` configuration file containing a setting for `_opt_collapse_mask`: If found the code makes sure that the option takes effect before the package is loaded. This means that inside projects using the *fastverse* and `options("collapse_mask")` to replace base R / *dplyr* functions, *collapse* cannot be loaded without the masking being applied, making it more secure to utilize this feature. For more information about function masking see `help("collapse-options")` and for `.fastverse` configuration files see the [fastverse vignette](https://fastverse.github.io/fastverse/articles/fastverse_intro.html#custom-fastverse-configurations-for-projects). 

* Added hidden `.list` methods for `fhdwithin/HDW` and `fhdbetween/HDB`. As for the other `.FAST_FUN` this is just a wrapper for the data frame method and meant to be used on unclassed data frames. 

* `ss()` supports unnamed lists / data frames. 

* The `t` and `w` arguments in 'grouped_df' methods (NSE) and where formula input is allowed, supports ad-hoc transformations. E.g. `wlddev %>% gby(iso3c) %>% flag(t = qG(date))` or `L(wlddev, 1, ~ iso3c, ~qG(date))`, similarly `qsu(wlddev, w = ~ log(POP))`, `wlddev %>% gby(iso3c) %>% collapg(w = log(POP))` or `wlddev %>% gby(iso3c) %>% nv() %>% fmean(w = log(POP))`.   

* Small improvements to `group()` algorithm, avoiding some cases where the hash function performed badly, particularly with integers. 

* Function `GRPnames` now has a `sep` argument to choose a separator other than `"."`.


# collapse 1.7.6

* Corrected a C-level bug in `gsplit` that could lead R to crash in some instances (`gsplit` is used internally in `fsummarise`, `fmutate`, `BY` and `collap` to perform computations with base R (non-optimized) functions). 

* Ensured that `BY.grouped_df` always (by default) returns grouping columns in aggregations i.e. `iris |> gby(Species) |> nv() |> BY(sum)` now gives the same as `iris |> gby(Species) |> nv() |> fsum()`.

* A `.` was added to the first argument of functions `fselect`, `fsubset`, `colorder` and `fgroup_by`, i.e. `fselect(x, ...) -> fselect(.x, ...)`. The reason for this is that over time I added the option to select-rename columns e.g. `fselect(mtcars, cylinders = cyl)`, which was not offered when these functions were created. This presents problems if columns should be renamed into `x`, e.g. `fselect(mtcars, x = cyl)` failed, see [#221](https://github.com/SebKrantz/collapse/issues/221). Renaming the first argument to `.x` somewhat guards against such situations. I think this change is worthwhile to implement, because it makes the package more robust going forward, and usually the first argument of these functions is never invoked explicitly. I really hope this breaks nobody's code. 

* Added a function `GRPN` to make it easy to add a column of group sizes e.g. `mtcars %>% fgroup_by(cyl,vs,am) %>% ftransform(Sizes = GRPN(.))` or `mtcars %>% ftransform(Sizes = GRPN(list(cyl, vs, am)))` or `GRPN(mtcars, by = ~cyl+vs+am)`. 

* Added `[.pwcor` and `[.pwcov`, to be able to subset correlation/covariance matrices without loosing the print formatting. 

# collapse 1.7.5

* Also ensuring tidyverse examples are in `\donttest{}` and building without the *dplyr* testing file to avoid issues with static code analysis on CRAN.

* 20-50% Speed improvement in `gsplit` (and therefore in `fsummarise`, `fmutate`, `collap` and `BY` *when invoked with base R functions*) when grouping with `GRP(..., sort = TRUE, return.order = TRUE)`. To enable this by default, the default for argument `return.order` in `GRP` was set to `sort`, which retains the ordering vector (needed for the optimization). Retaining the ordering vector uses up some memory which can possibly adversely affect computations with big data, but with big data `sort = FALSE` usually gives faster results anyway, and you can also always set `return.order = FALSE` (also in `fgroup_by`, `collap`), so this default gives the best of both worlds. 

<!-- also considering that including more information in the grouping object can (and will) lead to further optimizations in the future.  -->

* An ancient depreciated argument `sort.row` (replaced by `sort` in 2020) is now removed from `collap`. Also arguments `return.order` and `method` were added to `collap` providing full control of the grouping that happens internally.

# collapse 1.7.4

* Tests needed to be adjusted for the upcoming release of *dplyr* 1.0.8 which involves an API change in `mutate`. `fmutate` will not take over these changes i.e. `fmutate(..., .keep = "none")` will continue to work like `dplyr::transmute`. Furthermore, no more tests involving *dplyr* are run on CRAN, and I will also not follow along with any future *dplyr* API changes.  

* The C-API macro `installTrChar` (used in the new `massign` function) was replaced with `installChar` to maintain backwards compatibility with R versions prior to 3.6.0. Thanks @tedmoorman #213. 

* Minor improvements to `group()`, providing increased performance for doubles and also increased performance when the second grouping variable is integer, which turned out to be very slow in some instances. 

# collapse 1.7.3

* Removed tests involving the *weights* package (which is not available on R-devel CRAN checks).

* `fgroup_by` is more flexible, supporting computing columns e.g. `fgroup_by(GGDC10S, Variable, Decade = floor(Year / 10) * 10)` and various programming options e.g. `fgroup_by(GGDC10S, 1:3)`, `fgroup_by(GGDC10S, c("Variable", "Country"))`, or `fgroup_by(GGDC10S, is.character)`. You can also use column sequences e.g. `fgroup_by(GGDC10S, Country:Variable, Year)`, but this should not be mixed with computing columns. Compute expressions may also not include the `:` function. 

* More memory efficient attribute handling in C/C++ (using C-API macro `SHALLOW_DUPLICATE_ATTRIB` instead of `DUPLICATE_ATTRIB`) in most places.

<!--
* *plm* methods for 'pseries' were adjusted to allow fetching of 'index' attribute from an external pointer, to adjust for changes in *plm* which allow computations on 'pseries' to be conducted in data masking environments. 
-->
# collapse 1.7.2

* Ensured that the base pipe `|>` is not used in tests or examples, to avoid errors on CRAN checks with older versions of R. 

* Also adjusted `psacf` / `pspacf` / `psccf` to take advantage of the faster grouping by `group`. 

# collapse 1.7.1

* Fixed minor C/C++ issues flagged in CRAN checks. 

* Added option `ties = "last"` to `fmode`. 

* Added argument `stable.algo` to `qsu`. Setting `stable.algo = FALSE` toggles a faster calculation of the standard deviation, yielding 2x speedup on large datasets. 

* *Fast Statistical Functions* now internally use `group` for grouping data if both `g` and `TRA` arguments are used, yielding efficiency gains on unsorted data. 

* Ensured that `fmutate` and `fsummarise` can be called if *collapse* is not attached. 


# collapse 1.7.0

*collapse* 1.7.0, released mid January 2022, brings major improvements in the computational backend of the package, its data manipulation capabilities, and a whole set of new functions that enable more flexible and memory efficient R programming - significantly enhancing the language itself. For the vast majority of codes, updating to 1.7 should not cause any problems. 



<!--
### Changes to functionality

* `ffirst`, `flast` and `fnobs` don't have hidden `*.list` methods anymore. 
-->

### Changes to functionality

* `num_vars` is now implemented in C, yielding a massive performance increase over checking columns using `vapply(x, is.numeric, logical(1))`. It selects columns where `(is.double(x) || is.integer(x)) && !is.object(x)`. This provides the same results for most common classes found in data frames (e.g. factors and date columns are not numeric), however it is possible for users to define methods for `is.numeric` for other objects, which will not be respected by `num_vars` anymore. A prominent example are base R's 'ts' objects i.e. `is.numeric(AirPassengers)` returns `TRUE`, but `is.object(AirPassengers)` is also `TRUE` so the above yields `FALSE`, implying - if you happened to work with data frames of 'ts' columns - that `num_vars` will now not select those anymore. Please make me aware if there are other important classes that are found in data frames and where `is.numeric` returns `TRUE`. `num_vars` is also used internally in `collap` so this might affect your aggregations. 

* In `flag`, `fdiff` and `fgrowth`, if a plain numeric vector is passed to the `t` argument such that `is.double(t) && !is.object(t)`, it is coerced to integer using `as.integer(t)` and directly used as time variable, rather than applying ordered grouping first. This is to avoid the inefficiency of grouping, and owes to the fact that in most data imported into R with various packages, the time (year) variables are coded as double although they should be integer (I also don't know of any cases where time needs to be indexed by a non-date variable with decimal places). Note that the algorithm internally handles irregularity in the time variable so this is not a problem. Should this break any code, kindly raise an issue on GitHub.

* The function `setrename` now truly renames objects by reference (without creating a shallow copy). The same is true for `vlabels<-` (which was rewritten in C) and a new function `setrelabel`. Thus additional care needs to be taken (with use inside functions etc.) as the renaming will take global effects unless a shallow copy of the data was created by some prior operation inside the function. If in doubt, better use `frename` which creates a shallow copy. 

* Some improvements to the `BY` function, both in terms of performance and security. Performance is enhanced through a new C function `gsplit`, providing split-apply-combine computing speeds competitive with *dplyr* on a much broader range of R objects. Regarding Security: if the result of the computation has the same length as the original data, names / rownames and grouping columns (for grouped data) are only added to the result object if known to be valid, i.e. if the data was originally sorted by the grouping columns (information recorded by `GRP.default(..., sort = TRUE)`, which is called internally on non-factor/GRP/qG objects). This is because `BY` does not reorder data after the split-apply-combine step (unlike `dplyr::mutate`); data are simply recombined in the order of the groups. Because of this, in general, `BY` should be used to compute summary statistics (unless data are sorted before grouping). The added security makes this explicit. 

* Added a method `length.GRP` giving the length of a grouping object. This could break code calling `length` on a grouping object before (which just returned the length of the list).

* Functions renamed in collapse 1.6.0 will now print a message telling you to use the updated names. The functions under the old names will stay around for 1-3 more years. 

<!--
* A new function `gsplit` is used internally in `collap` and `BY`, and the new `fmutate` function: To perform faster split-apply-combine computing with functions other than those provided in this package. `gsplit` does not support 'POSIXlt' and other complex classes requiring simultaneous splitting of multiple atomic vectors, but works well with all standard classes (including 'factor', 'Date', 'POSIXct') and beyond that all classes consisting of a vector (including lists) with attributes attached. 
-->

* The passing of argument `order` instead of `sort` in function `GRP` (from a very early version of collapse), is now disabled.


### Bug Fixes

* Fixed a bug in some functions using Welfords Online Algorithm (`fvar`, `fsd`, `fscale` and `qsu`) to calculate variances, occurring when initial or final zero weights caused the running sum of weights in the algorithm to be zero, yielding a division by zero and `NA` as output although a value was expected. These functions now skip zero weights alongside missing weights, which also implies that you can pass a logical vector to the weights argument to very efficiently calculate statistics on a subset of data (e.g. using `qsu`). 

### Additions

#### Basic Computational Infrastructure

* Function `group` was added, providing a low-level interface to a new unordered grouping algorithm based on hashing in C and optimized for R's data structures. The algorithm was heavily inspired by the great `kit` package of Morgan Jacob, and now feeds into the package through multiple central functions (including `GRP` / `fgroup_by`, `funique` and `qF`) when invoked with argument `sort = FALSE`. It is also used in internal groupings performed in data transformation functions such as `fwithin` (when no factor or 'GRP' object is provided to the `g` argument). The speed of the algorithm is very promising (often superior to `radixorder`), and it could be used in more places still. I welcome any feedback on its performance on different datasets.    

* Function `gsplit` provides an efficient alternative to `split` based on grouping objects. It is used as a new backend to `rsplit` (which also supports data frame) as well as `BY`, `collap`, `fsummarise` and `fmutate` - for more efficient grouped operations with functions external to the package. 

* Added multiple functions to facilitate memory efficient programming (written in C). These include elementary mathematical operations by reference (`setop`, `%+=%`, `%-=%`, `%*=%`, `%/=%`), supporting computations involving integers and doubles on vectors, matrices and data frames (including row-wise operations via `setop`) with no copies at all. Furthermore a set of functions which check a single value against a vector without generating logical vectors: `whichv`, `whichNA` (operators `%==%` and `%!=%` which return indices and are significantly faster than `==`, especially inside functions like `fsubset`), `anyv` and `allv` (`allNA` was already added before). Finally, functions `setv` and `copyv` speed up operations involving the replacement of a value (`x[x == 5] <- 6`) or of a sequence of values from a equally sized object (`x[x == 5] <- y[x == 5]`, or `x[ind] <- y[ind]` where `ind` could be pre-computed vectors or indices) in vectors and data frames without generating any logical vectors or materializing vector subsets.

* Function `vlengths` was added as a more efficient alternative to `lengths` (without method dispatch, simply coded in C).

* Function `massign` provides a multivariate version of `assign` (written in C, and supporting all basic vector types). In addition the operator `%=%` was added as an efficient multiple assignment operator. (It is called `%=%` and not `%<-%` to facilitate the translation of Matlab or Python codes into R, and because the [zeallot](<https://cran.r-project.org/package=zeallot>) package already provides multiple-assignment operators (`%<-%` and `%->%`), which are significantly more versatile, but orders of magnitude slower than `%=%`)

#### High-Level Features

* Fully fledged `fmutate` function that provides functionality analogous to `dplyr::mutate` (sequential evaluation of arguments, including arbitrary tagged expressions and `across` statements). `fmutate` is optimized to work together with the packages *Fast Statistical and Data Transformation Functions*, yielding fast, vectorized execution, but also benefits from `gsplit` for other operations. 

* `across()` function implemented for use inside `fsummarise` and `fmutate`. It is also optimized for *Fast Statistical and Data Transformation Functions*, but performs well with other functions too. It has an additional arguments `.apply = FALSE` which will apply functions to the entire subset of the data instead of individual columns, and thus allows for nesting tibbles and estimating models or correlation matrices by groups etc.. `across()` also supports an arbitrary number of additional arguments which are split and evaluated by groups if necessary. Multiple `across()` statements can be combined with tagged vector expressions in a single call to `fsummarise` or `fmutate`. Thus the computational framework is pretty general and similar to *data.table*, although less efficient with big datasets. 

* Added functions `relabel` and `setrelabel` to make interactive dealing with variable labels a bit easier. Note that both functions operate by reference. (Through `vlabels<-` which is implemented in C. Taking a shallow copy of the data frame is useless in this case because variable labels are attributes of the columns, not of the frame). The only difference between the two is that `setrelabel` returns the result invisibly.  

* function shortcuts `rnm` and `mtt` added for `frename` and `fmutate`. `across` can also be abbreviated using `acr`. 

* Added two options that can be invoked before loading of the package to change the namespace: `options(collapse_mask = c(...))` can be set to export copies of selected (or all) functions in the package that start with `f` removing the leading `f` e.g. `fsubset` -> `subset` (both `fsubset` and `subset` will be exported). This allows masking base R and dplyr functions (even basic functions such as `sum`, `mean`, `unique` etc. if desired) with *collapse*'s fast functions, facilitating the optimization of existing codes and allowing you to work with *collapse* using a more natural namespace. The package has been internally insulated against such changes, but of course they might have major effects on existing codes. Also `options(collapse_F_to_FALSE = FALSE)` can be invoked to get rid of the lead operator `F`, which masks `base::F` (an issue raised by some people who like to use `T`/`F` instead of `TRUE`/`FALSE`). Read the help page `?collapse-options` for more information.     

### Improvements

* Package loads faster (because I don't fetch functions from some other C/C++ heavy packages in `.onLoad` anymore, which implied unnecessary loading of a lot of DLLs).

* `fsummarise` is now also fully featured supporting evaluation of arbitrary expressions and `across()` statements. Note that mixing *Fast Statistical Functions* with other functions in a single expression can yield unintended outcomes, read more at `?fsummarise`. 

* `funique` benefits from `group` in the default `sort = FALSE`, configuration, providing extra speed and unique values in first-appearance order in both the default and the data frame method, for all data types. 

* Function `ss` supports both empty `i` or `j`. 

* The printout of `fgroup_by` also shows minimum and maximum group size for unbalanced groupings. 
* In `ftransformv/settransformv` and `fcomputev`, the `vars` argument is also evaluated inside the data frame environment, allowing NSE specifications using column names e.g. `ftransformv(data, c(col1, col2:coln), FUN)`.

* `qF` with option `sort = FALSE` now generates factors with levels in first-appearance order (instead of a random order assigned by the hash function), and can also be called on an existing factor to recast the levels in first-appearance order. It is also faster with `sort = FALSE` (thanks to `group`). 

* `finteraction` has argument `sort = FALSE` to also take advantage of `group`. 

* `rsplit` has improved performance through `gsplit`, and an additional argument `use.names`, which can be used to return an unnamed list. 

* Speedup in `vtypes` and functions `num_vars`, `cat_vars`, `char_vars`, `logi_vars` and `fact_vars`. Note than `num_vars` behaves slightly differently as discussed above.

* `vlabels(<-)` / `setLabels` rewritten in C, giving a ~20x speed improvement. Note that they now operate by reference. 

* `vlabels`, `vclasses` and `vtypes` have a `use.names` argument. The default is `TRUE` (as before). 

* `colorder` can rename columns on the fly and also has a new mode `pos = "after"` to place all selected  columns after the first selected one, e.g.: `colorder(mtcars, cyl, vs_new = vs, am, pos = "after")`. The `pos = "after"` option was also added to `roworderv`. 

+ `add_stub` and `rm_stub` have an additional `cols` argument to apply a stub to certain columns only e.g. `add_stub(mtcars, "new_", cols = 6:9)`.

* `namlab` has additional arguments `N` and `Ndistinct`, allowing to display number of observations and distinct values next to variable names, labels and classes, to get a nice and quick overview of the variables in a large dataset. 

* `copyMostAttrib` only copies the `"row.names"` attribute when known to be valid. 

* `na_rm` can now be used to efficiently remove empty or `NULL` elements from a list. 

* `flag`, `fdiff` and `fgrowth` produce less messages (i.e. no message if you don't use a time variable in grouped operations, and messages about computations on highly irregular panel data only if data length exceeds 10 million obs.).

* The print methods of `pwcor` and `pwcov` now have a `return` argument, allowing users to obtain the formatted correlation matrix, for exporting purposes. 

* `replace_NA`, `recode_num` and `recode_char` have improved performance and an additional argument `set` to take advantage of `setv` to change (some) data by reference. For `replace_NA`, this feature is mature and setting `set = TRUE` will modify all selected columns in place and return the data invisibly. For `recode_num` and `recode_char` only a part of the transformations are done by reference, thus users will still have to assign the data to preserve changes. In the future, this will be improved so that `set = TRUE` toggles all transformations to be done by reference.   

<!-- 

* `TRA`, `varying`, `seqid` and `groupid` rewritten in C. Mainly to reduce the size of compiled code, allowing me to add functions while keeping the package size constant.  

* `TRA` has extra option `"replace_NA"` to only replace missing values with computed statistics. This is useful for example to do extremely fast imputation using `fmean(..., TRA = "replace_NA")` or `fmedian`. 

* C API (Dirk)
-->


# collapse 1.6.5

* Use of `VECTOR_PTR` in C API now gives an error on R-devel even if `USE_RINTERNALS` is defined. Thus this patch gets rid of all remaining usage of this macro to avoid errors on CRAN checks using the development version of R. 

* The print method for `qsu` now uses an apostrophe (') to designate million digits, instead of a comma (,). This is to avoid confusion with the decimal point, and the typical use of (,) for thousands (which I don't like). 

# collapse 1.6.4
Checks on the gcc11 compiler flagged an additional issue with a pointer pointing to element -1 of a C array (which I had done on purpose to index it with an R integer vector). 

# collapse 1.6.3
CRAN checks flagged a valgrind issue because of comparing an uninitialized value to something. 

# collapse 1.6.2
CRAN maintainers have asked me to remove a line in a Makevars file intended to reduce the size of Rcpp object files (which has been there since version 1.4). So the installed size of the package may now be larger.  

# collapse 1.6.1
A patch for 1.6.0 which fixes issues flagged by CRAN and adds a few handy extras. 

### Bug Fixes

* Puts examples using the new base pipe `|>` inside `\donttest{}` so that they don't fail CRAN tests on older R versions.

* Fixes a LTO issue caused by a small mistake in a header file (which does not have any implications to the user but was detected by CRAN checks).

### Additions

* Added a function `fcomputev`, which allows selecting columns and transforming them with a function in one go. The `keep` argument can be used to add columns to the selection that are not transformed. 

* Added a function `setLabels` as a wrapper around `vlabels<-` to facilitate setting variable labels inside pipes. 

* Function `rm_stub` now has an argument `regex = TRUE` which triggers a call to `gsub` and allows general removing of character sequences in column names on the fly. 

### Improvements

* `vlabels<-` and  `setLabels` now support list of variable labels or other attributes (i.e. the `value` is internally subset using `[[`, not `[`). Thus they are now general functions to attach a vector or list of attributes to columns in a list / data frame. 



# collapse 1.6.0
*collapse* 1.6.0, released end of June 2021, presents some significant improvements in the user-friendliness, compatibility and programmability of the package, as well as a few function additions. 

### Changes to Functionality

* `ffirst`, `flast`, `fnobs`, `fsum`, `fmin` and `fmax` were rewritten in C. The former three now also support list columns (where `NULL` or empty list elements are considered missing values when `na.rm = TRUE`), and are extremely fast for grouped aggregation if `na.rm = FALSE`. The latter three also support and return integers, with significant performance gains, even compared to base R. Code using these functions expecting an error for list-columns or expecting double output even if the input is integer should be adjusted. 

* *collapse* now directly supports *sf* data frames through functions like `fselect`, `fsubset`, `num_vars`, `qsu`, `descr`, `varying`, `funique`, `roworder`, `rsplit`, `fcompute` etc., which will take along the geometry column even if it is not explicitly selected (mirroring *dplyr* methods for *sf* data frames). This is mostly done internally at C-level, so functions remain simple and fast. Existing code that explicitly selects the geometry column is unaffected by the change, but code of the form `sf_data %>% num_vars %>% qDF %>% ...`, where columns excluding geometry were selected and the object later converted to a data frame, needs to be rewritten as `sf_data %>% qDF %>% num_vars %>% ...`. A short vignette was added describing the integration of *collapse* and *sf*. 

* I've received several requests for increased namespace consistency. *collapse* functions were named to be consistent with base R, *dplyr* and *data.table*, resulting in names like `is.Date`, `fgroup_by` or `settransformv`. To me this makes sense, but I've been convinced that a bit more consistency is advantageous. Towards that end I have decided to eliminate the '.' notation of base R and to remove some unexpected capitalizations in function names giving some people the impression I was using camel-case. The following functions are renamed:
`fNobs` -> `fnobs`, `fNdistinct` -> `fndistinct`, `pwNobs` -> `pwnobs`, `fHDwithin` -> `fhdwithin`,
`fHDbetween` -> `fhdbetween`, `as.factor_GRP` -> `as_factor_GRP`, `as.factor_qG` -> `as_factor_qG`, `is.GRP` -> `is_GRP`, `is.qG` -> `is_qG`, `is.unlistable` -> `is_unlistable`, `is.categorical` -> `is_categorical`, `is.Date` -> `is_date`, `as.numeric_factor` -> `as_numeric_factor`, `as.character_factor` -> `as_character_factor`, 
`Date_vars` -> `date_vars`. 
This is done in a very careful manner, the others will stick around for a long while (end of 2022), and the generics of `fNobs`, `fNdistinct`, `fHDbetween` and `fHDwithin` will be kept in the package for an indeterminate period, but their core methods will not be exported beyond 2022. I will start warning about these renamed functions in 2022. In the future I will undogmatically stick to a function naming style with lowercase function names and underslashes where words need to be split. Other function names will be kept. To say something about this: The quick-conversion functions `qDF` `qDT`, `qM`, `qF`, `qG` are consistent and in-line with *data.table* (`setDT` etc.), and similarly the operators `L`, `F`, `D`, `Dlog`, `G`, `B`, `W`, `HDB`, `HDW`. I'll keep `GRP`, `BY` and `TRA`, for lack of better names, parsimony and because they are central to the package. The camel case will be kept in helper functions `setDimnames` etc. because they work like *stats* `setNames` and do not modify the argument by reference (like `settransform` or `setrename` and various *data.table* functions). Functions `copyAttrib` and `copyMostAttrib` are exports of like-named functions in the C API and thus kept as they are. Finally, I want to keep `fFtest` the way it is because the F-distribution is widely recognized by a capital F. 

* I've updated the `wlddev` dataset with the latest data from the World Bank, and also added a variable giving the total population (which may be useful e.g. for population-weighted aggregations across regions). The extra column could invalidate codes used to demonstrate something (I had to adjust some examples, tests and code in vignettes).

### Additions

* Added a function `fcumsum` (written in C), permitting flexible (grouped, ordered) cumulative summations on matrix-like objects (integer or double typed) with extra methods for grouped data frames and panel series and data frames. Apart from the internal grouping, and an ordering argument allowing cumulative sums in a different order than data appear, `fcumsum` has 2 options to deal with missing values. The default (`na.rm = TRUE`) is to skip (preserve) missing values, whereas setting `fill = TRUE` allows missing values to be populated with the previous value of the cumulative sum (starting from 0). 

* Added a function `alloc` to efficiently generate vectors initialized with any value (faster than `rep_len`). 

* Added a function `pad` to efficiently pad vectors / matrices / data.frames with a value (default is `NA`). This function was mainly created to make it easy to expand results coming from a statistical model fitted on data with missing values to the original length. For example let `data <- na_insert(mtcars); mod <- lm(mpg ~ cyl, data)`, then we can do `settransform(data, resid = pad(resid(mod), mod$na.action))`, or we could do `pad(model.matrix(mod), mod$na.action)` or `pad(model.frame(mod), mod$na.action)` to receive matrices and data frames from model data matching the rows of `data`. `pad` is a general function that will also work with mixed-type data. It is also possible to pass a vector of indices matching the rows of the data to `pad`, in which case `pad` will fill gaps in those indices with a value/row in the data.  


### Improvements

* Full *data.table* support, including reference semantics (`set*`, `:=`)!! There is some complex C-level programming behind *data.table*'s operations by reference. Notably, additional (hidden) column pointers are allocated to be able to add columns without taking a shallow copy of the *data.table*, and an `".internal.selfref"` attribute containing an external pointer is used to check if any shallow copy was made using base R commands like `<-`. This is done to avoid even a shallow copy of the *data.table* in manipulations using `:=` (and is in my opinion not worth it as even large tables are shallow copied by base R (>=3.1.0) within microseconds and all of this complicates development immensely). Previously, *collapse* treated *data.table*'s like any other data frame, using shallow copies in manipulations and preserving the attributes (thus ignoring how *data.table* works internally). This produced a warning whenever you wanted to use *data.table* reference semantics (`set*`, `:=`) after passing the *data.table* through a *collapse* function such as `collap`, `fselect`, `fsubset`, `fgroup_by` etc. From v1.6.0, I have adopted essential C code from *data.table* to do the overallocation and generate the `".internal.selfref"` attribute, thus seamless workflows combining *collapse* and *data.table* are now possible. This comes at a cost of about 2-3 microseconds per function, as to do this I have to shallow copy the *data.table* again and add extra column pointers and an `".internal.selfref"` attribute telling *data.table* that this table was not copied (it seems to be the only way to do it for now). This integration encompasses all data manipulation functions in *collapse*, but not the *Fast Statistical Functions* themselves. Thus you can do `agDT <- DT %>% fselect(id, col1:coln) %>% collap(~id, fsum); agDT[, newcol := 1]`, but you would need to do add a `qDT` after a function like `fsum` if you want to use reference semantics without incurring a warning: `agDT <- DT %>% fselect(id, col1:coln) %>% fgroup_by(id) %>% fsum %>% qDT; agDT[, newcol := 1]`. *collapse* appears to be the first package that attempts to account for *data.table*'s internal working without importing *data.table*, and `qDT` is now the fastest way to create a fully functional *data.table* from any R object. A global option `"collapse_DT_alloccol"` was added to regulate how many columns *collapse* overallocates when creating *data.table*'s. The default is 100, which is lower than the *data.table* default of 1024. This was done to increase efficiency of the additional shallow copies, and may be changed by the user. 

* Programming enabled with `fselect` and `fgroup_by` (you can now pass vectors containing column names or indices). Note that instead of `fselect` you should use `get_vars` for standard eval programming.  

* `fselect` and `fsubset` support in-place renaming, e.g. `fselect(data, newname = var1, var3:varN)`,
`fsubset(data, vark > varp, newname = var1, var3:varN)`.

* `collap` supports renaming columns in the custom argument, e.g. `collap(data, ~ id, custom = list(fmean = c(newname = "var1", "var2"), fmode = c(newname = 3), flast = is_date))`. 

* Performance improvements: `fsubset` / `ss` return the data or perform a simple column subset without deep copying the data if all rows are selected through a logical expression. `fselect` and `get_vars`, `num_vars` etc. are slightly faster through data frame subsetting done fully in C. `ftransform` / `fcompute` use `alloc` instead of `base::rep` to replicate a scalar value which is slightly more efficient. 

* `fcompute` now has a `keep` argument, to preserve several existing columns when computing columns on a data frame. 

* `replace_NA` now has a `cols` argument, so we can do `replace_NA(data, cols = is.numeric)`, to replace `NA`'s in numeric columns. I note that for big numeric data `data.table::setnafill` is the most efficient solution. 

* `fhdbetween` and `fhdwithin` have an `effect` argument in *plm* methods, allowing centering on selected identifiers. The default is still to center on all panel identifiers. 

 <!-- * `flag`, `fdiff`, `fgrowth` methods for *tsibble* and *tibbletime* -->

* The plot method for panel series matrices and arrays `plot.psmat` was improved slightly. It now supports custom colours and drawing of a grid. 

* `settransform` and `settransformv` can now be called without attaching the package e.g. `collapse::settransform(data, ...)`. These errored before when *collapse* is not loaded because they are simply wrappers around `data <- ftransform(data, ...)`. I'd like to note from a [discussion](https://github.com/SebKrantz/collapse/issues/136) that avoiding shallow copies with `<-` (e.g. via `:=`) does not appear to yield noticeable performance gains. Where *data.table* is faster on big data this mostly has to do with parallelism and sometimes with algorithms, generally not memory efficiency.

* Functions `setAttrib`, `copyAttrib` and `copyMostAttrib` only make a shallow copy of lists, not of atomic vectors (which amounts to doing a full copy and is inefficient). Thus atomic objects are now modified in-place. 

* Small improvements: Calling `qF(x, ordered = FALSE)` on an ordered factor will remove the ordered class, the operators `L`, `F`, `D`, `Dlog`, `G`, `B`, `W`, `HDB`, `HDW` and functions like `pwcor` now work on unnamed matrices or data frames. 


# collapse 1.5.3

* A test that occasionally fails on Mac is removed, and all unit testing is now removed from CRAN. *collapse* has close to 10,000 unit tests covering all central pieces of code. Half of these tests depend on generated data, and for some reasons there is always a test or two that occasionally fail on some operating system (usually not Windows), requiring me to submit a patch. This is not constructive to either the development or the use of this package, therefore tests are now removed from CRAN. They are still run on codecov.io, and every new release is thoroughly tested on Windows. 

# collapse 1.5.2

### Changes to Functionality

* The first argument of `ftransform` was renamed to `.data` from `X`. This was done to enable the user to transform columns named "X". For the same reason the first argument of `frename` was renamed to `.x` from `x` (not `.data` to make it explicit that `.x` can be any R object with a "names" attribute). It is not possible to depreciate `X` and `x` without at the same time undoing the benefits of the argument renaming, thus this change is immediate and code breaking in rare cases where the first argument is explicitly set. 

* The function `is.regular` to check whether an R object is atomic or list-like is depreciated and will be removed before the end of the year. This was done to avoid a namespace clash with the *zoo* package (#127).

### Bug Fixes

<!-- 
* For reasons of efficiency, most statistical and transformation functions used the C macro `SHALLOW_DUPLICATE_ATTRIB` to copy column attributes in a data frame. Since this macro does not copy S4 object bits, it caused some problems with S4 object columns such as POSIXct (e.g. computing lags/leads, first and last values on these columns). This is now fixed, all statistical functions (apart from `fvar` and `fsd`) now use `DUPLICATE_ATTRIB` and thus preserve S4 object columns (#91). 
-->
<!-- Also `BY` now handles POSIXct properly. -->

* `unlist2d` produced a subsetting error if an empty list was present in the list-tree. This is now fixed, empty or `NULL` elements in the list-tree are simply ignored (#99).

### Additions

* A function `fsummarize` was added to facilitate translating *dplyr* / *data.table* code to *collapse*. Like `collap`, it is only very fast when used with the *Fast Statistical Functions*. 

* A function `t_list` is made available to efficiently transpose lists of lists. 

<!-- 
* A small set of row-wise statistical functions: `rowsums`, `rowmeans`, `rowNobs`, `rowmins` and `rowmaxs` is introduced for fast (grouped, weighted) row-wise computations on matrices and data frames. -->

### Improvements
<!-- 
* `ffirst` and `flast` were rewritten in C, with slightly better performance and reduced file size. -->

* C files are compiled -O3 on Windows, which gives a boost of around 20% for the grouping mechanism applied to character data.



# collapse 1.5.1
A small patch for 1.5.0 that:

* Fixes a numeric precision issue when grouping doubles (e.g. before `qF(wlddev$LIFEEX)` gave an error, now it works). 

* Fixes a minor issue with `fhdwithin` when applied to *pseries* and `fill = FALSE`.


# collapse 1.5.0
*collapse* 1.5.0, released early January 2021, presents important refinements and some additional functionality. 

### Back to CRAN

* I apologize for inconveniences caused by the temporal archival of *collapse* from December 19, 2020. This archival was caused by the archival of the important *lfe* package on the 4th of December. *collapse* depended on *lfe* for higher-dimensional centering, providing the `fhdbetween / fhdwithin` functions for generalized linear projecting / partialling out. To remedy the damage caused by the removal of *lfe*, I had to rewrite `fhdbetween / fhdwithin` to take advantage of the demeaning algorithm provided by *fixest*, which has some quite different mechanics. Beforehand, I made some significant changes to `fixest::demean` itself to make this integration happen. The CRAN deadline was the 18th of December, and I realized too late that I would not make this. A request to CRAN for extension was declined, so *collapse* got archived on the 19th. I have learned from this experience, and *collapse* is now sufficiently insulated that it will not be taken off CRAN even if all suggested packages were removed from CRAN. 

### Bug Fixes

* Segfaults in several *Fast Statistical Functions* when passed `numeric(0)` are fixed (thanks to @eshom and @acylam, [#101](https://github.com/SebKrantz/collapse/issues/101)). The default behavior is that all *collapse* functions return `numeric(0)` again, except for `fnobs`, `fndistinct` which return `0L`, and `fvar`, `fsd` which return `NA_real_`.

### Changes to Functionality

* Functions `fhdwithin / HDW` and `fhdbetween / HDB` have been reworked, delivering higher performance and greater functionality: For higher-dimensional centering and heterogeneous slopes, the `demean` function from the *fixest* package is imported (conditional on the availability of that package). The linear prediction  and partialling out functionality is now built around `flm` and also allows for weights and different fitting methods. 

* In `collap`, the default behavior of `give.names = "auto"` was altered when used together with the `custom` argument. Before the function name was always added to the column names. Now it is only added if a column is aggregated with two different functions. I apologize if this breaks any code dependent on the new names, but this behavior just better reflects most common use (applying only one function per column), as well as STATA's collapse. 

* For list processing functions like `get_elem`, `has_elem` etc. the default for the argument `DF.as.list` was changed from `TRUE` to `FALSE`. This means if a nested lists contains data frame's, these data frame's will not be searched for matching elements. This default also reflects the more common usage of these functions (extracting entire data frame's or computed quantities from nested lists rather than searching / subsetting lists of data frame's). The change also delivers a considerable performance gain. 


<!-- Missing value removal is still done using *data.table* source code, so these functions are now equipped for large and complex linear prediction and partialling-out problems. To fully utilize them users must install *fixest*. -->

* Vignettes were outsourced to the [website](<https://sebkrantz.github.io/collapse/articles/index.html>). This nearly halves the size of the source package, and should induce users to appreciate the built-in documentation. The website also makes for much more convenient reading and navigation of these book-style vignettes. 

<!-- , and also made available as PDF versions for download there -->

### Additions

* Added a set of 10 operators `%rr%`, `%r+%`, `%r-%`, `%r*%`, `%r/%`, `%cr%`, `%c+%`, `%c-%`, `%c*%`, `%c/%` to facilitate and speed up row- and column-wise arithmetic operations involving a vector and a matrix / data frame / list. For example `X %r*% v` efficiently multiplies every row of `X` with `v`. Note that more advanced functionality is already provided in `TRA()`, `dapply()` and the *Fast Statistical Functions*, but these operators are intuitive and very convenient to use in matrix or matrix-style code, or in piped expressions. 

* Added function `missing_cases` (opposite of `complete.cases` and faster for data frame's / lists).

* Added function `allNA` for atomic vectors. 

* New vignette about using *collapse* together with *data.table*, available [online](<https://sebkrantz.github.io/collapse/articles/index.html>).

### Improvements

* Time series functions and operators `flag / L / F`, `fdiff / D / Dlog` and `fgrowth / G` now natively support irregular time series and panels, and feature a 'complete approach' i.e. values are shifted around taking full account of the underlying time-dimension!

<!-- Thus *collapse* now provides a comprehensive, robust and extremely fast approach to working with time-dependent data in R. -->

* Functions `pwcor` and `pwcov` can now compute weighted correlations on the pairwise or complete observations, supported by C-code that is (conditionally) imported from the *weights* package. 

* `fFtest` now also supports weights.

* `collap` now provides an easy workaround to aggregate some columns using weights and others without. The user may simply append the names of *Fast Statistical Functions* with `_uw` to disable weights. Example: `collapse::collap(mtcars, ~ cyl, custom = list(fmean_uw = 3:4, fmean = 8:10), w = ~ wt)` aggregates columns 3 through 4 using a simple mean and columns 8 through 10 using the weighted mean. 

* The parallelism in `collap` using `parallel::mclapply` has been reworked to operate at the column-level, and not at the function level as before. It is still not available for Windows though. The default number of cores was set to `mc.cores = 2L`, which now gives an error on windows if `parallel = TRUE`.  

* function `recode_char` now has additional options `ignore.case` and `fixed` (passed to `grepl`), for enhanced recoding character data based on regular expressions. 

* `rapply2d` now has `classes` argument permitting more flexible use. 

* `na_rm` and some other internal functions were rewritten in C. `na_rm` is now 2x faster than `x[!is.na(x)]` with missing values and 10x faster without missing values. 


# collapse 1.4.2

* An improvement to the `[.GRP_df` method enabling the use of most *data.table* methods (such as `:=`) on a grouped *data.table* created with `fgroup_by`.

* Some documentation updates by Kevin Tappe.

# collapse 1.4.1
collapse 1.4.1 is a small patch for 1.4.0 that:

* fixes clang-UBSAN and rchk issues in 1.4.0 (minor bugs in compiled code resulting, in this case, from trying to coerce a `NaN` value to integer, and failing to protect a shallow copy of a variable).

* Adds a method `[.GRP_df` that allows robust subsetting of grouped objects created with `fgroup_by` (thanks to Patrice Kiener for flagging this).

# collapse 1.4.0
*collapse* 1.4.0, released early November 2020, presents some important refinements, particularly in the domain of attribute handling, as well as some additional functionality. The changes make *collapse* smarter, more broadly compatible and more secure, and should not break existing code.  <!-- , is a major update: -->

### Changes to Functionality

* *Deep Matrix Dispatch / Extended Time Series Support:* The default methods of all statistical and transformation functions dispatch to the matrix method if `is.matrix(x) && !inherits(x, "matrix")` evaluates to `TRUE`. This specification avoids invoking the default method on classed matrix-based objects (such as multivariate time series of the *xts* / *zoo* class) not inheriting a 'matrix' class, while still allowing the user to manually call the default method on matrices (objects with implicit or explicit 'matrix' class). The change implies that *collapse*'s generic statistical functions are now well suited to transform *xts* / *zoo* and many other time series and matrix-based classes. 

* *Fully Non-Destructive Piped Workflow:* `fgroup_by(x, ...)` now only adds a class *grouped_df*, not classes *table_df*, *tbl*, *grouped_df*, and preserves all classes of `x`. This implies that workflows such as `x %>% fgroup_by(...) %>% fmean` etc. yields an object `xAG` of the same class and attributes as `x`, not a tibble as before. *collapse* aims to be as broadly compatible, class-agnostic and attribute preserving as possible. 

<!-- Not a priority for now! Not really necessary at all, can always use base R converters! -->
* *Thorough and Controlled Object Conversions:* Quick conversion functions `qDF`, `qDT` and `qM` now have additional arguments `keep.attr` and `class` providing precise user control over object conversions in terms of classes and other attributes assigned / maintained. The default (`keep.attr = FALSE`) yields *hard* conversions removing all but essential attributes from the object. E.g. before `qM(EuStockMarkets)` would just have returned `EuStockMarkets` (because `is.matrix(EuStockMarkets)` is `TRUE`) whereas now the time series class and 'tsp' attribute are removed. `qM(EuStockMarkets, keep.attr = TRUE)` returns `EuStockMarkets` as before. 

<!--
In general `keep.attr = TRUE` gives a *soft* conversion were attributes necessary to establish the new data type ('dim', 'dimnames', 'names', 'row.names', 'class') are modified, but all other attributes are kept. 
This may be useful in some cases, for examples it is now possible to write something like `mtcars %>% fgroup_by(cyl, vs, am) %>% qM(TRUE) %>% fmean(attr(.,"groups"))`. 
 
 (ensuring that `qDF`, `qDT` and `qM` now truly behave like `as.data.frame`, `as.data.table` and `as.matrix`) -->

* *Smarter Attribute Handling:* Drawing on the guidance given in the R Internals manual, the following standards for optimal non-destructive attribute handling are formalized and communicated to the user: 

  + The default and matrix methods of the *Fast Statistical Functions* preserve attributes of the input in grouped aggregations ('names', 'dim' and 'dimnames' are suitably modified). If inputs are classed objects (e.g. factors, time series, checked by `is.object`), the class and other attributes are dropped. Simple (non-grouped) aggregations of vectors and matrices do not preserve attributes, unless `drop = FALSE` in the matrix method. An exemption is made in the default methods of functions `ffirst`, `flast` and `fmode`, which always preserve the attributes (as the input could well be a factor or date variable). 
  
  + The data frame methods are unaltered: All attributes of the data frame and columns in the data frame are preserved unless the computation result from each column is a scalar (not computing by groups) and `drop = TRUE` (the default). 
  
  + Transformations with functions like `flag`, `fwithin`, `fscale` etc. are also unaltered: All attributes of the input are preserved in the output (regardless of whether the input is a vector, matrix, data.frame or related classed object). The same holds for transformation options modifying the input ("-", "-+", "/", "+", "\*", "%%", "-%%") when using `TRA()` function or the `TRA = "..."` argument to the *Fast Statistical Functions*. 
  
  + For `TRA` 'replace' and 'replace_fill' options, the data type of the STATS is preserved, not of x. This provides better results particularly with functions like `fnobs` and `fndistinct`. E.g. previously `fnobs(letters, TRA = "replace")` would have returned the observation counts coerced to character, because `letters` is character. Now the result is integer typed. For attribute handling this means that the attributes of x are preserved unless x is a classed object and the data types of x and STATS do not match. An exemption to this rule is made if x is a factor and an integer (non-factor) replacement is offered to STATS. In that case the attributes of x are copied exempting the 'class' and 'levels' attribute, e.g. so that `fnobs(iris$Species, TRA = "replace")` gives an integer vector, not a (malformed) factor. In the unlikely event that STATS is a classed object, the attributes of STATS are preserved and the attributes of x discarded. 
  
<!-- This simple but thorough (and now formalized) system of attribute handling should be optimal in more than 90% of common applications. -->

  
 <!--  
  + The default methods of statistical functions returning numeric values (`fmean`, `fmedian`, `fsum`, `fprod`, `fvar`, `fsd`, `fmin`, `fmax` and `fnth`, `fnobs` and `fndistinct`) do not preserve the attributes of classed objects (e.g. univariate time series) in grouped aggregations. This is consistent with the matrix methods of these functions and avoids errors, particularly after computations on time series. -->
 
 <!--, random conversions to tibble are not part of its philosophy. -->

<!--
* All S3 generic functions with a `default` method for atomic vectors and a `matrix` method now have an additional internal dispatch from the `default` to the `matrix` method if a classed matrix object missing a 'matrix' class is passed to the generic. For example consider a matrix time series `x <- structure(matrix(1:9, ncol = 3), class = "ts", tsp = c(1, 3, 1))` inheriting only a 'ts' but not a 'matrix' class. In collapse 1.3.2 `fsum(x)` would invoke the default method and return a scalar value. Now `fsum(x)` returns the sum for each column in the matrix. The `matrix` method is only called from the `default` method if `is.matrix(x) && !inherits(x, "matrix")` evaluates to `TRUE`, thus it is still possible to manually invoke the default method on a matrix. As the example indicates, this change is warranted to improve the inherent compatibility of *collapse* with various time series and matrix based classes (such as *xts* / *zoo*). -->


* *Reduced Dependency Burden:* The dependency on the *lfe* package was made optional. Functions `fhdwithin` / `fhdbetween` can only perform higher-dimensional centering if *lfe* is available. Linear prediction and centering with a single factor (among a list of covariates) is still possible without installing *lfe*. This change means that *collapse* now only depends on base R and *Rcpp* and is supported down to R version 2.10. 

### Additions

* Added function `rsplit` for efficient (recursive) splitting of vectors and data frames. 

* Added function `fdroplevels` for very fast missing level removal + added argument `drop` to `qF` and `GRP.factor`, the default is `drop = FALSE`. The addition of `fdroplevels` also enhances the speed of the `fFtest` function.

* `fgrowth` supports annualizing / compounding growth rates through added `power` argument.

* A function `flm` was added for bare bones (weighted) linear regression fitting using different efficient methods: 4 from base R (`.lm.fit`, `solve`, `qr`, `chol`), using `fastLm` from *RcppArmadillo* (if installed), or `fastLm` from *RcppEigen* (if installed). 

* Added function `qTBL` to quickly convert R objects to tibble.

* helpers `setAttrib`, `copyAttrib` and `copyMostAttrib` exported for fast attribute handling in R (similar to `attributes<-()`, these functions return a shallow copy of the first argument with the set of attributes replaced, but do not perform checks for attribute validity like `attributes<-()`. This can yield large performance gains with big objects).  

* helper `cinv` added wrapping the expression `chol2inv(chol(x))` (efficient inverse of a symmetric, positive definite matrix via Choleski factorization). 

* A shortcut `gby` is now available to abbreviate the frequently used `fgroup_by` function. 

* A print method for grouped data frames of any class was added.

### Improvements

* Faster internal methods for factors for `funique`, `fmode` and `fndistinct`.

<!-- * `flag`, `fdiff`, `fgrowth` support *xts* / *zoo* via explicit methods for fast and secure computations on unordered data. --> 

* The *grouped_df* methods for `flag`, `fdiff`, `fgrowth` now also support multiple time variables to identify a panel e.g. `data %>% fgroup_by(region, person_id) %>% flag(1:2, list(month, day))`.

* More security features for `fsubset.data.frame` / `ss`, `ss` is now internal generic and also supports subsetting matrices. 

* In some functions (like `na_omit`), passing double values (e.g. `1` instead of integer `1L`) or negative indices to the `cols` argument produced an error or unexpected behavior. This is now fixed in all functions. 

* Fixed a bug in helper function `all_obj_equal` occurring if objects are not all equal. 

* Some performance improvements through increased use of pointers and C API functions.



# collapse 1.3.2
collapse 1.3.2, released mid September 2020: <!-- , is a minor update: -->

* Fixed a small bug in `fndistinct` for grouped distinct value counts on logical vectors. 

* Additional security for `ftransform`, which now efficiently checks the names of the data and replacement arguments for uniqueness, and also allows computing and transforming list-columns.  

* Added function `ftransformv` to facilitate transforming selected columns with function - a very efficient replacement for `dplyr::mutate_if` and `dplyr::mutate_at`. 

* `frename` now allows additional arguments to be passed to a renaming function.  

# collapse 1.3.1
collapse 1.3.1, released end of August 2020, is a patch for v1.3.0 that takes care of some unit test failures on certain operating systems (mostly because of numeric precision issues). It provides no changes to the code or functionality.

# collapse 1.3.0
collapse 1.3.0, released mid August 2020: <!-- , is another major update: -->

### Changes to Functionality

* `dapply` and `BY` now drop all unnecessary attributes if `return = "matrix"` or `return = "data.frame"` are explicitly requested (the default `return = "same"` still seeks to preserve the input data structure).

* `unlist2d` now saves integer rownames if `row.names = TRUE` and a list of matrices without rownames is passed, and `id.factor = TRUE` generates a normal factor not an ordered factor. It is however possible to write `id.factor = "ordered"` to get an ordered factor id.  

* `fdiff` argument `logdiff` renamed to `log`, and taking logs is now done in R (reduces size of C++ code and does not generate as many NaN's). `logdiff` may still be used, but it may be deactivated in the future. Also in the matrix and data.frame methods for `flag`, `fdiff` and `fgrowth`, columns are only stub-renamed if more than one lag/difference/growth rate is computed. 

### Additions

* Added `fnth` for fast (grouped, weighted) n'th element/quantile computations.

* Added `roworder(v)` and `colorder(v)` for fast row and column reordering.  

* Added `frename` and `setrename` for fast and flexible renaming (by reference).  

* Added function `fungroup`, as replacement for `dplyr::ungroup`, intended for use with `fgroup_by`. 

* `fmedian` now supports weights, computing a decently fast (grouped) weighted median based on radix ordering. 

* `fmode` now has the option to compute min and max mode, the default is still simply the first mode. 

* `fwithin` now supports quasi-demeaning (added argument `theta`) and can thus be used to manually estimate random-effects models. 

* `funique` is now generic with a default vector and data.frame method, providing fast unique values and rows of data. The default was changed to `sort = FALSE`.   

* The shortcut `gvr` was created for `get_vars(..., regex = TRUE)`. 

* A helper `.c` was introduced for non-standard concatenation (i.e. `.c(a, b) == c("a", "b")`). 

### Improvements

* `fmode` and `fndistinct` have become a bit faster.

* `fgroup_by` now preserves *data.table*'s.

* `ftransform` now also supports a data.frame as replacement argument, which automatically replaces matching columns and adds unmatched ones. Also `ftransform<-` was created as a more formal replacement method for this feature.

* `collap` columns selected through `cols` argument are returned in the order selected if `keep.col.order = FALSE`. Argument `sort.row` is depreciated, and replace by argument `sort`. In addition the `decreasing` and `na.last` arguments were added and handed down to `GRP.default`. 

* `radixorder` 'sorted' attribute is now always attached.

* `stats::D` which is masked when collapse is attached, is now preserved through methods `D.expression` and `D.call`. 

* `GRP` option `call = FALSE` to omit a call to `match.call` -> minor performance improvement.

* Several small performance improvements through rewriting some internal helper functions in C and reworking some R code. 

* Performance improvements for some helper functions, `setRownames` / `setColnames`, `na_insert` etc.  

* Increased scope of testing statistical functions. The functionality of the package is now secured by 7700 unit tests covering all central bits and pieces. 


# collapse 1.2.1
collapse 1.2.1, released end of May 2020: <!-- , is a patch for v1.2.0: -->

* Minor fixes for 1.2.0 issues that prevented correct installation on Mac OS X and a vignette rebuilding error on solaris.

* `fmode.grouped_df` with groups and weights now saves the sum of the weights instead of the max (this makes more sense as the max only applies if all elements are unique). 

# collapse 1.2.0
collapse 1.2.0, released mid May 2020: <!-- , is a major update of the package - changes and additions: -->

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

* `fdiff` now supports quasi-differences i.e. $x_t - \rho x_{t-1}$ and quasi-log differences i.e. $log(x_t) - \rho log(x_{t-1})$. an arbitrary $\rho$ can be supplied.

* Added a `Dlog` operator for faster access to log-differences. 

### Improvements
* Faster grouping with `GRP` and faster factor generation with added radix method + automatic dispatch between hash and radix method. `qF` is now ~ 5x faster than `as.factor` on character and around 30x faster on numeric data. Also `qG` was enhanced. 

* Further slight speed tweaks here and there. 

* `collap` now provides more control for weighted aggregations with additional arguments `w`, `keep.w` and `wFUN` to aggregate the weights as well. The defaults are `keep.w = TRUE` and `wFUN = fsum`. A specialty of `collap` remains that `keep.by` and `keep.w` also work for external objects passed, so code of the form `collap(data, by, FUN, catFUN, w = data$weights)` will now have an aggregated `weights` vector in the first column. 

<!-- In such cases use `keep.w = FALSE` to omit the weights or `collap(data, by, FUN, catFUN, w = ~ weights)` to keep the column order. -->

* `qsu` now also allows weights to be passed in formula i.e. `qsu(data, by = ~ group, pid = ~ panelid, w = ~ weights)`. 

* `fgrowth` has a `scale` argument, the default is `scale = 100` which provides growth rates in percentage terms (as before), but this may now be changed. 

* All statistical and transformation functions now have a hidden list method, so they can be applied to unclassed list-objects as well. An error is however provided in grouped operations with unequal-length columns. 



# collapse 1.1.0
collapse 1.1.0 released early April 2020: <!--  - some small fixes and additions: -->

* Fixed remaining gcc10, LTO and valgrind issues in C/C++ code, and added some more tests (there are now ~ 5300 tests ensuring that *collapse* statistical functions perform as expected).

* Fixed the issue that supplying an unnamed list to `GRP()`, i.e. `GRP(list(v1, v2))` would give an error. Unnamed lists are now automatically named 'Group.1', 'Group.2', etc...

* Fixed an issue where aggregating by a single id in `collap()` (i.e. `collap(data, ~ id1)`), the id would be coded as factor in the aggregated data.frame. All variables including id's now retain their class and attributes in the aggregated data.

* Added weights (`w`) argument to `fsum` and `fprod`. 

* Added an argument `mean = 0` to `fwithin / W`. This allows simple and grouped centering on an arbitrary mean, `0` being the default. For grouped centering `mean = "overall.mean"` can be specified, which will center data on the overall mean of the data. The logical argument `add.global.mean = TRUE` used to toggle this in *collapse* 1.0.0 is therefore depreciated. 

* Added arguments `mean = 0` (the default) and `sd = 1` (the default) to `fscale / STD`. These arguments now allow to (group) scale and center data to an arbitrary mean and standard deviation. Setting `mean = FALSE` will just scale data while preserving the mean(s). Special options for grouped scaling are `mean = "overall.mean"` (same as `fwithin / W`), and `sd = "within.sd"`, which will scale the data such that the standard deviation of each group is equal to the within- standard deviation (= the standard deviation computed on the group-centered data). Thus group scaling a panel-dataset with `mean = "overall.mean"` and `sd = "within.sd"` harmonizes the data across all groups in terms of both mean and variance. The fast algorithm for variance calculation toggled with `stable.algo = FALSE` was removed from `fscale`. Welford's numerically stable algorithm used by default is fast enough for all practical purposes. The fast algorithm is still available for `fvar` and `fsd`. 

* Added the modulus (`%%`) and subtract modulus (`-%%`) operations to `TRA()`. 

* Added the function `finteraction`, for fast interactions, and `as_character_factor` to coerce a factor, or all factors in a list, to character (analogous to `as_numeric_factor`). Also exported the function `ckmatch`, for matching with error message showing non-matched elements.


# collapse 1.0.0 and earlier

* First version of the package featuring only the functions `collap` and `qsu` based on code shared by Sebastian Krantz on R-devel, February 2019.

* Major rework of the package using Rcpp and data.table internals, introduction of fast statistical functions and operators and expansion of the scope of the package to a broad set of data transformation and exploration tasks. Several iterations of enhancing speed of R code. Seamless integration of *collapse* with *dplyr*, *plm* and *data.table*. CRAN release of *collapse* 1.0.0 on 19th March 2020. 

