# *collapse* 1.7 News Digest

## C-Level Development

* Multiple assignment `massign()`, and operator `%=%`: `c("T", "K") %=% dim(EuStockMarkets)` or `.c(T, K) %=% dim(EuStockMarkets)`, also `names(iris) %=% iris` etc.. Supports all vector types. 

* Mathematical operations by reference `setop()` and operators `%+=%`, `%-=%`, `%*=%`, `%/=%`, supporting vectors, matrices and data frames / lists. Operators work column-wise (like R's math): `mtcars %+=% mtcars`, `num_vars(iris) %/=% iris$Sepal.Length`, `Titanic %/=% sum(Titanic) * 100`. `setop()` also supports `rowwise` operations: `setop(EuStockMarkets, "-", fmean(EuStockMarkets), rowwise = TRUE)` (centering by reference). 

* Efficient comparisons of vectors with single values: `whichv` ( + operators `%==%` and `%!=%` and `whichNA`), `anyv`, `allv` (and `allNA`). Like `whichv`, `%==%` and `%!=%` return indices and are significantly faster than `==`, e.g. benchmark `fsubset(data, col %==% value)` against `fsubset(data, col == value)`. Can also use inside *data.table*. 

* Efficient replacement and copying of values in atomic objects and data frames (complementing `data.table::set`): `setv(mtcars, 0, 2)` (replace 0 with 2 everywhere), `setv(mtcars$carb, mtcars$cyl %==% 4, mtcars$mpg)` (same as `qDT(mtcars)[cyl == 4, carb := mpg]`, but no need for a *data.table*), `setv(mtcars[5:8], mtcars$cyl %==% 4, df[6:9])` (copy rows from equally sized data frames, here just shifting them by one column), etc. There is also `copyv` to do all this with a deep copy. These functions also benefit `replace_NA` and friends, which now have an additional `set` argument. 

* `group()` provides access to a new multivariate first-appearance-order grouping algorithm (thanks to Morgan Jacob from the great *kit* package for providing the hash function). The algorithm now feeds into the package through multiple central functions (including `GRP` / `fgroup_by`, `collap`, `funique` and `qF`) when invoked with argument `sort = FALSE`. More benchmarks need to be run, but in my experience more often than not the gain in grouping performance is striking. 

* `gsplit()` very efficiently splits vectors based on *collapse* grouping objects. It is used as a new backend to `rsplit` as well as `BY`, `collap`, `fsummarise` and `fmutate`. Overall this means that, assuming serial execution e.g. on M1 Mac, *collapse* is competitive with `dplyr` in operations using base R / user-defined functions (possibly even with *data.table*, but haven't explored much), competitive with *data.table* for aggregations and transformations with the *Fast Statistical Functions*, and remains fastest for weighted aggregation, computations with panel data and other specialized things. 

* `vlabel(<-)` was re-written in C, and user level functions `relabel` and `setrelabel` were added, providing consistent and fast support for dealing with variable labels (and further attributes to be distributed across data frame columns) in R. Functions `setLabels`, `namlab`, `descr` and `qsu` also benefit from this for retrieval of variable labels. 

* Functions `vlengths()` added as faster replacement for `lengths`, and function `vtypes` was also rewritten in C. The latter benefits column selection by data type using `num_vars`, `char_vars`, `fact_vars` etc. which is now fully C-based and super fast, especially for high-dimensional data. 

* `qsu` has an argument `stable.algo = FALSE` which gives 2x performance improvement for big data (sacrificing numeric accurancy of the standard deviation, but I use long doubles in C++ to still maximize numeric precision). 

* More efficient attribute handling in many functions (using shallow copies only). This also benefits the user-level functions `setAttrib`, `copyAttrib` and `copyMostAttrib`, which now take a shallow copy of lists and also a shallow copy of the attributes list. 

* `setrename` now truly renames objects by reference (without creating a shallow copy, except for *data.tables* where the external pointer needs to be renewed).

* `na_rm` can now also efficiently remove empty or `NULL` elements from a list.

* Added option `ties = "last"` to `fmode`. 

## R-Level Development

* Fully fledged `fmutate` function that provides functionality analogous to `dplyr::mutate` (sequential evaluation of arguments, including arbitrary tagged expressions and `across` statements). `fmutate` is optimized to work together with the packages *Fast Statistical and Data Transformation Functions*, yielding fast, vectorized execution, but also benefits from `gsplit` for other operations. 

* `across()` function implemented for use inside `fsummarise` and `fmutate`. It is also optimized for *Fast Statistical and Data Transformation Functions*, but performs well with other functions too. It has an additional arguments `.apply = FALSE` which will apply functions to the entire subset of the data instead of individual columns, and thus allows for nesting tibbles and estimating models or correlation matrices by groups etc.. `across()` also supports an arbitrary number of additional arguments which are split and evaluated by groups if necessary. Multiple `across()` statements can be combined with tagged vector expressions in a single call to `fsummarise` or `fmutate`. 

* `fsummarise` is now also fully featured supporting evaluation of arbitrary expressions and `across()` statements. 

* Added option `options("collapse_mask")`, which can be set before loading the package to export copies of selected (or all) functions in the package that start with `f` removing the leading `f` e.g. `options(collapse_mask = "fsubset")` means both `fsubset` and `subset` will be exported. This allows masking base R and *dplyr* functions (even basic functions such as `sum`, `mean`, `unique` etc. can be masked by *collapse* fast functions if desired). There are some shorthands e.g. `options(collapse_mask = "manip")` exports copies of all data manipulation functions, allowing translation of *dplyr* codes without much change of syntax, i.e. you can write `data %>% group_by(col) %>% summarise(...)` instead of `data %>% fgroup_by(col) %>% fsummarise(...)`. Masking basic statistical functions such as `options(collapse_mask = "all")` (mask everything) will also modify internal macros so that `sum`, `mean` etc. receive vectorized executions inside data manipulation functions. 
    
* `fgroup_by` is more flexible, supporting computing columns e.g. `fgroup_by(GGDC10S, Variable, Decade = floor(Year / 10) * 10)` and various programming options e.g. `fgroup_by(GGDC10S, 1:3)`, `fgroup_by(GGDC10S, c("Variable", "Country"))`, or `fgroup_by(GGDC10S, is.character)`. You can also use column sequences e.g. `fgroup_by(GGDC10S, Country:Variable, Year)`, but this should not be mixed with computing columns. Compute expressions may also not include the `:` function. 

* In `ftransformv/settransformv` and `fcomputev`, the `vars` argument is also evaluated inside the data frame environment, allowing NSE specifications using column names e.g. `ftransformv(data, c(col1, col2:coln), FUN)`.

* Function `ss` supports both empty `i` or `j`. 

* The printout of `fgroup_by` also shows minimum and maximum group size for unbalanced groupings. 

* `qF` with option `sort = FALSE` now generates factors with levels in first-appearance order (instead of a random order assigned by the hash function), and can also be called on an existing factor to recast the levels in first-appearance order. It is also faster with `sort = FALSE` (thanks to `group`). 

* `colorder` can rename columns on the fly and also has a new mode `pos = "after"` to place all selected  columns after the first selected one, e.g.: `colorder(mtcars, cyl, vs_new = vs, am, pos = "after")`. The `pos = "after"` option was also added to `roworderv`. 

+ `add_stub` and `rm_stub` have an additional `cols` argument to apply a stub to certain columns only e.g. `add_stub(mtcars, "new_", cols = 6:9)`.

* `namlab` has additional arguments `N` and `Ndistinct`, allowing to display number of observations and distinct values next to variable names, labels and classes, to get a nice and quick overview of the variables in a large dataset. 

* The print methods of `pwcor` and `pwcov` now have a `return` argument, allowing users to obtain the formatted correlation matrix, for exporting purposes. 

* Some improvements to the `BY` function, both in terms of performance and security. 

* Internally adjusted quite some functions e.g. `fwithin` etc. to take advantage of the faster grouping by `group`, wherever ordered grouping is not required. 


